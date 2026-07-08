from argparse import Namespace
from io import StringIO
import pickle as stdlib_pickle
import subprocess
import sys

import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

import phykit.services.tree.nearest_neighbor_interchange as nni_module
from phykit.services.tree.nearest_neighbor_interchange import (
    NearestNeighborInterchange,
    _split_branch_fields,
)


@pytest.fixture
def args():
    return Namespace(tree="/some/path/to/file.tre", output=None)


def test_module_import_does_not_import_biophylo_or_numpy():
    code = """
import sys
import phykit.services.tree.nearest_neighbor_interchange as module

assert hasattr(module.Phylo, "write")
assert "typing" not in sys.modules
assert "pickle" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_phylo_caches_resolved_write(monkeypatch):
    calls = []

    def cached_write(*args, **kwargs):
        calls.append((args, kwargs))
        return 1

    def uncached_write(*_args, **_kwargs):
        return 2

    lazy_phylo = nni_module._LazyPhylo()
    monkeypatch.setattr(Phylo, "write", cached_write)

    assert lazy_phylo.write([], "trees.nwk", "newick") == 1

    monkeypatch.setattr(Phylo, "write", uncached_write)

    assert lazy_phylo.write(["tree"], "more.nwk", "newick") == 1
    assert lazy_phylo.__dict__["write"] is cached_write
    assert calls == [
        (([], "trees.nwk", "newick"), {}),
        ((["tree"], "more.nwk", "newick"), {}),
    ]


def test_lazy_pickle_caches_resolved_copy_helpers(monkeypatch):
    calls = []

    def cached_dumps(*args, **kwargs):
        calls.append(("dumps", args, kwargs))
        return args[0]

    def cached_loads(*args, **kwargs):
        calls.append(("loads", args, kwargs))
        return args[0]

    def uncached_dumps(*_args, **_kwargs):
        return "uncached-dumps"

    def uncached_loads(*_args, **_kwargs):
        return "uncached-loads"

    lazy_pickle = nni_module._LazyPickle()
    monkeypatch.setattr(stdlib_pickle, "dumps", cached_dumps)
    monkeypatch.setattr(stdlib_pickle, "loads", cached_loads)

    protocol = lazy_pickle.HIGHEST_PROTOCOL
    assert lazy_pickle.loads(lazy_pickle.dumps("tree", protocol=protocol)) == "tree"

    monkeypatch.setattr(stdlib_pickle, "dumps", uncached_dumps)
    monkeypatch.setattr(stdlib_pickle, "loads", uncached_loads)

    assert lazy_pickle.loads(lazy_pickle.dumps("tree2", protocol=protocol)) == "tree2"
    assert lazy_pickle.__dict__["dumps"] is cached_dumps
    assert lazy_pickle.__dict__["loads"] is cached_loads
    assert lazy_pickle.__dict__["HIGHEST_PROTOCOL"] == stdlib_pickle.HIGHEST_PROTOCOL
    assert calls == [
        ("dumps", ("tree",), {"protocol": protocol}),
        ("loads", ("tree",), {}),
        ("dumps", ("tree2",), {"protocol": protocol}),
        ("loads", ("tree2",), {}),
    ]


class TestNearestNeighborInterchange:
    def test_init_sets_expected_attrs(self, args):
        service = NearestNeighborInterchange(args)
        assert service.tree_file_path == args.tree
        assert service.output_file_path == f"{args.tree}.nnis"
        assert service.json_output is False

    def test_process_args_honors_custom_output(self):
        args = Namespace(tree="/some/path/to/file.tre", output="/tmp/nnis.tre", json=True)
        service = NearestNeighborInterchange(args)
        assert service.output_file_path == "/tmp/nnis.tre"
        assert service.json_output is True

    def test_split_branch_fields_handles_commas_tabs_and_whitespace(self):
        assert _split_branch_fields("A,B") == ["A", "B"]
        assert _split_branch_fields("label\tA\tB") == ["label", "A", "B"]
        assert _split_branch_fields(" label , A\tB ") == ["label", "A", "B"]

    def test_resolve_branch_specs_parses_file_without_regex(self, args, tmp_path):
        branches = tmp_path / "branches.tsv"
        branches.write_text(
            "# ignored\n"
            "A,B\n"
            "label\tC\tD\n"
            " label2 , E\tF \n"
        )
        args.branches = str(branches)
        service = NearestNeighborInterchange(args)

        assert service._resolve_branch_specs() == [
            ("A|B", ["A", "B"]),
            ("label", ["C", "D"]),
            ("label2", ["E", "F"]),
        ]

    def test_fast_tree_copy_returns_distinct_object(self, args):
        service = NearestNeighborInterchange(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        copied = service._fast_tree_copy(tree)
        assert copied is not tree
        assert copied.format("newick") == tree.format("newick")

    def test_build_parent_map_uses_single_traversal(self, args, monkeypatch):
        service = NearestNeighborInterchange(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("standard parent map should not use find_clades")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        parents = service._build_parent_map(tree)

        assert len(parents) == 6
        assert parents[tree.root.clades[0]] is tree.root
        assert parents[tree.root.clades[0].clades[0]] is tree.root.clades[0]

    def test_build_parent_map_handles_mixed_child_counts(self, args, monkeypatch):
        service = NearestNeighborInterchange(args)
        tree = Phylo.read(StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"), "newick")

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("standard parent map should not use find_clades")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        parents = service._build_parent_map(tree)
        binary = tree.root.clades[1]
        trifurcation = tree.root.clades[2]

        assert parents[tree.root.clades[0]] is tree.root
        assert parents[binary] is tree.root
        assert parents[trifurcation] is tree.root
        assert parents[binary.clades[0]] is binary
        assert parents[binary.clades[1]] is binary
        assert parents[trifurcation.clades[0]] is trifurcation
        assert parents[trifurcation.clades[1]] is trifurcation
        assert parents[trifurcation.clades[2]] is trifurcation

    def test_nonterminals_level_order_handles_mixed_child_counts(self, args, monkeypatch):
        tree = Phylo.read(StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"), "newick")

        def fail_get_nonterminals(*args, **kwargs):
            raise AssertionError("standard nonterminal scan should be direct")

        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_get_nonterminals)

        nonterminals = NearestNeighborInterchange._nonterminals_level_order_direct(
            tree
        )

        assert nonterminals == [
            tree.root,
            tree.root.clades[1],
            tree.root.clades[2],
        ]

    def test_get_neighbors_uses_direct_level_order_nonterminals(
        self, args, monkeypatch
    ):
        service = NearestNeighborInterchange(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        service._fast_tree_copy = lambda tree: tree

        def fail_get_nonterminals(*args, **kwargs):
            raise AssertionError("standard NNI generation should not use get_nonterminals")

        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_get_nonterminals)

        neighbors = service.get_neighbors(tree)

        assert len(neighbors) == 2

    def test_get_neighbors_returns_list(self, args):
        service = NearestNeighborInterchange(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        neighbors = service.get_neighbors(tree)
        assert isinstance(neighbors, list)
        assert len(neighbors) >= 1

    def test_get_neighbors_exercises_right_child_branch(self, args):
        service = NearestNeighborInterchange(args)
        # This topology has a non-root internal clade as parent's right child,
        # which exercises the alternate branch in get_neighbors.
        tree = Phylo.read(StringIO("((A:1,(B:1,C:1):1):1,(D:1,E:1):1);"), "newick")
        neighbors = service.get_neighbors(tree)
        assert len(neighbors) >= 1
        assert all(hasattr(n, "root") for n in neighbors)

    def test_get_neighbors_exercises_left_child_branch(self, args):
        service = NearestNeighborInterchange(args)
        # This topology has a non-root internal clade as parent's left child.
        tree = Phylo.read(StringIO("((((A:1,B:1):1,C:1):1,D:1):1,(E:1,F:1):1);"), "newick")
        neighbors = service.get_neighbors(tree)
        assert len(neighbors) >= 1

    def test_generate_targeted_nnis_uses_direct_tip_name_lookup(self, args, monkeypatch):
        service = NearestNeighborInterchange(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("targeted NNI setup should use fast terminal names")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        output_trees, report = service._generate_targeted_nnis(
            tree,
            [("A|B", ["A", "B"])],
        )

        assert len(output_trees) == 2
        assert report == [{"label": "A|B", "taxa": ["A", "B"], "n_nnis": 2}]

    def test_run_json_payload(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", output="/tmp/nnis.tre", json=True)
        service = NearestNeighborInterchange(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        mocker.patch.object(NearestNeighborInterchange, "read_tree_file", return_value=tree)
        mocker.patch.object(NearestNeighborInterchange, "get_neighbors", return_value=[tree, tree])
        mocked_write = mocker.patch("phykit.services.tree.nearest_neighbor_interchange.Phylo.write")
        mocked_json = mocker.patch("phykit.services.tree.nearest_neighbor_interchange.print_json")

        service.run()

        mocked_write.assert_called_once()
        payload = mocked_json.call_args.args[0]
        assert payload == {
            "input_tree": "/some/path/to/file.tre",
            "total_trees": 3,
            "nni_neighbors": 2,
            "output_file": "/tmp/nnis.tre",
        }
