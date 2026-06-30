from argparse import Namespace
from io import StringIO
import subprocess
import sys

import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

from phykit.errors import PhykitUserError
from phykit.services.tree.base import Tree
import phykit.services.tree.consensus_tree as consensus_tree_module
from phykit.services.tree.consensus_tree import ConsensusTree


def _write(path, content):
    path.write_text(content)
    return str(path)


def test_module_import_does_not_import_biophylo_or_numpy():
    code = """
import sys
import phykit.services.tree.consensus_tree as module

assert hasattr(module.Phylo, "read")
assert hasattr(module.Consensus, "strict_consensus")
assert hasattr(module.Consensus, "majority_consensus")
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


class TestConsensusTree:
    def test_process_args_defaults_json_false(self):
        args = Namespace(trees="trees.txt", method="majority", missing_taxa="error")
        parsed = ConsensusTree(args).process_args(args)
        assert parsed["json_output"] is False
        assert parsed["method"] == "majority"

    def test_run_json_with_newick_lines(self, tmp_path, monkeypatch):
        tree_file = tmp_path / "trees.nwk"
        _write(tree_file, "((A,B),(C,D));\n((A,B),(C,D));\n")

        captured = {}
        svc = ConsensusTree(
            Namespace(
                trees=str(tree_file),
                method="majority",
                missing_taxa="error",
                json=True,
            )
        )
        monkeypatch.setattr(
            "phykit.services.tree.consensus_tree.print_json",
            lambda payload: captured.setdefault("payload", payload),
        )

        svc.run()

        payload = captured["payload"]
        assert payload["method"] == "majority"
        assert payload["taxa_in_consensus"] == 4
        assert payload["pruned_to_shared_taxa"] is False
        assert "A" in payload["tree_newick"] and "D" in payload["tree_newick"]

    def test_parse_trees_skips_comments_and_blank_lines(self, tmp_path):
        tree_file = tmp_path / "trees.nwk"
        _write(
            tree_file,
            "# ignored\n\n  ((A,B),(C,D));  \n# also ignored\n((A,B),(C,D));\n",
        )
        svc = ConsensusTree(
            Namespace(
                trees=str(tree_file),
                method="majority",
                missing_taxa="error",
                json=False,
            )
        )

        trees = svc._parse_trees_from_source(str(tree_file))

        assert len(trees) == 2
        assert [sorted(ConsensusTree._tips(tree)) for tree in trees] == [
            ["A", "B", "C", "D"],
            ["A", "B", "C", "D"],
        ]

    def test_parse_tree_path_list_avoids_per_row_path_objects(self, tmp_path, mocker):
        tree_a = tmp_path / "a.nwk"
        tree_b = tmp_path / "b.nwk"
        _write(tree_a, "((A,B),(C,D));\n")
        _write(tree_b, "((A,C),(B,D));\n")
        tree_list = tmp_path / "tree-list.txt"
        _write(tree_list, "a.nwk\nb.nwk\n")
        path_calls = 0
        original_path = __import__(
            "phykit.services.tree.consensus_tree",
            fromlist=["Path"],
        ).Path

        def counting_path(*args, **kwargs):
            nonlocal path_calls
            path_calls += 1
            return original_path(*args, **kwargs)

        mocker.patch(
            "phykit.services.tree.consensus_tree.Path",
            side_effect=counting_path,
        )
        svc = ConsensusTree(
            Namespace(
                trees=str(tree_list),
                method="majority",
                missing_taxa="error",
                json=False,
            )
        )

        trees = svc._parse_trees_from_source(str(tree_list))

        assert len(trees) == 2
        assert path_calls == 1

    def test_missing_taxa_error_mode_raises(self, tmp_path):
        tree_file = tmp_path / "trees.nwk"
        _write(tree_file, "((A,B),(C,D));\n((A,B),(C,E));\n")

        svc = ConsensusTree(
            Namespace(
                trees=str(tree_file),
                method="majority",
                missing_taxa="error",
                json=False,
            )
        )

        with pytest.raises(PhykitUserError):
            svc.run()

    def test_missing_taxa_shared_mode_prunes(self, tmp_path, monkeypatch):
        tree_file = tmp_path / "trees.nwk"
        _write(tree_file, "((A,B),(C,D));\n((A,B),(C,E));\n")

        captured = {}
        svc = ConsensusTree(
            Namespace(
                trees=str(tree_file),
                method="majority",
                missing_taxa="shared",
                json=True,
            )
        )
        monkeypatch.setattr(
            "phykit.services.tree.consensus_tree.print_json",
            lambda payload: captured.setdefault("payload", payload),
        )

        svc.run()

        payload = captured["payload"]
        assert payload["pruned_to_shared_taxa"] is True
        assert payload["taxa_in_consensus"] == 3
        assert "D" not in payload["tree_newick"]
        assert "E" not in payload["tree_newick"]

    def test_prune_to_taxa_avoids_generic_terminal_collection(self, monkeypatch):
        tree = Phylo.read(StringIO("((A,B),(C,D));"), "newick")
        original_get_terminals = type(tree).get_terminals
        calls = {"count": 0}

        def counted_get_terminals(self):
            calls["count"] += 1
            return original_get_terminals(self)

        monkeypatch.setattr(type(tree), "get_terminals", counted_get_terminals)

        ConsensusTree._prune_to_taxa(tree, {"A", "C"})

        assert calls["count"] == 0
        assert Tree.calculate_terminal_names_fast(tree) == ["A", "C"]

    def test_prune_to_taxa_uses_batch_prune_for_standard_tree(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:2):3,(C:4,D:5):6);"), "newick")

        def fail_prune(*_args, **_kwargs):
            raise AssertionError("per-target prune should not be used")

        monkeypatch.setattr(TreeMixin, "prune", fail_prune)

        ConsensusTree._prune_to_taxa(tree, {"B", "D"})

        assert [tip.name for tip in tree.get_terminals()] == ["B", "D"]
        assert [tip.branch_length for tip in tree.get_terminals()] == [5.0, 11.0]

    def test_prune_to_taxa_uses_direct_terminal_collection(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:2):3,(C:4,D:5):6);"), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("standard tree should use direct terminal traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        ConsensusTree._prune_to_taxa(tree, {"B", "D"})

        assert ConsensusTree._tips(tree) == {"B", "D"}

    def test_tips_uses_fast_terminal_names_for_parsed_tree(self, monkeypatch):
        tree = Phylo.read(StringIO("((A,B),(C,D));"), "newick")

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("get_terminals fallback should not be called")

        monkeypatch.setattr(tree, "get_terminals", fail_get_terminals)

        assert ConsensusTree._tips(tree) == {"A", "B", "C", "D"}

    def test_normalize_taxa_identical_sets_uses_fast_path(self, monkeypatch):
        svc = ConsensusTree(
            Namespace(
                trees="unused",
                method="majority",
                missing_taxa="error",
                json=False,
            )
        )
        trees = [
            Phylo.read(StringIO("((A,B),(C,D));"), "newick"),
            Phylo.read(StringIO("((A,C),(B,D));"), "newick"),
        ]

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("parsed trees should use fast terminal names")

        def fail_prune(*_args, **_kwargs):
            raise AssertionError("identical taxon sets should not be pruned")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(ConsensusTree, "_prune_to_taxa", fail_prune)

        normalized, pruned, taxa_count = svc._normalize_taxa(trees)

        assert normalized is trees
        assert pruned is False
        assert taxa_count == 4

    def test_all_tip_sets_identical_does_not_slice_rows(self):
        class NoSliceList(list):
            def __getitem__(self, key):
                if isinstance(key, slice):
                    raise AssertionError("tip-set scan should not slice")
                return super().__getitem__(key)

        tip_sets = NoSliceList([
            {"A", "B", "C"},
            {"A", "B", "C"},
            {"A", "B", "C"},
        ])

        assert consensus_tree_module._all_tip_sets_identical(tip_sets) is True

    def test_path_list_input(self, tmp_path, monkeypatch):
        t1 = tmp_path / "t1.tre"
        t2 = tmp_path / "t2.tre"
        list_file = tmp_path / "trees.txt"
        _write(t1, "((A,B),(C,D));\n")
        _write(t2, "((A,B),(C,D));\n")
        _write(list_file, f"{t1.name}\n{t2.name}\n")

        printed = {}
        svc = ConsensusTree(
            Namespace(
                trees=str(list_file),
                method="strict",
                missing_taxa="error",
                json=False,
            )
        )
        monkeypatch.setattr("builtins.print", lambda x: printed.setdefault("v", x))

        svc.run()

        assert "A" in printed["v"] and "D" in printed["v"]
