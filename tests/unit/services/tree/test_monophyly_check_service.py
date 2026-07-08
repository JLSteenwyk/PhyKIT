from argparse import Namespace
from io import StringIO
import subprocess
import sys

import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree, TreeMixin

from phykit.services.tree.monophyly_check import MonophylyCheck


def test_module_import_does_not_import_biophylo_or_numpy():
    code = """
import sys
import phykit.services.tree.monophyly_check as module

assert callable(module.print_json)
assert callable(module.read_single_column_file_to_list)
assert "json" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.stats_summary" not in sys.modules
assert "phykit.helpers.files" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def args():
    return Namespace(tree="/some/path/to/file.tre", list_of_taxa="/some/path/to/taxa.txt")


class _Node:
    def __init__(self, confidence):
        self.confidence = confidence


class _Clade:
    def get_nonterminals(self):
        return [_Node(90), _Node(80), _Node(None)]


class TestMonophylyCheck:
    def test_init_sets_expected_attrs(self, args):
        service = MonophylyCheck(args)
        assert service.tree_file_path == args.tree
        assert service.list_of_taxa == args.list_of_taxa
        assert service.json_output is False

    def test_get_bootstrap_statistics(self, args):
        service = MonophylyCheck(args)
        stats = service.get_bootstrap_statistics(_Clade())
        assert round(stats["mean"], 4) == 85.0
        assert stats["maximum"] == 90
        assert stats["minimum"] == 80

    def test_get_bootstrap_statistics_uses_direct_traversal(self, args, monkeypatch):
        service = MonophylyCheck(args)
        tree = Phylo.read(StringIO("((A:1,B:1)90:1,(C:1,D:1)80:1)95;"), "newick")

        def fail_get_nonterminals(*_args, **_kwargs):
            raise AssertionError("standard clades should use direct traversal")

        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_get_nonterminals)

        stats = service.get_bootstrap_statistics(tree.root)

        assert round(stats["mean"], 4) == 88.3333
        assert stats["maximum"] == 95
        assert stats["minimum"] == 80

    def test_get_bootstrap_statistics_handles_mixed_child_counts(
        self, args, monkeypatch
    ):
        service = MonophylyCheck(args)
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1)70:1,(D:1,E:1,F:1)90:1)80;"),
            "newick",
        )

        def fail_get_nonterminals(*_args, **_kwargs):
            raise AssertionError("standard clades should use direct traversal")

        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_get_nonterminals)

        stats = service.get_bootstrap_statistics(tree.root)

        assert round(stats["mean"], 4) == 80.0
        assert stats["maximum"] == 90
        assert stats["minimum"] == 70

    def test_populate_res_arr_monophyletic(self, args):
        service = MonophylyCheck(args)
        stats = {"mean": 90.0, "maximum": 100.0, "minimum": 80.0, "standard_deviation": 10.0}
        res = service.populate_res_arr([], stats, [])
        assert res[0][0] == "monophyletic"
        assert res[0][1:5] == [90.0, 100.0, 80.0, 10.0]

    def test_print_results_json(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", list_of_taxa="/some/path/to/taxa.txt", json=True)
        service = MonophylyCheck(args)
        mocked_json = mocker.patch("phykit.services.tree.monophyly_check.print_json")
        service.print_results(
            [
                ["not_monophyletic", 95.0, 100.0, 85.0, 7.0, ["z", "a"]],
                ["insufficient_taxon_representation"],
            ]
        )
        payload = mocked_json.call_args.args[0]
        assert payload["rows"] == payload["results"]
        assert payload["rows"][0]["status"] == "not_monophyletic"
        assert payload["rows"][0]["offending_taxa"] == ["a", "z"]
        assert payload["rows"][1] == {"status": "insufficient_taxon_representation"}

    def test_run_insufficient_taxa_exits(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", list_of_taxa="/some/path/to/taxa.txt")
        service = MonophylyCheck(args)
        mocker.patch.object(MonophylyCheck, "read_tree_file_unmodified", return_value=object())
        mocker.patch(
            "phykit.services.tree.monophyly_check.read_single_column_file_to_list",
            return_value=["only_one_taxon"],
        )
        mocker.patch.object(MonophylyCheck, "get_tip_names_from_tree", return_value=["only_one_taxon"])
        with pytest.raises(SystemExit) as excinfo:
            service.run()
        assert excinfo.value.code == 2

    def test_run_happy_path(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", list_of_taxa="/some/path/to/taxa.txt")
        service = MonophylyCheck(args)

        cached_tree = mocker.Mock()
        root_tree = mocker.Mock()
        clade_tree = mocker.Mock()
        root_tree.common_ancestor.return_value = clade_tree

        mocker.patch.object(service, "read_tree_file_unmodified", return_value=cached_tree)
        mocker.patch.object(service, "_fast_copy", return_value=root_tree)
        mocker.patch(
            "phykit.services.tree.monophyly_check.read_single_column_file_to_list",
            return_value=["a", "b", "c"],
        )
        mocker.patch.object(
            service,
            "get_tip_names_from_tree",
            side_effect=[["a", "b", "c", "d"], ["a", "b", "c"]],
        )
        mocker.patch.object(service, "shared_tips", return_value=["a", "b", "c"])
        mocker.patch.object(
            service,
            "get_bootstrap_statistics",
            return_value={"mean": 90.0, "maximum": 100.0, "minimum": 80.0, "standard_deviation": 10.0},
        )
        mocked_print_results = mocker.patch.object(service, "print_results")

        service.run()

        cached_tree.root_with_outgroup.assert_not_called()
        root_tree.root_with_outgroup.assert_called_once_with(["d"])
        mocked_print_results.assert_called_once()
        res_arr = mocked_print_results.call_args.args[0]
        assert res_arr[0][0] == "monophyletic"

    def test_run_reuses_unmodified_tree(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", list_of_taxa="/some/path/to/taxa.txt")
        service = MonophylyCheck(args)
        tree = object()
        clade = object()
        taxa = frozenset({"a", "b"})

        mocker.patch.object(
            service,
            "read_tree_file",
            side_effect=AssertionError("run should not copy the cached tree"),
        )
        read_unmodified = mocker.patch.object(
            service, "read_tree_file_unmodified", return_value=tree
        )
        mocker.patch(
            "phykit.services.tree.monophyly_check.read_single_column_file_to_list",
            return_value=["a", "b"],
        )
        mocker.patch.object(service, "get_tip_names_from_tree", return_value=["a", "b", "c"])
        resolve = mocker.patch.object(service, "_resolve_interest_clade", return_value=(clade, taxa))
        mocker.patch.object(
            service,
            "get_bootstrap_statistics",
            return_value={"mean": 90.0, "maximum": 100.0, "minimum": 80.0, "standard_deviation": 10.0},
        )
        mocker.patch.object(service, "print_results")

        service.run()

        read_unmodified.assert_called_once_with()
        resolve.assert_called_once_with(tree, taxa, frozenset({"a", "b", "c"}))

    def test_resolve_interest_clade_copies_before_reroot(self, mocker, args):
        service = MonophylyCheck(args)
        cached_tree = mocker.Mock()
        copied_tree = mocker.Mock()
        clade = mocker.Mock()
        copied_tree.common_ancestor.return_value = clade

        mocker.patch.object(service, "_find_exact_clade_by_taxa", return_value=None)
        mocker.patch.object(service, "shared_tips", return_value=["a", "b"])
        copy_tree = mocker.patch.object(service, "_fast_copy", return_value=copied_tree)
        mocker.patch.object(service, "get_tip_names_from_tree", return_value=["a", "b"])

        resolved, tips = service._resolve_interest_clade(
            cached_tree,
            frozenset({"a", "b"}),
            frozenset({"a", "b", "c"}),
        )

        copy_tree.assert_called_once_with(cached_tree)
        cached_tree.root_with_outgroup.assert_not_called()
        copied_tree.root_with_outgroup.assert_called_once_with(["c"])
        assert sorted(copied_tree.common_ancestor.call_args.args[0]) == ["a", "b"]
        assert resolved is clade
        assert tips == frozenset({"a", "b"})

    def test_resolve_interest_clade_skips_shared_tip_recalculation(
        self, mocker, args
    ):
        service = MonophylyCheck(args)
        cached_tree = mocker.Mock()
        copied_tree = mocker.Mock()
        clade = mocker.Mock()
        copied_tree.common_ancestor.return_value = clade

        mocker.patch.object(service, "_find_exact_clade_by_taxa", return_value=None)
        mocker.patch.object(
            service,
            "shared_tips",
            side_effect=AssertionError("taxa already intersect tree tips"),
        )
        mocker.patch.object(service, "_fast_copy", return_value=copied_tree)
        mocker.patch.object(service, "get_tip_names_from_tree", return_value=["a", "b"])

        resolved, tips = service._resolve_interest_clade(
            cached_tree,
            frozenset({"a", "b"}),
            frozenset({"a", "b", "c"}),
        )

        copied_tree.root_with_outgroup.assert_called_once_with(["c"])
        assert sorted(copied_tree.common_ancestor.call_args.args[0]) == ["a", "b"]
        assert resolved is clade
        assert tips == frozenset({"a", "b"})

    def test_resolve_interest_clade_skips_reroot_for_exact_clade(self, monkeypatch, args):
        service = MonophylyCheck(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_root_with_outgroup(*_args, **_kwargs):
            raise AssertionError("root_with_outgroup should not be called")

        def fail_common_ancestor(*_args, **_kwargs):
            raise AssertionError("common_ancestor should not be called")

        monkeypatch.setattr(tree, "root_with_outgroup", fail_root_with_outgroup)
        monkeypatch.setattr(tree, "common_ancestor", fail_common_ancestor)

        clade, tips = service._resolve_interest_clade(
            tree, frozenset(["A", "B"]), frozenset(["A", "B", "C", "D"])
        )

        assert sorted(t.name for t in clade.get_terminals()) == ["A", "B"]
        assert tips == frozenset(["A", "B"])

    def test_find_exact_clade_uses_count_path(self, monkeypatch, args):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        expected = tree.root.clades[0]

        def fail_is_terminal(*_args, **_kwargs):
            raise AssertionError("count path should not call is_terminal")

        monkeypatch.setattr(Clade, "is_terminal", fail_is_terminal)

        clade = MonophylyCheck._find_exact_clade_by_taxa(tree, frozenset(["A", "B"]))

        assert clade is expected

    def test_find_exact_clade_uses_direct_traversal(self, monkeypatch, args):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        expected = tree.root.clades[0]

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree path should not use find_clades")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        clade = MonophylyCheck._find_exact_clade_by_taxa(tree, frozenset(["A", "B"]))

        assert clade is expected

    def test_find_exact_clade_direct_handles_mixed_child_counts(self, monkeypatch, args):
        tree = Phylo.read(StringIO("(((A:1):1,(B:1,C:1,D:1):1):1,E:1);"), "newick")
        expected = tree.root.clades[0]

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree path should not use find_clades")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        clade = MonophylyCheck._find_exact_clade_by_taxa(
            tree,
            frozenset(["A", "B", "C", "D"]),
        )

        assert clade is expected

    def test_find_exact_clade_direct_handles_deep_ladder(self, monkeypatch, args):
        root = Clade()
        current = root
        taxa = []
        for index in range(1200):
            name = f"T{index}"
            taxa.append(name)
            next_internal = Clade()
            current.clades = [Clade(name=name), next_internal]
            current = next_internal
        current.clades = [Clade(name="T1200")]
        taxa.append("T1200")
        tree = Tree(root=root)

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree path should not use find_clades")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        clade = MonophylyCheck._find_exact_clade_by_taxa(tree, frozenset(taxa))

        assert clade is root

    def test_find_exact_clade_duplicate_names_use_set_fallback(self, args):
        tree = Phylo.read(StringIO("((A:1,A:1):1,B:1);"), "newick")

        clade = MonophylyCheck._find_exact_clade_by_taxa(tree, frozenset(["A", "B"]))

        assert clade is tree.root

    def test_print_results_text_paths(self, mocker, args):
        service = MonophylyCheck(args)
        service.json_output = False
        mocked_print = mocker.patch("builtins.print")

        service.print_results(
            [
                ["not_monophyletic", 95.0, 100.0, 85.0, 7.0, ["z", "a"]],
                ["monophyletic", 92.0, 100.0, 90.0, 3.0, []],
                ["insufficient_taxon_representation"],
            ]
        )

        mocked_print.assert_called_once_with(
            "not_monophyletic\t95.0\t100.0\t85.0\t7.0\ta;z\n"
            "monophyletic\t92.0\t100.0\t90.0\t3.0\n"
            "insufficient_taxon_representation"
        )
