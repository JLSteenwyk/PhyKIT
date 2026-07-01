from argparse import Namespace
from io import StringIO
import subprocess
import sys

import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

import phykit.services.tree.prune_tree as prune_module
from phykit.services.tree.prune_tree import PruneTree, _strip_branch_label


def test_module_import_does_not_import_json_helpers_or_typing():
    code = """
import sys
before = set(sys.modules)
import phykit.services.tree.prune_tree as module
imported = set(sys.modules) - before

assert callable(module.print_json)
assert callable(module.read_single_column_file_to_list)
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.files" not in sys.modules
assert "typing" not in sys.modules
assert "re" not in imported
assert module._BRANCH_LABEL_RE is None
assert module._BRANCH_LABEL_SUB is None
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def args():
    return Namespace(
        tree="/some/path/to/file.tre",
        list_of_taxa="/some/path/to/taxa.txt",
        output=None,
        keep=None,
    )


class _Tip:
    def __init__(self, name):
        self.name = name


class _Tree:
    def __init__(self, tips):
        self._tips = tips

    def get_terminals(self):
        return [_Tip(tip) for tip in self._tips]

    def count_terminals(self):
        return len(self._tips)


class TestPruneTree:
    def test_init_sets_expected_attrs_and_defaults(self, args):
        service = PruneTree(args)
        assert service.tree_file_path == args.tree
        assert service.list_of_taxa == args.list_of_taxa
        assert service.output_file_path == f"{args.tree}.pruned"
        assert service.keep is True
        assert service.json_output is False

    def test_process_args_honors_explicit_keep_false(self):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/out.tre",
            keep=False,
            json=True,
        )
        service = PruneTree(args)
        assert service.output_file_path == "/tmp/out.tre"
        assert service.keep is False
        assert service.json_output is True

    def test_strip_branch_label_skips_regex_for_plain_names(self, monkeypatch):
        class FailRegex:
            def sub(self, *_args, **_kwargs):
                raise AssertionError("plain tip names should not run branch-label regex")

        monkeypatch.setattr(prune_module, "_BRANCH_LABEL_RE", FailRegex())

        assert _strip_branch_label("plain_tip") == "plain_tip"
        assert _strip_branch_label("") == ""
        assert _strip_branch_label(None) is None

    def test_strip_branch_label_removes_labeled_segments(self):
        assert _strip_branch_label("tip{FG}") == "tip"
        assert _strip_branch_label("tip{FG}{BG}") == "tip"

    def test_run_prunes_given_taxa_and_emits_json(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/out.tre",
            keep=False,
            json=True,
        )
        tree = _Tree(["a", "b", "c"])

        mocker.patch.object(PruneTree, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(PruneTree, "_fast_copy", return_value=tree)
        mocker.patch(
            "phykit.services.tree.prune_tree.read_single_column_file_to_list",
            return_value=["a", "b"],
        )
        mocked_prune = mocker.patch.object(PruneTree, "prune_tree_using_taxa_list", return_value=tree)
        mocked_write = mocker.patch.object(PruneTree, "write_tree_file")
        mocked_json = mocker.patch("phykit.services.tree.prune_tree.print_json")

        service = PruneTree(args)
        service.run()

        mocked_prune.assert_called_once_with(tree, ["a", "b"])
        mocked_write.assert_called_once_with(tree, "/tmp/out.tre")
        payload = mocked_json.call_args.args[0]
        assert payload["keep_input_taxa"] is False
        assert payload["taxa_pruned"] == ["a", "b"]
        assert payload["pruned_count"] == 2
        assert payload["remaining_tips"] == 3

    def test_run_json_uses_direct_remaining_tip_count(self, mocker, monkeypatch):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/out.tre",
            keep=False,
            json=True,
        )
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        mocker.patch.object(PruneTree, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(PruneTree, "_fast_copy", return_value=tree)
        mocker.patch(
            "phykit.services.tree.prune_tree.read_single_column_file_to_list",
            return_value=["A"],
        )
        mocker.patch.object(
            PruneTree,
            "prune_tree_using_taxa_list",
            return_value=tree,
        )
        mocker.patch.object(PruneTree, "write_tree_file")
        mocked_json = mocker.patch("phykit.services.tree.prune_tree.print_json")

        def fail_count_terminals(*_args, **_kwargs):
            raise AssertionError("JSON count should use direct terminal count")

        monkeypatch.setattr(TreeMixin, "count_terminals", fail_count_terminals)

        service = PruneTree(args)
        service.run()

        assert mocked_json.call_args.args[0]["remaining_tips"] == 4

    def test_run_keep_mode_prunes_complement_of_given_taxa(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/out.tre",
            keep=True,
            json=False,
        )
        tree = _Tree(["a", "b", "c"])

        mocker.patch.object(PruneTree, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(PruneTree, "_fast_copy", return_value=tree)
        mocker.patch(
            "phykit.services.tree.prune_tree.read_single_column_file_to_list",
            return_value=["a", "b"],
        )
        mocked_prune = mocker.patch.object(PruneTree, "prune_tree_using_taxa_list", return_value=tree)
        mocked_write = mocker.patch.object(PruneTree, "write_tree_file")

        service = PruneTree(args)
        service.run()

        mocked_prune.assert_called_once_with(tree, ["c"])
        mocked_write.assert_called_once_with(tree, "/tmp/out.tre")

    def test_run_keep_mode_uses_set_membership(self, mocker):
        class NoContainsList(list):
            def __contains__(self, item):
                raise AssertionError("keep-mode membership should use a set")

        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/out.tre",
            keep=True,
            json=False,
        )
        tree = _Tree(["a", "b", "c"])

        mocker.patch.object(PruneTree, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(PruneTree, "_fast_copy", return_value=tree)
        mocker.patch(
            "phykit.services.tree.prune_tree.read_single_column_file_to_list",
            return_value=NoContainsList(["a", "b"]),
        )
        mocked_prune = mocker.patch.object(PruneTree, "prune_tree_using_taxa_list", return_value=tree)
        mocker.patch.object(PruneTree, "write_tree_file")

        service = PruneTree(args)
        service.run()

        mocked_prune.assert_called_once_with(tree, ["c"])

    def test_run_keep_mode_uses_fast_tip_names_for_parsed_tree(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/out.tre",
            keep=True,
            json=False,
        )
        tree = Phylo.read(StringIO("((a:1,b:1):1,c:2);"), "newick")

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("tip-name selection should not call get_terminals")

        mocker.patch.object(PruneTree, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(PruneTree, "_fast_copy", return_value=tree)
        mocker.patch(
            "phykit.services.tree.prune_tree.read_single_column_file_to_list",
            return_value=["a", "b"],
        )
        mocker.patch.object(tree, "get_terminals", side_effect=fail_get_terminals)
        mocked_prune = mocker.patch.object(PruneTree, "prune_tree_using_taxa_list", return_value=tree)
        mocker.patch.object(PruneTree, "write_tree_file")

        service = PruneTree(args)
        service.run()

        mocked_prune.assert_called_once_with(tree, ["c"])

    def test_run_ignore_branch_labels_uses_fast_tip_names_for_parsed_tree(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/out.tre",
            keep=True,
            json=False,
            ignore_branch_labels=True,
        )
        tree = Phylo.read(StringIO("((a{FG}:1,b:1):1,c{BG}:2);"), "newick")

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("tip-name selection should not call get_terminals")

        mocker.patch.object(PruneTree, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(PruneTree, "_fast_copy", return_value=tree)
        mocker.patch(
            "phykit.services.tree.prune_tree.read_single_column_file_to_list",
            return_value=["a", "b"],
        )
        mocker.patch.object(tree, "get_terminals", side_effect=fail_get_terminals)
        mocked_prune = mocker.patch.object(PruneTree, "prune_tree_using_taxa_list", return_value=tree)
        mocker.patch.object(PruneTree, "write_tree_file")

        service = PruneTree(args)
        service.run()

        mocked_prune.assert_called_once_with(tree, ["c{BG}"])

    def test_run_keep_all_uses_read_only_tree_without_copy_or_prune(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/out.tre",
            keep=True,
            json=True,
        )
        tree = Phylo.read(StringIO("((a:1,b:1):1,c:2);"), "newick")

        mocker.patch.object(PruneTree, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(
            PruneTree,
            "read_tree_file",
            side_effect=AssertionError("empty-prune path should not copy via reader"),
        )
        mocker.patch(
            "phykit.services.tree.prune_tree.read_single_column_file_to_list",
            return_value=["a", "b", "c"],
        )
        mocked_prune = mocker.patch.object(
            PruneTree,
            "prune_tree_using_taxa_list",
            side_effect=AssertionError("empty-prune path should not prune"),
        )
        mocked_write = mocker.patch.object(PruneTree, "write_tree_file")
        mocked_json = mocker.patch("phykit.services.tree.prune_tree.print_json")

        service = PruneTree(args)
        mocked_copy = mocker.spy(service, "_fast_copy")
        service.run()

        mocked_copy.assert_not_called()
        mocked_prune.assert_not_called()
        mocked_write.assert_called_once_with(tree, "/tmp/out.tre")
        payload = mocked_json.call_args.args[0]
        assert payload["taxa_pruned"] == []
        assert payload["pruned_count"] == 0
        assert payload["remaining_tips"] == 3

    def test_run_copies_before_pruning_cached_tree(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/out.tre",
            keep=False,
            json=False,
        )
        tree = Phylo.read(StringIO("((a:1,b:1):1,c:2);"), "newick")

        mocker.patch.object(PruneTree, "read_tree_file_unmodified", return_value=tree)
        mocker.patch(
            "phykit.services.tree.prune_tree.read_single_column_file_to_list",
            return_value=["a"],
        )
        mocked_write = mocker.patch.object(PruneTree, "write_tree_file")

        service = PruneTree(args)
        mocked_copy = mocker.spy(service, "_fast_copy")
        service.run()

        mocked_copy.assert_called_once_with(tree)
        written_tree = mocked_write.call_args.args[0]
        assert written_tree is not tree
        assert [tip.name for tip in tree.get_terminals()] == ["a", "b", "c"]
        assert [tip.name for tip in written_tree.get_terminals()] == ["b", "c"]
