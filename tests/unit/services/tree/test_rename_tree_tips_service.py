from argparse import Namespace
from io import StringIO
import subprocess
import sys

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
import pytest

from phykit.services.tree.rename_tree_tips import RenameTreeTips


def test_module_import_does_not_import_json_or_typing():
    code = """
import sys
import phykit.services.tree.rename_tree_tips as module
assert callable(module.print_json)
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def args():
    return Namespace(
        tree="/some/path/to/file.tre",
        idmap="/some/path/to/idmap.txt",
        output=None,
    )


class _Tip:
    def __init__(self, name):
        self.name = name


class _Tree:
    def __init__(self, tips):
        self._tips = [_Tip(tip) for tip in tips]

    def get_terminals(self):
        return self._tips


class TestRenameTreeTips:
    def test_init_sets_expected_attrs(self, args):
        service = RenameTreeTips(args)
        assert service.tree_file_path == args.tree
        assert service.idmap == args.idmap
        assert service.output_file_path == f"{args.tree}.renamed"
        assert service.json_output is False

    def test_process_args_honors_custom_output(self):
        args = Namespace(
            tree="/some/path/to/file.tre",
            idmap="/some/path/to/idmap.txt",
            output="/tmp/out.tre",
            json=True,
        )
        service = RenameTreeTips(args)
        assert service.output_file_path == "/tmp/out.tre"
        assert service.json_output is True

    def test_read_id_map_parses_two_column_file(self, tmp_path):
        idmap_file = tmp_path / "idmap.txt"
        idmap_file.write_text("a A\nb B\n")
        args = Namespace(tree="/x.tre", idmap=str(idmap_file), output=None)
        service = RenameTreeTips(args)
        assert service.read_id_map() == {"a": "A", "b": "B"}

    def test_read_id_map_parses_whitespace_tokens_from_bytes(self, tmp_path):
        idmap_file = tmp_path / "idmap.txt"
        idmap_file.write_bytes(b"a\tA\nb   B\na A2\n")
        args = Namespace(tree="/x.tre", idmap=str(idmap_file), output=None)
        service = RenameTreeTips(args)

        assert service.read_id_map() == {"a": "A2", "b": "B"}

    def test_read_id_map_missing_file_exits(self):
        args = Namespace(tree="/x.tre", idmap="/does/not/exist.txt", output=None)
        service = RenameTreeTips(args)
        with pytest.raises(SystemExit) as excinfo:
            service.read_id_map()
        assert excinfo.value.code == 2

    def test_replace_tip_names_renames_and_counts(self, args):
        service = RenameTreeTips(args)
        tree = _Tree(["a", "b", "c"])
        renamed_tree, renamed_count = service.replace_tip_names(tree, {"a": "A", "c": "C"})
        assert renamed_count == 2
        assert [tip.name for tip in renamed_tree.get_terminals()] == ["A", "b", "C"]

    def test_replace_tip_names_standard_tree_uses_direct_traversal(self, args, monkeypatch):
        service = RenameTreeTips(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,c:1);"), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("generic terminal traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        renamed_tree, renamed_count = service.replace_tip_names(
            tree,
            {"a": "A", "c": "C"},
        )

        assert renamed_tree is tree
        assert renamed_count == 2
        assert tree.root.clades[0].clades[0].name == "A"
        assert tree.root.clades[0].clades[1].name == "b"
        assert tree.root.clades[1].name == "C"

    def test_replace_tip_names_standard_tree_preserves_falsey_mapped_names(self, args):
        service = RenameTreeTips(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,c:1);"), "newick")

        renamed_tree, renamed_count = service.replace_tip_names(
            tree,
            {"a": "", "c": None},
        )

        assert renamed_tree is tree
        assert renamed_count == 2
        assert tree.root.clades[0].clades[0].name == ""
        assert tree.root.clades[1].name is None

    def test_count_matching_tip_names_standard_tree_uses_direct_traversal(
        self, args, monkeypatch
    ):
        service = RenameTreeTips(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,c:1);"), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("generic terminal traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        assert service.count_matching_tip_names(tree, {"a": "A", "c": "C"}) == 2

    def test_standard_tree_helpers_handle_mixed_child_counts(self, args, monkeypatch):
        service = RenameTreeTips(args)
        newick = "(a:1,(b:1,c:1):1,(d:1,e:1,f:1):1);"
        idmap = {"a": "A", "c": "C", "f": ""}

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("standard tree helper should use direct traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        count_tree = Phylo.read(StringIO(newick), "newick")
        assert service.count_matching_tip_names(count_tree, idmap) == 3
        assert service.has_matching_tip_name(count_tree, idmap) is True
        assert service.has_matching_tip_name(count_tree, {"missing": "M"}) is False

        rename_tree = Phylo.read(StringIO(newick), "newick")
        renamed_tree, renamed_count = service.replace_tip_names(rename_tree, idmap)

        assert renamed_tree is rename_tree
        assert renamed_count == 3
        assert rename_tree.root.clades[0].name == "A"
        assert rename_tree.root.clades[1].clades[0].name == "b"
        assert rename_tree.root.clades[1].clades[1].name == "C"
        assert rename_tree.root.clades[2].clades[2].name == ""

    def test_empty_idmap_skips_tree_traversal(self, args, monkeypatch):
        service = RenameTreeTips(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,c:1);"), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("empty id map should not traverse terminals")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        assert service.has_matching_tip_name(tree, {}) is False
        assert service.count_matching_tip_names(tree, {}) == 0

        renamed_tree, renamed_count = service.replace_tip_names(tree, {})

        assert renamed_tree is tree
        assert renamed_count == 0
        assert tree.root.clades[0].clades[0].name == "a"
        assert tree.root.clades[0].clades[1].name == "b"
        assert tree.root.clades[1].name == "c"

    def test_replace_standard_tree_tip_names_does_not_require_reversed_children(
        self, args
    ):
        service = RenameTreeTips(args)

        class NoReverseList(list):
            def __reversed__(self):
                raise AssertionError("renaming does not require child order")

        class Clade:
            def __init__(self, name=None, clades=None):
                self.name = name
                self.clades = NoReverseList(clades or [])

        tree = type(
            "Tree",
            (),
            {
                "root": Clade(
                    clades=[
                        Clade("a"),
                        Clade(clades=[Clade("b"), Clade("c")]),
                    ],
                )
            },
        )()

        renamed_tree, renamed_count = service._replace_standard_tree_tip_names(
            tree,
            {"a": "A", "c": "C"},
        )

        assert renamed_tree is tree
        assert renamed_count == 2
        assert tree.root.clades[0].name == "A"
        assert tree.root.clades[1].clades[0].name == "b"
        assert tree.root.clades[1].clades[1].name == "C"

    def test_run_no_matching_tips_uses_read_only_tree_without_copy(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            idmap="/some/path/to/idmap.txt",
            output="/tmp/out.tre",
            json=True,
        )
        service = RenameTreeTips(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,c:1);"), "newick")
        mocker.patch.object(service, "read_id_map", return_value={"x": "X"})
        mocker.patch.object(service, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(
            service,
            "read_tree_file",
            side_effect=AssertionError("no-match path should not copy via reader"),
        )
        mocked_copy = mocker.spy(service, "_fast_copy")
        mocked_write = mocker.patch.object(RenameTreeTips, "write_tree_file")
        mocked_json = mocker.patch("phykit.services.tree.rename_tree_tips.print_json")

        service.run()

        mocked_copy.assert_not_called()
        mocked_write.assert_called_once_with(tree, "/tmp/out.tre")
        assert [term.name for term in tree.get_terminals()] == ["a", "b", "c"]
        assert mocked_json.call_args.args[0]["renamed_tips"] == 0

    def test_run_matching_tips_copies_before_renaming(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            idmap="/some/path/to/idmap.txt",
            output="/tmp/out.tre",
            json=False,
        )
        service = RenameTreeTips(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,c:1);"), "newick")
        mocker.patch.object(service, "read_id_map", return_value={"a": "A"})
        mocker.patch.object(service, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(
            service,
            "count_matching_tip_names",
            side_effect=AssertionError("non-JSON run should only check for any match"),
        )
        mocked_copy = mocker.spy(service, "_fast_copy")
        mocked_write = mocker.patch.object(RenameTreeTips, "write_tree_file")

        service.run()

        mocked_copy.assert_called_once_with(tree)
        written_tree = mocked_write.call_args.args[0]
        assert written_tree is not tree
        assert [term.name for term in tree.get_terminals()] == ["a", "b", "c"]
        assert [term.name for term in written_tree.get_terminals()] == ["A", "b", "c"]

    def test_run_json_matching_tips_uses_boolean_precheck(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            idmap="/some/path/to/idmap.txt",
            output="/tmp/out.tre",
            json=True,
        )
        service = RenameTreeTips(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,c:1);"), "newick")
        mocker.patch.object(service, "read_id_map", return_value={"a": "A", "c": "C"})
        mocker.patch.object(service, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(
            service,
            "count_matching_tip_names",
            side_effect=AssertionError("JSON precheck only needs a boolean match"),
        )
        mocked_write = mocker.patch.object(RenameTreeTips, "write_tree_file")
        mocked_json = mocker.patch("phykit.services.tree.rename_tree_tips.print_json")

        service.run()

        written_tree = mocked_write.call_args.args[0]
        assert [term.name for term in written_tree.get_terminals()] == ["A", "b", "C"]
        assert mocked_json.call_args.args[0]["renamed_tips"] == 2

    def test_run_writes_and_emits_json(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            idmap="/some/path/to/idmap.txt",
            output="/tmp/out.tre",
            json=True,
        )
        tree = _Tree(["a", "b"])
        mocker.patch.object(RenameTreeTips, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(RenameTreeTips, "read_id_map", return_value={"a": "A"})
        mocker.patch.object(RenameTreeTips, "_fast_copy", return_value=tree)
        mocked_write = mocker.patch.object(RenameTreeTips, "write_tree_file")
        mocked_json = mocker.patch("phykit.services.tree.rename_tree_tips.print_json")

        service = RenameTreeTips(args)
        service.run()

        assert [tip.name for tip in tree.get_terminals()] == ["A", "b"]
        mocked_write.assert_called_once_with(tree, "/tmp/out.tre")
        payload = mocked_json.call_args.args[0]
        assert payload["input_tree"] == "/some/path/to/file.tre"
        assert payload["idmap"] == "/some/path/to/idmap.txt"
        assert payload["output_file"] == "/tmp/out.tre"
        assert payload["renamed_tips"] == 1
