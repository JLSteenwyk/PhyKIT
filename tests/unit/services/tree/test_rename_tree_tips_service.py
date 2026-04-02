from argparse import Namespace

import pytest

from phykit.services.tree.rename_tree_tips import RenameTreeTips


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

    def test_run_writes_and_emits_json(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            idmap="/some/path/to/idmap.txt",
            output="/tmp/out.tre",
            json=True,
        )
        tree = _Tree(["a", "b"])
        mocker.patch.object(RenameTreeTips, "read_tree_file", return_value=tree)
        mocker.patch("phykit.services.tree.rename_tree_tips.pickle.loads", return_value=tree)
        mocker.patch.object(RenameTreeTips, "read_id_map", return_value={"a": "A"})
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
