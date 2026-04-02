from argparse import Namespace

import pytest

from phykit.services.tree.prune_tree import PruneTree


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

    def test_run_prunes_given_taxa_and_emits_json(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/out.tre",
            keep=False,
            json=True,
        )
        tree = _Tree(["a", "b", "c"])

        mocker.patch.object(PruneTree, "read_tree_file", return_value=tree)
        mocker.patch("phykit.services.tree.prune_tree.pickle.loads", return_value=tree)
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

    def test_run_keep_mode_prunes_complement_of_given_taxa(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/out.tre",
            keep=True,
            json=False,
        )
        tree = _Tree(["a", "b", "c"])

        mocker.patch.object(PruneTree, "read_tree_file", return_value=tree)
        mocker.patch("phykit.services.tree.prune_tree.pickle.loads", return_value=tree)
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
