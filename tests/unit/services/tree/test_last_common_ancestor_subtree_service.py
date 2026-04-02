from argparse import Namespace

import pytest

from phykit.services.tree.last_common_ancestor_subtree import LastCommonAncestorSubtree


@pytest.fixture
def args():
    return Namespace(
        tree="/some/path/to/file.tre",
        list_of_taxa="/some/path/to/taxa.txt",
        output=None,
    )


class _Subtree:
    def count_terminals(self):
        return 3


class _Tree:
    def common_ancestor(self, taxa):
        assert taxa == ["a", "b"]
        return _Subtree()


class TestLastCommonAncestorSubtree:
    def test_init_sets_expected_attrs(self, args):
        service = LastCommonAncestorSubtree(args)
        assert service.tree_file_path == args.tree
        assert service.list_of_taxa == args.list_of_taxa
        assert service.output_file_path == f"{args.tree}.subtree.tre"
        assert service.json_output is False

    def test_process_args_honors_output_and_json(self):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/subtree.tre",
            json=True,
        )
        service = LastCommonAncestorSubtree(args)
        assert service.output_file_path == "/tmp/subtree.tre"
        assert service.json_output is True

    def test_run_writes_subtree_and_json(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            list_of_taxa="/some/path/to/taxa.txt",
            output="/tmp/subtree.tre",
            json=True,
        )
        service = LastCommonAncestorSubtree(args)
        tree = _Tree()
        mocker.patch.object(LastCommonAncestorSubtree, "read_tree_file", return_value=tree)
        mocker.patch("phykit.services.tree.last_common_ancestor_subtree.pickle.loads", return_value=tree)
        mocker.patch(
            "phykit.services.tree.last_common_ancestor_subtree.read_single_column_file_to_list",
            return_value=["a", "b"],
        )
        mocked_write = mocker.patch.object(LastCommonAncestorSubtree, "write_tree_file")
        mocked_json = mocker.patch("phykit.services.tree.last_common_ancestor_subtree.print_json")

        service.run()

        mocked_write.assert_called_once()
        payload = mocked_json.call_args.args[0]
        assert payload["taxa_count"] == 2
        assert payload["subtree_tips"] == 3
        assert payload["output_file"] == "/tmp/subtree.tre"
