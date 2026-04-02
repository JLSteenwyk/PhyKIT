from argparse import Namespace

import pytest

from phykit.services.tree.internode_labeler import InternodeLabeler


@pytest.fixture
def args():
    return Namespace(tree="/some/path/to/file.tre", output=None)


class _Node:
    def __init__(self):
        self.confidence = None


class _Tree:
    def __init__(self):
        self.nonterms = [_Node(), _Node(), _Node()]

    def get_nonterminals(self):
        return self.nonterms


class TestInternodeLabeler:
    def test_init_sets_expected_attrs(self, args):
        service = InternodeLabeler(args)
        assert service.tree_file_path == args.tree
        assert service.output_file_path == "/some/path/to/file.tre.internode_labels.tre"
        assert service.json_output is False

    def test_process_args_honors_custom_output(self):
        args = Namespace(tree="/some/path/to/file.tre", output="/tmp/out.tre", json=True)
        service = InternodeLabeler(args)
        assert service.output_file_path == "/tmp/out.tre"
        assert service.json_output is True

    def test_add_labels_to_tree(self, args):
        service = InternodeLabeler(args)
        tree = _Tree()
        count = service.add_labels_to_tree(tree)
        assert count == 3
        assert [node.confidence for node in tree.get_nonterminals()] == [1, 2, 3]

    def test_run_json_payload(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", output="/tmp/out.tre", json=True)
        service = InternodeLabeler(args)
        tree = _Tree()
        mocker.patch.object(InternodeLabeler, "read_tree_file", return_value=tree)
        mocker.patch("phykit.services.tree.internode_labeler.pickle.loads", return_value=tree)
        mocked_write = mocker.patch.object(InternodeLabeler, "write_tree_file")
        mocked_json = mocker.patch("phykit.services.tree.internode_labeler.print_json")

        service.run()

        mocked_write.assert_called_once_with(tree, "/tmp/out.tre")
        payload = mocked_json.call_args.args[0]
        assert payload == {
            "input_tree": "/some/path/to/file.tre",
            "labeled_internodes": 3,
            "output_file": "/tmp/out.tre",
        }
