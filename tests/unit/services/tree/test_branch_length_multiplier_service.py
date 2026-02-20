from argparse import Namespace

import pytest

from phykit.services.tree.branch_length_multiplier import BranchLengthMultiplier


@pytest.fixture
def args():
    return Namespace(tree="/some/path/to/file.tre", factor=2.0, output=None)


class _Node:
    def __init__(self, branch_length):
        self.branch_length = branch_length


class _Tree:
    def __init__(self):
        self.nonterms = [_Node(1.0), _Node(None)]
        self.terms = [_Node(2.0)]

    def get_nonterminals(self):
        return self.nonterms

    def get_terminals(self):
        return self.terms


class TestBranchLengthMultiplier:
    def test_init_sets_expected_attrs(self, args):
        service = BranchLengthMultiplier(args)
        assert service.tree_file_path == args.tree
        assert service.factor == 2.0
        assert service.output_file_path == "/some/path/to/file.tre.factor_2.0.tre"
        assert service.json_output is False

    def test_process_args_honors_custom_output(self):
        args = Namespace(tree="/some/path/to/file.tre", factor=3.0, output="/tmp/out.tre", json=True)
        service = BranchLengthMultiplier(args)
        assert service.output_file_path == "/tmp/out.tre"
        assert service.json_output is True

    def test_multiply_branch_lengths_by_factor(self, args):
        service = BranchLengthMultiplier(args)
        tree = _Tree()
        scaled_count = service.multiply_branch_lengths_by_factor(tree, 3.0)
        assert scaled_count == 2
        assert tree.nonterms[0].branch_length == 3.0
        assert tree.nonterms[1].branch_length is None
        assert tree.terms[0].branch_length == 6.0

    def test_run_json_payload(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=2.0, output="/tmp/out.tre", json=True)
        service = BranchLengthMultiplier(args)
        tree = _Tree()
        mocker.patch.object(BranchLengthMultiplier, "read_tree_file", return_value=tree)
        mocker.patch("phykit.services.tree.branch_length_multiplier.copy.deepcopy", return_value=tree)
        mocker.patch.object(BranchLengthMultiplier, "multiply_branch_lengths_by_factor", return_value=2)
        mocked_write = mocker.patch.object(BranchLengthMultiplier, "write_tree_file")
        mocked_json = mocker.patch("phykit.services.tree.branch_length_multiplier.print_json")

        service.run()
        mocked_write.assert_called_once_with(tree, "/tmp/out.tre")
        payload = mocked_json.call_args.args[0]
        assert payload == {
            "input_tree": "/some/path/to/file.tre",
            "factor": 2.0,
            "scaled_branches": 2,
            "output_file": "/tmp/out.tre",
        }
