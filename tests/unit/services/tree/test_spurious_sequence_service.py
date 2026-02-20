from argparse import Namespace

import pytest

from phykit.services.tree.spurious_sequence import SpuriousSequence


@pytest.fixture
def args():
    return Namespace(tree="/some/path/to/file.tre", factor=None)


class _Terminal:
    def __init__(self, name, branch_length):
        self.name = name
        self.branch_length = branch_length


class _Tree:
    def __init__(self, terminals):
        self._terminals = terminals

    def get_terminals(self):
        return self._terminals


class TestSpuriousSequence:
    def test_init_sets_expected_attrs(self, args):
        service = SpuriousSequence(args)
        assert service.tree_file_path == args.tree
        assert service.factor == 20
        assert service.json_output is False

    def test_get_branch_lengths_and_names_uses_terminal_branches_only(self, args):
        service = SpuriousSequence(args)
        tree = _Tree(
            [
                _Terminal("a", 1.0),
                _Terminal("b", None),
                _Terminal("c", 3.0),
            ]
        )
        branch_lengths, name_map = service.get_branch_lengths_and_their_names(tree)
        assert branch_lengths == [1.0, 3.0]
        assert name_map == {"a": 1.0, "c": 3.0}

    def test_identify_spurious_sequence(self, args):
        service = SpuriousSequence(args)
        tree = _Tree([_Terminal("a", 1.0), _Terminal("b", 3.0), _Terminal("c", 5.0)])
        name_map, threshold, median = service.identify_spurious_sequence(tree, factor=2)
        assert median == 3.0
        assert threshold == 6.0
        assert name_map["c"] == 5.0

    def test_run_prints_none_when_no_spurious_sequences(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=20, json=False)
        service = SpuriousSequence(args)
        mocker.patch.object(SpuriousSequence, "read_tree_file", return_value=_Tree([]))
        mocker.patch.object(
            SpuriousSequence,
            "identify_spurious_sequence",
            return_value=({"a": 1.0, "b": 2.0}, 100.0, 1.5),
        )
        mocked_print = mocker.patch("builtins.print")
        service.run()
        mocked_print.assert_called_once_with("None")

    def test_run_prints_spurious_sequences(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=2, json=False)
        service = SpuriousSequence(args)
        mocker.patch.object(SpuriousSequence, "read_tree_file", return_value=_Tree([]))
        mocker.patch.object(
            SpuriousSequence,
            "identify_spurious_sequence",
            return_value=({"a": 10.0, "b": 1.0}, 5.0, 2.5),
        )
        mocked_print = mocker.patch("builtins.print")
        service.run()
        mocked_print.assert_called_once_with("a\t10.0\t5.0\t2.5")

    def test_run_json_payload(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=2, json=True)
        service = SpuriousSequence(args)
        mocker.patch.object(SpuriousSequence, "read_tree_file", return_value=_Tree([]))
        mocker.patch.object(
            SpuriousSequence,
            "identify_spurious_sequence",
            return_value=({"a": 10.12345, "b": 1.0}, 5.0, 2.5),
        )
        mocked_json = mocker.patch("phykit.services.tree.spurious_sequence.print_json")
        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload["rows"] == payload["spurious_sequences"]
        assert payload["rows"][0] == {
            "taxon": "a",
            "branch_length": 10.1235,
            "threshold": 5.0,
            "median": 2.5,
        }
