from argparse import Namespace

import pytest

from phykit.services.tree.collapse_branches import CollapseBranches


@pytest.fixture
def args():
    return Namespace(tree="/some/path/to/file.tre", support=90.0, output=None)


class _Tree:
    def __init__(self):
        self._nonterm_count = 5
        self.called = False
        self.predicate_result = None

    def get_nonterminals(self):
        return [object()] * self._nonterm_count

    def collapse_all(self, predicate):
        self.called = True
        class _Clade:
            def __init__(self, confidence):
                self.confidence = confidence
        self.predicate_result = (
            predicate(_Clade(80.0)),
            predicate(_Clade(95.0)),
            predicate(_Clade(None)),
        )
        self._nonterm_count = 3


class TestCollapseBranches:
    def test_init_sets_expected_attrs(self, args):
        service = CollapseBranches(args)
        assert service.tree_file_path == args.tree
        assert service.support == 90.0
        assert service.output_file_path == "/some/path/to/file.tre.collapsed_90.0.tre"
        assert service.json_output is False

    def test_process_args_honors_custom_output(self):
        args = Namespace(tree="/some/path/to/file.tre", support=80.0, output="/tmp/out.tre", json=True)
        service = CollapseBranches(args)
        assert service.output_file_path == "/tmp/out.tre"
        assert service.json_output is True

    def test_run_collapses_and_writes(self, mocker, args):
        service = CollapseBranches(args)
        tree = _Tree()
        mocker.patch.object(CollapseBranches, "read_tree_file", return_value=tree)
        mocker.patch("phykit.services.tree.collapse_branches.copy.deepcopy", return_value=tree)
        mocked_write = mocker.patch.object(CollapseBranches, "write_tree_file")

        service.run()

        assert tree.called is True
        assert tree.predicate_result == (True, False, None)
        mocked_write.assert_called_once_with(tree, "/some/path/to/file.tre.collapsed_90.0.tre")

    def test_run_json_payload(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", support=90.0, output="/tmp/out.tre", json=True)
        service = CollapseBranches(args)
        tree = _Tree()
        mocker.patch.object(CollapseBranches, "read_tree_file", return_value=tree)
        mocker.patch("phykit.services.tree.collapse_branches.copy.deepcopy", return_value=tree)
        mocker.patch.object(CollapseBranches, "write_tree_file")
        mocked_json = mocker.patch("phykit.services.tree.collapse_branches.print_json")

        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload == {
            "input_tree": "/some/path/to/file.tre",
            "support_threshold": 90.0,
            "collapsed_branches": 2,
            "output_file": "/tmp/out.tre",
        }
