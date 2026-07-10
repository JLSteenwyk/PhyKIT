from argparse import Namespace
from io import StringIO
import subprocess
import sys

from Bio import Phylo
from Bio.Phylo.Newick import Clade, Tree
from Bio.Phylo.BaseTree import TreeMixin
import pytest

from phykit.services.tree.collapse_branches import CollapseBranches


def test_module_import_defers_typing_and_json_helper():
    code = (
        "import sys; "
        "import phykit.services.tree.collapse_branches as module; "
        "assert callable(module.print_json); "
        "assert 'typing' not in sys.modules; "
        "assert 'json' not in sys.modules; "
        "assert 'phykit.helpers.json_output' not in sys.modules"
    )

    subprocess.run([sys.executable, "-c", code], check=True)


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


class _TreeWithoutCounts(_Tree):
    def get_nonterminals(self):
        raise AssertionError("non-JSON run should not count internodes")


class _TreeWithoutStandardTraversal(_Tree):
    @property
    def root(self):
        raise AttributeError


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
        mocker.patch.object(CollapseBranches, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(service, "_fast_copy", return_value=tree)
        mocked_write = mocker.patch.object(CollapseBranches, "write_tree_file")

        service.run()

        assert tree.called is True
        assert tree.predicate_result == (True, False, None)
        mocked_write.assert_called_once_with(tree, "/some/path/to/file.tre.collapsed_90.0.tre")

    def test_run_without_json_skips_internode_counts(self, mocker, args):
        service = CollapseBranches(args)
        tree = _TreeWithoutCounts()
        mocker.patch.object(CollapseBranches, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(service, "_fast_copy", return_value=tree)
        mocked_write = mocker.patch.object(CollapseBranches, "write_tree_file")

        service.run()

        assert tree.called is True
        mocked_write.assert_called_once_with(tree, "/some/path/to/file.tre.collapsed_90.0.tre")

    def test_run_skips_collapse_all_when_standard_tree_has_no_weak_support(
        self, mocker, args
    ):
        service = CollapseBranches(args)
        tree = Tree(
            root=Clade(
                confidence=100.0,
                clades=[
                    Clade(confidence=95.0, clades=[Clade(name="A"), Clade(name="B")]),
                    Clade(confidence=97.0, clades=[Clade(name="C"), Clade(name="D")]),
                ],
            )
        )

        def fail_collapse_all(*_args, **_kwargs):
            raise AssertionError("collapse_all should not be called")

        tree.collapse_all = fail_collapse_all
        mocker.patch.object(CollapseBranches, "read_tree_file_unmodified", return_value=tree)
        mocked_copy = mocker.spy(service, "_fast_copy")
        mocked_write = mocker.patch.object(CollapseBranches, "write_tree_file")

        service.run()

        mocked_copy.assert_not_called()
        mocked_write.assert_called_once_with(tree, "/some/path/to/file.tre.collapsed_90.0.tre")

    def test_run_copies_standard_tree_before_collapsing_weak_support(
        self, mocker, args
    ):
        service = CollapseBranches(args)
        tree = Tree(
            root=Clade(
                confidence=100.0,
                clades=[
                    Clade(
                        branch_length=1.0,
                        confidence=80.0,
                        clades=[
                            Clade(branch_length=1.0, name="A"),
                            Clade(branch_length=1.0, name="B"),
                        ],
                    ),
                    Clade(
                        branch_length=1.0,
                        confidence=97.0,
                        clades=[
                            Clade(branch_length=1.0, name="C"),
                            Clade(branch_length=1.0, name="D"),
                        ],
                    ),
                ],
            )
        )
        mocker.patch.object(CollapseBranches, "read_tree_file_unmodified", return_value=tree)
        mocked_copy = mocker.spy(service, "_fast_copy")
        mocked_write = mocker.patch.object(CollapseBranches, "write_tree_file")

        service.run()

        mocked_copy.assert_called_once_with(tree)
        written_tree = mocked_write.call_args.args[0]
        assert written_tree is not tree
        assert len(tree.root.clades) == 2

    def test_run_non_json_standard_tree_uses_early_boolean_scan(
        self, mocker, args
    ):
        service = CollapseBranches(args)
        tree = Tree(
            root=Clade(
                confidence=80.0,
                clades=[
                    Clade(confidence=95.0, clades=[Clade(name="A"), Clade(name="B")]),
                    Clade(confidence=97.0, clades=[Clade(name="C"), Clade(name="D")]),
                ],
            )
        )

        def fail_combined_scan(*_args, **_kwargs):
            raise AssertionError("non-JSON runs only need an early boolean scan")

        mocker.patch.object(CollapseBranches, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(
            service,
            "_scan_standard_tree_for_collapse",
            side_effect=fail_combined_scan,
        )
        mocked_copy = mocker.spy(service, "_fast_copy")
        mocked_write = mocker.patch.object(CollapseBranches, "write_tree_file")

        service.run()

        mocked_copy.assert_called_once_with(tree)
        written_tree = mocked_write.call_args.args[0]
        assert written_tree is not tree
        assert len(tree.root.clades) == 2

    def test_run_json_noop_uses_combined_standard_tree_scan(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", support=90.0, output="/tmp/out.tre", json=True)
        service = CollapseBranches(args)
        tree = Tree(
            root=Clade(
                confidence=100.0,
                clades=[
                    Clade(confidence=95.0, clades=[Clade(name="A"), Clade(name="B")]),
                    Clade(confidence=97.0, clades=[Clade(name="C"), Clade(name="D")]),
                ],
            )
        )

        def fail_count(*_args, **_kwargs):
            raise AssertionError("combined scan should provide the initial count")

        def fail_collapse_all(*_args, **_kwargs):
            raise AssertionError("collapse_all should not be called")

        tree.collapse_all = fail_collapse_all
        mocker.patch.object(CollapseBranches, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(CollapseBranches, "write_tree_file")
        mocker.patch.object(service, "count_internal_nodes", side_effect=fail_count)
        mocked_copy = mocker.spy(service, "_fast_copy")
        mocked_json = mocker.patch("phykit.services.tree.collapse_branches.print_json")

        service.run()

        mocked_copy.assert_not_called()
        payload = mocked_json.call_args.args[0]
        assert payload["collapsed_branches"] == 0

    def test_run_falls_back_to_collapse_all_for_nonstandard_tree(self, mocker, args):
        service = CollapseBranches(args)
        tree = _TreeWithoutStandardTraversal()
        mocker.patch.object(CollapseBranches, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(service, "_fast_copy", return_value=tree)
        mocker.patch.object(CollapseBranches, "write_tree_file")

        service.run()

        assert tree.called is True

    def test_count_internal_nodes_uses_direct_traversal(self, args, monkeypatch):
        service = CollapseBranches(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_get_nonterminals(*_args, **_kwargs):
            raise AssertionError("generic internal-node traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_get_nonterminals)

        assert service.count_internal_nodes(tree) == 3

    def test_count_internal_nodes_direct_traversal_handles_multifurcations(
        self, args, monkeypatch
    ):
        service = CollapseBranches(args)
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(name="A"),
                    Clade(clades=[Clade(name="B"), Clade(name="C")]),
                    Clade(clades=[Clade(name="D"), Clade(name="E")]),
                ],
            )
        )

        def fail_get_nonterminals(*_args, **_kwargs):
            raise AssertionError("generic internal-node traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_get_nonterminals)

        assert service.count_internal_nodes(tree) == 3

    def test_scan_standard_tree_for_collapse_counts_and_flags_weak_support(self):
        tree = Tree(
            root=Clade(
                confidence=100.0,
                clades=[
                    Clade(confidence=95.0, clades=[Clade(name="A"), Clade(name="B")]),
                    Clade(confidence=75.0, clades=[Clade(name="C"), Clade(name="D")]),
                ],
            )
        )

        assert CollapseBranches._scan_standard_tree_for_collapse(tree, 90.0) == (
            3,
            True,
        )

    def test_run_json_payload(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", support=90.0, output="/tmp/out.tre", json=True)
        service = CollapseBranches(args)
        tree = _Tree()
        mocker.patch.object(CollapseBranches, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(service, "_fast_copy", return_value=tree)
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
