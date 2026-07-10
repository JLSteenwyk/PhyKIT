from argparse import Namespace
from io import StringIO
import subprocess
import sys

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
import pytest

from phykit.services.tree.branch_length_multiplier import BranchLengthMultiplier


def test_lightweight_tree_output_modules_defer_heavy_imports():
    modules = [
        "phykit.services.tree.branch_length_multiplier",
        "phykit.services.tree.total_tree_length",
        "phykit.services.tree.dvmc",
        "phykit.services.tree.internode_labeler",
    ]
    code = (
        "import sys; "
        f"modules = {modules!r}; "
        "__import__(modules[0]); "
        "assert 'typing' not in sys.modules; "
        "[__import__(module) for module in modules]; "
        "assert 'json' not in sys.modules; "
        "assert 'phykit.helpers.json_output' not in sys.modules; "
        "assert 'Bio.Phylo' not in sys.modules; "
        "assert 'numpy' not in sys.modules"
    )

    subprocess.run([sys.executable, "-c", code], check=True)


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

    def find_clades(self, order=None):
        return iter(self.nonterms + self.terms)


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

    def test_multiply_standard_tree_uses_direct_traversal(self, args, monkeypatch):
        service = BranchLengthMultiplier(args)
        tree = Phylo.read(StringIO("((A:2,B:3):5,C:7);"), "newick")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        scaled_count = service.multiply_branch_lengths_by_factor(tree, 2.0)

        assert scaled_count == 4
        assert tree.root.branch_length is None
        assert tree.root.clades[0].branch_length == 10.0
        assert tree.root.clades[0].clades[0].branch_length == 4.0
        assert tree.root.clades[0].clades[1].branch_length == 6.0
        assert tree.root.clades[1].branch_length == 14.0

    def test_standard_tree_helper_scales_root_branch_when_present(self, args):
        service = BranchLengthMultiplier(args)
        tree = Phylo.read(StringIO("((A:2,B:3):5,C:7):11;"), "newick")

        scaled_count = service._multiply_standard_tree_branch_lengths(tree, 3.0)

        assert scaled_count == 5
        assert tree.root.branch_length == 33.0
        assert tree.root.clades[0].branch_length == 15.0
        assert tree.root.clades[0].clades[0].branch_length == 6.0
        assert tree.root.clades[0].clades[1].branch_length == 9.0
        assert tree.root.clades[1].branch_length == 21.0

    def test_standard_tree_helper_scales_mixed_child_counts(self, args):
        service = BranchLengthMultiplier(args)

        class Clade:
            def __init__(self, branch_length=None, clades=None):
                self.branch_length = branch_length
                self.clades = clades or []

        tree = type(
            "Tree",
            (),
            {
                "root": Clade(
                    None,
                    [
                        Clade(1.0),
                        Clade(2.0, [Clade(3.0)]),
                        Clade(None, [Clade(4.0), Clade(5.0), Clade(None)]),
                    ],
                )
            },
        )()

        scaled_count = service._multiply_standard_tree_branch_lengths(tree, 2.0)

        assert scaled_count == 5
        assert tree.root.clades[0].branch_length == 2.0
        assert tree.root.clades[1].branch_length == 4.0
        assert tree.root.clades[1].clades[0].branch_length == 6.0
        assert tree.root.clades[2].branch_length is None
        assert tree.root.clades[2].clades[0].branch_length == 8.0
        assert tree.root.clades[2].clades[1].branch_length == 10.0
        assert tree.root.clades[2].clades[2].branch_length is None

    def test_count_branch_lengths_uses_direct_traversal(self, args, monkeypatch):
        service = BranchLengthMultiplier(args)
        tree = Phylo.read(StringIO("((A:2,B:3):5,C:7);"), "newick")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        assert service.count_branch_lengths(tree) == 4

    def test_standard_tree_count_helper_counts_direct_clade_lists(self, args):
        service = BranchLengthMultiplier(args)

        class Clade:
            def __init__(self, branch_length=None, clades=None):
                self.branch_length = branch_length
                self.clades = clades or []

        tree = type(
            "Tree",
            (),
            {
                "root": Clade(
                    None,
                    [
                        Clade(1.0),
                        Clade(2.0, [Clade(None), Clade(3.0)]),
                    ],
                )
            },
        )()

        count = service._count_standard_tree_branch_lengths(tree)

        assert count == 3

    def test_run_factor_one_uses_read_only_tree_without_scaling(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=1.0, output="/tmp/out.tre")
        service = BranchLengthMultiplier(args)
        tree = _Tree()
        read_tree = mocker.patch.object(
            service,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(
            service,
            "read_tree_file",
            side_effect=AssertionError("factor-one path should not copy the tree"),
        )
        scale = mocker.patch.object(
            service,
            "multiply_branch_lengths_by_factor",
            side_effect=AssertionError("factor-one path should not mutate branches"),
        )
        mocked_write = mocker.patch.object(BranchLengthMultiplier, "write_tree_file")

        service.run()

        read_tree.assert_called_once_with()
        scale.assert_not_called()
        mocked_write.assert_called_once_with(tree, "/tmp/out.tre")
        assert tree.nonterms[0].branch_length == 1.0
        assert tree.terms[0].branch_length == 2.0

    def test_run_json_factor_one_counts_branch_lengths_without_scaling(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            factor=1.0,
            output="/tmp/out.tre",
            json=True,
        )
        service = BranchLengthMultiplier(args)
        tree = _Tree()
        mocker.patch.object(service, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(
            service,
            "multiply_branch_lengths_by_factor",
            side_effect=AssertionError("factor-one path should not mutate branches"),
        )
        mocker.patch.object(BranchLengthMultiplier, "write_tree_file")
        mocked_json = mocker.patch("phykit.services.tree.branch_length_multiplier.print_json")

        service.run()

        payload = mocked_json.call_args.args[0]
        assert payload["scaled_branches"] == 2
        assert tree.nonterms[0].branch_length == 1.0
        assert tree.terms[0].branch_length == 2.0

    def test_run_json_payload(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", factor=2.0, output="/tmp/out.tre", json=True)
        service = BranchLengthMultiplier(args)
        tree = _Tree()
        mocker.patch.object(BranchLengthMultiplier, "read_tree_file", return_value=tree)
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
