from argparse import Namespace
import subprocess
import sys

import pytest

from phykit.services.tree.evolutionary_rate import EvolutionaryRate


def test_module_import_does_not_import_json_or_typing():
    code = """
import sys
import phykit.services.tree.evolutionary_rate as module
assert callable(module.print_json)
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def args():
    return Namespace(tree="/some/path/to/file.tre")


class _Tree:
    def total_branch_length(self):
        return 12.3456

    def count_terminals(self):
        return 4


class TestEvolutionaryRate:
    def test_init_sets_expected_attrs(self, args):
        service = EvolutionaryRate(args)
        assert service.tree_file_path == args.tree
        assert service.json_output is False

    def test_run_prints_rate(self, mocker, args):
        service = EvolutionaryRate(args)
        mocker.patch.object(service, "_get_simple_newick_summary", return_value=None)
        mocker.patch.object(
            EvolutionaryRate,
            "read_tree_file_unmodified",
            return_value=_Tree(),
        )
        mocked_print = mocker.patch("builtins.print")
        service.run()
        mocked_print.assert_called_once_with(3.0864)

    def test_run_uses_unmodified_tree_read(self, mocker, args):
        tree = _Tree()
        service = EvolutionaryRate(args)
        mocker.patch.object(service, "_get_simple_newick_summary", return_value=None)
        read_tree = mocker.patch.object(
            service,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocked_print = mocker.patch("builtins.print")

        service.run()

        read_tree.assert_called_once_with()
        mocked_print.assert_called_once_with(3.0864)

    def test_run_uses_simple_newick_summary_before_tree_read(self, mocker, args):
        service = EvolutionaryRate(args)
        mocker.patch.object(
            service,
            "_get_simple_newick_summary",
            return_value=(("a", "b", "c", "d"), 12.3456, 3.0),
        )
        read_tree = mocker.patch.object(
            service,
            "read_tree_file_unmodified",
            side_effect=AssertionError("simple Newick should skip tree parsing"),
        )
        mocked_print = mocker.patch("builtins.print")

        service.run()

        read_tree.assert_not_called()
        mocked_print.assert_called_once_with(3.0864)

    def test_run_json_output(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", json=True)
        service = EvolutionaryRate(args)
        mocker.patch.object(service, "_get_simple_newick_summary", return_value=None)
        mocker.patch.object(
            EvolutionaryRate,
            "read_tree_file_unmodified",
            return_value=_Tree(),
        )
        mocked_json = mocker.patch("phykit.services.tree.evolutionary_rate.print_json")
        service.run()
        mocked_json.assert_called_once_with({"evolutionary_rate": 3.0864})
