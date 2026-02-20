from argparse import Namespace

import pytest

from phykit.services.tree.evolutionary_rate import EvolutionaryRate


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
        mocker.patch.object(EvolutionaryRate, "read_tree_file", return_value=_Tree())
        mocked_print = mocker.patch("builtins.print")
        service.run()
        mocked_print.assert_called_once_with(3.0864)

    def test_run_json_output(self, mocker):
        args = Namespace(tree="/some/path/to/file.tre", json=True)
        service = EvolutionaryRate(args)
        mocker.patch.object(EvolutionaryRate, "read_tree_file", return_value=_Tree())
        mocked_json = mocker.patch("phykit.services.tree.evolutionary_rate.print_json")
        service.run()
        mocked_json.assert_called_once_with({"evolutionary_rate": 3.0864})
