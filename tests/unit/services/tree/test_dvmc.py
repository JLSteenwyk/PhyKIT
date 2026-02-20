import pytest
from argparse import Namespace
import numpy as np

from phykit.services.tree.dvmc import DVMC
import phykit.services.tree.dvmc as dvmc_module


@pytest.fixture
def args():
    kwargs = dict(tree="/some/path/to/file.tre")
    return Namespace(**kwargs)


class TestDVMC(object):
    def test_init_sets_tree_file_path(self, args):
        d = DVMC(args)
        assert d.output_file_path is None

    def test_read_file_reads_tree_file_path(self, mocker, args):
        # Mock the cached tree read method instead of Phylo.read
        mock_cached_read = mocker.patch("phykit.services.tree.base.Tree._cached_tree_read")
        mock_get_hash = mocker.patch("phykit.services.tree.base.Tree._get_file_hash", return_value="test_hash")

        d = DVMC(args)
        d.read_tree_file()

        # Verify the cached read was called with the correct parameters
        mock_get_hash.assert_called_with(args.tree)
        mock_cached_read.assert_called_with(args.tree, "newick", "test_hash")

    def test_process_args_defaults_json_false(self):
        parsed = DVMC(Namespace(tree="x.tre")).process_args(Namespace(tree="x.tre"))
        assert parsed["json_output"] is False

    def test_determine_dvmc_matches_numpy_formula(self, tree_simple):
        d = DVMC(Namespace(tree="x.tre"))
        result = d.determine_dvmc(tree_simple)
        distances = np.array([tree_simple.distance(term) for term in tree_simple.get_terminals()])
        expected = np.std(distances, ddof=1)
        assert result == pytest.approx(expected)

    def test_run_json_output(self, tree_simple, monkeypatch):
        captured = {}
        d = DVMC(Namespace(tree="x.tre", json=True))
        monkeypatch.setattr(d, "read_tree_file", lambda: tree_simple)
        monkeypatch.setattr(dvmc_module, "print_json", lambda payload: captured.setdefault("payload", payload))

        d.run()

        assert "dvmc" in captured["payload"]
        assert isinstance(captured["payload"]["dvmc"], float)

    def test_run_terminal_output(self, tree_simple, monkeypatch, capsys):
        d = DVMC(Namespace(tree="x.tre", json=False))
        monkeypatch.setattr(d, "read_tree_file", lambda: tree_simple)

        d.run()

        out, _ = capsys.readouterr()
        assert out.strip() != ""
