import pytest
import builtins
import subprocess
import sys
from argparse import Namespace
from io import StringIO
import numpy as np
from Bio import Phylo

from phykit.services.tree.dvmc import DVMC
import phykit.services.tree.dvmc as dvmc_module


def test_module_import_defers_typing_numpy_json_and_biophylo():
    code = """
import sys
import phykit.services.tree.dvmc as module

assert callable(module.print_json)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Phylo" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


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

    def test_determine_dvmc_fast_path_does_not_call_distance(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:2):0.5;"), "newick")
        d = DVMC(Namespace(tree="x.tre"))

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("terminal fallback should not be called")

        def fail_depths(*args, **kwargs):
            raise AssertionError("depths fallback should not be called")

        def fail_distance(*args, **kwargs):
            raise AssertionError("distance fallback should not be called")

        monkeypatch.setattr(tree, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(tree, "depths", fail_depths)
        monkeypatch.setattr(tree, "distance", fail_distance)
        result = d.determine_dvmc(tree)

        expected = np.std(np.array([2.0, 2.0, 2.0]), ddof=1)
        assert result == pytest.approx(expected)

    def test_determine_dvmc_standard_tree_does_not_import_numpy(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:2):1,C:4):0.5;"), "newick")
        d = DVMC(Namespace(tree="x.tre"))
        original_import = builtins.__import__

        def fail_numpy_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "numpy" or name.startswith("numpy."):
                raise AssertionError("standard DVMC path should not import NumPy")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fail_numpy_import)

        result = d.determine_dvmc(tree)

        expected = 1.0
        assert result == pytest.approx(expected)

    def test_determine_dvmc_standard_tree_handles_mixed_child_counts(self):
        tree = Phylo.read(StringIO("(A:1,(B:2,C:4):1,D:5):0.5;"), "newick")

        result = DVMC._determine_dvmc_standard_tree(tree)

        expected = np.std(np.array([1.0, 3.0, 5.0, 5.0]), ddof=1)
        assert result == pytest.approx(expected)

    def test_determine_dvmc_fallback_does_not_import_numpy(self, monkeypatch):
        class FallbackTree:
            root = object()

            def get_terminals(self):
                return ["a", "b", "c"]

            def depths(self):
                raise AttributeError

            def distance(self, terminal):
                return {"a": 1.0, "b": 2.0, "c": 4.0}[terminal]

        tree = FallbackTree()
        d = DVMC(Namespace(tree="x.tre"))
        monkeypatch.setattr(
            d,
            "calculate_terminal_root_distances_fast",
            lambda _tree: None,
        )
        previous_numpy = sys.modules.pop("numpy", None)
        original_import = builtins.__import__

        def fail_numpy_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "numpy" or name.startswith("numpy."):
                raise AssertionError("fallback DVMC path should not import NumPy")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fail_numpy_import)
        try:
            result = d.determine_dvmc(tree)
        finally:
            if previous_numpy is not None:
                sys.modules["numpy"] = previous_numpy

        expected = np.std(np.array([1.0, 2.0, 4.0]), ddof=1)
        assert result == pytest.approx(expected)

    def test_determine_dvmc_reuses_terminal_list_for_count(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:2):0.5;"), "newick")
        d = DVMC(Namespace(tree="x.tre"))

        def fail_count_terminals(*args, **kwargs):
            raise AssertionError("terminal count should come from get_terminals")

        monkeypatch.setattr(tree, "count_terminals", fail_count_terminals)
        result = d.determine_dvmc(tree)

        expected = np.std(np.array([2.0, 2.0, 2.0]), ddof=1)
        assert result == pytest.approx(expected)

    def test_run_json_output(self, tree_simple, monkeypatch):
        captured = {}
        d = DVMC(Namespace(tree="x.tre", json=True))
        monkeypatch.setattr(
            d,
            "_get_simple_newick_terminal_distance_stats",
            lambda *_args: None,
        )
        monkeypatch.setattr(d, "read_tree_file_unmodified", lambda: tree_simple)
        monkeypatch.setattr(dvmc_module, "print_json", lambda payload: captured.setdefault("payload", payload))

        d.run()

        assert "dvmc" in captured["payload"]
        assert isinstance(captured["payload"]["dvmc"], float)

    def test_run_uses_unmodified_tree_read(self, tree_simple, mocker):
        d = DVMC(Namespace(tree="x.tre", json=False))
        mocker.patch.object(
            d,
            "_get_simple_newick_terminal_distance_stats",
            return_value=None,
        )
        read_tree = mocker.patch.object(
            d,
            "read_tree_file_unmodified",
            return_value=tree_simple,
        )
        mocked_print = mocker.patch("builtins.print")

        d.run()

        read_tree.assert_called_once_with()
        mocked_print.assert_called_once()

    def test_run_uses_simple_newick_terminal_distance_summary(
        self,
        mocker,
        tmp_path,
    ):
        tree_path = tmp_path / "tree.tre"
        tree_path.write_text("((A:1,B:2):1,C:4):0.5;")
        d = DVMC(Namespace(tree=str(tree_path), json=False))
        read_tree = mocker.patch.object(
            d,
            "read_tree_file_unmodified",
            side_effect=AssertionError("simple Newick should use summary path"),
        )
        mocked_print = mocker.patch("builtins.print")

        d.run()

        read_tree.assert_not_called()
        expected = np.std(np.array([2.0, 3.0, 4.0]), ddof=1)
        mocked_print.assert_called_once_with(round(expected, 4))

    def test_run_terminal_output(self, tree_simple, monkeypatch, capsys):
        d = DVMC(Namespace(tree="x.tre", json=False))
        monkeypatch.setattr(
            d,
            "_get_simple_newick_terminal_distance_stats",
            lambda *_args: None,
        )
        monkeypatch.setattr(d, "read_tree_file_unmodified", lambda: tree_simple)

        d.run()

        out, _ = capsys.readouterr()
        assert out.strip() != ""
