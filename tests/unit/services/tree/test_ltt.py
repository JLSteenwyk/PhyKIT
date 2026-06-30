import pytest
import json
import math
import os
import subprocess
import sys
from io import StringIO
from argparse import Namespace
from pathlib import Path
from unittest.mock import Mock

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

from phykit.services.tree.ltt import LTT, _ltt_from_internal_depths, _max_tip_height
from phykit.errors import PhykitUserError


def test_max_tip_height_scans_tips_once():
    class SinglePassTips:
        def __init__(self, tips):
            self.tips = tips
            self.iterations = 0

        def __iter__(self):
            self.iterations += 1
            if self.iterations > 1:
                raise AssertionError("tips should be scanned once")
            return iter(self.tips)

    tips = [object(), object(), object()]
    depths = {tips[0]: 1.25, tips[1]: 4.5, tips[2]: 2.0}
    iterable = SinglePassTips(tips)

    assert _max_tip_height(iterable, depths, 0.5) == 4.0
    assert iterable.iterations == 1


def test_ltt_from_internal_depths_does_not_slice_rows():
    class NoSliceList(list):
        def __getitem__(self, key):
            if isinstance(key, slice):
                raise AssertionError("LTT row construction should not slice")
            return super().__getitem__(key)

    internal_depths = NoSliceList([0.0, 1.25, 2.5])

    assert _ltt_from_internal_depths(internal_depths, 4.0) == [
        (0.0, 2),
        (1.25, 3),
        (2.5, 4),
        (4.0, 4),
    ]


def test_module_import_does_not_import_json_plot_numpy_or_matplotlib():
    code = """
import sys
import phykit.services.tree.ltt as module
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "numpy" not in sys.modules
assert "matplotlib" not in sys.modules
assert "matplotlib.pyplot" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
BALANCED_TREE = str(SAMPLE_FILES / "ltt_test_balanced.tre")
LADDER_TREE = str(SAMPLE_FILES / "ltt_test_ladder.tre")
RECENT_BURST_TREE = str(SAMPLE_FILES / "ltt_test_recent_burst.tre")
EARLY_BURST_TREE = str(SAMPLE_FILES / "ltt_test_early_burst.tre")


class TestGammaStat:
    def test_terminal_clades_preserves_order_with_mixed_child_counts(self):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        tips = LTT._terminal_clades(tree)

        assert [tip.name for tip in tips] == ["A", "B", "C", "D", "E", "F"]

    def test_balanced_tree_matches_r(self):
        """R: gammaStat(balanced_8) = -1.414214"""
        tree = Phylo.read(
            StringIO("(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);"),
            "newick",
        )
        gamma, p_val, _, _ = LTT._compute_gamma(tree)
        assert abs(gamma - (-1.414214)) < 1e-4

    def test_ladder_tree_matches_r(self):
        """R: gammaStat(ladder_5) = -0.7142857"""
        tree = Phylo.read(
            StringIO("((((A:1,B:1):1,C:2):1,D:3):1,E:4);"), "newick"
        )
        gamma, p_val, _, _ = LTT._compute_gamma(tree)
        assert abs(gamma - (-0.7142857)) < 1e-4

    def test_recent_burst_matches_r(self):
        """R: gammaStat(recent_burst) = 2.282479"""
        tree = Phylo.read(
            StringIO(
                "(((((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1):0.1,"
                "((E:0.1,F:0.1):0.1,(G:0.1,H:0.1):0.1):0.1):3,"
                "I:3.4):0.6,J:4);"
            ),
            "newick",
        )
        gamma, p_val, _, _ = LTT._compute_gamma(tree)
        assert abs(gamma - 2.282479) < 1e-4

    def test_early_burst_matches_r(self):
        """R: gammaStat(early_burst) = -3.536202"""
        tree = Phylo.read(
            StringIO(
                "((((((A:3,B:3):0.1,C:3.1):0.1,D:3.2):0.1,E:3.3):0.1,F:3.4):0.1,G:3.5);"
            ),
            "newick",
        )
        gamma, p_val, _, _ = LTT._compute_gamma(tree)
        assert abs(gamma - (-3.536202)) < 1e-4

    def test_positive_gamma_significant(self):
        """Recent burst should have significant positive gamma."""
        tree = Phylo.read(
            StringIO(
                "(((((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1):0.1,"
                "((E:0.1,F:0.1):0.1,(G:0.1,H:0.1):0.1):0.1):3,"
                "I:3.4):0.6,J:4);"
            ),
            "newick",
        )
        gamma, p_val, _, _ = LTT._compute_gamma(tree)
        assert gamma > 0
        assert p_val < 0.05

    def test_negative_gamma_early_burst(self):
        """Early burst should have negative gamma."""
        tree = Phylo.read(
            StringIO(
                "((((((A:3,B:3):0.1,C:3.1):0.1,D:3.2):0.1,E:3.3):0.1,F:3.4):0.1,G:3.5);"
            ),
            "newick",
        )
        gamma, p_val, _, _ = LTT._compute_gamma(tree)
        assert gamma < 0

    def test_gamma_p_value_uses_standard_library_erfc(self, monkeypatch):
        tree = Phylo.read(
            StringIO("(((A:1,B:1):1,(C:1,D:1):1):1,E:3);"), "newick"
        )
        real_import = __import__

        def fail_scipy_import(name, *args, **kwargs):
            if name == "scipy" or name.startswith("scipy."):
                raise AssertionError("gamma p-value should not import scipy")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr("builtins.__import__", fail_scipy_import)

        gamma, p_val, _, _ = LTT._compute_gamma(tree)

        assert p_val == pytest.approx(math.erfc(abs(gamma) / math.sqrt(2.0)))

    def test_too_few_tips_raises(self):
        """Need at least 3 tips."""
        tree = Phylo.read(StringIO("(A:1,B:1);"), "newick")
        with pytest.raises(SystemExit):
            LTT._compute_gamma(tree)

    def test_gamma_fast_path_does_not_call_distance(self, monkeypatch):
        newick = "(((A:1,B:1):1,(C:1,D:1):1):1,E:3):0.5;"
        expected_gamma, _, expected_bt, expected_g = LTT._compute_gamma(
            Phylo.read(StringIO(newick), "newick")
        )
        tree = Phylo.read(StringIO(newick), "newick")

        def fail_distance(*args, **kwargs):
            raise AssertionError("distance fallback should not be called")

        monkeypatch.setattr(tree, "distance", fail_distance)
        gamma, _, bt, g = LTT._compute_gamma(tree)

        assert gamma == pytest.approx(expected_gamma)
        assert bt == expected_bt
        assert g == expected_g

    def test_gamma_default_tips_use_direct_terminal_traversal(self, monkeypatch):
        newick = "(((A:1,B:1):1,(C:1,D:1):1):1,E:3):0.5;"
        expected_gamma, _, expected_bt, expected_g = LTT._compute_gamma(
            Phylo.read(StringIO(newick), "newick")
        )
        tree = Phylo.read(StringIO(newick), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("generic terminal traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        gamma, _, bt, g = LTT._compute_gamma(tree)

        assert gamma == pytest.approx(expected_gamma)
        assert bt == expected_bt
        assert g == expected_g

    def test_depths_from_root_uses_direct_standard_tree_traversal(self, monkeypatch):
        tree = Phylo.read(
            StringIO("(((A:1,B:1):1,(C:1,D:1):1):1,E:3):0.5;"),
            "newick",
        )

        def fail_depths(*_args, **_kwargs):
            raise AssertionError("generic tree.depths() should not be used")

        monkeypatch.setattr(tree, "depths", fail_depths)

        root, depths, root_depth = LTT._depths_from_root(tree)

        assert root is tree.root
        assert root_depth == 0.0
        assert max(depths[tip] for tip in tree.get_terminals()) == pytest.approx(3.0)

    def test_internal_depths_from_root_handles_mixed_child_counts(self, monkeypatch):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):2,(D:1,E:1,F:1):3);"),
            "newick",
        )
        depth_data = LTT._depths_from_root(tree)

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("generic clade traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        with_root = LTT._internal_depths_from_root(
            tree,
            depth_data,
            include_root=True,
        )
        without_root = LTT._internal_depths_from_root(
            tree,
            depth_data,
            include_root=False,
        )

        assert with_root == pytest.approx([0.0, 2.0, 3.0])
        assert without_root == pytest.approx([2.0, 3.0])

    def test_gamma_uses_provided_tree_context(self, monkeypatch):
        newick = "(((A:1,B:1):1,(C:1,D:1):1):1,E:3):0.5;"
        expected_gamma, _, expected_bt, expected_g = LTT._compute_gamma(
            Phylo.read(StringIO(newick), "newick")
        )
        tree = Phylo.read(StringIO(newick), "newick")
        tips = list(tree.get_terminals())
        depth_data = LTT._depths_from_root(tree)

        def fail_traversal(*args, **kwargs):
            raise AssertionError("provided context should avoid setup traversal")

        monkeypatch.setattr(tree, "get_terminals", fail_traversal)
        monkeypatch.setattr(tree, "depths", fail_traversal)
        monkeypatch.setattr(tree, "distance", fail_traversal)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)

        gamma, _, bt, g = LTT._compute_gamma(
            tree, tips=tips, depth_data=depth_data
        )

        assert gamma == pytest.approx(expected_gamma)
        assert bt == expected_bt
        assert g == expected_g

    def test_combined_gamma_and_ltt_matches_standalone_helpers(self):
        newick = "(((A:1,B:1):1,(C:1,D:1):1):1,E:3):0.5;"
        tree = Phylo.read(StringIO(newick), "newick")
        tips = LTT._terminal_clades(tree)
        depth_data = LTT._depths_from_root(tree)

        expected_gamma, expected_p, expected_bt, expected_g = LTT._compute_gamma(
            tree, tips=tips, depth_data=depth_data
        )
        expected_ltt = LTT._compute_ltt(tree, tips=tips, depth_data=depth_data)

        gamma, p_value, bt, g, ltt_data = LTT._compute_gamma_and_ltt(
            tree, tips=tips, depth_data=depth_data
        )

        assert gamma == pytest.approx(expected_gamma)
        assert p_value == pytest.approx(expected_p)
        assert bt == expected_bt
        assert g == expected_g
        assert ltt_data == expected_ltt


class TestLTTData:
    def test_ltt_starts_at_two(self):
        """LTT should start with 2 lineages at time 0."""
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:2);"), "newick")
        ltt = LTT._compute_ltt(tree)
        assert ltt[0] == (0.0, 2)

    def test_ltt_ends_at_n_tips(self):
        """LTT should end at N lineages."""
        tree = Phylo.read(
            StringIO("(((A:1,B:1):1,(C:1,D:1):1):1,E:3);"), "newick"
        )
        ltt = LTT._compute_ltt(tree)
        assert ltt[-1][1] == 5

    def test_ltt_monotonically_increasing(self):
        """Lineage count should only increase."""
        tree = Phylo.read(
            StringIO(
                "(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);"
            ),
            "newick",
        )
        ltt = LTT._compute_ltt(tree)
        lineages = [pt[1] for pt in ltt]
        for i in range(1, len(lineages)):
            assert lineages[i] >= lineages[i - 1]

    def test_ltt_correct_steps(self):
        """3-taxon tree: 2 lineages at root, 3 at one split."""
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:2);"), "newick")
        ltt = LTT._compute_ltt(tree)
        # At root: 2 lineages
        # At the AB split (1.0 from root): 3 lineages
        assert ltt[0] == (0.0, 2)
        assert any(abs(t - 1.0) < 1e-6 and n == 3 for t, n in ltt)

    def test_ltt_fast_path_does_not_call_distance(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:2):0.5;"), "newick")

        def fail_distance(*args, **kwargs):
            raise AssertionError("distance fallback should not be called")

        monkeypatch.setattr(tree, "distance", fail_distance)
        ltt = LTT._compute_ltt(tree)

        assert ltt == [(0.0, 2), (1.0, 3), (2.0, 3)]

    def test_ltt_default_tips_use_direct_terminal_traversal(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:2):0.5;"), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("generic terminal traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        assert LTT._compute_ltt(tree) == [(0.0, 2), (1.0, 3), (2.0, 3)]

    def test_ltt_uses_provided_tree_context(self, monkeypatch):
        newick = "(((A:1,B:1):1,(C:1,D:1):1):1,E:3):0.5;"
        expected = LTT._compute_ltt(Phylo.read(StringIO(newick), "newick"))
        tree = Phylo.read(StringIO(newick), "newick")
        tips = list(tree.get_terminals())
        depth_data = LTT._depths_from_root(tree)

        def fail_traversal(*args, **kwargs):
            raise AssertionError("provided context should avoid setup traversal")

        monkeypatch.setattr(tree, "get_terminals", fail_traversal)
        monkeypatch.setattr(tree, "depths", fail_traversal)
        monkeypatch.setattr(tree, "distance", fail_traversal)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)

        assert LTT._compute_ltt(tree, tips=tips, depth_data=depth_data) == expected


class TestLTTPlot:
    def test_plot_creates_file(self, tmp_path):
        """Plot should create a PNG file."""
        tree = Phylo.read(
            StringIO("(((A:1,B:1):1,(C:1,D:1):1):1,E:3);"), "newick"
        )
        ltt_data = LTT._compute_ltt(tree)
        plot_path = str(tmp_path / "ltt.png")
        instance = LTT(Namespace(tree="dummy.tre"))
        instance._plot_ltt(ltt_data, plot_path, gamma=-0.5, p_value=0.6)
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(tree="t.tre")
        svc = LTT(args)
        assert svc.tree_file_path == "t.tre"
        assert svc.verbose is False
        assert svc.json_output is False
        assert svc.plot_output is None

    def test_all_options(self):
        args = Namespace(
            tree="t.tre", verbose=True, json=True, plot_output="out.png"
        )
        svc = LTT(args)
        assert svc.verbose is True
        assert svc.json_output is True
        assert svc.plot_output == "out.png"


class TestRun:
    def test_run_uses_unmodified_tree_read(self, monkeypatch):
        tree = object()
        tips = ["a", "b", "c"]
        depth_data = object()
        args = Namespace(tree="dummy.tre", verbose=False, json=False, plot_output=None)
        svc = LTT(args)

        read_tree = Mock(return_value=tree)
        monkeypatch.setattr(svc, "read_tree_file_unmodified", read_tree)
        monkeypatch.setattr(
            svc,
            "read_tree_file",
            Mock(side_effect=AssertionError("copying tree reader should not be used")),
        )
        monkeypatch.setattr(svc, "validate_tree", Mock())
        monkeypatch.setattr(svc, "_terminal_clades", Mock(return_value=tips))
        monkeypatch.setattr(svc, "_depths_from_root", Mock(return_value=depth_data))
        monkeypatch.setattr(
            svc,
            "_compute_gamma_and_ltt",
            Mock(return_value=(0.0, 1.0, [], [], [(0.0, 2), (1.0, 3)])),
        )
        monkeypatch.setattr(svc, "_output_text", Mock())

        svc.run()

        read_tree.assert_called_once_with()

    def test_text_output(self, capsys):
        args = Namespace(
            tree=BALANCED_TREE, verbose=False, json=False, plot_output=None
        )
        svc = LTT(args)
        svc.run()
        out, _ = capsys.readouterr()
        parts = out.strip().split("\t")
        assert len(parts) == 2
        gamma_val = float(parts[0])
        p_val = float(parts[1])
        assert abs(gamma_val - (-1.4142)) < 1e-2

    def test_verbose_output(self, capsys):
        args = Namespace(
            tree=BALANCED_TREE, verbose=True, json=False, plot_output=None
        )
        svc = LTT(args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Branching times" in out
        assert "Lineage-through-time" in out

    def test_output_text_batches_verbose_report(self, monkeypatch):
        args = Namespace(tree=BALANCED_TREE, verbose=True, json=False, plot_output=None)
        svc = LTT(args)
        printed = Mock()
        monkeypatch.setattr("builtins.print", printed)

        svc._output_text(
            -1.23456,
            0.04567,
            [(0.0, 2), (0.5, 3)],
            [0.25, 0.75],
            [],
        )

        printed.assert_called_once_with(
            "-1.2346\t0.0457\n"
            "\n"
            "Branching times (node ages):\n"
            "  1\t0.250000\n"
            "  2\t0.750000\n"
            "\n"
            "Lineage-through-time:\n"
            "  time_from_root\tn_lineages\n"
            "  0.000000\t2\n"
            "  0.500000\t3"
        )

    def test_json_output(self, capsys):
        args = Namespace(
            tree=BALANCED_TREE, verbose=False, json=True, plot_output=None
        )
        svc = LTT(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "gamma" in data
        assert "p_value" in data
        assert "branching_times" in data
        assert "ltt" in data
        assert abs(data["gamma"] - (-1.414214)) < 1e-4

    def test_plot_output(self, tmp_path):
        plot_path = str(tmp_path / "test_ltt.png")
        args = Namespace(
            tree=BALANCED_TREE, verbose=False, json=False, plot_output=plot_path
        )
        svc = LTT(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0

    def test_too_few_tips_raises(self, tmp_path):
        tree_file = tmp_path / "two_tips.tre"
        tree_file.write_text("(A:1,B:1);\n")
        args = Namespace(
            tree=str(tree_file), verbose=False, json=False, plot_output=None
        )
        svc = LTT(args)
        with pytest.raises(SystemExit):
            svc.run()
