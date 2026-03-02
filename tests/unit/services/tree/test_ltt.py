import pytest
import json
import os
from io import StringIO
from argparse import Namespace
from pathlib import Path

from Bio import Phylo

from phykit.services.tree.ltt import LTT
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
BALANCED_TREE = str(SAMPLE_FILES / "ltt_test_balanced.tre")
LADDER_TREE = str(SAMPLE_FILES / "ltt_test_ladder.tre")
RECENT_BURST_TREE = str(SAMPLE_FILES / "ltt_test_recent_burst.tre")
EARLY_BURST_TREE = str(SAMPLE_FILES / "ltt_test_early_burst.tre")


class TestGammaStat:
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

    def test_too_few_tips_raises(self):
        """Need at least 3 tips."""
        tree = Phylo.read(StringIO("(A:1,B:1);"), "newick")
        with pytest.raises(SystemExit):
            LTT._compute_gamma(tree)


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


class TestLTTPlot:
    def test_plot_creates_file(self, tmp_path):
        """Plot should create a PNG file."""
        tree = Phylo.read(
            StringIO("(((A:1,B:1):1,(C:1,D:1):1):1,E:3);"), "newick"
        )
        ltt = LTT._compute_ltt(tree)
        plot_path = str(tmp_path / "ltt.png")
        LTT._plot_ltt(ltt, plot_path, gamma=-0.5, p_value=0.6)
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
