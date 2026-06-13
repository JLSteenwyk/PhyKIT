import sys
from pathlib import Path

import pytest
from mock import patch

from phykit.phykit import Phykit

here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
INVARIANT_TRAIT = str(SAMPLE_FILES / "tree_simple_invariant_trait.tsv")


@pytest.mark.integration
class TestPhylogeneticSignalInvariant:
    """An invariant (zero-variance) trait has undefined phylogenetic signal;
    the command should report NA rather than emit nan/inf."""

    def test_blombergs_k_invariant_trait(self, capsys):
        testargs = [
            "phykit",
            "phylogenetic_signal",
            "-t", TREE_SIMPLE,
            "-d", INVARIANT_TRAIT,
            "-m", "blombergs_k",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        captured = capsys.readouterr()
        assert captured.out.strip() == "NA\tNA\tNA"
        assert "zero variance" in captured.err

    def test_lambda_invariant_trait(self, capsys):
        testargs = [
            "phykit",
            "phylogenetic_signal",
            "-t", TREE_SIMPLE,
            "-d", INVARIANT_TRAIT,
            "-m", "lambda",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        captured = capsys.readouterr()
        assert captured.out.strip() == "NA\tNA\tNA\tNA"
        assert "zero variance" in captured.err
