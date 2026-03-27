import pytest
import sys
from mock import patch
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestPhyloAnovaIntegration(object):
    @patch("builtins.print")
    def test_basic_anova(self, mocked_print):
        testargs = [
            "phykit", "phylo_anova",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--traits", f"{here.parent.parent.parent}/sample_files/tree_simple_anova_single.tsv",
            "--permutations", "99",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_alias_panova(self, mocked_print):
        testargs = [
            "phykit", "panova",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--traits", f"{here.parent.parent.parent}/sample_files/tree_simple_anova_single.tsv",
            "--permutations", "99",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_alias_phylo_manova(self, mocked_print):
        testargs = [
            "phykit", "phylo_manova",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--traits", f"{here.parent.parent.parent}/sample_files/tree_simple_anova_traits.tsv",
            "--permutations", "99",
            "--seed", "42",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_with_pairwise(self, mocked_print):
        testargs = [
            "phykit", "panova",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--traits", f"{here.parent.parent.parent}/sample_files/tree_simple_anova_single.tsv",
            "--permutations", "99",
            "--seed", "42",
            "--pairwise",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1
