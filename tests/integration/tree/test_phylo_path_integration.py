import pytest
import sys
from mock import patch
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestPhyloPathIntegration(object):
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print):
        testargs = [
            "phykit", "phylo_path",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--traits", f"{here.parent.parent.parent}/sample_files/tree_simple_multi_traits.tsv",
            "--models", f"{here.parent.parent.parent}/sample_files/phylo_path_models.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_alias_ppath(self, mocked_print):
        testargs = [
            "phykit", "ppath",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--traits", f"{here.parent.parent.parent}/sample_files/tree_simple_multi_traits.tsv",
            "--models", f"{here.parent.parent.parent}/sample_files/phylo_path_models.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_alias_phylopath(self, mocked_print):
        testargs = [
            "phykit", "phylopath",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--traits", f"{here.parent.parent.parent}/sample_files/tree_simple_multi_traits.tsv",
            "--models", f"{here.parent.parent.parent}/sample_files/phylo_path_models.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_best_only(self, mocked_print):
        testargs = [
            "phykit", "ppath",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--traits", f"{here.parent.parent.parent}/sample_files/tree_simple_multi_traits.tsv",
            "--models", f"{here.parent.parent.parent}/sample_files/phylo_path_models.txt",
            "--best-only",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1
