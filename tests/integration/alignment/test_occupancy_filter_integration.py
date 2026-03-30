import pytest
import sys
from mock import patch
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)
SAMPLE = f"{here.parent.parent.parent}/sample_files/occupancy_test"


@pytest.mark.integration
class TestOccupancyFilterIntegration(object):
    @patch("builtins.print")
    def test_basic_fasta(self, mocked_print, tmp_path):
        testargs = [
            "phykit", "occupancy_filter",
            "-l", f"{SAMPLE}/aln_list.txt",
            "-f", "fasta", "-t", "2",
            "-o", str(tmp_path),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_alias_occ_filter(self, mocked_print, tmp_path):
        testargs = [
            "phykit", "occ_filter",
            "-l", f"{SAMPLE}/aln_list.txt",
            "-f", "fasta", "-t", "2",
            "-o", str(tmp_path),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_trees_format(self, mocked_print, tmp_path):
        testargs = [
            "phykit", "filter_occupancy",
            "-l", f"{SAMPLE}/tree_list.txt",
            "-f", "trees", "-t", "2",
            "-o", str(tmp_path),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1
