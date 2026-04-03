import pytest
import sys
from mock import patch
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)
SAMPLE = here.parent.parent.parent / "sample_files"


@pytest.mark.integration
class TestChronogramIntegration(object):
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print, tmp_path):
        testargs = [
            "phykit", "chronogram",
            "-t", str(SAMPLE / "ultrametric_tree.tre"),
            "--root-age", "70",
            "--plot-output", str(tmp_path / "chrono.png"),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_alias_chrono(self, mocked_print, tmp_path):
        testargs = [
            "phykit", "chrono",
            "-t", str(SAMPLE / "ultrametric_tree.tre"),
            "--root-age", "70",
            "--plot-output", str(tmp_path / "chrono.png"),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_alias_time_tree(self, mocked_print, tmp_path):
        testargs = [
            "phykit", "time_tree",
            "-t", str(SAMPLE / "ultrametric_tree.tre"),
            "--root-age", "70",
            "--plot-output", str(tmp_path / "chrono.png"),
            "--circular",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1
