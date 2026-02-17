import sys
from pathlib import Path

import pytest
from mock import patch

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestPlotAlignmentQC:
    @patch("phykit.services.alignment.plot_alignment_qc.PlotAlignmentQC.run")
    def test_plot_alignment_qc(self, mocked_run):
        testargs = [
            "phykit",
            "plot_alignment_qc",
            f"{here.parent.parent.parent}/sample_files/alignment_outlier_taxa.fa",
            "-o",
            "dummy.png",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_run.called

    @patch("phykit.services.alignment.plot_alignment_qc.PlotAlignmentQC.run")
    def test_plot_alignment_qc_alias(self, mocked_run):
        testargs = [
            "phykit",
            "paqc",
            f"{here.parent.parent.parent}/sample_files/alignment_outlier_taxa.fa",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_run.called
