import pytest
import sys
from mock import patch
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)
SAMPLE = here.parent.parent.parent / "sample_files"


@pytest.mark.integration
class TestTransferAnnotationsIntegration(object):
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print, tmp_path):
        testargs = [
            "phykit", "transfer_annotations",
            "--source", str(SAMPLE / "wastral_annotated.tre"),
            "--target", str(SAMPLE / "raxml_unannotated.tre"),
            "-o", str(tmp_path / "out.tre"),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_alias_transfer_annot(self, mocked_print, tmp_path):
        testargs = [
            "phykit", "transfer_annot",
            "--source", str(SAMPLE / "wastral_annotated.tre"),
            "--target", str(SAMPLE / "raxml_unannotated.tre"),
            "-o", str(tmp_path / "out.tre"),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1

    @patch("builtins.print")
    def test_alias_annotate_tree(self, mocked_print, tmp_path):
        testargs = [
            "phykit", "annotate_tree",
            "--source", str(SAMPLE / "wastral_annotated.tre"),
            "--target", str(SAMPLE / "raxml_unannotated.tre"),
            "-o", str(tmp_path / "out.tre"),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_count >= 1
