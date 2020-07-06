import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestGCContent(object):
    @patch("builtins.print")
    def test_gc_content_clustal(self, mocked_print):
        expected_result = 0.2273
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/simple_aln.clustal",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]
    
    @patch("builtins.print")
    def test_gc_content_maf(self, mocked_print):
        expected_result = 0.2273
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/simple_aln.maf",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gc_content_mauve(self, mocked_print):
        expected_result = 0.2273
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/simple_aln.mauve",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gc_content_phylip(self, mocked_print):
        expected_result = 0.2273
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/simple_aln.phylip",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gc_content_stockholm(self, mocked_print):
        expected_result = 0.2273
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/simple_aln.stockholm",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]