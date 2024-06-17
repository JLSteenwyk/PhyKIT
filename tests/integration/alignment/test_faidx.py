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
    def test_faidx(self, mocked_print):
        expected_result = ">1\nA-GTAT"
        testargs = [
            "phykit",
            "faidx",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            '-e',
            '1'
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_faidx_alias0(self, mocked_print):
        expected_result = ">1\nA-GTAT"
        testargs = [
            "phykit",
            "get_entry",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            '-e',
            '1'
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_faidx_alias1(self, mocked_print):
        expected_result = ">1\nA-GTAT"
        testargs = [
            "phykit",
            "ge",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            '-e',
            '1'
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_faidx_alias1(self, mocked_print):
        expected_result = ">4\nML*"
        testargs = [
            "phykit",
            "faidx",
            f"{here.parent.parent.parent}/sample_files/test_alignment.prot.faa",
            '-e',
            '4'
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_faidx_multiple_entries(self, mocked_print):
        expected_result0 = ">1\nA-GTAT"
        expected_result1 = ">2\nA-G-AT"
        testargs = [
            "phykit",
            "faidx",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            '-e',
            '1,2'
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result0),
            call(expected_result1)
        ]
