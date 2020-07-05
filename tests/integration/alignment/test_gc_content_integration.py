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
    def test_gc_content0(self, mocked_print):
        expected_result = 0.2273
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gc_content1(self, mocked_print):
        expected_result = 0.2941
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/test_alignment_0.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gc_content2(self, mocked_print):
        expected_result = 0.25
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gc_content2(self, mocked_print):
        expected_result = 0.3
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/test_alignment_2.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gc_content_incorrect_input_file(self, mocked_print):
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/test_trees.txt",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2