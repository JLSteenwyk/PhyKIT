import pytest
import sys
import json
from mock import patch, call
from pathlib import Path

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
    def test_gc_content3(self, mocked_print):
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
    def test_gc_content_alias(self, mocked_print):
        expected_result = 0.3
        testargs = [
            "phykit",
            "gc",
            f"{here.parent.parent.parent}/sample_files/test_alignment_2.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_gc_content_incorrect_input_file(self, mocked_print):

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 1

    @patch("builtins.print")
    def test_gc_content_verbose(self, mocked_print):
        expected_result0 = "1\t0.5"
        expected_result1 = "2\t0.3333"
        expected_result2 = "3\t0.3333"
        expected_result3 = "4\t0.0"
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/test_alignment_2.fa",
            "-v",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result0),
            call(expected_result1),
            call(expected_result2),
            call(expected_result3),
        ]

    @patch("builtins.print")
    def test_gc_content_zero_division(self, mocked_print):
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/zero_length_alignment.fa",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit):
                Phykit()
        
        assert mocked_print.mock_calls == [
            call("Input file has an unacceptable format. Please check input file argument."),
        ]

    @patch("builtins.print")
    def test_gc_content_json(self, mocked_print):
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload == {"gc_content": 0.2273, "verbose": False}

    @patch("builtins.print")
    def test_gc_content_json_verbose(self, mocked_print):
        testargs = [
            "phykit",
            "gc_content",
            f"{here.parent.parent.parent}/sample_files/test_alignment_2.fa",
            "-v",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["verbose"] is True
        assert payload["rows"][0] == payload["sequences"][0]
        assert payload["sequences"][0] == {"taxon": "1", "gc_content": 0.5}
