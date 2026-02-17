import pytest
import sys
import json
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestRCV(object):
    @patch("builtins.print")
    def test_rcv0(self, mocked_print):
        expected_result = 0.292
        testargs = [
            "phykit",
            "rcv",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_rcv1(self, mocked_print):
        expected_result = 0.2708
        testargs = [
            "phykit",
            "rcv",
            f"{here.parent.parent.parent}/sample_files/test_alignment_0.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_rcv2(self, mocked_print):
        expected_result = 0.25
        testargs = [
            "phykit",
            "rcv",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_rcv_alias(self, mocked_print):
        expected_result = 0.25
        testargs = [
            "phykit",
            "relative_composition_variability",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_rcv_incorrect_input_file(self, mocked_print):
        testargs = [
            "phykit",
            "rcv",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.f",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_rcv_json(self, mocked_print):
        testargs = [
            "phykit",
            "rcv",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload == {"rcv": 0.292}
