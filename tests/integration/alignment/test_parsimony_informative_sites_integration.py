import pytest
import sys
import json
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestParsimonyInformativeSites(object):
    @patch("builtins.print")
    def test_parsimony_informative_sites0(self, mocked_print):
        expected_result = "3\t6\t50.0"
        testargs = [
            "phykit",
            "parsimony_informative_sites",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_parsimony_informative_sites2(self, mocked_print):
        expected_result = "1\t9\t11.1111"
        testargs = [
            "phykit",
            "parsimony_informative_sites",
            f"{here.parent.parent.parent}/sample_files/test_alignment_0.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_parsimony_informative_sites1(self, mocked_print):
        expected_result = "0\t3\t0.0"
        testargs = [
            "phykit",
            "parsimony_informative_sites",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_parsimony_informative_sites_alias(self, mocked_print):
        expected_result = "0\t3\t0.0"
        testargs = [
            "phykit",
            "pis",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_parsimony_informative_sites_incorrect_input_file(self, mocked_print):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 1

    @patch("builtins.print")
    def test_parsimony_informative_sites_json(self, mocked_print):
        testargs = [
            "phykit",
            "parsimony_informative_sites",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload == {
            "parsimony_informative_sites": 3,
            "alignment_length": 6,
            "percent_parsimony_informative_sites": 50.0,
        }
