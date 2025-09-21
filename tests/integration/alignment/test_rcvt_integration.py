import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestRCVT(object):
    @patch("builtins.print")
    def test_rcvt0(self, mocked_print):
        expected_result_0 = dedent(
            """1\t0.0667"""
        )
        expected_result_1 = dedent(
            """2\t0.04"""
        )
        expected_result_2 = dedent(
            """3\t0.04"""  
        )
        expected_result_3 = dedent(
            """4\t0.08"""
        )
        expected_result_4 = dedent(
            """5\t0.1333"""  
        )
        testargs = [
            "phykit",
            "rcvt",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
            call(expected_result_4),
        ]

    @patch("builtins.print")
    def test_rcvt1(self, mocked_print):
        expected_result_0 = dedent(
            """1\t0.0833"""
        )
        expected_result_1 = dedent(
            """2\t0.0694"""
        )
        expected_result_2 = dedent(
            """3\t0.0417"""  
        )
        expected_result_3 = dedent(
            """4\t0.1111"""
        )
        testargs = [
            "phykit",
            "rcvt",
            f"{here.parent.parent.parent}/sample_files/test_alignment_0.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
        ]

    @patch("builtins.print")
    def test_rcvt2(self, mocked_print):
        expected_result_0 = dedent(
            """1\t0.0417"""
        )
        expected_result_1 = dedent(
            """2\t0.0417"""
        )
        expected_result_2 = dedent(
            """3\t0.0417"""  
        )
        expected_result_3 = dedent(
            """4\t0.125"""
        )
        testargs = [
            "phykit",
            "rcvt",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
        ]

    @patch("builtins.print")
    def test_rcvt_alias(self, mocked_print):
        expected_result_0 = dedent(
            """1\t0.0417"""
        )
        expected_result_1 = dedent(
            """2\t0.0417"""
        )
        expected_result_2 = dedent(
            """3\t0.0417"""  
        )
        expected_result_3 = dedent(
            """4\t0.125"""
        )
        testargs = [
            "phykit",
            "relative_composition_variability_taxon",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(expected_result_0),
            call(expected_result_1),
            call(expected_result_2),
            call(expected_result_3),
        ]

    @patch("builtins.print")
    def test_rcvt_incorrect_input_file(self, mocked_print):
        testargs = [
            "phykit",
            "rcvt",
            f"{here.parent.parent.parent}/sample_files/does_not_exist",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
