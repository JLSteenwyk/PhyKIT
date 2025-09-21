import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestLBScore(object):
    @patch("builtins.print")
    def test_lb_score0(self, mocked_print):
        testargs = [
            "phykit",
            "lb_score",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        expected = [
            dict(label="mean", value=-12.5),
            dict(label="median", value=-27.806),
            dict(label="25th percentile", value=-31.0492),
            dict(label="75th percentile", value=-12.9039),
            dict(label="minimum", value=-39.2834),
            dict(label="maximum", value=65.6709),
            dict(label="standard deviation", value=35.2669),
            dict(label="variance", value=1243.7527),
        ]

        for print_call, expected_call in zip(mocked_print.call_args_list, expected):
            print_call_args, _ = print_call
            [print_label, print_value] = print_call_args[0].split(': ')
            assert print_label == expected_call["label"]
            assert isclose(float(print_value), expected_call["value"], rel_tol = 0.0001)

    @patch("builtins.print")
    def test_lb_score1(self, mocked_print):
        testargs = [
            "phykit",
            "lb_score",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("mean: -10.0"),
            call("median: -12.7812"),
            call("25th percentile: -21.7398"),
            call("75th percentile: 1.082"),
            call("minimum: -22.4599"),
            call("maximum: 6.964"),
            call("standard deviation: 12.7705"),
            call("variance: 163.086")
        ]

    @patch("builtins.print")
    def test_lb_score_alias0(self, mocked_print):
        testargs = [
            "phykit",
            "long_branch_score",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("mean: -10.0"),
            call("median: -12.7812"),
            call("25th percentile: -21.7398"),
            call("75th percentile: 1.082"),
            call("minimum: -22.4599"),
            call("maximum: 6.964"),
            call("standard deviation: 12.7705"),
            call("variance: 163.086")
        ]

    @patch("builtins.print")
    def test_lb_score_alias1(self, mocked_print):
        testargs = [
            "phykit",
            "lbs",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("mean: -10.0"),
            call("median: -12.7812"),
            call("25th percentile: -21.7398"),
            call("75th percentile: 1.082"),
            call("minimum: -22.4599"),
            call("maximum: 6.964"),
            call("standard deviation: 12.7705"),
            call("variance: 163.086")
        ]

    @patch("builtins.print")
    def test_lb_score_verbose(self, mocked_print):
        expected_result0 = "Aspergillus_fischeri_IBT_3003\t-21.8103"
        expected_result1 = "Aspergillus_fischeri_IBT_3007\t-21.7868"
        expected_result2 = "Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1\t-22.4599"
        expected_result3 = "Aspergillus_fischeri_NRRL4585\t-21.4924"
        expected_result4 = "Aspergillus_fumigatus_Af293\t4.1569"
        expected_result5 = "Aspergillus_fumigatus_CEA10\t0.9815"
        expected_result6 = "Aspergillus_fumigatus_HMR_AF_270\t6.964"
        expected_result7 = "Aspergillus_fumigatus_Z5\t1.1156"
        expected_result8 = "Aspergillus_oerlinghausenensis_CBS139183\t-4.0699"
        expected_result9 = "Aspergillus_fischeri_NRRL4161\t-21.5985"

        testargs = [
            "phykit",
            "lb_score",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-v"
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call(expected_result0),
            call(expected_result1),
            call(expected_result2),
            call(expected_result3),
            call(expected_result4),
            call(expected_result5),
            call(expected_result6),
            call(expected_result7),
            call(expected_result8),
            call(expected_result9),
        ]

    @patch("builtins.print")
    def test_lb_score_incorrect_input(self, mocked_print):
        testargs = [
            "phykit",
            "lbs",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tr",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_lb_score_zero_division_error(self, mocked_print):
        testargs = [
            "phykit",
            "lbs",
            f"{here.parent.parent.parent}/sample_files/tree_simple_blm_5.tre.factor_0.0.tre",
        ]

        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()
        
        assert mocked_print.mock_calls == [
            call("Invalid tree. Tree should contain branch lengths"),
        ]
