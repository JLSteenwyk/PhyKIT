import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestCovaryingEvolutionaryRates(object):
    @patch("builtins.print")
    def test_covarying_evolutionary_rates(self, mocked_print):
        expected_result = "0.6769\t0.065228"
        testargs = [
            "phykit",
            "covarying_evolutionary_rates",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_1.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple_2.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_alias(self, mocked_print):
        expected_result = "0.6769\t0.065228"
        testargs = [
            "phykit",
            "cover",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_1.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple_2.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_verbose(self, mocked_print):
        testargs = [
            "phykit",
            "cover",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_1.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple_2.tre",
            "-v"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call(f"{-0.1297}\t{-1.2712}"),
            call(f"{0.3588}\t{1.9021}"),
            call(f"{0.3588}\t{0.179}"),
            call(f"{1.2555}\t{0.7115}"),
            call(f"{0.3588}\t{-0.1114}"),
            call(f"{0.3588}\t{0.179}"),
            call(f"{-0.1375}\t{-0.1157}"),
            call(f"{-2.4235}\t{-1.4731}")
        ]

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_incorrect_tree0(self, mocked_print):
        expected_result = "0.6769\t0.065228"
        testargs = [
            "phykit",
            "cover",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tr",
            f"{here.parent.parent.parent}/sample_files/tree_simple_1.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple_2.tre",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_incorrect_tree1(self, mocked_print):
        expected_result = "0.6769\t0.065228"
        testargs = [
            "phykit",
            "cover",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_1.tr",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple_2.tre",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_incorrect_reference(self, mocked_print):
        expected_result = "0.6769\t0.065228"
        testargs = [
            "phykit",
            "cover",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_1.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple_2.tr",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_incorrect_tree_topology(self, mocked_print):
        expected_result = "0.6769\t0.065228"
        testargs = [
            "phykit",
            "cover",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_1.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple_3_incorrect_topology.tre ",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_no_common_tips(self, mocked_print):
        testargs = [
            "phykit",
            "cover",
            f"{here.parent.parent.parent}/sample_files/tree_simple_no_match_0.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_no_match_1.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple_1.tre",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_tree_topologies_do_not_match(self, mocked_print):
        expected_result = "0.6769\t0.065228"
        testargs = [
            "phykit",
            "cover",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_wrong_topology.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple_3_incorrect_topology.tre ",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
