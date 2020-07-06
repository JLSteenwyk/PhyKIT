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