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
        expected_result = "0.5755\t0.03959"
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
        expected_result = "0.5755\t0.03959"
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
            call(f"{-0.333}\t{-0.2246}"),
            call(f"{0.2747}\t{0.1191}"),
            call(f"{0.2747}\t{0.1191}"),
            call(f"{1.3905}\t{0.7501}"),
            call(f"{0.2747}\t{0.1191}"),
            call(f"{0.2747}\t{0.1191}"),
            call(f"{-0.3428}\t{-0.2302}"),
            call(f"{-3.1873}\t{-1.8388}"),
            call(f"{0.2747}\t{0.1191}"),
            call(f"{0.2747}\t{-1.3055}"),
            call(f"{0.2747}\t{2.6302}"),
            call(f"{0.2747}\t{0.1191}"),
            call(f"{0.2747}\t{-0.4956}")
        ]

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_incorrect_tree0(self, mocked_print):
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

    @patch("builtins.print")
    def test_covarying_evolutionary_rates(self, mocked_print):
        expected_result = "0.5755\t0.03959"
        testargs = [
            "phykit",
            "covarying_evolutionary_rates",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_outlier_branch.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple_2.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]
