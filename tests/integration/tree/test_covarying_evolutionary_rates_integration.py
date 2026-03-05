import pytest
import sys
import json
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestCovaryingEvolutionaryRates(object):
    @patch("builtins.print")
    def test_covarying_evolutionary_rates(self, mocked_print):
        expected_result = "0.5436\t0.054828"
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
        expected_result = "0.5436\t0.054828"
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
            call(f"{-0.333}\t{-1.3686}\traccoon"),
            call(f"{0.2747}\t{1.7774}\tbear"),
            call(f"{0.2747}\t{0.0692}\tsea_lion"),
            call(f"{1.3905}\t{0.597}\tseal"),
            call(f"{0.2747}\t{-0.2188}\tmonkey"),
            call(f"{0.2747}\t{0.0692}\tcat"),
            call(f"{-0.3428}\t{-0.223}\tweasel"),
            call(f"{-3.1873}\t{-1.5688}\tdog"),
            call(f"{0.2747}\t{0.0692}\traccoon;bear"),
            call(f"{0.2747}\t{1.5686}\tsea_lion;seal;monkey;cat;weasel"),
            call(f"{0.2747}\t{-1.4737}\tsea_lion;seal"),
            call(f"{0.2747}\t{0.0692}\tmonkey;cat;weasel"),
            call(f"{0.2747}\t{0.6333}\tmonkey;cat")
        ]

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_incorrect_tree0(self, mocked_print):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_incorrect_tree1(self, mocked_print):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_incorrect_reference(self, mocked_print):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_incorrect_tree_topology(self, mocked_print):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_no_common_tips(self, mocked_print):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_tree_topologies_do_not_match(self, mocked_print):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_outlier(self, mocked_print):
        expected_result = "0.5874\t0.044626"
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

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_json(self, mocked_print):
        testargs = [
            "phykit",
            "covarying_evolutionary_rates",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_1.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple_2.tre",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload == {
            "correlation": 0.5436,
            "p_value": 0.054828,
            "verbose": False,
        }

    @patch("builtins.print")
    def test_covarying_evolutionary_rates_json_verbose(self, mocked_print):
        testargs = [
            "phykit",
            "cover",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_1.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple_2.tre",
            "-v",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["verbose"] is True
        assert payload["rows"][0] == payload["branches"][0]
        assert payload["branches"][0] == {
            "branch": "raccoon",
            "tree_one_rate": -1.3686,
            "tree_zero_rate": -0.333,
        }

    @patch("phykit.services.tree.covarying_evolutionary_rates.CovaryingEvolutionaryRates._plot_covarying_rates_scatter")
    @patch("builtins.print")
    def test_covarying_evolutionary_rates_plot(self, mocked_print, mocked_plot):
        testargs = [
            "phykit",
            "covarying_evolutionary_rates",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_1.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple_2.tre",
            "--plot",
            "--plot-output",
            "cover_test_plot.png",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_plot.called
        assert any(
            call.args[0] == "Saved covarying rates plot: cover_test_plot.png"
            for call in mocked_print.mock_calls
        )

    @patch("phykit.services.tree.covarying_evolutionary_rates.CovaryingEvolutionaryRates._plot_covarying_rates_scatter")
    @patch("builtins.print")
    def test_covarying_evolutionary_rates_plot_json(self, mocked_print, mocked_plot):
        testargs = [
            "phykit",
            "covarying_evolutionary_rates",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_1.tre",
            "-r",
            f"{here.parent.parent.parent}/sample_files/tree_simple_2.tre",
            "--plot",
            "--plot-output",
            "cover_test_plot_json.png",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_plot.called
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["plot_output"] == "cover_test_plot_json.png"
