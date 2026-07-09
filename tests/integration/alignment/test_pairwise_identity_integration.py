import pytest
import sys
import json
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


def summary_call(
    mean,
    median,
    twenty_fifth,
    seventy_fifth,
    minimum,
    maximum,
    standard_deviation,
    variance,
):
    return call(
        (
            f"mean: {mean}\n"
            f"median: {median}\n"
            f"25th percentile: {twenty_fifth}\n"
            f"75th percentile: {seventy_fifth}\n"
            f"minimum: {minimum}\n"
            f"maximum: {maximum}\n"
            f"standard deviation: {standard_deviation}\n"
            f"variance: {variance}"
        )
    )


@pytest.mark.integration
class TestPairwiseIdentity(object):
    @patch("builtins.print")
    def test_pairwise_identity0(self, mocked_print):
        testargs = [
            "phykit",
            "pairwise_identity",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            summary_call(
                "0.4833", "0.5", "0.3333", "0.6667",
                "0.1667", "0.8333", "0.2284", "0.0522",
            )
        ]

    @patch("builtins.print")
    def test_pairwise_identity1(self, mocked_print):
        testargs = [
            "phykit",
            "pairwise_identity",
            f"{here.parent.parent.parent}/sample_files/test_alignment_0.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            summary_call(
                "0.7593", "0.7778", "0.6944", "0.8611",
                "0.5556", "0.8889", "0.1299", "0.0169",
            )
        ]

    @patch("builtins.print")
    def test_pairwise_identity2(self, mocked_print):
        testargs = [
            "phykit",
            "pairwise_identity",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            summary_call(
                "0.8333", "0.8333", "0.6667", "1.0",
                "0.6667", "1.0", "0.1826", "0.0333",
            )
        ]

    @patch("builtins.print")
    def test_pairwise_identity_alias(self, mocked_print):
        testargs = [
            "phykit",
            "pi",
            f"{here.parent.parent.parent}/sample_files/test_alignment_1.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            summary_call(
                "0.8333", "0.8333", "0.6667", "1.0",
                "0.6667", "1.0", "0.1826", "0.0333",
            )
        ]

    @patch("builtins.print")
    def test_pairwise_identity_incorrect_input_file(self, mocked_print):

        testargs = ["phykit", "pairwise_identity", "missing.fa"]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            with patch.object(sys, "argv", testargs):
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2
        assert mocked_print.mock_calls == [
            call("missing.fa corresponds to no such file."),
            call("Please check file name and pathing"),
        ]

    @patch("builtins.print")
    def test_pairwise_identity_verbose(self, mocked_print):
        expected_result = (
            "1\t2\t0.8333\n"
            "1\t3\t0.5\n"
            "1\t4\t0.1667\n"
            "1\t5\t0.1667\n"
            "2\t3\t0.6667\n"
            "2\t4\t0.3333\n"
            "2\t5\t0.3333\n"
            "3\t4\t0.6667\n"
            "3\t5\t0.5\n"
            "4\t5\t0.6667"
        )
        testargs = [
            "phykit",
            "pairwise_identity",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-v",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_pairwise_identity_10009at7524_aa_exclude_gaps(self, mocked_print):
        testargs = [
            "phykit",
            "pairwise_identity",
            f"{here.parent.parent.parent}/sample_files/10009at7524_aa.aln",
            "-e",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            summary_call(
                "0.8136", "0.8423", "0.8096", "0.8692",
                "0.6192", "0.9269", "0.0831", "0.0069",
            )
        ]

    @patch("builtins.print")
    def test_pairwise_identity_10009at7524_aa_exclude_gaps_long_arg(self, mocked_print):
        testargs = [
            "phykit",
            "pairwise_identity",
            f"{here.parent.parent.parent}/sample_files/10009at7524_aa.aln",
            "--exclude_gaps",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            summary_call(
                "0.8136", "0.8423", "0.8096", "0.8692",
                "0.6192", "0.9269", "0.0831", "0.0069",
            )
        ]

    @patch("builtins.print")
    def test_pairwise_identity_10009at7524_aa(self, mocked_print):
        testargs = [
            "phykit",
            "pairwise_identity",
            f"{here.parent.parent.parent}/sample_files/10009at7524_aa.aln",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            summary_call(
                "0.8789", "0.9154", "0.8462", "0.95",
                "0.6462", "1.0", "0.0968", "0.0094",
            )
        ]

    @patch("builtins.print")
    def test_pairwise_identity_json_summary(self, mocked_print):
        testargs = [
            "phykit",
            "pairwise_identity",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["verbose"] is False
        assert payload["exclude_gaps"] is False
        assert round(payload["summary"]["mean"], 4) == 0.4833

    @patch("builtins.print")
    def test_pairwise_identity_json_verbose(self, mocked_print):
        testargs = [
            "phykit",
            "pairwise_identity",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-v",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["verbose"] is True
        assert payload["exclude_gaps"] is False
        assert payload["rows"][0] == payload["pairs"][0]
        assert payload["pairs"][0] == {"taxon_a": "1", "taxon_b": "2", "identity": 0.8333}

    @patch("phykit.services.alignment.pairwise_identity.PairwiseIdentity._plot_pairwise_identity_heatmap")
    @patch("builtins.print")
    def test_pairwise_identity_plot(self, mocked_print, mocked_plot):
        testargs = [
            "phykit",
            "pairwise_identity",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--plot",
            "--plot-output",
            "pairwise_identity_heatmap_test.png",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_plot.called
        assert any(
            call.args[0] == "Saved pairwise identity heatmap: pairwise_identity_heatmap_test.png"
            for call in mocked_print.mock_calls
        )

    @patch("phykit.services.alignment.pairwise_identity.PairwiseIdentity._plot_pairwise_identity_heatmap")
    @patch("builtins.print")
    def test_pairwise_identity_plot_json(self, mocked_print, mocked_plot):
        testargs = [
            "phykit",
            "pairwise_identity",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "--plot",
            "--plot-output",
            "pairwise_identity_heatmap_test_json.png",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_plot.called
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["plot_output"] == "pairwise_identity_heatmap_test_json.png"
