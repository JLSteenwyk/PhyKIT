import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


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
            call("mean: 0.4833"),
            call("median: 0.5"),
            call("25th percentile: 0.3333"),
            call("75th percentile: 0.6667"),
            call("minimum: 0.1667"),
            call("maximum: 0.8333"),
            call("standard deviation: 0.2284"),
            call("variance: 0.0522")
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
            call("mean: 0.7593"),
            call("median: 0.7778"),
            call("25th percentile: 0.6944"),
            call("75th percentile: 0.8611"),
            call("minimum: 0.5556"),
            call("maximum: 0.8889"),
            call("standard deviation: 0.1299"),
            call("variance: 0.0169")
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
            call("mean: 0.8333"),
            call("median: 0.8333"),
            call("25th percentile: 0.6667"),
            call("75th percentile: 1.0"),
            call("minimum: 0.6667"),
            call("maximum: 1.0"),
            call("standard deviation: 0.1826"),
            call("variance: 0.0333")
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
            call("mean: 0.8333"),
            call("median: 0.8333"),
            call("25th percentile: 0.6667"),
            call("75th percentile: 1.0"),
            call("minimum: 0.6667"),
            call("maximum: 1.0"),
            call("standard deviation: 0.1826"),
            call("variance: 0.0333")
        ]

    @patch("builtins.print")
    def test_pairwise_identity_incorrect_input_file(self, mocked_print):
        testargs = [
            "phykit",
            "pairwise_identity",
            f"{here.parent.parent.parent}/sample_files/test_trees.txt",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_pairwise_identity_verbose(self, mocked_print):
        testargs = [
            "phykit",
            "pairwise_identity",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-v",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("1\t2\t0.8333"),
            call("1\t3\t0.5"),
            call("1\t4\t0.1667"),
            call("1\t5\t0.1667"),
            call("2\t3\t0.6667"),
            call("2\t4\t0.3333"),
            call("2\t5\t0.3333"),
            call("3\t4\t0.6667"),
            call("3\t5\t0.5"),
            call("4\t5\t0.6667"),
        ]

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
            call("mean: 0.8136"),
            call("median: 0.8423"),
            call("25th percentile: 0.8096"),
            call("75th percentile: 0.8692"),
            call("minimum: 0.6192"),
            call("maximum: 0.9269"),
            call("standard deviation: 0.0831"),
            call("variance: 0.0069")
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
            call("mean: 0.8136"),
            call("median: 0.8423"),
            call("25th percentile: 0.8096"),
            call("75th percentile: 0.8692"),
            call("minimum: 0.6192"),
            call("maximum: 0.9269"),
            call("standard deviation: 0.0831"),
            call("variance: 0.0069")
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
            call("mean: 0.8789"),
            call("median: 0.9154"),
            call("25th percentile: 0.8462"),
            call("75th percentile: 0.95"),
            call("minimum: 0.6462"),
            call("maximum: 1.0"),
            call("standard deviation: 0.0968"),
            call("variance: 0.0094")
        ]
