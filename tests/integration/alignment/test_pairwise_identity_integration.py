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
            call("mean: 0.4667"),
            call("median: 0.5"),
            call("25th percentile: 0.3333"),
            call("75th percentile: 0.625"),
            call("minimum: 0.1667"),
            call("maximum: 0.8333"),
            call("standard deviation: 0.2194"),
            call("variance: 0.0481")
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
            call("1-2\t0.8333333333333334"),
            call("1-3\t0.5"),
            call("1-4\t0.16666666666666666"),
            call("1-5\t0.16666666666666666"),
            call("2-3\t0.6666666666666666"),
            call("2-4\t0.3333333333333333"),
            call("2-5\t0.3333333333333333"),
            call("3-4\t0.6666666666666666"),
            call("3-5\t0.5"),
            call("4-5\t0.5"),
        ]