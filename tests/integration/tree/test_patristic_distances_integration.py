import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestPatristicDistances(object):
    @patch("builtins.print")
    def test_patristic_distances(self, mocked_print):
        testargs = [
            "phykit",
            "patristic_distances",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("mean: 76.1974"),
            call("median: 49.5888"),
            call("25th percentile: 40.5054"),
            call("75th percentile: 108.1385"),
            call("minimum: 24.0"),
            call("maximum: 152.8813"),
            call("standard deviation: 45.4698"),
            call("variance: 2067.502")
        ]

    @patch("builtins.print")
    def test_patristic_distances_alias(self, mocked_print):
        testargs = [
            "phykit",
            "pd",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("mean: 76.1974"),
            call("median: 49.5888"),
            call("25th percentile: 40.5054"),
            call("75th percentile: 108.1385"),
            call("minimum: 24.0"),
            call("maximum: 152.8813"),
            call("standard deviation: 45.4698"),
            call("variance: 2067.502")
        ]

    @patch("builtins.print")
    def test_patristic_distances_verbose(self, mocked_print):
        testargs = [
            "phykit",
            "patristic_distances",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-v"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("raccoon-bear\t26.0"),
            call("raccoon-sea_lion\t43.4461"),
            call("raccoon-seal\t43.4521"),
            call("raccoon-monkey\t147.4653"),
            call("raccoon-cat\t93.7467"),
            call("raccoon-weasel\t44.8935"),
            call("raccoon-dog\t45.5071"),
            call("bear-sea_lion\t31.047"),
            call("bear-seal\t31.053"),
            call("bear-monkey\t135.0661"),
            call("bear-cat\t81.3475"),
            call("bear-weasel\t32.4944"),
            call("bear-dog\t33.108"),
            call("sea_lion-seal\t24.0"),
            call("sea_lion-monkey\t143.0726"),
            call("sea_lion-cat\t89.354"),
            call("sea_lion-weasel\t40.5009"),
            call("sea_lion-dog\t48.8621"),
            call("seal-monkey\t143.0786"),
            call("seal-cat\t89.36"),
            call("seal-weasel\t40.5069"),
            call("seal-dog\t48.8681"),
            call("monkey-cat\t148.0"),
            call("monkey-weasel\t140.3308"),
            call("monkey-dog\t152.8813"),
            call("cat-weasel\t86.6122"),
            call("cat-dog\t99.1627"),
            call("weasel-dog\t50.3095"),
        ]

    @patch("builtins.print")
    def test_patristic_distances_wrong_input(self, mocked_print):
        testargs = [
            "phykit",
            "pd",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tr",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
