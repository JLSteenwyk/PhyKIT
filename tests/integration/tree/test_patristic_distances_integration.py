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
            call("raccoon\tbear\t26.0"),
            call("raccoon\tsea_lion\t43.4461"),
            call("raccoon\tseal\t43.4521"),
            call("raccoon\tmonkey\t147.4653"),
            call("raccoon\tcat\t93.7467"),
            call("raccoon\tweasel\t44.8935"),
            call("raccoon\tdog\t45.5071"),
            call("bear\tsea_lion\t31.047"),
            call("bear\tseal\t31.053"),
            call("bear\tmonkey\t135.0661"),
            call("bear\tcat\t81.3475"),
            call("bear\tweasel\t32.4944"),
            call("bear\tdog\t33.108"),
            call("sea_lion\tseal\t24.0"),
            call("sea_lion\tmonkey\t143.0726"),
            call("sea_lion\tcat\t89.354"),
            call("sea_lion\tweasel\t40.5009"),
            call("sea_lion\tdog\t48.8621"),
            call("seal\tmonkey\t143.0786"),
            call("seal\tcat\t89.36"),
            call("seal\tweasel\t40.5069"),
            call("seal\tdog\t48.8681"),
            call("monkey\tcat\t148.0"),
            call("monkey\tweasel\t140.3308"),
            call("monkey\tdog\t152.8813"),
            call("cat\tweasel\t86.6122"),
            call("cat\tdog\t99.1627"),
            call("weasel\tdog\t50.3095"),
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
