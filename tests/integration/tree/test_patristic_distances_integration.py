import pytest
import sys
import json
from mock import patch, call
from pathlib import Path

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

        assert mocked_print.call_args.args[0].splitlines() == [
            "raccoon\tbear\t26.0",
            "raccoon\tsea_lion\t43.4461",
            "raccoon\tseal\t43.4521",
            "raccoon\tmonkey\t147.4653",
            "raccoon\tcat\t93.7467",
            "raccoon\tweasel\t44.8935",
            "raccoon\tdog\t45.5071",
            "bear\tsea_lion\t31.047",
            "bear\tseal\t31.053",
            "bear\tmonkey\t135.0661",
            "bear\tcat\t81.3475",
            "bear\tweasel\t32.4944",
            "bear\tdog\t33.108",
            "sea_lion\tseal\t24.0",
            "sea_lion\tmonkey\t143.0726",
            "sea_lion\tcat\t89.354",
            "sea_lion\tweasel\t40.5009",
            "sea_lion\tdog\t48.8621",
            "seal\tmonkey\t143.0786",
            "seal\tcat\t89.36",
            "seal\tweasel\t40.5069",
            "seal\tdog\t48.8681",
            "monkey\tcat\t148.0",
            "monkey\tweasel\t140.3308",
            "monkey\tdog\t152.8813",
            "cat\tweasel\t86.6122",
            "cat\tdog\t99.1627",
            "weasel\tdog\t50.3095",
        ]

    @patch("builtins.print")
    def test_patristic_distances_wrong_input(self, mocked_print):
        testargs = [
            "phykit",
            "patristic_distances",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tr",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_patristic_distances_json_summary(self, mocked_print):
        testargs = [
            "phykit",
            "patristic_distances",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["verbose"] is False
        assert round(payload["summary"]["mean"], 4) == 76.1974

    @patch("builtins.print")
    def test_patristic_distances_json_verbose(self, mocked_print):
        testargs = [
            "phykit",
            "patristic_distances",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-v",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["verbose"] is True
        assert payload["rows"][0] == payload["pairs"][0]
        assert payload["pairs"][0] == {
            "taxon_a": "raccoon",
            "taxon_b": "bear",
            "patristic_distance": 26.0,
        }
