from mock import patch, call
from pathlib import Path
import pytest
import sys
import json

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestTipLabels(object):
    @patch("builtins.print")
    def test_tip_labels(self, mocked_print):
        testargs = [
            "phykit",
            "tip_labels",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call('raccoon\nbear\nsea_lion\nseal\nmonkey\ncat\nweasel\ndog')
        ]

    @patch("builtins.print")
    def test_tip_labels_alias0(self, mocked_print):
        testargs = [
            "phykit",
            "labels",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call('raccoon\nbear\nsea_lion\nseal\nmonkey\ncat\nweasel\ndog')
        ]

    @patch("builtins.print")
    def test_tip_labels_alias1(self, mocked_print):
        testargs = [
            "phykit",
            "tree_labels",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call('raccoon\nbear\nsea_lion\nseal\nmonkey\ncat\nweasel\ndog')
        ]

    @patch("builtins.print")
    def test_tip_labels_alias2(self, mocked_print):
        testargs = [
            "phykit",
            "tl",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [
            call('raccoon\nbear\nsea_lion\nseal\nmonkey\ncat\nweasel\ndog')
        ]

    @patch("builtins.print")
    def test_tip_labels_incorrect_file_path(self, mocked_print):
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 1

    @patch("builtins.print")
    def test_tip_labels_json(self, mocked_print):
        testargs = [
            "phykit",
            "tip_labels",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["rows"][0] == {"taxon": "raccoon"}
        assert payload["tips"] == [
            "raccoon", "bear", "sea_lion", "seal", "monkey", "cat", "weasel", "dog"
        ]
