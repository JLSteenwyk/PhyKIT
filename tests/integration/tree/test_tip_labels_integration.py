from mock import patch, call
from pathlib import Path
import pytest
import sys

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestTipLabels(object):
    @patch("builtins.print")
    def test_tip_labels(self, mocked_print):
        expected_result = "None"
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
        expected_result = "None"
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
        expected_result = "None"
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
        expected_result = "None"
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
        expected_result = "None"
        testargs = [
            "phykit",
            "tip_labels",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tr",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
