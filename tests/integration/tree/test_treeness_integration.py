from mock import patch, call
from pathlib import Path
import pytest
import sys

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestTreeness(object):
    @patch("builtins.print")
    def test_treeness0(self, mocked_print):
        expected_result = 0.1899
        testargs = [
            "phykit",
            "treeness",
            f"{here.parent.parent.parent}/sample_files/Yeasts_2832_eMRC_reference_renamed.tree",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_treeness1(self, mocked_print):
        expected_result = 0.7695
        testargs = [
            "phykit",
            "treeness",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_treeness_alias(self, mocked_print):
        expected_result = 0.7695
        testargs = [
            "phykit",
            "tness",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.mock_calls == [call(expected_result)]

    @patch("builtins.print")
    def test_treeness_incorrect_file_path(self, mocked_print):
        testargs = [
            "phykit",
            "tness",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tr",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
