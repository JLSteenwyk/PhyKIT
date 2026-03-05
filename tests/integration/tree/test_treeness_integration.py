from mock import patch, call
from pathlib import Path
import pytest
import sys
import json

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

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 1

    @patch("builtins.print")
    def test_treeness_json(self, mocked_print):
        testargs = [
            "phykit",
            "treeness",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload == {"treeness": 0.7695}
