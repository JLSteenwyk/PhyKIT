import pytest
import sys
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestHiddenParalogyCheck(object):
    @patch("builtins.print")
    def test_hidden_paralogy_check(self, mocked_print):
        testargs = [
            "phykit",
            "hidden_paralogy_check",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-c",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.hidden_paralogy_check.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("monophyletic"),
            call("not_monophyletic"),
            call("insufficient_taxon_representation"),
        ]

    @patch("builtins.print")
    def test_hidden_paralogy_check_long(self, mocked_print):
        testargs = [
            "phykit",
            "hidden_paralogy_check",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "--clade",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.hidden_paralogy_check.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("monophyletic"),
            call("not_monophyletic"),
            call("insufficient_taxon_representation"),
        ]

    @patch("builtins.print")
    def test_hidden_paralogy_check_alias(self, mocked_print):
        testargs = [
            "phykit",
            "clan_check",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-c",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.hidden_paralogy_check.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("monophyletic"),
            call("not_monophyletic"),
            call("insufficient_taxon_representation"),
        ]

    @patch("builtins.print")
    def test_hidden_paralogy_check_file_error0(self, mocked_print):
        testargs = [
            "phykit",
            "clan_check",
            "file doesn't exist",
            "-c",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.hidden_paralogy_check.txt",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
    
    @patch("builtins.print")
    def test_hidden_paralogy_check_file_error1(self, mocked_print):
        testargs = [
            "phykit",
            "clan_check",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-c",
            "file doesn't exist",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
