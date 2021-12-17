import pytest
import sys
from math import isclose
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestTree(object):
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
            call("monophyletic\t100\t100\t100\t0.0\tAspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5"),
            call("not_monophyletic\t95.7143\t100\t85\t7.3193\tAspergillus_fumigatus_Af293;Aspergillus_fischeri_IBT_3007"),
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
            call("monophyletic\t100\t100\t100\t0.0\tAspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5"),
            call("not_monophyletic\t95.7143\t100\t85\t7.3193\tAspergillus_fumigatus_Af293;Aspergillus_fischeri_IBT_3007"),
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
            call("monophyletic\t100\t100\t100\t0.0\tAspergillus_fumigatus_Af293;Aspergillus_fumigatus_CEA10;Aspergillus_fumigatus_HMR_AF_270;Aspergillus_fumigatus_Z5"),
            call("not_monophyletic\t95.7143\t100\t85\t7.3193\tAspergillus_fumigatus_Af293;Aspergillus_fischeri_IBT_3007"),
            call("insufficient_taxon_representation"),
        ]

    @patch("builtins.print")
    def test_hidden_paralogy_check_file_error0(self, mocked_print):
        testargs = [
            "phykit",
            "clan_check",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.t",
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
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.hidden_paralogy_check.t",
        ]
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2