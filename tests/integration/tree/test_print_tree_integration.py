from pathlib import Path
import pytest
from mock import patch
import sys

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestPrintTree(object):
    @patch("builtins.print")
    def test_print_tree0(self, mocked_print):
        expected_result = """
             _________ raccoon
            |
            |___ bear
            |
            |      _____ sea_lion
            |  ___|
            | |   |_____ seal
            _|_|
            | |            _____________________________________________________ monkey
            | | __________|
            | ||          |________________________ cat
            |  |
            |  |_________ weasel
            |
            |____________ dog
        """
        testargs = [
            "phykit",
            "print_tree",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()


    @patch("builtins.print")
    def test_print_tree1(self, mocked_print):
        expected_result = """
            , Aspergillus_fischeri_IBT_3003
            |
            , Aspergillus_fischeri_IBT_3007
            |
            , Aspergillus_fischeri_NRRL181.GCF_0001...
            |
            , Aspergillus_fischeri_NRRL4585
            |
            |                                , Aspergillus_fumigatus_Af293
            _|                               ,|
            |                               || Aspergillus_fumigatus_CEA10
            |          _____________________|
            |         |                     | _ Aspergillus_fumigatus_HMR_AF_270
            |_________|                     ||
            |         |                      | Aspergillus_fumigatus_Z5
            |         |
            |         |____ Aspergillus_oerlinghausenensis_CBS139183
            |
            | Aspergillus_fischeri_NRRL4161
        """
        testargs = [
            "phykit",
            "print_tree",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()


    @patch("builtins.print")
    def test_print_tree_wrong_input(self, mocked_print):
        testargs = [
            "phykit",
            "print_tree",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tr",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_print_tree_alias0(self, mocked_print):
        expected_result = """
            , Aspergillus_fischeri_IBT_3003
            |
            , Aspergillus_fischeri_IBT_3007
            |
            , Aspergillus_fischeri_NRRL181.GCF_0001...
            |
            , Aspergillus_fischeri_NRRL4585
            |
            |                                , Aspergillus_fumigatus_Af293
            _|                               ,|
            |                               || Aspergillus_fumigatus_CEA10
            |          _____________________|
            |         |                     | _ Aspergillus_fumigatus_HMR_AF_270
            |_________|                     ||
            |         |                      | Aspergillus_fumigatus_Z5
            |         |
            |         |____ Aspergillus_oerlinghausenensis_CBS139183
            |
            | Aspergillus_fischeri_NRRL4161
        """

        testargs = [
            "phykit",
            "print",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
        

    @patch("builtins.print")
    def test_print_tree_alias1(self, mocked_print):
        expected_result = """
            , Aspergillus_fischeri_IBT_3003
            |
            , Aspergillus_fischeri_IBT_3007
            |
            , Aspergillus_fischeri_NRRL181.GCF_0001...
            |
            , Aspergillus_fischeri_NRRL4585
            |
            |                                , Aspergillus_fumigatus_Af293
            _|                               ,|
            |                               || Aspergillus_fumigatus_CEA10
            |          _____________________|
            |         |                     | _ Aspergillus_fumigatus_HMR_AF_270
            |_________|                     ||
            |         |                      | Aspergillus_fumigatus_Z5
            |         |
            |         |____ Aspergillus_oerlinghausenensis_CBS139183
            |
            | Aspergillus_fischeri_NRRL4161
        """
        
        testargs = [
            "phykit",
            "pt",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_print_tree_remove_branch_lengths_short(self, mocked_print):
        expected_result = """
             ____ Aspergillus_fischeri_IBT_3003
            |
            |     ____ Aspergillus_fischeri_IBT_3007
            |    |
            |____|     ____ Aspergillus_fischeri_NRRL181.GCF_0001...
            |    |    |
            |    |____|     ____ Aspergillus_fischeri_NRRL4585
            |         |    |
            |         |    |               ____ Aspergillus_fumigatus_Af293
            _|         |____|          ____|
            |              |         |    |____ Aspergillus_fumigatus_CEA10
            |              |     ____|
            |              |    |    |     ____ Aspergillus_fumigatus_HMR_AF_270
            |              |____|    |____|
            |                   |         |____ Aspergillus_fumigatus_Z5
            |                   |
            |                   |____ Aspergillus_oerlinghausenensis_CBS139183
            |
            |____ Aspergillus_fischeri_NRRL4161
        """
        
        testargs = [
            "phykit",
            "print",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "-r"
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_print_tree_remove_branch_lengths_long(self, mocked_print):
        expected_result = """
             ____ Aspergillus_fischeri_IBT_3003
            |
            |     ____ Aspergillus_fischeri_IBT_3007
            |    |
            |____|     ____ Aspergillus_fischeri_NRRL181.GCF_0001...
            |    |    |
            |    |____|     ____ Aspergillus_fischeri_NRRL4585
            |         |    |
            |         |    |               ____ Aspergillus_fumigatus_Af293
            _|         |____|          ____|
            |              |         |    |____ Aspergillus_fumigatus_CEA10
            |              |     ____|
            |              |    |    |     ____ Aspergillus_fumigatus_HMR_AF_270
            |              |____|    |____|
            |                   |         |____ Aspergillus_fumigatus_Z5
            |                   |
            |                   |____ Aspergillus_oerlinghausenensis_CBS139183
            |
            |____ Aspergillus_fischeri_NRRL4161
        """
        
        testargs = [
            "phykit",
            "print",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            "--remove"
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
