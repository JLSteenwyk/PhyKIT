import pytest
import sys
from mock import patch, call
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestMonophylyCheck(object):
    @patch("builtins.print")
    def test_monophyly_check(self, mocked_print):
        testargs = [
            "phykit",
            "monophyly_check",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.monophyly_check.true.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("monophyletic\t100\t100\t100\t0.0"),
        ]

    @patch("builtins.print")
    def test_monophyly_checkmissing_tip_name(self, mocked_print):
        testargs = [
            "phykit",
            "monophyly_check",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.monophyly_check.true.absent_tip_name.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("monophyletic\t100\t100\t100\t0.0"),
        ]

    @patch("builtins.print")
    def test_monophyly_check_not_true(self, mocked_print):
        testargs = [
            "phykit",
            "monophyly_check",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.monophyly_check.false.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("not_monophyletic\t95.7143\t100\t85\t7.3193\tAspergillus_fischeri_IBT_3003;Aspergillus_fischeri_IBT_3007;Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1;Aspergillus_fischeri_NRRL4585;Aspergillus_oerlinghausenensis_CBS139183"),
        ]

    @patch("builtins.print")
    def test_monophyly_check_check_alias(self, mocked_print):
        testargs = [
            "phykit",
            "is_monophyletic",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.monophyly_check.true.absent_tip_name.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert mocked_print.mock_calls == [
            call("monophyletic\t100\t100\t100\t0.0"),
        ]

    @patch("builtins.print")
    def test_monophyly_check_file_error0(self, mocked_print):
        testargs = [
            "phykit",
            "is_monophyletic",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.t",
            f"{here.parent.parent.parent}/sample_files/wrong_file",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
    
    @patch("builtins.print")
    def test_monophyly_check_file_error1(self, mocked_print):
        testargs = [
            "phykit",
            "is_monophyletic",
            f"{here.parent.parent.parent}/sample_files/wrong_file",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.monophyly_check.true.absent_tip_name.txt",
        ]
        with patch.object(sys, "argv", testargs):
            with pytest.raises(SystemExit) as pytest_wrapped_e:
                Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2
