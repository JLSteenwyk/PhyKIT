import pytest
import sys
import json
from mock import patch
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestMonophylyCheck(object):
    def test_monophyly_check(self, capsys):
        testargs = [
            "phykit",
            "monophyly_check",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.monophyly_check.true.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert capsys.readouterr().out == "monophyletic\t100\t100\t100\t0.0\n"

    def test_monophyly_checkmissing_tip_name(self, capsys):
        testargs = [
            "phykit",
            "monophyly_check",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.monophyly_check.true.absent_tip_name.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert capsys.readouterr().out == "monophyletic\t100\t100\t100\t0.0\n"

    def test_monophyly_check_not_true(self, capsys):
        testargs = [
            "phykit",
            "monophyly_check",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.monophyly_check.false.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert capsys.readouterr().out == (
            "not_monophyletic\t95.7143\t100\t85\t7.3193\t"
            "Aspergillus_fischeri_IBT_3003;Aspergillus_fischeri_IBT_3007;"
            "Aspergillus_fischeri_NRRL181.GCF_000149645.1_ASM14964v1;"
            "Aspergillus_fischeri_NRRL4585;Aspergillus_oerlinghausenensis_CBS139183\n"
        )

    def test_monophyly_check_check_alias(self, capsys):
        testargs = [
            "phykit",
            "is_monophyletic",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.monophyly_check.true.absent_tip_name.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert capsys.readouterr().out == "monophyletic\t100\t100\t100\t0.0\n"

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

        assert pytest_wrapped_e.type is SystemExit
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

        assert pytest_wrapped_e.type is SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_monophyly_check_json(self, mocked_print):
        testargs = [
            "phykit",
            "monophyly_check",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.tre",
            f"{here.parent.parent.parent}/sample_files/small_Aspergillus_tree.monophyly_check.false.txt",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["rows"][0] == payload["results"][0]
        assert payload["results"][0]["status"] == "not_monophyletic"
        assert round(payload["results"][0]["mean_support"], 4) == 95.7143

    def test_monophyly_check_no_support_values(self, capsys):
        # a tree without internal support values (e.g., a time tree) should
        # still report monophyly, with "NA" in place of support statistics
        testargs = [
            "phykit",
            "monophyly_check",
            f"{here.parent.parent.parent}/sample_files/small_no_support_tree.tre",
            f"{here.parent.parent.parent}/sample_files/small_no_support_tree.monophyly_check.true.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert capsys.readouterr().out == "monophyletic\tNA\tNA\tNA\tNA\n"

    def test_monophyly_check_no_support_values_not_true(self, capsys):
        testargs = [
            "phykit",
            "monophyly_check",
            f"{here.parent.parent.parent}/sample_files/small_no_support_tree.tre",
            f"{here.parent.parent.parent}/sample_files/small_no_support_tree.monophyly_check.false.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert capsys.readouterr().out == (
            "not_monophyletic\tNA\tNA\tNA\tNA\t"
            "Aspergillus_b;Aspergillus_c;Penicillium_e\n"
        )

    @patch("builtins.print")
    def test_monophyly_check_no_support_values_json(self, mocked_print):
        testargs = [
            "phykit",
            "monophyly_check",
            f"{here.parent.parent.parent}/sample_files/small_no_support_tree.tre",
            f"{here.parent.parent.parent}/sample_files/small_no_support_tree.monophyly_check.true.txt",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["results"][0]["status"] == "monophyletic"
        assert payload["results"][0]["mean_support"] is None
        assert payload["results"][0]["offending_taxa"] == []
