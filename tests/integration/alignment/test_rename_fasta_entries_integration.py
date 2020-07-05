import pytest
import sys
from mock import patch, call
from pathlib import Path
from textwrap import dedent

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestRenameFastaEntries(object):
    @patch("builtins.print")
    def test_rename_fasta_entries(self, mocked_print):
        testargs = [
            "phykit",
            "rename_fasta_entries",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-i",
            f"{here.parent.parent.parent}/sample_files/simple_fasta_idmap.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/simple.fa.renamed.fa", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{here.parent.parent.parent}/sample_files/simple.fa.renamed.fa", "r") as out_renamed:
            out_renamed_content = out_renamed.read()

        assert expected_fa_content == out_renamed_content

    @patch("builtins.print")
    def test_rename_fasta_entries_outputting(self, mocked_print):
        testargs = [
            "phykit",
            "rename_fasta_entries",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-i",
            f"{here.parent.parent.parent}/sample_files/simple_fasta_idmap.txt",
            "-o",
            f"{here.parent.parent.parent}/sample_files/simple_customized_output_name.fa"
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/simple_customized_output_name.fa", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{here.parent.parent.parent}/sample_files/simple.fa.renamed.fa", "r") as out_renamed:
            out_renamed_content = out_renamed.read()

        assert expected_fa_content == out_renamed_content

    @patch("builtins.print")
    def test_rename_fasta_entries_alias(self, mocked_print):
        testargs = [
            "phykit",
            "rename_fasta",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-i",
            f"{here.parent.parent.parent}/sample_files/simple_fasta_idmap.txt",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        with open(f"{here.parent.parent}/expected/simple.fa.renamed.fa", "r") as expected_fa:
            expected_fa_content = expected_fa.read()

        with open(f"{here.parent.parent.parent}/sample_files/simple.fa.renamed.fa", "r") as out_renamed:
            out_renamed_content = out_renamed.read()

        assert expected_fa_content == out_renamed_content

    @patch("builtins.print")
    def test_rename_fasta_entries_incorrect_input_alignment_file(self, mocked_print):
        testargs = [
            "phykit",
            "rename_fasta_entries",
            f"{here.parent.parent.parent}/sample_files/simple.f",
            "-i",
            f"{here.parent.parent.parent}/sample_files/simple_fasta_idmap.txt",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    @patch("builtins.print")
    def test_rename_fasta_entries_incorrect_input_idmap_file(self, mocked_print):
        testargs = [
            "phykit",
            "rename_fasta_entries",
            f"{here.parent.parent.parent}/sample_files/simple.fa",
            "-i",
            f"{here.parent.parent.parent}/sample_files/simple_fasta_idmap.tx",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            Phykit()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2