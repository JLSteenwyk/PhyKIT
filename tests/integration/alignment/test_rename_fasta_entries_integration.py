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