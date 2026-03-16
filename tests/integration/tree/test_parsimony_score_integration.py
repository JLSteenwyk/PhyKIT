"""
Integration tests for parsimony_score (Fitch parsimony).

Cross-validated against R's phangorn::parsimony() which gives a
Fitch score of 4 for tree_simple.tre + tree_simple_alignment.fa.
R validation script: tests/r_validation/validate_parsimony.R

To reproduce in R:
    library(ape); library(phangorn)
    tree <- multi2di(read.tree("tests/sample_files/tree_simple.tre"))
    aln <- read.phyDat("tests/sample_files/tree_simple_alignment.fa",
                        format="fasta")
    parsimony(tree, aln, method="fitch")  # 4
"""
from mock import patch
from pathlib import Path
import json
import sys

import pytest

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestParsimonyScore:
    @patch("builtins.print")
    def test_parsimony_score_matches_r(self, mocked_print):
        # R's phangorn::parsimony() gives 4 for this tree+alignment
        testargs = [
            "phykit",
            "parsimony_score",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-a", f"{here.parent.parent.parent}/sample_files/tree_simple_alignment.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_args_list[0].args[0] == 4

    @patch("builtins.print")
    def test_parsimony_alias_pars(self, mocked_print):
        testargs = [
            "phykit",
            "pars",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-a", f"{here.parent.parent.parent}/sample_files/tree_simple_alignment.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_args_list[0].args[0] == 4

    @patch("builtins.print")
    def test_parsimony_alias_parsimony(self, mocked_print):
        testargs = [
            "phykit",
            "parsimony",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-a", f"{here.parent.parent.parent}/sample_files/tree_simple_alignment.fa",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert mocked_print.call_args_list[0].args[0] == 4

    @patch("builtins.print")
    def test_parsimony_json(self, mocked_print):
        testargs = [
            "phykit",
            "parsimony_score",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-a", f"{here.parent.parent.parent}/sample_files/tree_simple_alignment.fa",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        # Cross-validated against R: phangorn::parsimony() = 4
        assert payload["parsimony_score"] == 4
        assert payload["alignment_length"] == 16
        assert payload["n_taxa"] == 8
