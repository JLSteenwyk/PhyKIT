"""
Integration tests for independent_contrasts (Felsenstein's PIC).

Cross-validated against R's ape::pic(). The sum of squared contrasts
is invariant to polytomy resolution order and matches R exactly.
R validation script: tests/r_validation/validate_pic.R

To reproduce in R:
    library(ape)
    tree <- multi2di(read.tree("tests/sample_files/tree_simple.tre"))
    traits <- read.delim("tests/sample_files/tree_simple_traits.tsv",
                         comment.char="#", header=FALSE)
    pic(setNames(traits$V2, traits$V1), tree)
"""
from mock import patch
from pathlib import Path
import json
import sys

import pytest

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestIndependentContrasts:
    @patch("builtins.print")
    def test_pic_text_output(self, mocked_print):
        testargs = [
            "phykit",
            "independent_contrasts",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_traits.tsv",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        output = "\n".join(
            call.args[0] if call.args else ""
            for call in mocked_print.call_args_list
        )
        assert "Number of contrasts: 7" in output

    @patch("builtins.print")
    def test_pic_json_sum_of_squares_matches_r(self, mocked_print):
        # R's ape::pic() gives sum(contrasts^2) = 0.307253
        testargs = [
            "phykit",
            "independent_contrasts",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_traits.tsv",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_contrasts"] == 7
        ss = sum(c["contrast"] ** 2 for c in payload["contrasts"])
        assert abs(ss - 0.307253) < 0.001

    @patch("builtins.print")
    def test_pic_alias_pic(self, mocked_print):
        testargs = [
            "phykit",
            "pic",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_traits.tsv",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_contrasts"] == 7

    @patch("builtins.print")
    def test_pic_alias_phylo_contrasts(self, mocked_print):
        testargs = [
            "phykit",
            "phylo_contrasts",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/tree_simple_traits.tsv",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_contrasts"] == 7
