"""
Integration tests for Kuhner-Felsenstein (branch score) distance.

Expected values cross-validated against R's phangorn::KF.dist().
R validation script: tests/r_validation/validate_kf_distance.R

To reproduce the expected value in R:
    library(phangorn)
    t0 <- read.tree("tests/sample_files/tree_simple.tre")
    t1 <- read.tree("tests/sample_files/tree_simple_other_topology.tre")
    KF.dist(t0, t1)  # 51.9041
"""
from mock import patch, call
from pathlib import Path
import pytest
import sys
import json

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestKFDistance(object):
    @patch("builtins.print")
    def test_kf_distance(self, mocked_print):
        # Cross-validated against R's phangorn::KF.dist() = 51.9041
        testargs = [
            "phykit",
            "kf_distance",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_other_topology.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        output = mocked_print.call_args.args[0]
        plain_kf, normalized_kf = output.split("\t")
        assert float(plain_kf) == 51.9041

    @patch("builtins.print")
    def test_kf_distance_alias_kuhner_felsenstein(self, mocked_print):
        testargs = [
            "phykit",
            "kuhner_felsenstein_distance",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_other_topology.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        output = mocked_print.call_args.args[0]
        plain_kf = float(output.split("\t")[0])
        assert plain_kf == 51.9041

    @patch("builtins.print")
    def test_kf_distance_alias_kf_dist(self, mocked_print):
        testargs = [
            "phykit",
            "kf_dist",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_other_topology.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        output = mocked_print.call_args.args[0]
        plain_kf = float(output.split("\t")[0])
        assert plain_kf == 51.9041

    @patch("builtins.print")
    def test_kf_distance_alias_kf(self, mocked_print):
        testargs = [
            "phykit",
            "kf",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_other_topology.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        output = mocked_print.call_args.args[0]
        plain_kf = float(output.split("\t")[0])
        assert plain_kf == 51.9041

    @patch("builtins.print")
    def test_kf_distance_json(self, mocked_print):
        testargs = [
            "phykit",
            "kf_distance",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_other_topology.tre",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["plain_kf"] == 51.9041
        assert "normalized_kf" in payload

    @patch("builtins.print")
    def test_kf_distance_identical_trees(self, mocked_print):
        testargs = [
            "phykit",
            "kf_distance",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        output = mocked_print.call_args.args[0]
        plain_kf = float(output.split("\t")[0])
        assert plain_kf == 0.0

    @patch("builtins.print")
    def test_kf_distance_different_taxa(self, mocked_print):
        # Trees with different tip sets — should prune to shared taxa.
        # Cross-validated against R's phangorn::KF.dist() after pruning
        # to shared taxa: KF = 51.5134
        testargs = [
            "phykit",
            "kf_distance",
            f"{here.parent.parent.parent}/sample_files/tree_simple.tre",
            f"{here.parent.parent.parent}/sample_files/tree_simple_other_topology_incomplete_taxon_representation.tre",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        output = mocked_print.call_args.args[0]
        plain_kf = float(output.split("\t")[0])
        assert plain_kf == 51.5134
