import json
import os
import sys

import pytest

TREE_SIMPLE = "tests/sample_files/tree_simple.tre"
GENE_TREES = "tests/sample_files/gene_trees_simple.nwk"


class TestBasicInvocation:
    def test_basic_invocation(self, capsys):
        from phykit.phykit import Phykit
        sys.argv = [
            "phykit", "hybridization",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
        ]
        Phykit()
        captured = capsys.readouterr()
        assert "Hybridization Analysis" in captured.out
        assert "Estimated reticulation events:" in captured.out
        assert "Non-significant branches:" in captured.out


class TestAliasHybrid:
    def test_alias_hybrid(self, capsys):
        from phykit.phykit import Phykit
        sys.argv = [
            "phykit", "hybrid",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
        ]
        Phykit()
        captured = capsys.readouterr()
        assert "Hybridization Analysis" in captured.out


class TestAliasReticulation:
    def test_alias_reticulation(self, capsys):
        from phykit.phykit import Phykit
        sys.argv = [
            "phykit", "reticulation",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
        ]
        Phykit()
        captured = capsys.readouterr()
        assert "Hybridization Analysis" in captured.out


class TestJsonOutput:
    def test_json_output(self, capsys):
        from phykit.phykit import Phykit
        sys.argv = [
            "phykit", "hybridization",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
            "--json",
        ]
        Phykit()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert "n_branches" in data
        assert "n_gene_trees" in data
        assert "n_reticulations" in data
        assert "branches" in data
        assert isinstance(data["branches"], list)
        assert data["n_gene_trees"] == 10


class TestWithSupportThreshold:
    def test_with_support_threshold(self, capsys):
        from phykit.phykit import Phykit
        sys.argv = [
            "phykit", "hybridization",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
            "--support", "70",
        ]
        Phykit()
        captured = capsys.readouterr()
        assert "Hybridization Analysis" in captured.out


class TestWithPlot:
    def test_with_plot(self, tmp_path, capsys):
        output = str(tmp_path / "hybrid_integration.png")
        from phykit.phykit import Phykit
        sys.argv = [
            "phykit", "hybridization",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
            "--plot", output,
        ]
        Phykit()
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0


class TestWithAlphaArg:
    def test_with_alpha(self, capsys):
        from phykit.phykit import Phykit
        sys.argv = [
            "phykit", "hybridization",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
            "--alpha", "0.01",
        ]
        Phykit()
        captured = capsys.readouterr()
        assert "Hybridization Analysis" in captured.out
