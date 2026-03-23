import json
import os
import sys
import tempfile

import pytest
from argparse import Namespace

from phykit.errors import PhykitUserError

TREE_SIMPLE = "tests/sample_files/tree_simple.tre"
GENE_TREES = "tests/sample_files/gene_trees_simple.nwk"


def _make_args(**kwargs):
    defaults = dict(
        tree=TREE_SIMPLE,
        gene_trees=GENE_TREES,
        support=None,
        alpha=0.05,
        json=False,
        plot_output=None,
    )
    defaults.update(kwargs)
    return Namespace(**defaults)


def _make_svc(**kwargs):
    from phykit.services.tree.hybridization import Hybridization
    return Hybridization(_make_args(**kwargs))


class TestPerBranchCounts:
    def test_per_branch_counts(self):
        svc = _make_svc()
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result, _ = svc._count_topologies(species_tree, gene_trees)
        assert isinstance(result, dict)
        assert len(result) > 0
        for branch_key, data in result.items():
            assert "split" in data
            assert "n_concordant" in data
            assert "n_alt1" in data
            assert "n_alt2" in data
            total = data["n_concordant"] + data["n_alt1"] + data["n_alt2"]
            assert total <= 10


class TestAsymmetryRatioRange:
    def test_asymmetry_ratio_range(self):
        svc = _make_svc()
        # Test with various inputs
        for n1, n2 in [(5, 5), (9, 1), (1, 9), (3, 7), (0, 10), (10, 0)]:
            result = svc._test_asymmetry(n1, n2)
            if result["asymmetry_ratio"] is not None:
                assert 0.5 <= result["asymmetry_ratio"] <= 1.0

    def test_asymmetry_ratio_none_for_zero(self):
        svc = _make_svc()
        result = svc._test_asymmetry(0, 0)
        assert result["asymmetry_ratio"] is None


class TestReticulationCountNonneg:
    def test_reticulation_count_nonneg(self, capsys):
        svc = _make_svc()
        svc.run()
        captured = capsys.readouterr()
        # Parse the reticulation count from text output
        for line in captured.out.splitlines():
            if "Estimated reticulation events:" in line:
                count = int(line.split(":")[-1].strip())
                assert count >= 0
                break


class TestSupportThreshold:
    def test_support_threshold(self):
        svc_no_thresh = _make_svc()
        svc_with_thresh = _make_svc(support=99999.0)

        species_tree = svc_no_thresh.read_tree_file()
        gene_trees = svc_no_thresh._parse_gene_trees(GENE_TREES)

        result_no, _ = svc_no_thresh._count_topologies(species_tree, gene_trees)
        result_with, _ = svc_with_thresh._count_topologies(species_tree, gene_trees,
                                                           support_threshold=99999.0)

        # With a very high threshold, all gene tree bipartitions should be
        # collapsed, so concordant counts should differ or be zero
        # (gene trees in sample don't have support values, so confidence is None
        #  and the threshold check is skipped; counts should be equal)
        # This verifies the threshold logic doesn't crash.
        assert isinstance(result_with, dict)

    def test_support_threshold_with_supported_trees(self, tmp_path):
        """Test with gene trees that have support values."""
        from Bio import Phylo
        from io import StringIO

        sp_file = tmp_path / "species.tre"
        sp_file.write_text("((a,b),(c,d));")

        # Gene trees with support values
        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = [
            "((a,b)90,(c,d)90);",
            "((a,b)50,(c,d)50);",
            "((a,c)90,(b,d)90);",
        ]
        gt_file.write_text("\n".join(gt_lines))

        svc_low = _make_svc(tree=str(sp_file), gene_trees=str(gt_file), support=40.0)
        svc_high = _make_svc(tree=str(sp_file), gene_trees=str(gt_file), support=80.0)

        species_tree = svc_low.read_tree_file()
        gt_low = svc_low._parse_gene_trees(str(gt_file))
        gt_high = svc_high._parse_gene_trees(str(gt_file))

        result_low, _ = svc_low._count_topologies(species_tree, gt_low, support_threshold=40.0)
        result_high, _ = svc_high._count_topologies(species_tree, gt_high, support_threshold=80.0)

        # Both should produce results without error
        assert isinstance(result_low, dict)
        assert isinstance(result_high, dict)


class TestTextOutput:
    def test_text_output(self, capsys):
        svc = _make_svc()
        svc.run()
        captured = capsys.readouterr()
        assert "Hybridization Analysis" in captured.out
        assert "======================" in captured.out
        assert "Species tree branches:" in captured.out
        assert "Gene trees:" in captured.out
        assert "Estimated reticulation events:" in captured.out
        assert "Non-significant branches:" in captured.out


class TestJsonOutput:
    def test_json_output(self, capsys):
        svc = _make_svc(json=True)
        svc.run()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert "n_branches" in data
        assert "n_gene_trees" in data
        assert "alpha" in data
        assert "n_reticulations" in data
        assert "branches" in data
        assert isinstance(data["branches"], list)
        assert len(data["branches"]) > 0

        # Check branch entry has expected keys
        branch = data["branches"][0]
        for key in ["taxa", "n_concordant", "n_alt1", "n_alt2",
                     "gcf", "asymmetry_ratio", "p_value", "fdr_p",
                     "significant", "hybrid_score", "favored_alt"]:
            assert key in branch, f"Missing key: {key}"

    def test_json_n_gene_trees(self, capsys):
        svc = _make_svc(json=True)
        svc.run()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert data["n_gene_trees"] == 10


class TestPlotCreated:
    def test_plot_created(self, tmp_path):
        output = str(tmp_path / "test_hybrid.png")
        svc = _make_svc(plot_output=output)
        svc.run()
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0


class TestCircularPlot:
    def test_circular_plot(self, tmp_path):
        output = str(tmp_path / "test_hybrid_circular.png")
        args = _make_args(plot_output=output)
        args.circular = True
        from phykit.services.tree.hybridization import Hybridization
        svc = Hybridization(args)
        svc.run()
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0


class TestNoSignificantBranches:
    def test_no_significant_branches(self, capsys):
        # With a very low alpha, nothing should be significant
        svc = _make_svc(alpha=1e-100)
        svc.run()
        captured = capsys.readouterr()
        assert "Estimated reticulation events: 0" in captured.out

    def test_no_significant_json(self, capsys):
        svc = _make_svc(alpha=1e-100, json=True)
        svc.run()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert data["n_reticulations"] == 0
        for branch in data["branches"]:
            assert branch["significant"] is False
            assert branch["hybrid_score"] == 0.0


class TestProcessArgs:
    def test_default_args(self):
        svc = _make_svc()
        assert svc.gene_trees_path == GENE_TREES
        assert svc.support_threshold is None
        assert svc.alpha == 0.05
        assert svc.json_output is False
        assert svc.plot_output is None

    def test_custom_args(self):
        svc = _make_svc(support=70.0, alpha=0.01)
        assert svc.support_threshold == 70.0
        assert svc.alpha == 0.01


class TestFDR:
    def test_empty_list(self):
        from phykit.services.tree.hybridization import Hybridization
        assert Hybridization._fdr([]) == []

    def test_single_pvalue(self):
        from phykit.services.tree.hybridization import Hybridization
        result = Hybridization._fdr([0.03])
        assert result == [0.03]

    def test_known_correction(self):
        from phykit.services.tree.hybridization import Hybridization
        pvals = [0.01, 0.04, 0.03]
        result = Hybridization._fdr(pvals)
        assert result[0] == pytest.approx(0.03)
        assert result[1] == pytest.approx(0.04)
        assert result[2] == pytest.approx(0.04)


class TestCLI:
    def test_hybridization_alias(self, capsys):
        from phykit.phykit import Phykit
        sys.argv = [
            "phykit", "hybrid",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
        ]
        Phykit()
        captured = capsys.readouterr()
        assert "Hybridization Analysis" in captured.out

    def test_reticulation_alias(self, capsys):
        from phykit.phykit import Phykit
        sys.argv = [
            "phykit", "reticulation",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
        ]
        Phykit()
        captured = capsys.readouterr()
        assert "Hybridization Analysis" in captured.out
