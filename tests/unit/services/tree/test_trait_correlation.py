import json
import os
import tempfile

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.trait_correlation import TraitCorrelation
from phykit.helpers.trait_parsing import parse_multi_trait_file
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


def _make_args(output_path, **overrides):
    defaults = dict(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        output=output_path,
        alpha=0.05,
        cluster=False,
        gene_trees=None,
        json=False,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


def _build_service(output_path, **overrides):
    args = _make_args(output_path, **overrides)
    return TraitCorrelation(args)


class TestTraitCorrelation:
    def test_correlation_matrix_symmetric(self, tmp_path):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        svc.run()

        # Re-run internal logic to inspect the matrix
        corr, pmat, trait_names, n = self._compute_correlation(svc)
        for i in range(len(trait_names)):
            for j in range(len(trait_names)):
                assert abs(corr[i, j] - corr[j, i]) < 1e-10, \
                    f"corr[{i},{j}] != corr[{j},{i}]"

    def test_correlation_diagonal_is_one(self, tmp_path):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        svc.run()

        corr, pmat, trait_names, n = self._compute_correlation(svc)
        for i in range(len(trait_names)):
            assert abs(corr[i, i] - 1.0) < 1e-10

    def test_correlation_range(self, tmp_path):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        svc.run()

        corr, pmat, trait_names, n = self._compute_correlation(svc)
        assert np.all(corr >= -1.0 - 1e-10)
        assert np.all(corr <= 1.0 + 1e-10)

    def test_p_values_diagonal_is_one(self, tmp_path):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        svc.run()

        corr, pmat, trait_names, n = self._compute_correlation(svc)
        for i in range(len(trait_names)):
            assert abs(pmat[i, i] - 1.0) < 1e-10

    def test_p_values_nonnegative(self, tmp_path):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        svc.run()

        corr, pmat, trait_names, n = self._compute_correlation(svc)
        assert np.all(pmat >= 0.0)

    def test_creates_png(self, tmp_path):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        svc.run()
        assert os.path.exists(out)

    def test_cluster_mode(self, tmp_path):
        out = str(tmp_path / "corr_cluster.png")
        svc = _build_service(out, cluster=True)
        svc.run()
        assert os.path.exists(out)

    def test_json_output(self, tmp_path, capsys):
        out = str(tmp_path / "corr.png")
        args = _make_args(out, json=True)
        svc = TraitCorrelation(args)
        svc.run()

        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert "n_taxa" in payload
        assert "n_traits" in payload
        assert "trait_names" in payload
        assert "correlation_matrix" in payload
        assert "p_value_matrix" in payload
        assert "significant_pairs" in payload
        assert "output_file" in payload
        assert payload["n_traits"] == 3
        assert len(payload["correlation_matrix"]) == 3
        assert len(payload["correlation_matrix"][0]) == 3

    def test_significance_stars(self, tmp_path):
        out = str(tmp_path / "corr.png")
        svc = _build_service(out)
        # Test the stars function
        assert svc._significance_stars(0.0001) == "***"
        assert svc._significance_stars(0.005) == "**"
        assert svc._significance_stars(0.02) == "*"
        assert svc._significance_stars(0.1) == ""

    def test_text_output_format(self, tmp_path, capsys):
        out = str(tmp_path / "corr.png")
        args = _make_args(out, json=False)
        svc = TraitCorrelation(args)
        svc.run()

        captured = capsys.readouterr()
        text = captured.out
        assert "Phylogenetic Trait Correlation" in text
        assert "body_mass" in text
        assert "brain_size" in text
        assert "longevity" in text
        assert "Taxa:" in text
        assert "Traits:" in text
        assert "Significant pairs" in text

    def _compute_correlation(self, svc):
        """Helper: re-compute the correlation matrix from a built service."""
        from phykit.services.tree.vcv_utils import build_vcv_matrix
        from scipy.stats import t as t_dist

        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(
            svc.trait_data_path, tree_tips, min_columns=3
        )
        ordered_names = sorted(traits.keys())
        vcv = build_vcv_matrix(tree, ordered_names)

        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])

        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Y_centered = Y - a_hat

        phylo_cov = (Y_centered.T @ C_inv @ Y_centered) / (n - 1)
        std_devs = np.sqrt(np.diag(phylo_cov))
        corr_matrix = phylo_cov / np.outer(std_devs, std_devs)
        np.fill_diagonal(corr_matrix, 1.0)

        p_matrix = np.ones((p, p))
        for i in range(p):
            for j in range(i + 1, p):
                r = corr_matrix[i, j]
                if abs(r) >= 1.0 - 1e-10:
                    pval = 0.0
                else:
                    t_stat = r * np.sqrt((n - 2) / (1 - r ** 2))
                    pval = 2.0 * t_dist.sf(abs(t_stat), df=n - 2)
                p_matrix[i, j] = pval
                p_matrix[j, i] = pval

        return corr_matrix, p_matrix, trait_names, n
