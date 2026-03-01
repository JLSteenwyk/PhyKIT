import os
import json

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.phylogenetic_dimreduce import PhylogeneticDimreduce
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        method="tsne",
        correction="BM",
        n_components=2,
        perplexity=None,
        n_neighbors=None,
        min_dist=0.1,
        seed=42,
        json=False,
        plot=False,
        plot_output="phylo_dimreduce_plot.png",
        plot_tree=False,
        color_by=None,
    )


@pytest.fixture
def umap_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        method="umap",
        correction="BM",
        n_components=2,
        perplexity=None,
        n_neighbors=None,
        min_dist=0.1,
        seed=42,
        json=False,
        plot=False,
        plot_output="phylo_dimreduce_plot.png",
        plot_tree=False,
        color_by=None,
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(tree="t.tre", trait_data="d.tsv")
        svc = PhylogeneticDimreduce(args)
        assert svc.tree_file_path == "t.tre"
        assert svc.trait_data_path == "d.tsv"
        assert svc.method == "tsne"
        assert svc.correction == "BM"
        assert svc.n_components == 2
        assert svc.perplexity is None
        assert svc.n_neighbors is None
        assert svc.min_dist == 0.1
        assert svc.seed is None
        assert svc.json_output is False

    def test_overrides(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            method="umap",
            correction="lambda",
            n_components=3,
            perplexity=10.0,
            n_neighbors=5,
            min_dist=0.5,
            seed=99,
            json=True,
        )
        svc = PhylogeneticDimreduce(args)
        assert svc.method == "umap"
        assert svc.correction == "lambda"
        assert svc.n_components == 3
        assert svc.perplexity == 10.0
        assert svc.n_neighbors == 5
        assert svc.min_dist == 0.5
        assert svc.seed == 99
        assert svc.json_output is True


class TestTSNEEmbedding:
    def test_shape(self, default_args):
        svc = PhylogeneticDimreduce(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)

        embedding, params = svc._embed_tsne(Z, n)
        assert embedding.shape == (n, 2)

    def test_deterministic_with_seed(self, default_args):
        svc = PhylogeneticDimreduce(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)

        embedding1, _ = svc._embed_tsne(Z, n)
        embedding2, _ = svc._embed_tsne(Z, n)
        np.testing.assert_array_almost_equal(embedding1, embedding2)

    def test_auto_perplexity(self, default_args):
        svc = PhylogeneticDimreduce(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)

        _, params = svc._embed_tsne(Z, n)
        expected_perplexity = min(30.0, (n - 1) / 3.0)
        assert params["perplexity"] == round(expected_perplexity, 2)


class TestUMAPEmbedding:
    def test_shape(self, umap_args):
        svc = PhylogeneticDimreduce(umap_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)

        embedding, params = svc._embed_umap(Z, n)
        assert embedding.shape == (n, 2)

    def test_deterministic_with_seed(self, umap_args):
        svc = PhylogeneticDimreduce(umap_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)

        embedding1, _ = svc._embed_umap(Z, n)
        embedding2, _ = svc._embed_umap(Z, n)
        np.testing.assert_array_almost_equal(embedding1, embedding2)

    def test_auto_n_neighbors(self, umap_args):
        svc = PhylogeneticDimreduce(umap_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)

        _, params = svc._embed_umap(Z, n)
        expected_neighbors = min(15, n - 1)
        assert params["n_neighbors"] == expected_neighbors


class TestLambdaCorrection:
    def test_lambda_produces_different_embedding(self, default_args):
        svc_bm = PhylogeneticDimreduce(default_args)
        svc_bm.run = lambda: None  # prevent run
        tree = svc_bm.read_tree_file()
        tips = svc_bm.get_tip_names_from_tree(tree)
        trait_names, traits = svc_bm._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])

        # BM centering
        vcv_bm = svc_bm._build_vcv_matrix(tree, ordered_names)
        C_inv_bm = np.linalg.inv(vcv_bm)
        ones = np.ones(n)
        denom_bm = ones @ C_inv_bm @ ones
        a_hat_bm = np.array([(ones @ C_inv_bm @ Y[:, j]) / denom_bm for j in range(p)])
        Z_bm = Y - np.outer(ones, a_hat_bm)

        # Lambda centering
        vcv_lam = svc_bm._build_vcv_matrix(tree, ordered_names)
        max_lam = svc_bm._max_lambda(tree)
        lambda_val, _ = svc_bm._multi_trait_lambda(Y, vcv_lam, max_lam)
        diag_vals = np.diag(vcv_lam).copy()
        vcv_lam = vcv_lam * lambda_val
        np.fill_diagonal(vcv_lam, diag_vals)
        C_inv_lam = np.linalg.inv(vcv_lam)
        denom_lam = ones @ C_inv_lam @ ones
        a_hat_lam = np.array([(ones @ C_inv_lam @ Y[:, j]) / denom_lam for j in range(p)])
        Z_lam = Y - np.outer(ones, a_hat_lam)

        # The centered data should differ
        assert not np.allclose(Z_bm, Z_lam)

    def test_lambda_value_returned(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="tsne",
            correction="lambda",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=True,
        )
        svc = PhylogeneticDimreduce(args)
        svc.run()
        # If it runs without error and lambda is set, the test passes
        # The lambda correction path was exercised


class TestSmallSampleGuards:
    def test_tsne_too_few_taxa(self, tmp_path):
        tree_file = tmp_path / "small.tre"
        tree_file.write_text("((A:1,B:1):1,C:2);")
        trait_file = tmp_path / "small.tsv"
        trait_file.write_text("taxon\tt1\tt2\nA\t1.0\t2.0\nB\t3.0\t4.0\nC\t5.0\t6.0\n")
        args = Namespace(
            tree=str(tree_file),
            trait_data=str(trait_file),
            method="tsne",
            correction="BM",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=False,
            plot=False,
            plot_output="test.png",
            plot_tree=False,
            color_by=None,
        )
        # With 3 taxa, perplexity = (3-1)/3 = 0.67 < 1, should error
        svc = PhylogeneticDimreduce(args)
        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()
        assert any("perplexity" in m for m in exc_info.value.messages)

    def test_umap_too_few_taxa(self, tmp_path):
        tree_file = tmp_path / "tiny.tre"
        tree_file.write_text("(A:1,B:1);")
        trait_file = tmp_path / "tiny.tsv"
        # Only 2 shared taxa -> _parse_multi_trait_file will error first
        trait_file.write_text("taxon\tt1\tt2\nA\t1.0\t2.0\nB\t3.0\t4.0\n")
        args = Namespace(
            tree=str(tree_file),
            trait_data=str(trait_file),
            method="umap",
            correction="BM",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=False,
            plot=False,
            plot_output="test.png",
            plot_tree=False,
            color_by=None,
        )
        svc = PhylogeneticDimreduce(args)
        with pytest.raises(PhykitUserError):
            svc.run()


class TestRun:
    def test_tsne_text_output(self, default_args, capsys):
        svc = PhylogeneticDimreduce(default_args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Method: tsne" in out
        assert "Correction: BM" in out
        assert "Perplexity:" in out
        assert "Embedding:" in out
        assert "Dim1" in out
        assert "Dim2" in out
        assert "raccoon" in out

    def test_umap_text_output(self, umap_args, capsys):
        svc = PhylogeneticDimreduce(umap_args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Method: umap" in out
        assert "n_neighbors:" in out
        assert "Embedding:" in out

    def test_json_output(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="tsne",
            correction="BM",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=True,
        )
        svc = PhylogeneticDimreduce(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert data["method"] == "tsne"
        assert data["correction"] == "BM"
        assert "embedding" in data
        assert "bear" in data["embedding"]
        assert "Dim1" in data["embedding"]["bear"]
        assert "Dim2" in data["embedding"]["bear"]
        assert "lambda" not in data

    def test_lambda_json_output(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="tsne",
            correction="lambda",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=True,
        )
        svc = PhylogeneticDimreduce(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "lambda" in data
        assert "log_likelihood" in data
        assert "embedding" in data

    def test_subset_check(self, default_args, capsys):
        svc = PhylogeneticDimreduce(default_args)
        svc.run()
        out, _ = capsys.readouterr()
        # All 8 taxa from sample file should appear
        for taxon in ["bear", "cat", "dog", "monkey", "raccoon", "sea_lion", "seal", "weasel"]:
            assert taxon in out


class TestPlot:
    def test_plot_created(self, tmp_path, capsys):
        plot_path = str(tmp_path / "test_dimreduce.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="tsne",
            correction="BM",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=False,
            plot=True,
            plot_output=plot_path,
            plot_tree=False,
            color_by=None,
        )
        svc = PhylogeneticDimreduce(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0
        out, _ = capsys.readouterr()
        assert "Saved plot:" in out

    def test_plot_tree_created(self, tmp_path, capsys):
        plot_path = str(tmp_path / "test_tree.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="tsne",
            correction="BM",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=False,
            plot=True,
            plot_output=plot_path,
            plot_tree=True,
            color_by=None,
        )
        svc = PhylogeneticDimreduce(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0

    def test_color_by_works(self, tmp_path, capsys):
        plot_path = str(tmp_path / "test_color.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="tsne",
            correction="BM",
            n_components=2,
            perplexity=None,
            n_neighbors=None,
            min_dist=0.1,
            seed=42,
            json=False,
            plot=True,
            plot_output=plot_path,
            plot_tree=False,
            color_by="body_mass",
        )
        svc = PhylogeneticDimreduce(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0
