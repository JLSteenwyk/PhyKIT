import os

import pytest
import json
import numpy as np
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.phylogenetic_pca import PhylogeneticPCA
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
        method="BM",
        mode="cov",
        json=False,
    )


@pytest.fixture
def corr_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        method="BM",
        mode="corr",
        json=False,
    )


@pytest.fixture
def lambda_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        method="lambda",
        mode="cov",
        json=False,
    )


@pytest.fixture
def lambda_corr_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        method="lambda",
        mode="corr",
        json=False,
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(tree="t.tre", trait_data="d.tsv")
        svc = PhylogeneticPCA(args)
        assert svc.tree_file_path == "t.tre"
        assert svc.trait_data_path == "d.tsv"
        assert svc.method == "BM"
        assert svc.mode == "cov"
        assert svc.json_output is False

    def test_overrides(self):
        args = Namespace(
            tree="t.tre", trait_data="d.tsv",
            json=True, method="lambda", mode="corr",
        )
        svc = PhylogeneticPCA(args)
        assert svc.json_output is True
        assert svc.method == "lambda"
        assert svc.mode == "corr"


class TestParseMultiTraitFile:
    def test_valid_file(self, default_args):
        svc = PhylogeneticPCA(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        assert len(traits) == 8
        assert trait_names == ["body_mass", "brain_size", "longevity"]
        assert traits["raccoon"] == pytest.approx([1.04, 1.60, 2.71])
        assert traits["weasel"] == pytest.approx([-0.30, 0.85, 1.79])

    def test_missing_file(self, default_args):
        svc = PhylogeneticPCA(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_multi_trait_file("/nonexistent/path.tsv", ["a", "b", "c"])

    def test_non_numeric_value(self, default_args, tmp_path):
        trait_file = tmp_path / "bad.tsv"
        trait_file.write_text("taxon\ttrait1\ntaxon1\tabc\n")
        svc = PhylogeneticPCA(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_multi_trait_file(str(trait_file), ["taxon1"])

    def test_wrong_column_count(self, default_args, tmp_path):
        trait_file = tmp_path / "bad.tsv"
        trait_file.write_text("taxon\ttrait1\ttrait2\ntaxon1\t1.0\n")
        svc = PhylogeneticPCA(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_multi_trait_file(str(trait_file), ["taxon1"])

    def test_header_only(self, default_args, tmp_path):
        trait_file = tmp_path / "header_only.tsv"
        trait_file.write_text("taxon\ttrait1\ttrait2\n")
        svc = PhylogeneticPCA(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_multi_trait_file(str(trait_file), ["taxon1", "taxon2", "taxon3"])

    def test_comments_and_blanks(self, default_args, tmp_path):
        trait_file = tmp_path / "good.tsv"
        trait_file.write_text(
            "# comment\n\ntaxon\tt1\tt2\n"
            "taxon1\t1.0\t2.0\ntaxon2\t3.0\t4.0\ntaxon3\t5.0\t6.0\n"
        )
        svc = PhylogeneticPCA(default_args)
        trait_names, traits = svc._parse_multi_trait_file(
            str(trait_file), ["taxon1", "taxon2", "taxon3"]
        )
        assert len(traits) == 3
        assert trait_names == ["t1", "t2"]

    def test_taxon_mismatch_warns(self, default_args, tmp_path, capsys):
        trait_file = tmp_path / "partial.tsv"
        trait_file.write_text(
            "taxon\tt1\ntaxon1\t1.0\ntaxon2\t2.0\ntaxon3\t3.0\nextra\t4.0\n"
        )
        svc = PhylogeneticPCA(default_args)
        trait_names, traits = svc._parse_multi_trait_file(
            str(trait_file), ["taxon1", "taxon2", "taxon3", "taxon4"]
        )
        _, err = capsys.readouterr()
        assert "Warning" in err
        assert len(traits) == 3

    def test_too_few_shared_taxa(self, default_args, tmp_path):
        trait_file = tmp_path / "few.tsv"
        trait_file.write_text("taxon\tt1\ntaxon1\t1.0\ntaxon2\t2.0\n")
        svc = PhylogeneticPCA(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_multi_trait_file(str(trait_file), ["taxon1", "taxon2"])


class TestBuildVCVMatrix:
    def test_symmetric(self, default_args):
        svc = PhylogeneticPCA(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        vcv = svc._build_vcv_matrix(tree, tips)
        np.testing.assert_array_almost_equal(vcv, vcv.T)

    def test_diagonal_is_root_to_tip(self, default_args):
        svc = PhylogeneticPCA(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        vcv = svc._build_vcv_matrix(tree, tips)
        for i, name in enumerate(tips):
            expected = tree.distance(tree.root, name)
            assert vcv[i, i] == pytest.approx(expected, rel=1e-6)

    def test_correct_shape(self, default_args):
        svc = PhylogeneticPCA(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        vcv = svc._build_vcv_matrix(tree, tips)
        assert vcv.shape == (8, 8)


class TestPhylogeneticPCACov:
    """Test BM + cov mode against R phytools::phyl.pca reference values."""

    # R reference values from phyl.pca(tree, traits, method="BM", mode="cov")
    REF_EIGENVALUES = [0.0801798810, 0.0029237493, 0.0003075388]
    REF_LOADINGS = np.array([
        [0.7347392318, -0.4226132689, 0.5306187767],
        [0.5274002700, -0.1360637256, -0.8386510703],
        [0.4266230379, 0.8960383293, 0.1229149949],
    ])
    REF_SCORES = {
        "bear":     [ 0.9781912622, -0.02928980193,  0.02678735213],
        "cat":      [-1.4420092510,  0.17646695023,  0.09307185066],
        "dog":      [-0.6517229303, -0.09142653369, -0.04614241827],
        "monkey":   [-0.6404655004,  1.09725068252, -0.20806754274],
        "raccoon":  [-0.8671208830,  0.06719921291,  0.11461079162],
        "sea_lion": [ 0.8603996846, -0.19926813823, -0.11510172228],
        "seal":     [ 0.5225838100,  0.27753903060, -0.05910744959],
        "weasel":   [-2.6397148510, -0.08880647561, -0.08051186175],
    }

    def test_eigenvalues(self, default_args):
        svc = PhylogeneticPCA(default_args)
        tree = svc.read_tree_file()
        svc._validate_tree(tree)
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
        R = (Z.T @ C_inv @ Z) / (n - 1)
        eigenvalues, _ = np.linalg.eigh(R)
        eigenvalues = np.sort(eigenvalues)[::-1]

        np.testing.assert_allclose(eigenvalues, self.REF_EIGENVALUES, atol=1e-4)

    def test_loadings(self, default_args):
        svc = PhylogeneticPCA(default_args)
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
        R = (Z.T @ C_inv @ Z) / (n - 1)
        eigenvalues, eigenvectors = np.linalg.eigh(R)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]

        # Compare abs values (eigenvectors defined up to sign)
        np.testing.assert_allclose(
            np.abs(eigenvectors), np.abs(self.REF_LOADINGS), atol=1e-4
        )

    def test_scores(self, default_args):
        svc = PhylogeneticPCA(default_args)
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
        R = (Z.T @ C_inv @ Z) / (n - 1)
        eigenvalues, eigenvectors = np.linalg.eigh(R)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]
        scores = Z @ eigenvectors

        # Compare each column up to global sign flip
        for k, name in enumerate(ordered_names):
            ref = np.array(self.REF_SCORES[name])
            computed = scores[k]
            # Either signs match or all signs are flipped per PC
            for pc in range(p):
                assert abs(abs(computed[pc]) - abs(ref[pc])) < 1e-4, (
                    f"Score mismatch for {name}, PC{pc+1}: "
                    f"computed={computed[pc]}, ref={ref[pc]}"
                )


class TestPhylogeneticPCACorr:
    """Test BM + corr mode against R phytools::phyl.pca reference values."""

    REF_EIGENVALUES = [2.8490059295, 0.1402470393, 0.0107470312]
    REF_LOADINGS = np.array([
        [-0.5820001250, -0.4638777874,  0.66790212815],
        [-0.5861515834, -0.3299925945, -0.73995351804],
        [-0.5636507569,  0.8221449300,  0.07984696836],
    ])
    REF_SCORES = {
        "bear":     [-5.718644630, -0.2520427517,  0.1680858429],
        "cat":      [ 8.308065874,  1.4817516853,  0.5169817433],
        "dog":      [ 4.015546433, -0.6377122605, -0.2647078216],
        "monkey":   [ 1.158323429,  7.3219170395, -1.3711746799],
        "raccoon":  [ 5.154737886,  0.7073985615,  0.6616876057],
        "sea_lion": [-4.828148632, -1.6179405137, -0.6477830531],
        "seal":     [-3.785695422,  1.7805146765, -0.3805726389],
        "weasel":   [15.786800820, -0.4888037720, -0.4806870806],
    }

    def test_eigenvalues(self, corr_args):
        svc = PhylogeneticPCA(corr_args)
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
        R = (Z.T @ C_inv @ Z) / (n - 1)

        D = np.sqrt(np.diag(R))
        D_inv = np.diag(1.0 / D)
        R_corr = D_inv @ R @ D_inv
        eigenvalues, _ = np.linalg.eigh(R_corr)
        eigenvalues = np.sort(eigenvalues)[::-1]

        np.testing.assert_allclose(eigenvalues, self.REF_EIGENVALUES, atol=1e-4)

    def test_loadings(self, corr_args):
        svc = PhylogeneticPCA(corr_args)
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
        R = (Z.T @ C_inv @ Z) / (n - 1)

        D = np.sqrt(np.diag(R))
        D_inv = np.diag(1.0 / D)
        R_corr = D_inv @ R @ D_inv
        eigenvalues, eigenvectors = np.linalg.eigh(R_corr)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]

        np.testing.assert_allclose(
            np.abs(eigenvectors), np.abs(self.REF_LOADINGS), atol=1e-4
        )

    def test_scores(self, corr_args):
        svc = PhylogeneticPCA(corr_args)
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
        R = (Z.T @ C_inv @ Z) / (n - 1)

        D = np.sqrt(np.diag(R))
        D_inv = np.diag(1.0 / D)
        R_corr = D_inv @ R @ D_inv
        eigenvalues, eigenvectors = np.linalg.eigh(R_corr)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]

        Z_std = Z @ D_inv
        scores = Z_std @ eigenvectors

        for k, name in enumerate(ordered_names):
            ref = np.array(self.REF_SCORES[name])
            computed = scores[k]
            for pc in range(p):
                assert abs(abs(computed[pc]) - abs(ref[pc])) < 1e-4, (
                    f"Score mismatch for {name}, PC{pc+1}: "
                    f"computed={computed[pc]}, ref={ref[pc]}"
                )


class TestPhylogeneticPCALambda:
    """Test lambda + cov mode against R phytools::phyl.pca reference values."""

    REF_LAMBDA = 6.610696e-05
    REF_LOG_LIKELIHOOD = -6.016939
    REF_EIGENVALUES = [0.076301214816, 0.002241847818, 0.000315397635]

    def test_lambda_value(self, lambda_args):
        svc = PhylogeneticPCA(lambda_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)

        ordered_names = sorted(traits.keys())
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        max_lam = svc._max_lambda(tree)

        lambda_val, log_likelihood = svc._multi_trait_lambda(Y, vcv, max_lam)

        assert lambda_val == pytest.approx(self.REF_LAMBDA, abs=1e-3)
        # Log-likelihood may differ slightly from R due to optimizer internals;
        # both converge to lambda ≈ 0 and the LL surface is flat there.
        assert log_likelihood == pytest.approx(self.REF_LOG_LIKELIHOOD, abs=0.2)

    def test_eigenvalues(self, lambda_args):
        svc = PhylogeneticPCA(lambda_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        max_lam = svc._max_lambda(tree)

        lambda_val, _ = svc._multi_trait_lambda(Y, vcv, max_lam)

        # Transform VCV
        diag_vals = np.diag(vcv).copy()
        vcv_t = vcv * lambda_val
        np.fill_diagonal(vcv_t, diag_vals)

        C_inv = np.linalg.inv(vcv_t)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Z = Y - np.outer(ones, a_hat)
        R = (Z.T @ C_inv @ Z) / (n - 1)

        eigenvalues, _ = np.linalg.eigh(R)
        eigenvalues = np.sort(eigenvalues)[::-1]

        np.testing.assert_allclose(eigenvalues, self.REF_EIGENVALUES, atol=1e-4)


class TestRun:
    def test_bm_cov_text_output(self, default_args, capsys):
        svc = PhylogeneticPCA(default_args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Eigenvalues:" in out
        assert "Loadings:" in out
        assert "Scores:" in out
        assert "PC1" in out
        assert "body_mass" in out
        assert "raccoon" in out

    def test_bm_corr_text_output(self, corr_args, capsys):
        svc = PhylogeneticPCA(corr_args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Eigenvalues:" in out
        assert "Loadings:" in out
        assert "Scores:" in out

    def test_lambda_text_output(self, lambda_args, capsys):
        svc = PhylogeneticPCA(lambda_args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Lambda:" in out
        assert "Log-likelihood:" in out

    def test_json_output_bm_cov(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="BM",
            mode="cov",
            json=True,
        )
        svc = PhylogeneticPCA(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "eigenvalues" in data
        assert "proportion_of_variance" in data
        assert "loadings" in data
        assert "scores" in data
        assert "lambda" not in data

    def test_json_output_lambda(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="lambda",
            mode="cov",
            json=True,
        )
        svc = PhylogeneticPCA(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "eigenvalues" in data
        assert "lambda" in data
        assert "log_likelihood" in data

    def test_run_eigenvalues_match_reference(self, capsys):
        """End-to-end: run() eigenvalues match R reference."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="BM",
            mode="cov",
            json=True,
        )
        svc = PhylogeneticPCA(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        ev = data["eigenvalues"]
        assert ev["PC1"] == pytest.approx(0.0801798810, abs=1e-4)
        assert ev["PC2"] == pytest.approx(0.0029237493, abs=1e-4)
        assert ev["PC3"] == pytest.approx(0.0003075388, abs=1e-4)

    def test_plot_creates_file(self, tmp_path, capsys):
        pytest.importorskip("matplotlib")
        plot_path = str(tmp_path / "test_pca.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="BM",
            mode="cov",
            json=False,
            plot=True,
            plot_output=plot_path,
            plot_tree=False,
            color_by=None,
        )
        svc = PhylogeneticPCA(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0
        out, _ = capsys.readouterr()
        assert "Saved PCA plot:" in out

    def test_plot_json_includes_plot_output(self, tmp_path, capsys):
        pytest.importorskip("matplotlib")
        plot_path = str(tmp_path / "test_pca.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="BM",
            mode="cov",
            json=True,
            plot=True,
            plot_output=plot_path,
            plot_tree=False,
            color_by=None,
        )
        svc = PhylogeneticPCA(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "plot_output" in data
        assert data["plot_output"] == plot_path

    def test_no_plot_by_default(self, tmp_path, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="BM",
            mode="cov",
            json=False,
        )
        svc = PhylogeneticPCA(args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Saved PCA plot:" not in out

    def test_plot_tree_creates_file(self, tmp_path, capsys):
        pytest.importorskip("matplotlib")
        plot_path = str(tmp_path / "test_tree.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="BM",
            mode="cov",
            json=False,
            plot=True,
            plot_output=plot_path,
            plot_tree=True,
            color_by=None,
        )
        svc = PhylogeneticPCA(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0

    def test_plot_color_by_column(self, tmp_path, capsys):
        pytest.importorskip("matplotlib")
        plot_path = str(tmp_path / "test_color.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="BM",
            mode="cov",
            json=False,
            plot=True,
            plot_output=plot_path,
            plot_tree=False,
            color_by="body_mass",
        )
        svc = PhylogeneticPCA(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0

    def test_plot_tree_and_color_by_together(self, tmp_path, capsys):
        pytest.importorskip("matplotlib")
        plot_path = str(tmp_path / "test_both.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="BM",
            mode="cov",
            json=False,
            plot=True,
            plot_output=plot_path,
            plot_tree=True,
            color_by="body_mass",
        )
        svc = PhylogeneticPCA(args)
        svc.run()
        assert os.path.exists(plot_path)
        assert os.path.getsize(plot_path) > 0


class TestReconstructAncestralScores:
    def test_all_nodes_get_scores(self, default_args):
        svc = PhylogeneticPCA(default_args)
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
        R = (Z.T @ C_inv @ Z) / (n - 1)
        eigenvalues, eigenvectors = np.linalg.eigh(R)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]
        scores = Z @ eigenvectors

        node_estimates, node_distances, tree_pruned = \
            svc._reconstruct_ancestral_scores(tree, scores, ordered_names)

        # All clades should have estimates
        all_clades = list(tree_pruned.find_clades())
        for clade in all_clades:
            assert id(clade) in node_estimates

    def test_tips_match_original_scores(self, default_args):
        svc = PhylogeneticPCA(default_args)
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
        R = (Z.T @ C_inv @ Z) / (n - 1)
        eigenvalues, eigenvectors = np.linalg.eigh(R)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]
        scores = Z @ eigenvectors

        node_estimates, _, tree_pruned = \
            svc._reconstruct_ancestral_scores(tree, scores, ordered_names)

        name_to_idx = {name: i for i, name in enumerate(ordered_names)}
        for tip in tree_pruned.get_terminals():
            tip_scores = node_estimates[id(tip)]
            expected = scores[name_to_idx[tip.name]]
            np.testing.assert_array_almost_equal(tip_scores, expected)

    def test_root_within_tip_range(self, default_args):
        svc = PhylogeneticPCA(default_args)
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
        R = (Z.T @ C_inv @ Z) / (n - 1)
        eigenvalues, eigenvectors = np.linalg.eigh(R)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]
        scores = Z @ eigenvectors

        node_estimates, _, tree_pruned = \
            svc._reconstruct_ancestral_scores(tree, scores, ordered_names)

        root_scores = node_estimates[id(tree_pruned.root)]
        # Root PC1 should be within range of tip PC1 values
        tip_pc1 = scores[:, 0]
        assert root_scores[0] >= min(tip_pc1) - 0.5
        assert root_scores[0] <= max(tip_pc1) + 0.5


class TestParseColorBy:
    def test_column_name_match(self, default_args):
        svc = PhylogeneticPCA(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])

        values, categories, kind = svc._parse_color_by(
            "body_mass", trait_names, Y, ordered_names
        )
        assert kind == "continuous"
        assert len(values) == len(ordered_names)
        assert len(categories) == 0

    def test_continuous_file(self, default_args, tmp_path):
        svc = PhylogeneticPCA(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])

        color_file = tmp_path / "colors.tsv"
        lines = [f"{name}\t{i * 1.5}\n" for i, name in enumerate(ordered_names)]
        color_file.write_text("".join(lines))

        values, categories, kind = svc._parse_color_by(
            str(color_file), trait_names, Y, ordered_names
        )
        assert kind == "continuous"
        assert len(values) == len(ordered_names)

    def test_discrete_file(self, default_args, tmp_path):
        svc = PhylogeneticPCA(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])

        color_file = tmp_path / "groups.tsv"
        lines = [f"{name}\t{'A' if i % 2 == 0 else 'B'}\n"
                 for i, name in enumerate(ordered_names)]
        color_file.write_text("".join(lines))

        values, categories, kind = svc._parse_color_by(
            str(color_file), trait_names, Y, ordered_names
        )
        assert kind == "discrete"
        assert set(categories) == {"A", "B"}
        assert len(values) == len(ordered_names)

    def test_invalid_input(self, default_args):
        svc = PhylogeneticPCA(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered_names = sorted(traits.keys())
        p = len(trait_names)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])

        with pytest.raises(PhykitUserError):
            svc._parse_color_by(
                "nonexistent_column_and_not_a_file", trait_names, Y, ordered_names
            )
