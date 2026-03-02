import json
import sys

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.phylogenetic_regression import PhylogeneticRegression
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
        response="brain_size",
        predictors=["body_mass"],
        method="BM",
        json=False,
    )


@pytest.fixture
def lambda_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        response="brain_size",
        predictors=["body_mass"],
        method="lambda",
        json=False,
    )


@pytest.fixture
def multi_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        response="brain_size",
        predictors=["body_mass", "longevity"],
        method="BM",
        json=False,
    )


@pytest.fixture
def json_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        response="brain_size",
        predictors=["body_mass"],
        method="BM",
        json=True,
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            response="y",
            predictors=["x1"],
        )
        svc = PhylogeneticRegression.__new__(PhylogeneticRegression)
        parsed = svc.process_args(args)
        assert parsed["tree_file_path"] == "t.tre"
        assert parsed["trait_data_path"] == "d.tsv"
        assert parsed["response"] == "y"
        assert parsed["predictors"] == ["x1"]
        assert parsed["method"] == "BM"
        assert parsed["json_output"] is False

    def test_overrides(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            response="y",
            predictors=["x1", "x2"],
            method="lambda",
            json=True,
        )
        svc = PhylogeneticRegression.__new__(PhylogeneticRegression)
        parsed = svc.process_args(args)
        assert parsed["method"] == "lambda"
        assert parsed["json_output"] is True
        assert parsed["predictors"] == ["x1", "x2"]


class TestValidation:
    def test_missing_response_column(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="nonexistent",
            predictors=["body_mass"],
            method="BM",
            json=False,
        )
        with pytest.raises(SystemExit):
            svc = PhylogeneticRegression(args)
            svc.run()

    def test_missing_predictor_column(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["nonexistent"],
            method="BM",
            json=False,
        )
        with pytest.raises(SystemExit):
            svc = PhylogeneticRegression(args)
            svc.run()

    def test_response_in_predictors(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["brain_size"],
            method="BM",
            json=False,
        )
        with pytest.raises(SystemExit):
            svc = PhylogeneticRegression(args)
            svc.run()


class TestPGLSBM:
    """Test BM regression against R reference values.

    Reference values computed in R 4.4.0 using manual GLS with ape::vcv()
    (raw VCV, matching caper::pgls behavior):

        library(ape)
        tree <- read.tree("tree_simple.tre")
        traits <- read.delim("tree_simple_multi_traits.tsv", row.names=1)
        ord <- sort(rownames(traits)); traits_s <- traits[ord, ]
        C <- vcv(tree)[ord, ord]; C_inv <- solve(C)
        y <- traits_s$brain_size; X <- cbind(1, traits_s$body_mass)
        beta <- solve(t(X)%*%C_inv%*%X) %*% t(X)%*%C_inv%*%y
        # Intercept=0.9972138  body_mass=0.7085628
    """

    def test_coefficients_match_r(self, default_args):
        """Coefficients must match R manual GLS with raw VCV."""
        svc = PhylogeneticRegression(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(svc.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("brain_size")
        pred_idx = trait_names.index("body_mass")

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)

        beta_hat, residuals, sigma2, var_beta = svc._fit_gls(y, X, C_inv)
        se = np.sqrt(np.diag(var_beta))

        # R reference: Intercept=0.9972138 SE=0.0871177 t=11.4467407
        assert beta_hat[0] == pytest.approx(0.9972139, abs=1e-4)
        assert se[0] == pytest.approx(0.0871177, abs=1e-4)

        # R reference: body_mass=0.7085628 SE=0.0450911 t=15.7140257
        assert beta_hat[1] == pytest.approx(0.7085628, abs=1e-4)
        assert se[1] == pytest.approx(0.0450911, abs=1e-4)

    def test_model_stats_match_r(self, default_args):
        """R², F-stat, log-lik, AIC must match R manual GLS."""
        svc = PhylogeneticRegression(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(svc.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("brain_size")
        pred_idx = trait_names.index("body_mass")

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)

        beta_hat, residuals, sigma2, var_beta = svc._fit_gls(y, X, C_inv)
        fitted = X @ beta_hat

        r2, adj_r2, f_stat, f_p = svc._compute_model_stats(
            y, fitted, residuals, C_inv, k=1, n=n
        )

        # R reference
        assert r2 == pytest.approx(0.9762781, abs=1e-4)
        assert adj_r2 == pytest.approx(0.9723244, abs=1e-4)
        assert f_stat == pytest.approx(246.9306, abs=0.01)
        assert f_p == pytest.approx(4.2092e-06, rel=0.01)

        # Residual SE (REML): 0.0249942
        rse = np.sqrt(sigma2)
        assert rse == pytest.approx(0.0249942, abs=1e-4)

    def test_log_likelihood_and_aic_match_r(self, default_args):
        """ML log-likelihood and AIC must match R."""
        svc = PhylogeneticRegression(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(svc.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("brain_size")
        pred_idx = trait_names.index("body_mass")

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)
        beta_hat, residuals, _, _ = svc._fit_gls(y, X, C_inv)

        # Compute ML log-likelihood (same formula as in run())
        sign, logdet_C = np.linalg.slogdet(vcv)
        sigma2_ml = float(residuals @ C_inv @ residuals) / n
        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sigma2_ml) + logdet_C + n)
        aic = -2.0 * ll + 2.0 * (1 + 2)

        # R reference: ML log-likelihood=6.0558469, AIC=-6.1116939
        assert ll == pytest.approx(6.0558, abs=0.001)
        assert aic == pytest.approx(-6.1117, abs=0.001)

    def test_t_stats_and_pvalues(self, default_args):
        svc = PhylogeneticRegression(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(svc.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("brain_size")
        pred_idx = trait_names.index("body_mass")

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)

        beta_hat, residuals, sigma2, var_beta = svc._fit_gls(y, X, C_inv)
        se = np.sqrt(np.diag(var_beta))
        t_stats = beta_hat / se

        from scipy.stats import t as t_distribution
        df_resid = n - 2
        p_values = 2.0 * t_distribution.sf(np.abs(t_stats), df=df_resid)

        # R reference: t-values 11.4467407, 15.7140257
        assert t_stats[0] == pytest.approx(11.4467, abs=0.01)
        assert t_stats[1] == pytest.approx(15.7140, abs=0.01)

        # R reference: p-values 0.0000266772, 0.0000042092
        assert p_values[0] == pytest.approx(2.6677e-05, rel=0.01)
        assert p_values[1] == pytest.approx(4.2092e-06, rel=0.01)


class TestPGLSLambda:
    def test_lambda_estimation(self, lambda_args):
        svc = PhylogeneticRegression(lambda_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(svc.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("brain_size")
        pred_idx = trait_names.index("body_mass")

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        max_lam = svc._max_lambda(tree)

        lambda_val, ll = svc._estimate_lambda(y, X, vcv, max_lam)

        # Lambda should be between 0 and max_lambda
        assert 0 <= lambda_val <= max_lam
        # Log-likelihood should be a finite number
        assert np.isfinite(ll)

    def test_lambda_coefficients(self, lambda_args):
        svc = PhylogeneticRegression(lambda_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(svc.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("brain_size")
        pred_idx = trait_names.index("body_mass")

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        max_lam = svc._max_lambda(tree)

        lambda_val, _ = svc._estimate_lambda(y, X, vcv, max_lam)

        # Transform VCV
        diag_vals = np.diag(vcv).copy()
        vcv_t = vcv * lambda_val
        np.fill_diagonal(vcv_t, diag_vals)

        C_inv = np.linalg.inv(vcv_t)
        beta_hat, _, _, _ = svc._fit_gls(y, X, C_inv)

        # Coefficients should still be reasonable
        assert 0.5 < beta_hat[0] < 1.5
        assert 0.4 < beta_hat[1] < 1.0


class TestMultipleRegression:
    """Test multiple regression against R reference values.

    R reference (manual GLS with raw VCV):
        Intercept=0.6544307 SE=0.3311317 t=1.9763454 p=0.1050679
        body_mass=0.6193379 SE=0.0943992 t=6.5608407 p=0.0012333
        longevity=0.1628710 SE=0.1519291 t=1.0720198 p=0.3327044
        Residual SE (REML)=0.0246890
        ML log-likelihood=6.8834005  AIC=-5.7668009
    """

    def test_multiple_predictors_match_r(self, multi_args):
        svc = PhylogeneticRegression(multi_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(svc.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("brain_size")
        pred_indices = [trait_names.index("body_mass"), trait_names.index("longevity")]

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([[traits[name][j] for j in pred_indices] for name in ordered_names]),
        ])
        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C_inv = np.linalg.inv(vcv)

        beta_hat, residuals, sigma2, var_beta = svc._fit_gls(y, X, C_inv)
        se = np.sqrt(np.diag(var_beta))

        # Should have 3 coefficients: intercept + 2 predictors
        assert len(beta_hat) == 3
        assert var_beta.shape == (3, 3)

        # R reference values
        assert beta_hat[0] == pytest.approx(0.6544307, abs=1e-4)
        assert beta_hat[1] == pytest.approx(0.6193379, abs=1e-4)
        assert beta_hat[2] == pytest.approx(0.1628710, abs=1e-4)

        assert se[0] == pytest.approx(0.3311317, abs=1e-4)
        assert se[1] == pytest.approx(0.0943992, abs=1e-4)
        assert se[2] == pytest.approx(0.1519291, abs=1e-4)

        rse = np.sqrt(sigma2)
        assert rse == pytest.approx(0.024689, abs=1e-4)


class TestRun:
    def test_text_output(self, default_args, capsys):
        svc = PhylogeneticRegression(default_args)
        svc.run()
        captured = capsys.readouterr()
        assert "Phylogenetic Generalized Least Squares (PGLS)" in captured.out
        assert "Formula: brain_size ~ body_mass" in captured.out
        assert "(Intercept)" in captured.out
        assert "body_mass" in captured.out
        assert "R-squared" in captured.out
        assert "F-statistic" in captured.out
        assert "Log-likelihood" in captured.out
        assert "AIC" in captured.out

    def test_json_output(self, json_args, capsys):
        svc = PhylogeneticRegression(json_args)
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert "formula" in payload
        assert "coefficients" in payload
        assert "(Intercept)" in payload["coefficients"]
        assert "body_mass" in payload["coefficients"]
        assert "r_squared" in payload
        assert "adj_r_squared" in payload
        assert "f_statistic" in payload
        assert "log_likelihood" in payload
        assert "aic" in payload
        assert "residuals" in payload
        assert "fitted_values" in payload

    def test_json_lambda_output(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="lambda",
            json=True,
        )
        svc = PhylogeneticRegression(args)
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert "lambda" in payload

    def test_text_lambda_output(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="lambda",
            json=False,
        )
        svc = PhylogeneticRegression(args)
        svc.run()
        captured = capsys.readouterr()
        assert "Estimated lambda" in captured.out

    def test_signif_codes_displayed(self, default_args, capsys):
        svc = PhylogeneticRegression(default_args)
        svc.run()
        captured = capsys.readouterr()
        assert "Signif. codes:" in captured.out
        assert "***" in captured.out

    def test_multiple_predictors_output(self, multi_args, capsys):
        svc = PhylogeneticRegression(multi_args)
        svc.run()
        captured = capsys.readouterr()
        assert "body_mass" in captured.out
        assert "longevity" in captured.out
        assert "brain_size ~ body_mass + longevity" in captured.out


GENE_TREES_FILE = str(SAMPLE_FILES / "gene_trees_simple.nwk")


class TestDiscordanceVCV:
    def test_run_with_gene_trees_json(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="BM",
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticRegression(args)
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert "vcv_metadata" in payload
        assert payload["vcv_metadata"]["n_gene_trees"] == 10
        assert "coefficients" in payload

    def test_discordance_changes_coefficients(self, capsys):
        """With discordant gene trees, PGLS coefficients should differ."""
        args_no_gt = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="BM",
            json=True,
        )
        svc = PhylogeneticRegression(args_no_gt)
        svc.run()
        out, _ = capsys.readouterr()
        coefs_no_gt = json.loads(out)["coefficients"]

        args_gt = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="BM",
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticRegression(args_gt)
        svc.run()
        out, _ = capsys.readouterr()
        coefs_gt = json.loads(out)["coefficients"]

        # Coefficients should differ (not exactly equal)
        assert coefs_gt["body_mass"]["estimate"] != pytest.approx(
            coefs_no_gt["body_mass"]["estimate"], abs=1e-6
        )

    def test_text_output_with_gene_trees(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            response="brain_size",
            predictors=["body_mass"],
            method="BM",
            json=False,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticRegression(args)
        svc.run()
        captured = capsys.readouterr()
        assert "PGLS" in captured.out
        assert "body_mass" in captured.out
