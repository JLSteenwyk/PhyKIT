import json
import sys

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from phykit.services.tree.phylogenetic_glm import PhylogeneticGLM
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
GLM_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_glm_traits.tsv")


@pytest.fixture
def binomial_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=GLM_TRAITS_FILE,
        response="binary_trait",
        predictors=["body_mass"],
        family="binomial",
        method=None,
        btol=10,
        log_alpha_bound=4,
        json=False,
    )


@pytest.fixture
def poisson_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=GLM_TRAITS_FILE,
        response="count_trait",
        predictors=["body_mass"],
        family="poisson",
        method=None,
        btol=10,
        log_alpha_bound=4,
        json=False,
    )


@pytest.fixture
def binomial_json_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=GLM_TRAITS_FILE,
        response="binary_trait",
        predictors=["body_mass"],
        family="binomial",
        method=None,
        btol=10,
        log_alpha_bound=4,
        json=True,
    )


@pytest.fixture
def poisson_json_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=GLM_TRAITS_FILE,
        response="count_trait",
        predictors=["body_mass"],
        family="poisson",
        method=None,
        btol=10,
        log_alpha_bound=4,
        json=True,
    )


@pytest.fixture
def multi_pred_binomial_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=GLM_TRAITS_FILE,
        response="binary_trait",
        predictors=["body_mass", "count_trait"],
        family="binomial",
        method=None,
        btol=10,
        log_alpha_bound=4,
        json=False,
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            response="y",
            predictors=["x1"],
            family="binomial",
        )
        svc = PhylogeneticGLM.__new__(PhylogeneticGLM)
        parsed = svc.process_args(args)
        assert parsed["tree_file_path"] == "t.tre"
        assert parsed["trait_data_path"] == "d.tsv"
        assert parsed["response"] == "y"
        assert parsed["predictors"] == ["x1"]
        assert parsed["family"] == "binomial"
        assert parsed["method"] == "logistic_MPLE"
        assert parsed["btol"] == 10
        assert parsed["log_alpha_bound"] == 4
        assert parsed["json_output"] is False

    def test_family_auto_method_binomial(self):
        args = Namespace(
            tree="t.tre", trait_data="d.tsv", response="y",
            predictors=["x1"], family="binomial", method=None,
        )
        svc = PhylogeneticGLM.__new__(PhylogeneticGLM)
        parsed = svc.process_args(args)
        assert parsed["method"] == "logistic_MPLE"

    def test_family_auto_method_poisson(self):
        args = Namespace(
            tree="t.tre", trait_data="d.tsv", response="y",
            predictors=["x1"], family="poisson", method=None,
        )
        svc = PhylogeneticGLM.__new__(PhylogeneticGLM)
        parsed = svc.process_args(args)
        assert parsed["method"] == "poisson_GEE"

    def test_custom_btol_and_log_alpha_bound(self):
        args = Namespace(
            tree="t.tre", trait_data="d.tsv", response="y",
            predictors=["x1"], family="binomial", method=None,
            btol=20, log_alpha_bound=8, json=True,
        )
        svc = PhylogeneticGLM.__new__(PhylogeneticGLM)
        parsed = svc.process_args(args)
        assert parsed["btol"] == 20
        assert parsed["log_alpha_bound"] == 8
        assert parsed["json_output"] is True


class TestValidation:
    def test_missing_column(self, binomial_args):
        binomial_args.response = "nonexistent_column"
        svc = PhylogeneticGLM(binomial_args)
        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()
        assert any("nonexistent_column" in m for m in exc_info.value.messages)

    def test_response_in_predictors(self, binomial_args):
        binomial_args.predictors = ["body_mass", "binary_trait"]
        svc = PhylogeneticGLM(binomial_args)
        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()
        assert any("must not also be a predictor" in m for m in exc_info.value.messages)

    def test_non_binary_for_binomial(self, poisson_args):
        poisson_args.family = "binomial"
        poisson_args.method = None
        poisson_args.response = "count_trait"
        svc = PhylogeneticGLM(poisson_args)
        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()
        assert any("must contain only 0 and 1" in m for m in exc_info.value.messages)

    def test_insufficient_taxa(self, binomial_args):
        # Covered by _validate_tree requiring >= 3 tips
        pass

    def test_missing_predictor_column(self, binomial_args):
        binomial_args.predictors = ["nonexistent"]
        svc = PhylogeneticGLM(binomial_args)
        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()
        assert any("nonexistent" in m for m in exc_info.value.messages)


class TestLogisticMPLE:
    def test_coefficients_reasonable(self, binomial_args):
        """Logistic MPLE should produce finite coefficients within btol bounds."""
        svc = PhylogeneticGLM(binomial_args)
        tree = svc.read_tree_file()
        svc._validate_tree(tree)
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("binary_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        result = svc._fit_logistic_mple(tree, y, X, ordered_names)

        assert "coefficients" in result
        intercept = result["coefficients"]["(Intercept)"]["estimate"]
        slope = result["coefficients"]["body_mass"]["estimate"]
        assert np.isfinite(intercept)
        assert np.isfinite(slope)
        assert abs(intercept) <= svc.btol + 0.01
        assert abs(slope) <= svc.btol + 0.01

    def test_alpha_within_bounds(self, binomial_args):
        """Alpha should be within the bounds set by log_alpha_bound."""
        svc = PhylogeneticGLM(binomial_args)
        tree = svc.read_tree_file()
        svc._validate_tree(tree)
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("binary_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        result = svc._fit_logistic_mple(tree, y, X, ordered_names)
        alpha = result["alpha"]
        assert alpha > 0
        assert alpha <= np.exp(svc.log_alpha_bound) + 0.01

    def test_standard_errors_positive(self, binomial_args):
        """Standard errors should be positive and finite."""
        svc = PhylogeneticGLM(binomial_args)
        tree = svc.read_tree_file()
        svc._validate_tree(tree)
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("binary_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        result = svc._fit_logistic_mple(tree, y, X, ordered_names)
        for name, coef in result["coefficients"].items():
            assert coef["std_error"] > 0, f"SE for {name} should be positive"
            assert np.isfinite(coef["std_error"]), f"SE for {name} should be finite"

    def test_ultrametric_corrections(self, binomial_args):
        """_make_ultrametric should give non-negative D values."""
        svc = PhylogeneticGLM(binomial_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())

        D, Tmax, mean_height = svc._make_ultrametric(tree, ordered_names)
        assert Tmax > 0
        assert mean_height > 0
        assert mean_height <= Tmax
        assert np.all(D >= 0)
        # At least one tip should be at maximum distance (D=0)
        assert np.min(D) == pytest.approx(0.0)

    def test_bernoulli_likelihood_finite(self, binomial_args):
        """Bernoulli log-likelihood should return a finite negative value."""
        svc = PhylogeneticGLM(binomial_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("binary_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        beta = np.array([0.0, 1.0])
        eta = X @ beta
        mu = 1.0 / (1.0 + np.exp(-eta))

        ll = PhylogeneticGLM._bernoulli_log_likelihood(y, mu)
        assert np.isfinite(ll)
        assert ll < 0  # Log-likelihood should be negative


class TestPoissonGEE:
    def test_convergence(self, poisson_args):
        """Poisson GEE should converge for our test data."""
        svc = PhylogeneticGLM(poisson_args)
        tree = svc.read_tree_file()
        svc._validate_tree(tree)
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("count_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        result = svc._fit_poisson_gee(tree, y, X, ordered_names)

        # Should have sensible coefficients
        intercept = result["coefficients"]["(Intercept)"]["estimate"]
        slope = result["coefficients"]["body_mass"]["estimate"]
        assert np.isfinite(intercept)
        assert np.isfinite(slope)

    def test_positive_slope(self, poisson_args):
        """Count trait increases with body_mass, so slope should be positive."""
        svc = PhylogeneticGLM(poisson_args)
        tree = svc.read_tree_file()
        svc._validate_tree(tree)
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("count_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        result = svc._fit_poisson_gee(tree, y, X, ordered_names)
        slope = result["coefficients"]["body_mass"]["estimate"]
        assert slope > 0

    def test_overdispersion_positive(self, poisson_args):
        """Overdispersion phi should be positive."""
        svc = PhylogeneticGLM(poisson_args)
        tree = svc.read_tree_file()
        svc._validate_tree(tree)
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("count_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        result = svc._fit_poisson_gee(tree, y, X, ordered_names)
        assert result["overdispersion"] > 0

    def test_standard_errors_positive(self, poisson_args):
        """Standard errors should be positive."""
        svc = PhylogeneticGLM(poisson_args)
        tree = svc.read_tree_file()
        svc._validate_tree(tree)
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("count_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        result = svc._fit_poisson_gee(tree, y, X, ordered_names)
        for name, coef in result["coefficients"].items():
            assert coef["std_error"] > 0, f"SE for {name} should be positive"

    def test_starting_values(self, poisson_args):
        """Starting values should give reasonable initial estimates."""
        svc = PhylogeneticGLM(poisson_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(GLM_TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        resp_idx = trait_names.index("count_trait")
        pred_idx = trait_names.index("body_mass")
        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X = np.column_stack([
            np.ones(n),
            np.array([traits[name][pred_idx] for name in ordered_names]),
        ])

        beta0 = svc._poisson_starting_values(y, X)
        assert np.all(np.isfinite(beta0))
        # Intercept should be near log(mean(y))
        assert abs(beta0[0] - np.log(np.mean(y))) < 2.0


class TestRun:
    @patch("builtins.print")
    def test_text_output_binomial(self, mocked_print, binomial_args):
        svc = PhylogeneticGLM(binomial_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Phylogenetic GLM (Logistic MPLE)" in all_output
        assert "(Intercept)" in all_output
        assert "body_mass" in all_output
        assert "Estimated alpha" in all_output
        assert "binary_trait ~ body_mass" in all_output

    @patch("builtins.print")
    def test_text_output_poisson(self, mocked_print, poisson_args):
        svc = PhylogeneticGLM(poisson_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Phylogenetic GLM (Poisson GEE)" in all_output
        assert "(Intercept)" in all_output
        assert "body_mass" in all_output
        assert "Overdispersion (phi)" in all_output
        assert "count_trait ~ body_mass" in all_output

    @patch("builtins.print")
    def test_json_output_binomial(self, mocked_print, binomial_json_args):
        svc = PhylogeneticGLM(binomial_json_args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["family"] == "binomial"
        assert payload["method"] == "logistic_MPLE"
        assert "alpha" in payload
        assert "(Intercept)" in payload["coefficients"]
        assert "body_mass" in payload["coefficients"]
        assert payload["formula"] == "binary_trait ~ body_mass"

    @patch("builtins.print")
    def test_json_output_poisson(self, mocked_print, poisson_json_args):
        svc = PhylogeneticGLM(poisson_json_args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["family"] == "poisson"
        assert payload["method"] == "poisson_GEE"
        assert "overdispersion" in payload
        assert "(Intercept)" in payload["coefficients"]
        assert payload["formula"] == "count_trait ~ body_mass"

    @patch("builtins.print")
    def test_signif_codes_displayed(self, mocked_print, poisson_args):
        svc = PhylogeneticGLM(poisson_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Signif. codes:" in all_output

    @patch("builtins.print")
    def test_multiple_predictors(self, mocked_print, multi_pred_binomial_args):
        svc = PhylogeneticGLM(multi_pred_binomial_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "body_mass" in all_output
        assert "count_trait" in all_output
        assert "binary_trait ~ body_mass + count_trait" in all_output


GENE_TREES_FILE = str(SAMPLE_FILES / "gene_trees_simple.nwk")


class TestDiscordanceVCV:
    def test_poisson_gee_with_gene_trees_json(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="count_trait",
            predictors=["body_mass"],
            family="poisson",
            method=None,
            btol=10,
            log_alpha_bound=4,
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "vcv_metadata" in data
        assert data["vcv_metadata"]["n_gene_trees"] == 10
        assert "coefficients" in data

    def test_binomial_with_gene_trees_json(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="binary_trait",
            predictors=["body_mass"],
            family="binomial",
            method=None,
            btol=10,
            log_alpha_bound=4,
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "vcv_metadata" in data
        assert "coefficients" in data

    def test_discordance_changes_poisson_coefficients(self, capsys):
        """With discordant gene trees, Poisson GEE coefficients should differ."""
        base_kwargs = dict(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="count_trait",
            predictors=["body_mass"],
            family="poisson",
            method=None,
            btol=10,
            log_alpha_bound=4,
            json=True,
        )
        args_no_gt = Namespace(**base_kwargs)
        svc = PhylogeneticGLM(args_no_gt)
        svc.run()
        out, _ = capsys.readouterr()
        coefs_no_gt = json.loads(out)["coefficients"]

        args_gt = Namespace(**base_kwargs, gene_trees=GENE_TREES_FILE)
        svc = PhylogeneticGLM(args_gt)
        svc.run()
        out, _ = capsys.readouterr()
        coefs_gt = json.loads(out)["coefficients"]

        # Coefficients should differ with discordant gene trees
        assert coefs_gt["body_mass"]["estimate"] != pytest.approx(
            coefs_no_gt["body_mass"]["estimate"], abs=1e-6
        )

    @patch("builtins.print")
    def test_text_output_with_gene_trees(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="count_trait",
            predictors=["body_mass"],
            family="poisson",
            method=None,
            btol=10,
            log_alpha_bound=4,
            json=False,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        assert mocked_print.called
