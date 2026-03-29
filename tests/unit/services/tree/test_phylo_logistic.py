import json
import sys

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from phykit.services.tree.phylo_logistic import PhyloLogistic
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
BINARY_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_binary_traits.tsv")
GLM_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_glm_traits.tsv")


@pytest.fixture
def basic_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=BINARY_TRAITS_FILE,
        response="has_wings",
        predictor="body_mass",
        method="logistic_MPLE",
        json=False,
    )


@pytest.fixture
def json_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=BINARY_TRAITS_FILE,
        response="has_wings",
        predictor="body_mass",
        method="logistic_MPLE",
        json=True,
    )


@pytest.fixture
def glm_traits_args():
    """Uses the GLM traits file which has binary_trait column."""
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=GLM_TRAITS_FILE,
        response="binary_trait",
        predictor="body_mass",
        method="logistic_MPLE",
        json=False,
    )


@pytest.fixture
def multi_predictor_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=GLM_TRAITS_FILE,
        response="binary_trait",
        predictor="body_mass,count_trait",
        method="logistic_MPLE",
        json=False,
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            response="y",
            predictor="x1",
            method="logistic_MPLE",
        )
        svc = PhyloLogistic.__new__(PhyloLogistic)
        parsed = svc.process_args(args)
        assert parsed["tree_file_path"] == "t.tre"
        assert parsed["trait_data_path"] == "d.tsv"
        assert parsed["response"] == "y"
        assert parsed["predictors"] == ["x1"]
        assert parsed["method"] == "logistic_MPLE"
        assert parsed["json_output"] is False

    def test_comma_separated_predictors(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            response="y",
            predictor="x1,x2,x3",
            method="logistic_MPLE",
        )
        svc = PhyloLogistic.__new__(PhyloLogistic)
        parsed = svc.process_args(args)
        assert parsed["predictors"] == ["x1", "x2", "x3"]


class TestTransformedVCV:
    def test_transformed_vcv_alpha_zero(self, basic_args):
        """When alpha ~ 0, the transformed VCV should be the standard BM VCV."""
        from phykit.services.tree.vcv_utils import build_vcv_matrix

        svc = PhyloLogistic(basic_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(BINARY_TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())

        vcv_bm = build_vcv_matrix(tree, ordered_names)
        vcv_ou, diag_corr = svc._build_logistic_vcv(
            tree, 1e-12, ordered_names, build_vcv_matrix
        )

        # diag_corr should be all ones for alpha ~ 0
        np.testing.assert_allclose(diag_corr, np.ones(len(ordered_names)), atol=1e-6)
        # VCV should be the standard BM VCV
        np.testing.assert_allclose(vcv_ou, vcv_bm, atol=1e-4)

    def test_transformed_vcv_positive_alpha(self, basic_args):
        """Positive alpha should give a valid (PD) VCV matrix."""
        from phykit.services.tree.vcv_utils import build_vcv_matrix

        svc = PhyloLogistic(basic_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = svc._parse_multi_trait_file(BINARY_TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())

        vcv, diag_corr = svc._build_logistic_vcv(
            tree, 0.05, ordered_names, build_vcv_matrix
        )

        # VCV should be symmetric
        np.testing.assert_allclose(vcv, vcv.T, atol=1e-10)

        # VCV should be positive definite
        eigenvalues = np.linalg.eigvalsh(vcv)
        assert np.all(eigenvalues > -1e-10), f"VCV has negative eigenvalue: {min(eigenvalues)}"

        # Diagonal correction should be positive and < 1
        assert np.all(diag_corr > 0)
        assert np.all(diag_corr <= 1.0 + 1e-10)


class TestFit:
    def test_fit_returns_coefficients(self, glm_traits_args):
        """Fitting should return finite beta, SE, z, p for each coefficient."""
        from phykit.services.tree.vcv_utils import build_vcv_matrix

        svc = PhyloLogistic(glm_traits_args)
        tree = svc.read_tree_file()
        svc.validate_tree(tree, min_tips=3, require_branch_lengths=True)
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

        result = svc._fit(tree, y, X, ordered_names, build_vcv_matrix)

        assert "coefficients" in result
        intercept = result["coefficients"]["(Intercept)"]
        assert np.isfinite(intercept["estimate"])
        assert np.isfinite(intercept["std_error"])
        assert intercept["std_error"] > 0
        assert np.isfinite(intercept["z_value"])
        assert np.isfinite(intercept["p_value"])
        assert 0 <= intercept["p_value"] <= 1

    def test_alpha_positive(self, glm_traits_args):
        """Estimated alpha should be positive."""
        from phykit.services.tree.vcv_utils import build_vcv_matrix

        svc = PhyloLogistic(glm_traits_args)
        tree = svc.read_tree_file()
        svc.validate_tree(tree, min_tips=3, require_branch_lengths=True)
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

        result = svc._fit(tree, y, X, ordered_names, build_vcv_matrix)
        assert result["alpha"] > 0

    def test_convergence(self, glm_traits_args):
        """Optimizer should converge (code 0 or 2 for L-BFGS-B)."""
        from phykit.services.tree.vcv_utils import build_vcv_matrix

        svc = PhyloLogistic(glm_traits_args)
        tree = svc.read_tree_file()
        svc.validate_tree(tree, min_tips=3, require_branch_lengths=True)
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

        result = svc._fit(tree, y, X, ordered_names, build_vcv_matrix)
        # 0 = converged, 1 = max iterations, 2 = abnormal termination in line search
        # For small datasets, code 2 is acceptable (numerical precision at boundary)
        assert result["convergence"] in (0, 2)


class TestValidation:
    def test_missing_column(self, basic_args):
        basic_args.response = "nonexistent_column"
        svc = PhyloLogistic(basic_args)
        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()
        assert any("nonexistent_column" in m for m in exc_info.value.messages)

    def test_response_in_predictors(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=BINARY_TRAITS_FILE,
            response="has_wings",
            predictor="has_wings",
            method="logistic_MPLE",
            json=False,
        )
        svc = PhyloLogistic(args)
        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()
        assert any("must not also be a predictor" in m for m in exc_info.value.messages)


class TestRun:
    @patch("builtins.print")
    def test_creates_results(self, mocked_print, glm_traits_args):
        """Text output should contain key sections."""
        svc = PhyloLogistic(glm_traits_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Phylogenetic Logistic Regression" in all_output
        assert "(Intercept)" in all_output
        assert "body_mass" in all_output
        assert "Alpha:" in all_output
        assert "Log-likelihood:" in all_output
        assert "Penalized log-likelihood:" in all_output
        assert "AIC:" in all_output
        assert "Signif. codes:" in all_output

    def test_json_output(self, capsys):
        """JSON should have all expected fields."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=GLM_TRAITS_FILE,
            response="binary_trait",
            predictor="body_mass",
            method="logistic_MPLE",
            json=True,
        )
        svc = PhyloLogistic(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert data["method"] == "logistic_MPLE"
        assert data["response"] == "binary_trait"
        assert data["predictors"] == ["body_mass"]
        assert "alpha" in data
        assert data["alpha"] > 0
        assert "(Intercept)" in data["coefficients"]
        assert "body_mass" in data["coefficients"]
        assert "log_likelihood" in data
        assert "penalized_log_likelihood" in data
        assert "aic" in data
        assert "convergence" in data
        assert "n_taxa" in data

    @patch("builtins.print")
    def test_multiple_predictors(self, mocked_print, multi_predictor_args):
        """Comma-separated predictors should work."""
        svc = PhyloLogistic(multi_predictor_args)
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "body_mass" in all_output
        assert "count_trait" in all_output
        assert "binary_trait ~ body_mass + count_trait" in all_output
