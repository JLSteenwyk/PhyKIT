import builtins
import json
import os
import subprocess
import sys
import tempfile

import pytest
import numpy as np
from argparse import Namespace
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, TreeMixin

import phykit.services.tree.ancestral_reconstruction as ancestral_module
from phykit.services.tree.ancestral_reconstruction import (
    AncestralReconstruction,
    _value_range,
)
from phykit.errors import PhykitUserError


def test_module_import_does_not_import_numpy():
    code = """
import sys
import phykit.services.tree.ancestral_reconstruction as module
assert callable(module.print_json)
assert hasattr(module.np, "__getattr__")
assert hasattr(module.pickle, "dumps")
assert hasattr(module.heapq, "merge")
assert "heapq" not in sys.modules
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "numpy" not in sys.modules
assert "pickle" not in sys.modules
assert "phykit.services.tree.vcv_utils" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.discrete_models" not in sys.modules
assert "phykit.helpers.circular_layout" not in sys.modules
assert "phykit.helpers.color_annotations" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = ancestral_module._LazyNumpy()

    ones_attr = lazy_np.ones

    assert lazy_np.__dict__["ones"] is ones_attr
    assert lazy_np.ones is ones_attr
    assert lazy_np._module is not None


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")
DISCRETE_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_discrete_traits.tsv")


@pytest.fixture(autouse=True)
def clear_discrete_asr_plot_cache():
    previous_cache = ancestral_module._DISCRETE_ASR_PLOT_CACHE.copy()
    ancestral_module._DISCRETE_ASR_PLOT_CACHE.clear()
    try:
        yield
    finally:
        ancestral_module._DISCRETE_ASR_PLOT_CACHE.clear()
        ancestral_module._DISCRETE_ASR_PLOT_CACHE.update(previous_cache)


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        trait=None,
        method="fast",
        ci=False,
        plot=None,
        json=False,
    )


@pytest.fixture
def ci_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        trait=None,
        method="fast",
        ci=True,
        plot=None,
        json=False,
    )


@pytest.fixture
def ml_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        trait=None,
        method="ml",
        ci=True,
        plot=None,
        json=False,
    )


@pytest.fixture
def json_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        trait=None,
        method="fast",
        ci=False,
        plot=None,
        json=True,
    )


@pytest.fixture
def multi_trait_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        trait="body_mass",
        method="fast",
        ci=False,
        plot=None,
        json=False,
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            trait=None,
            method="fast",
            ci=False,
            plot=None,
            json=False,
        )
        svc = AncestralReconstruction(args)
        assert svc.method == "fast"
        assert svc.ci is False
        assert svc.plot_output is None
        assert svc.json_output is False
        assert svc.trait_column is None

    def test_overrides(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            trait="body_mass",
            method="ml",
            ci=True,
            plot="out.png",
            json=True,
        )
        svc = AncestralReconstruction(args)
        assert svc.method == "ml"
        assert svc.ci is True
        assert svc.plot_output == "out.png"
        assert svc.json_output is True
        assert svc.trait_column == "body_mass"


class TestTraitParsing:
    def test_single_trait_format(self, default_args):
        svc = AncestralReconstruction(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        traits = svc._parse_single_trait_data(TRAITS_FILE, tree_tips)
        assert len(traits) == 8
        assert pytest.approx(traits["raccoon"]) == 1.04
        assert pytest.approx(traits["bear"]) == 2.39
        assert pytest.approx(traits["weasel"]) == -0.30

    def test_single_trait_skips_comments_and_blanks(self, default_args, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "   # ignored\n"
            "\n"
            "raccoon\t1.04\n"
            "bear\t2.39\n"
            "weasel\t-0.30\n"
        )
        svc = AncestralReconstruction(default_args)

        traits = svc._parse_single_trait_data(
            str(trait_file), ["raccoon", "bear", "weasel"]
        )

        assert traits == {"raccoon": 1.04, "bear": 2.39, "weasel": -0.30}

    def test_single_trait_all_shared_emits_no_warnings(
        self, default_args, tmp_path, capsys
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("raccoon\t1.04\nbear\t2.39\nweasel\t-0.30\n")
        svc = AncestralReconstruction(default_args)

        traits = svc._parse_single_trait_data(
            str(trait_file), ["raccoon", "bear", "weasel"]
        )

        assert traits == {"raccoon": 1.04, "bear": 2.39, "weasel": -0.30}
        assert capsys.readouterr().err == ""

    def test_single_trait_wrong_column_count(self, default_args, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "raccoon\t1.04\n"
            "bear\t2.39\textra\n"
            "weasel\t-0.30\n"
        )
        svc = AncestralReconstruction(default_args)

        with pytest.raises(PhykitUserError):
            svc._parse_single_trait_data(
                str(trait_file), ["raccoon", "bear", "weasel"]
            )

    def test_multi_trait_format(self, default_args):
        svc = AncestralReconstruction(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        traits = svc._parse_multi_trait_data(
            MULTI_TRAITS_FILE, tree_tips, "body_mass"
        )
        assert len(traits) == 8
        assert pytest.approx(traits["raccoon"]) == 1.04
        assert pytest.approx(traits["bear"]) == 2.39

    def test_multi_trait_missing_column(self, default_args):
        svc = AncestralReconstruction(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        with pytest.raises(PhykitUserError):
            svc._parse_multi_trait_data(
                MULTI_TRAITS_FILE, tree_tips, "nonexistent_col"
            )

    def test_multi_trait_skips_comments_while_preserving_data_row_numbers(
        self, default_args, tmp_path
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "\n".join(
                [
                    "   # ignored before header",
                    "taxon\tbody_mass\tlength",
                    "\t# ignored between rows",
                    "raccoon\t1.04\t10",
                    "bear\t2.39\t20",
                    "weasel\t-0.30\t30",
                ]
            )
            + "\n"
        )
        svc = AncestralReconstruction(default_args)

        traits = svc._parse_multi_trait_data(
            str(trait_file), ["raccoon", "bear", "weasel"], "body_mass"
        )

        assert traits == {"raccoon": 1.04, "bear": 2.39, "weasel": -0.30}

        trait_file.write_text(
            "\n".join(
                [
                    "taxon\tbody_mass\tlength",
                    "raccoon\t1.04\t10",
                    "# ignored between rows",
                    "bear\t2.39",
                ]
            )
            + "\n"
        )

        with pytest.raises(PhykitUserError) as excinfo:
            svc._parse_multi_trait_data(
                str(trait_file), ["raccoon", "bear", "weasel"], "body_mass"
            )
        assert "Line 3 has 2 columns; expected 3." in excinfo.value.messages

    def test_multi_trait_all_shared_emits_no_warnings(
        self, default_args, tmp_path, capsys
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "taxon\tbody_mass\tlength\n"
            "raccoon\t1.04\t10\n"
            "bear\t2.39\t20\n"
            "weasel\t-0.30\t30\n"
        )
        svc = AncestralReconstruction(default_args)

        traits = svc._parse_multi_trait_data(
            str(trait_file), ["raccoon", "bear", "weasel"], "body_mass"
        )

        assert traits == {"raccoon": 1.04, "bear": 2.39, "weasel": -0.30}
        assert capsys.readouterr().err == ""

    def test_file_not_found(self, default_args):
        svc = AncestralReconstruction(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_single_trait_data("nonexistent.tsv", [])


class TestFastAnc:
    def test_point_estimates(self, ci_args):
        svc = AncestralReconstruction(ci_args)
        tree = svc.read_tree_file()
        import copy
        tree_copy = copy.deepcopy(tree)
        tree_tips = svc.get_tip_names_from_tree(tree_copy)
        traits = svc._parse_single_trait_data(TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())
        x = np.array([traits[name] for name in ordered_names])
        node_labels = svc._label_internal_nodes(tree_copy)

        estimates, cis, sigma2, ll = svc._fast_anc(
            tree_copy, x, ordered_names, node_labels
        )

        # Root estimate should be present
        root_label = None
        for clade in tree_copy.find_clades(order="preorder"):
            if not clade.is_terminal():
                root_label = node_labels[id(clade)]
                break
        assert root_label in estimates
        # Sigma^2 should be positive
        assert sigma2 > 0
        # Log-likelihood should be finite
        assert np.isfinite(ll)

    def test_sigma2_positive(self, ci_args):
        svc = AncestralReconstruction(ci_args)
        tree = svc.read_tree_file()
        import copy
        tree_copy = copy.deepcopy(tree)
        tree_tips = svc.get_tip_names_from_tree(tree_copy)
        traits = svc._parse_single_trait_data(TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())
        x = np.array([traits[name] for name in ordered_names])
        node_labels = svc._label_internal_nodes(tree_copy)

        _, _, sigma2, _ = svc._fast_anc(
            tree_copy, x, ordered_names, node_labels
        )
        assert sigma2 > 0

    def test_cis_present(self, ci_args):
        svc = AncestralReconstruction(ci_args)
        tree = svc.read_tree_file()
        import copy
        tree_copy = copy.deepcopy(tree)
        tree_tips = svc.get_tip_names_from_tree(tree_copy)
        traits = svc._parse_single_trait_data(TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())
        x = np.array([traits[name] for name in ordered_names])
        node_labels = svc._label_internal_nodes(tree_copy)

        estimates, cis, _, _ = svc._fast_anc(
            tree_copy, x, ordered_names, node_labels
        )
        # CIs should be present since ci=True
        assert len(cis) > 0
        for label, (lo, hi) in cis.items():
            assert lo < estimates[label] < hi

    def test_fast_anc_uses_direct_tree_traversal(self, ci_args, monkeypatch):
        svc = AncestralReconstruction(ci_args)
        tree = Phylo.read(
            StringIO("((A:1,B:2):1,(C:1,D:2):1):0;"),
            "newick",
        )
        ordered_names = ["A", "B", "C", "D"]
        x = np.array([0.0, 1.0, 2.0, 3.0])
        node_labels = svc._label_internal_nodes(tree)
        expected = svc._fast_anc(tree, x, ordered_names, node_labels)

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError(
                "fast ancestral reconstruction should use direct traversal"
            )

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        observed = svc._fast_anc(tree, x, ordered_names, node_labels)

        assert observed[0] == pytest.approx(expected[0])
        assert observed[1] == pytest.approx(expected[1])
        assert observed[2] == pytest.approx(expected[2])
        assert observed[3] == pytest.approx(expected[3])

    def test_fast_anc_reuses_preorder_for_parent_map(self, ci_args, monkeypatch):
        svc = AncestralReconstruction(ci_args)
        tree = Phylo.read(
            StringIO("((A:1,B:2):1,(C:1,D:2):1):0;"),
            "newick",
        )
        ordered_names = ["A", "B", "C", "D"]
        x = np.array([0.0, 1.0, 2.0, 3.0])
        node_labels = svc._label_internal_nodes(tree)

        def fail_parent_map(*_args, **_kwargs):
            raise AssertionError("parent map should reuse preorder clades")

        monkeypatch.setattr(svc, "_build_parent_map", fail_parent_map)

        estimates, cis, sigma2, ll = svc._fast_anc(
            tree, x, ordered_names, node_labels
        )

        assert set(estimates) == {"N1", "N2", "N3"}
        assert set(cis) == {"N1", "N2", "N3"}
        assert sigma2 > 0
        assert np.isfinite(ll)

    def test_sigma2_from_contrasts_streams_polytomy_children(self, ci_args):
        svc = AncestralReconstruction(ci_args)
        tree = Phylo.read(StringIO("(A:1,B:1,C:1);"), "newick")
        values = {
            tip.name: float(index + 1)
            for index, tip in enumerate(tree.root.clades)
        }
        node_estimates = {id(tip): values[tip.name] for tip in tree.root.clades}
        node_variances = {id(tip): 0.0 for tip in tree.root.clades}

        sigma2 = svc._compute_sigma2_from_contrasts(
            tree,
            np.array([1.0, 2.0, 3.0]),
            node_estimates,
            node_variances,
            ["A", "B", "C"],
        )

        first = ((1.0 - 2.0) ** 2) / 2.0
        combined_est = (1.0 + 2.0) / 2.0
        combined_var = 0.5
        second = ((3.0 - combined_est) ** 2) / (1.0 + combined_var)
        assert sigma2 == pytest.approx((first + second) / 2.0)


class TestAncML:
    def test_point_estimates(self, ml_args):
        svc = AncestralReconstruction(ml_args)
        tree = svc.read_tree_file()
        import copy
        tree_copy = copy.deepcopy(tree)
        tree_tips = svc.get_tip_names_from_tree(tree_copy)
        traits = svc._parse_single_trait_data(TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())
        x = np.array([traits[name] for name in ordered_names])
        node_labels = svc._label_internal_nodes(tree_copy)

        estimates, cis, sigma2, ll = svc._anc_ml(
            tree_copy, x, ordered_names, node_labels
        )

        assert len(estimates) > 0
        assert sigma2 > 0
        assert np.isfinite(ll)

    def test_cis_present(self, ml_args):
        svc = AncestralReconstruction(ml_args)
        tree = svc.read_tree_file()
        import copy
        tree_copy = copy.deepcopy(tree)
        tree_tips = svc.get_tip_names_from_tree(tree_copy)
        traits = svc._parse_single_trait_data(TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())
        x = np.array([traits[name] for name in ordered_names])
        node_labels = svc._label_internal_nodes(tree_copy)

        estimates, cis, _, _ = svc._anc_ml(
            tree_copy, x, ordered_names, node_labels
        )
        assert len(cis) > 0
        # For non-root nodes, CI lower should be < estimate < CI upper
        for label, (lo, hi) in cis.items():
            assert lo <= estimates[label] <= hi

    def test_log_likelihood(self, ml_args):
        svc = AncestralReconstruction(ml_args)
        tree = svc.read_tree_file()
        import copy
        tree_copy = copy.deepcopy(tree)
        tree_tips = svc.get_tip_names_from_tree(tree_copy)
        traits = svc._parse_single_trait_data(TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())
        x = np.array([traits[name] for name in ordered_names])
        node_labels = svc._label_internal_nodes(tree_copy)

        _, _, _, ll = svc._anc_ml(
            tree_copy, x, ordered_names, node_labels
        )
        assert ll < 0  # Log-likelihood should be negative for real data

    def test_cross_covariance_fast_path_does_not_call_distance_get_path_or_mrca(
        self, ml_args, monkeypatch
    ):
        svc = AncestralReconstruction(ml_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:2):0.5;"), "newick")
        ordered_names = ["A", "B", "C"]
        x = np.array([0.0, 2.0, 4.0])
        node_labels = svc._label_internal_nodes(tree)
        expected = svc._anc_ml(tree, x, ordered_names, node_labels)

        tree = Phylo.read(StringIO("((A:1,B:1):1,C:2):0.5;"), "newick")
        node_labels = svc._label_internal_nodes(tree)

        def fail_distance(self, *args, **kwargs):
            raise AssertionError("distance fallback should not be called")

        def fail_common_ancestor(self, *args, **kwargs):
            raise AssertionError("common_ancestor fallback should not be called")

        def fail_get_path(self, *args, **kwargs):
            raise AssertionError("get_path should not be called")

        monkeypatch.setattr(type(tree), "distance", fail_distance)
        monkeypatch.setattr(type(tree), "common_ancestor", fail_common_ancestor)
        monkeypatch.setattr(type(tree), "get_path", fail_get_path)
        observed = svc._anc_ml(tree, x, ordered_names, node_labels)

        assert observed[0] == pytest.approx(expected[0])
        assert observed[1] == pytest.approx(expected[1])
        assert observed[2] == pytest.approx(expected[2])
        assert observed[3] == pytest.approx(expected[3])

    def test_ml_cross_covariance_uses_direct_tree_traversal(
        self, ml_args, monkeypatch
    ):
        svc = AncestralReconstruction(ml_args)
        tree = Phylo.read(StringIO("((A:1,B:2):1,(C:3,D:4):2):0;"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        x = np.array([0.0, 1.0, 3.0, 4.0])
        node_labels = svc._label_internal_nodes(tree)
        expected = svc._anc_ml(tree, x, ordered_names, node_labels)

        def fail_traversal(*_args, **_kwargs):
            raise AssertionError("ML cross-covariance should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        observed = svc._anc_ml(tree, x, ordered_names, node_labels)

        assert observed[0] == pytest.approx(expected[0])
        assert observed[1] == pytest.approx(expected[1])
        assert observed[2] == pytest.approx(expected[2])
        assert observed[3] == pytest.approx(expected[3])

    def test_batched_ml_estimates_and_cis_match_scalar_formulas(self, ml_args):
        svc = AncestralReconstruction(ml_args)
        tree = Phylo.read(StringIO("((A:1,B:2):1,(C:3,D:4):2):0;"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        x = np.array([0.0, 1.0, 3.0, 4.0])
        node_labels = svc._label_internal_nodes(tree)

        estimates, cis, sigma2_reml, _ = svc._anc_ml(
            tree, x, ordered_names, node_labels
        )

        vcv = svc._build_vcv_matrix(tree, ordered_names)
        c_inv = np.linalg.inv(vcv)
        ones = np.ones(len(ordered_names))
        a_hat = float(ones @ c_inv @ x) / float(ones @ c_inv @ ones)
        residuals = x - a_hat
        expected_sigma2 = float(residuals @ c_inv @ residuals) / (
            len(ordered_names) - 1
        )
        assert sigma2_reml == pytest.approx(expected_sigma2)

        depths = tree.depths()
        root = tree.root
        root_depth = depths[root]
        terminals = {terminal.name: terminal for terminal in tree.get_terminals()}

        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            label = node_labels[id(clade)]
            c_i = np.zeros(len(ordered_names))
            for idx, tip_name in enumerate(ordered_names):
                mrca = tree.common_ancestor(clade, terminals[tip_name])
                c_i[idx] = depths[mrca] - root_depth

            expected_mean = a_hat + float(c_i @ c_inv @ residuals)
            assert estimates[label] == pytest.approx(expected_mean)

            if clade is root:
                cond_var = expected_sigma2 / float(ones @ c_inv @ ones)
            else:
                d_i = depths[clade] - root_depth
                cond_var = expected_sigma2 * (d_i - float(c_i @ c_inv @ c_i))
            cond_var = max(cond_var, 0.0)
            se = np.sqrt(cond_var)
            assert cis[label] == pytest.approx(
                (expected_mean - 1.96 * se, expected_mean + 1.96 * se)
            )

    def test_compute_log_likelihood_cholesky_matches_inverse_formula(self, ml_args):
        svc = AncestralReconstruction(ml_args)
        tree = Phylo.read(StringIO("((A:1,B:2):1,(C:3,D:4):2):0;"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        x = np.array([0.0, 1.0, 3.0, 4.0])
        sigma2 = 0.75

        observed = svc._compute_log_likelihood(tree, x, ordered_names, sigma2)

        vcv = svc._build_vcv_matrix(tree, ordered_names)
        C = sigma2 * vcv
        C_inv = np.linalg.inv(C)
        sign, logdet = np.linalg.slogdet(C)
        ones = np.ones(len(ordered_names))
        a_hat = float(ones @ C_inv @ x) / float(ones @ C_inv @ ones)
        residuals = x - a_hat
        expected = -0.5 * (
            len(ordered_names) * np.log(2 * np.pi)
            + logdet
            + float(residuals @ C_inv @ residuals)
        )

        assert sign > 0
        assert observed == pytest.approx(expected)

    def test_compute_log_likelihood_cholesky_avoids_inverse(
        self, ml_args, monkeypatch
    ):
        svc = AncestralReconstruction(ml_args)
        tree = Phylo.read(StringIO("((A:1,B:2):1,(C:3,D:4):2):0;"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        x = np.array([0.0, 1.0, 3.0, 4.0])

        def fail_inv(*_args, **_kwargs):
            raise AssertionError("positive-definite log-likelihood should use Cholesky")

        monkeypatch.setattr(ancestral_module.np.linalg, "inv", fail_inv)

        observed = svc._compute_log_likelihood(
            tree, x, ordered_names, sigma2=0.75
        )

        assert np.isfinite(observed)

    def test_ml_reuses_inverse_weighted_residuals(self, ml_args, monkeypatch):
        args = Namespace(**vars(ml_args))
        args.ci = False
        svc = AncestralReconstruction(args)
        tree = Phylo.read(StringIO("((A:1,B:2):1,(C:3,D:4):2):0;"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        x = np.array([0.0, 1.0, 3.0, 4.0])
        node_labels = svc._label_internal_nodes(tree)
        expected = svc._anc_ml(tree, x, ordered_names, node_labels)
        original_inv = np.linalg.inv
        vector_products = 0

        class CountingInverse:
            def __init__(self, matrix):
                self.matrix = matrix

            def __matmul__(self, other):
                nonlocal vector_products
                if getattr(other, "ndim", 0) == 1:
                    vector_products += 1
                    if vector_products > 2:
                        raise AssertionError(
                            "C_inv @ residuals should reuse existing products"
                        )
                return self.matrix @ other

        def counting_inv(matrix):
            return CountingInverse(original_inv(matrix))

        monkeypatch.setattr(ancestral_module.np.linalg, "inv", counting_inv)

        observed = svc._anc_ml(tree, x, ordered_names, node_labels)

        assert vector_products == 2
        assert observed[0] == pytest.approx(expected[0])
        assert observed[1] == pytest.approx(expected[1])
        assert observed[2] == pytest.approx(expected[2])
        assert observed[3] == pytest.approx(expected[3])


class TestMethodAgreement:
    """fast and ml should produce identical point estimates."""

    def test_estimates_match(self):
        import copy

        fast_args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            trait=None,
            method="fast",
            ci=False,
            plot=None,
            json=False,
        )
        ml_args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            trait=None,
            method="ml",
            ci=False,
            plot=None,
            json=False,
        )

        svc_fast = AncestralReconstruction(fast_args)
        tree_fast = svc_fast.read_tree_file()
        tree_fast = copy.deepcopy(tree_fast)
        tree_tips = svc_fast.get_tip_names_from_tree(tree_fast)
        traits = svc_fast._parse_single_trait_data(TRAITS_FILE, tree_tips)
        ordered_names = sorted(traits.keys())
        x = np.array([traits[name] for name in ordered_names])
        labels_fast = svc_fast._label_internal_nodes(tree_fast)
        est_fast, _, _, _ = svc_fast._fast_anc(
            tree_fast, x, ordered_names, labels_fast
        )

        svc_ml = AncestralReconstruction(ml_args)
        tree_ml = svc_ml.read_tree_file()
        tree_ml = copy.deepcopy(tree_ml)
        labels_ml = svc_ml._label_internal_nodes(tree_ml)
        est_ml, _, _, _ = svc_ml._anc_ml(
            tree_ml, x, ordered_names, labels_ml
        )

        # Both should have same labels and identical estimates
        common_labels = set(est_fast.keys()) & set(est_ml.keys())
        assert len(common_labels) == len(est_fast)
        for label in common_labels:
            assert est_fast[label] == pytest.approx(est_ml[label], abs=1e-6)


class TestRValidation:
    """Validate against R 4.4.0 phytools::fastAnc() and ape::ace().

    R validation code:
        library(phytools); library(ape)
        tree <- read.tree("tree_simple.tre")
        x <- c(raccoon=1.04, bear=2.39, sea_lion=2.30, seal=1.88,
               monkey=0.60, cat=0.56, weasel=-0.30, dog=1.18)
        fa <- fastAnc(tree, x, vars=TRUE, CI=TRUE)
    """

    # R fastAnc results keyed by descendant sets:
    # Node 9  (all 8 tips):     est=1.6446924, CI=[0.8936713, 2.3957134]
    # Node 10 (bear, raccoon):  est=1.7012405, CI=[0.9696667, 2.4328142]
    # Node 11 (5 taxa):         est=1.4564597, CI=[0.6386751, 2.2742444]
    # Node 12 (sea_lion, seal): est=1.8090745, CI=[0.9756607, 2.6424883]
    # Node 13 (cat,monkey,weasel): est=1.2565917, CI=[0.3554623, 2.1577210]
    # Node 14 (cat, monkey):    est=0.9894725, CI=[-0.5653569, 2.5443018]

    R_ESTIMATES = {
        frozenset(["bear", "cat", "dog", "monkey", "raccoon",
                    "sea_lion", "seal", "weasel"]): 1.6446924,
        frozenset(["bear", "raccoon"]): 1.7012405,
        frozenset(["cat", "monkey", "sea_lion", "seal", "weasel"]): 1.4564597,
        frozenset(["sea_lion", "seal"]): 1.8090745,
        frozenset(["cat", "monkey", "weasel"]): 1.2565917,
        frozenset(["cat", "monkey"]): 0.9894725,
    }

    R_CIS = {
        frozenset(["bear", "cat", "dog", "monkey", "raccoon",
                    "sea_lion", "seal", "weasel"]): (0.8936713, 2.3957134),
        frozenset(["bear", "raccoon"]): (0.9696667, 2.4328142),
        frozenset(["cat", "monkey", "sea_lion", "seal", "weasel"]): (0.6386751, 2.2742444),
        frozenset(["sea_lion", "seal"]): (0.9756607, 2.6424883),
        frozenset(["cat", "monkey", "weasel"]): (0.3554623, 2.1577210),
        frozenset(["cat", "monkey"]): (-0.5653569, 2.5443018),
    }

    def _get_estimates_by_descendants(self, method):
        import copy
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE, trait=None,
            method=method, ci=True, plot=None, json=False,
        )
        svc = AncestralReconstruction(args)
        tree = svc.read_tree_file()
        tc = copy.deepcopy(tree)
        tips = svc.get_tip_names_from_tree(tc)
        traits = svc._parse_single_trait_data(TRAITS_FILE, tips)
        ordered = sorted(traits.keys())
        x = np.array([traits[n] for n in ordered])
        labels = svc._label_internal_nodes(tc)

        if method == "fast":
            est, cis, sigma2, ll = svc._fast_anc(tc, x, ordered, labels)
        else:
            est, cis, sigma2, ll = svc._anc_ml(tc, x, ordered, labels)

        # Map labels to descendant sets
        result = {}
        ci_result = {}
        for clade in tc.find_clades(order="preorder"):
            if not clade.is_terminal() and id(clade) in labels:
                label = labels[id(clade)]
                descs = frozenset(svc._get_descendant_tips(tc, clade))
                if label in est:
                    result[descs] = est[label]
                if label in cis:
                    ci_result[descs] = cis[label]
        return result, ci_result, sigma2

    def test_fast_estimates_match_r_fastAnc(self):
        est, _, _ = self._get_estimates_by_descendants("fast")
        for descs, r_val in self.R_ESTIMATES.items():
            assert est[descs] == pytest.approx(r_val, abs=1e-4), (
                f"fast estimate for {descs} = {est[descs]}, R = {r_val}"
            )

    def test_fast_cis_match_r_fastAnc(self):
        _, cis, _ = self._get_estimates_by_descendants("fast")
        for descs, (r_lo, r_hi) in self.R_CIS.items():
            pk_lo, pk_hi = cis[descs]
            assert pk_lo == pytest.approx(r_lo, abs=1e-4), (
                f"fast CI lower for {descs}: PK={pk_lo}, R={r_lo}"
            )
            assert pk_hi == pytest.approx(r_hi, abs=1e-4), (
                f"fast CI upper for {descs}: PK={pk_hi}, R={r_hi}"
            )

    def test_ml_estimates_match_r_fastAnc(self):
        est, _, _ = self._get_estimates_by_descendants("ml")
        for descs, r_val in self.R_ESTIMATES.items():
            assert est[descs] == pytest.approx(r_val, abs=1e-4), (
                f"ml estimate for {descs} = {est[descs]}, R = {r_val}"
            )

    def test_fast_sigma2_matches_r(self):
        """R's pic-based sigma2 = 0.04389322."""
        _, _, sigma2 = self._get_estimates_by_descendants("fast")
        assert sigma2 == pytest.approx(0.04389322, abs=1e-3)


class TestNodeLabeling:
    def test_internal_nodes_labeled(self, default_args):
        svc = AncestralReconstruction(default_args)
        tree = svc.read_tree_file()
        labels = svc._label_internal_nodes(tree)
        # All internal nodes should have labels
        for clade in tree.find_clades(order="preorder"):
            if not clade.is_terminal():
                assert id(clade) in labels
                assert labels[id(clade)].startswith("N")

    def test_labels_unique(self, default_args):
        svc = AncestralReconstruction(default_args)
        tree = svc.read_tree_file()
        labels = svc._label_internal_nodes(tree)
        label_values = list(labels.values())
        assert len(label_values) == len(set(label_values))

    def test_label_internal_nodes_uses_direct_preorder_traversal(
        self, default_args, monkeypatch
    ):
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(
            StringIO("((A:1,B:1)left:1,(C:1,D:1):1);"),
            "newick",
        )

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("node labeling should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        labels = svc._label_internal_nodes(tree)

        left, right = tree.root.clades
        assert labels[id(tree.root)] == "N1"
        assert labels[id(left)] == "left"
        assert labels[id(right)] == "N2"

    def test_build_parent_map_uses_direct_preorder_traversal(
        self, default_args, monkeypatch
    ):
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("parent map should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        parent_map = svc._build_parent_map(tree)

        left, right = tree.root.clades
        assert parent_map[id(left)] is tree.root
        assert parent_map[id(right)] is tree.root
        assert parent_map[id(left.clades[0])] is left
        assert parent_map[id(right.clades[1])] is right

    def test_build_parent_map_handles_mixed_child_counts(
        self, default_args, monkeypatch
    ):
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("parent map should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        parent_map = svc._build_parent_map(tree)

        terminal, binary, trifurcating = tree.root.clades
        assert id(tree.root) not in parent_map
        assert parent_map[id(terminal)] is tree.root
        assert parent_map[id(binary)] is tree.root
        assert parent_map[id(trifurcating)] is tree.root
        assert all(parent_map[id(child)] is binary for child in binary.clades)
        assert all(
            parent_map[id(child)] is trifurcating for child in trifurcating.clades
        )

    def test_get_descendant_tips_uses_direct_clade_traversal(
        self, default_args, monkeypatch
    ):
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("descendant tips should use direct traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        assert svc._get_descendant_tips(tree, tree.root) == ["A", "B", "C", "D"]

    def test_collect_descendant_tip_counts_uses_direct_traversal(
        self, default_args, monkeypatch
    ):
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_traversal(*_args, **_kwargs):
            raise AssertionError("descendant counts should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        counts = svc._collect_descendant_tip_counts(tree)

        left, right = tree.root.clades
        assert counts[id(tree.root)] == 4
        assert counts[id(left)] == 2
        assert counts[id(right)] == 2
        assert counts[id(left.clades[0])] == 1

    def test_collect_descendant_tip_names_uses_direct_traversal(
        self, default_args, monkeypatch
    ):
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((B:1,A:1):1,(D:1,C:1):1);"), "newick")

        def fail_traversal(*_args, **_kwargs):
            raise AssertionError("descendant names should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        names = svc._collect_descendant_tip_names(tree)

        left, right = tree.root.clades
        assert list(names[id(tree.root)]) == ["A", "B", "C", "D"]
        assert list(names[id(left)]) == ["A", "B"]
        assert list(names[id(right)]) == ["C", "D"]

    def test_collect_descendant_tip_names_concatenates_ordered_ranges(
        self, default_args, monkeypatch
    ):
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_merge(*_args, **_kwargs):
            raise AssertionError("ordered descendant ranges should concatenate")

        monkeypatch.setattr(ancestral_module.heapq, "merge", fail_merge)

        names = svc._collect_descendant_tip_names(tree)

        assert list(names[id(tree.root)]) == ["A", "B", "C", "D"]


class TestRun:
    @patch("builtins.print")
    def test_text_output(self, mocked_print, default_args):
        svc = AncestralReconstruction(default_args)
        svc.run()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Ancestral State Reconstruction" in all_output
        assert "fast" in all_output
        assert "Sigma-squared" in all_output
        assert "N1 (root)" in all_output

    @patch("builtins.print")
    def test_json_output(self, mocked_print, json_args):
        svc = AncestralReconstruction(json_args)
        svc.run()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "method" in payload
        assert payload["method"] == "fast"
        assert "ancestral_estimates" in payload
        assert "tip_values" in payload
        assert "sigma2" in payload
        assert "log_likelihood" in payload
        assert "N1" in payload["ancestral_estimates"]
        assert payload["ancestral_estimates"]["N1"]["is_root"] is True

    def test_run_uses_fast_tip_name_helper_for_continuous_prune_setup(
        self, default_args, mocker
    ):
        svc = AncestralReconstruction(default_args)
        spy = mocker.spy(svc, "get_tip_names_from_tree")
        mocker.patch("builtins.print")

        svc.run()

        assert spy.call_count == 1

    def test_run_uses_unmodified_tree_read_for_continuous(self, default_args, mocker):
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tip_names = ["A", "B", "C", "D"]

        read_unmodified = mocker.patch.object(
            svc, "read_tree_file_unmodified", return_value=tree
        )
        mocker.patch.object(
            svc,
            "read_tree_file",
            side_effect=AssertionError("run should use read_tree_file_unmodified"),
        )
        validate = mocker.patch.object(svc, "validate_tree")
        run_continuous = mocker.patch.object(svc, "_run_continuous")

        svc.run()

        read_unmodified.assert_called_once_with()
        validate.assert_called_once()
        run_continuous.assert_called_once_with(tree, tip_names)

    def test_continuous_run_skips_copy_for_all_shared_tree(
        self, default_args, mocker
    ):
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tip_names = ["A", "B", "C", "D"]
        traits = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}

        mocker.patch.object(svc, "_parse_single_trait_data", return_value=traits)
        fast_copy = mocker.patch.object(
            svc,
            "_fast_copy",
            side_effect=AssertionError("all-shared analysis should not copy tree"),
        )
        mocker.patch.object(svc, "_label_internal_nodes", return_value={"N1": "N1"})
        fast_anc = mocker.patch.object(
            svc, "_fast_anc", return_value=({"N1": 1.0}, None, 1.0, -1.0)
        )
        mocker.patch.object(svc, "_format_result", return_value={"ok": True})
        mocker.patch.object(svc, "_print_text_output")

        svc._run_continuous(tree, tip_names)

        fast_copy.assert_not_called()
        assert fast_anc.call_args.args[0] is tree

    def test_continuous_run_copies_before_pruning_missing_tree_tips(
        self, default_args, monkeypatch, mocker
    ):
        class OrderedTraitValues(dict):
            def __contains__(self, key):
                raise AssertionError("ordered prune path should not scan membership")

        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree_copy = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        pruned_tree = Phylo.read(StringIO("((A:1,B:1):1,C:1);"), "newick")
        tip_names = ["A", "B", "C", "D"]
        traits = OrderedTraitValues({"A": 1.0, "B": 2.0, "C": 3.0})

        monkeypatch.setattr(ancestral_module.Tree, "_ORDERED_MAPPING_PRUNE_MIN_SIZE", 0)
        mocker.patch.object(svc, "_parse_single_trait_data", return_value=traits)
        fast_copy = mocker.patch.object(svc, "_fast_copy", return_value=tree_copy)
        prune = mocker.patch.object(
            svc, "prune_tree_using_taxa_list", return_value=pruned_tree
        )
        mocker.patch.object(svc, "_label_internal_nodes", return_value={"N1": "N1"})
        fast_anc = mocker.patch.object(
            svc, "_fast_anc", return_value=({"N1": 1.0}, None, 1.0, -1.0)
        )
        mocker.patch.object(svc, "_format_result", return_value={"ok": True})
        mocker.patch.object(svc, "_print_text_output")

        svc._run_continuous(tree, tip_names)

        fast_copy.assert_called_once_with(tree)
        prune.assert_called_once_with(tree_copy, ["D"])
        assert fast_anc.call_args.args[0] is pruned_tree

    def test_text_output_uses_cached_descendant_counts(
        self, default_args, monkeypatch, mocker
    ):
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        node_labels = svc._label_internal_nodes(tree)
        node_estimates = {label: 1.0 for label in node_labels.values()}

        def fail_traversal(*_args, **_kwargs):
            raise AssertionError("text output should use cached descendant counts")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        mocked_print = mocker.patch("builtins.print")

        svc._print_text_output(
            method="fast",
            trait_name="trait",
            n_tips=4,
            sigma2=1.0,
            log_likelihood=-1.0,
            node_estimates=node_estimates,
            node_cis={},
            node_labels=node_labels,
            tree=tree,
        )

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "N1 (root)" in all_output

    def test_text_output_batches_continuous_report(self, default_args, mocker):
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        node_labels = svc._label_internal_nodes(tree)
        node_estimates = {
            label: index + 1.25
            for index, label in enumerate(node_labels.values())
        }
        mocked_print = mocker.patch("builtins.print")

        svc._print_text_output(
            method="fast",
            trait_name="trait",
            n_tips=4,
            sigma2=1.0,
            log_likelihood=-2.5,
            node_estimates=node_estimates,
            node_cis={},
            node_labels=node_labels,
            tree=tree,
        )

        mocked_print.assert_called_once_with(
            "Ancestral State Reconstruction\n"
            "\n"
            "Method: fast (Felsenstein's contrasts)\n"
            "Trait: trait\n"
            "Number of tips: 4\n"
            "\n"
            "Log-likelihood: -2.5000\n"
            "Sigma-squared (BM rate): 1.000000\n"
            "\n"
            "Ancestral estimates:\n"
            "  Node         Descendants    Estimate\n"
            "  N1 (root)              4      1.2500\n"
            "  N2                     2      2.2500\n"
            "  N3                     2      3.2500"
        )

    def test_format_result_uses_cached_descendant_names(
        self, default_args, monkeypatch
    ):
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((B:1,A:1):1,(D:1,C:1):1);"), "newick")
        node_labels = svc._label_internal_nodes(tree)
        node_estimates = {label: 1.0 for label in node_labels.values()}

        def fail_traversal(*_args, **_kwargs):
            raise AssertionError("result formatting should use cached descendants")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        result = svc._format_result(
            method="fast",
            trait_name="trait",
            n_tips=4,
            sigma2=1.0,
            log_likelihood=-1.0,
            node_estimates=node_estimates,
            node_cis={},
            node_labels=node_labels,
            tree=tree,
            trait_values={"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0},
        )

        assert result["ancestral_estimates"]["N1"]["descendants"] == [
            "A",
            "B",
            "C",
            "D",
        ]

    @patch("builtins.print")
    def test_multi_trait_run(self, mocked_print, multi_trait_args):
        svc = AncestralReconstruction(multi_trait_args)
        svc.run()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "body_mass" in all_output

    @patch("builtins.print")
    def test_ml_run(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            trait=None,
            method="ml",
            ci=True,
            plot=None,
            json=False,
        )
        svc = AncestralReconstruction(args)
        svc.run()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "ml" in all_output
        assert "95% CI" in all_output


class TestPlot:
    def test_value_range_scans_values_once(self):
        class SinglePassValues:
            def __init__(self, values):
                self.values = values
                self.iterations = 0

            def __iter__(self):
                self.iterations += 1
                if self.iterations > 1:
                    raise AssertionError("values should be scanned once")
                return iter(self.values)

        values = SinglePassValues([3.0, 1.0, 4.0, 2.0])

        assert _value_range(values) == (1.0, 4.0)
        assert values.iterations == 1
        assert _value_range([]) is None
        assert _value_range([2.0, 2.0]) == (2.0, 2.0)

    @pytest.mark.parametrize("circular", [False, True])
    def test_contmap_plot_uses_direct_tree_traversal(
        self, default_args, monkeypatch, tmp_path, circular
    ):
        default_args.circular = circular
        default_args.plot_ci = True
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        node_labels = svc._label_internal_nodes(tree)
        node_estimates = {
            label: float(index)
            for index, label in enumerate(node_labels.values(), start=1)
        }
        node_cis = {
            label: (value - 0.1, value + 0.1)
            for label, value in node_estimates.items()
        }
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}

        def fail_traversal(*args, **kwargs):
            raise AssertionError("plot setup should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        output_path = tmp_path / f"asr_contmap_direct_{circular}.png"
        svc._plot_contmap(
            tree, node_estimates, node_labels, trait_values,
            "trait", str(output_path), node_cis=node_cis,
        )

        assert output_path.exists()

    def test_iter_preorder_preserves_order_without_reversed(self, default_args):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("_iter_preorder should push children directly")

        svc = AncestralReconstruction(default_args)
        root = Clade(name="root")
        left = Clade(name="left")
        middle = Clade(name="middle")
        right = Clade(name="right")
        left.clades = NoReversedList([Clade(name="left_a"), Clade(name="left_b")])
        middle.clades = NoReversedList([Clade(name="middle_a")])
        right.clades = NoReversedList()
        root.clades = NoReversedList([left, middle, right])

        order = [clade.name for clade in svc._iter_preorder(root)]

        assert order == [
            "root",
            "left",
            "left_a",
            "left_b",
            "middle",
            "middle_a",
            "right",
        ]

    def test_contmap_plot_reuses_clade_lists_for_layout_helpers(
        self, default_args, monkeypatch, tmp_path
    ):
        default_args.circular = True
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        node_labels = svc._label_internal_nodes(tree)
        node_estimates = {
            label: float(index)
            for index, label in enumerate(node_labels.values(), start=1)
        }
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        layout_calls = {}
        original_node_positions = ancestral_module.compute_node_positions
        original_circular_coords = ancestral_module.compute_circular_coords

        def spy_node_positions(*args, **kwargs):
            layout_calls["node_preorder"] = kwargs.get("preorder_clades")
            return original_node_positions(*args, **kwargs)

        def spy_circular_coords(*args, **kwargs):
            layout_calls["circular_preorder"] = kwargs.get("preorder_clades")
            layout_calls["circular_tips"] = kwargs.get("terminal_clades")
            return original_circular_coords(*args, **kwargs)

        monkeypatch.setattr(
            ancestral_module, "compute_node_positions", spy_node_positions
        )
        monkeypatch.setattr(
            ancestral_module, "compute_circular_coords", spy_circular_coords
        )

        output_path = tmp_path / "asr_contmap_layout_lists.png"
        svc._plot_contmap(
            tree, node_estimates, node_labels, trait_values,
            "trait", str(output_path), node_cis=None,
        )

        expected_preorder = list(svc._iter_preorder(tree.root))
        expected_tips = [clade for clade in expected_preorder if not clade.clades]
        assert layout_calls["circular_preorder"] is layout_calls["node_preorder"]
        assert len(layout_calls["node_preorder"]) == len(expected_preorder)
        assert len(layout_calls["circular_tips"]) == len(expected_tips)
        assert all(
            actual is expected
            for actual, expected in zip(
                layout_calls["node_preorder"], expected_preorder
            )
        )
        assert all(
            actual is expected
            for actual, expected in zip(layout_calls["circular_tips"], expected_tips)
        )
        assert output_path.exists()

    def test_rectangular_contmap_batches_gradient_branches(
        self, default_args, monkeypatch, tmp_path
    ):
        import matplotlib
        import matplotlib.axes

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection

        default_args.circular = False
        default_args.no_title = True
        default_args.ylabel_fontsize = 0
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        node_labels = svc._label_internal_nodes(tree)
        node_estimates = {
            label: float(index)
            for index, label in enumerate(node_labels.values(), start=1)
        }
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        line_collections = []
        original_add_collection = matplotlib.axes.Axes.add_collection

        def fail_plot(*args, **kwargs):
            raise AssertionError("contMap branches should use LineCollection")

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(matplotlib.axes.Axes, "add_collection", capture_collection)

        output_path = tmp_path / "asr_contmap_batched.png"
        try:
            svc._plot_contmap(
                tree, node_estimates, node_labels, trait_values,
                "trait", str(output_path), node_cis=None,
            )
        finally:
            plt.close("all")

        assert len(line_collections) >= 2
        assert output_path.exists()

    def test_rectangular_contmap_uses_collection_scalar_colormap(
        self, default_args, monkeypatch, tmp_path
    ):
        import matplotlib
        import matplotlib.axes

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection

        default_args.circular = False
        default_args.no_title = True
        default_args.ylabel_fontsize = 0
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        node_labels = svc._label_internal_nodes(tree)
        node_estimates = {
            label: float(index)
            for index, label in enumerate(node_labels.values(), start=1)
        }
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        collection_arrays = []
        original_add_collection = matplotlib.axes.Axes.add_collection

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                array = collection.get_array()
                if array is not None:
                    collection_arrays.append(array)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "add_collection", capture_collection)

        output_path = tmp_path / "asr_contmap_scalar_colormap.png"
        try:
            svc._plot_contmap(
                tree,
                node_estimates,
                node_labels,
                trait_values,
                "trait",
                str(output_path),
                node_cis=None,
            )
        finally:
            plt.close("all")

        array_lengths = {len(array) for array in collection_arrays}
        assert 300 in array_lengths
        assert 6 in array_lengths
        assert output_path.exists()

    def test_circular_contmap_batches_gradient_branches(
        self, default_args, monkeypatch, tmp_path
    ):
        import matplotlib
        import matplotlib.axes

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection

        default_args.circular = True
        default_args.no_title = True
        default_args.ylabel_fontsize = 0
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        node_labels = svc._label_internal_nodes(tree)
        node_estimates = {
            label: float(index)
            for index, label in enumerate(node_labels.values(), start=1)
        }
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        line_collections = []
        original_add_collection = matplotlib.axes.Axes.add_collection

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "add_collection", capture_collection)

        output_path = tmp_path / "asr_contmap_circular_batched.png"
        try:
            svc._plot_contmap(
                tree, node_estimates, node_labels, trait_values,
                "trait", str(output_path), node_cis=None,
            )
        finally:
            plt.close("all")

        segment_counts = {len(collection.get_segments()) for collection in line_collections}
        assert 180 in segment_counts
        assert 3 in segment_counts
        assert output_path.exists()

    @pytest.mark.parametrize("circular", [False, True])
    def test_contmap_batches_ci_overlays(
        self, default_args, monkeypatch, tmp_path, circular
    ):
        import matplotlib
        import matplotlib.axes

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection

        default_args.circular = circular
        default_args.no_title = True
        default_args.ylabel_fontsize = 0
        svc = AncestralReconstruction(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        node_labels = svc._label_internal_nodes(tree)
        node_estimates = {
            label: float(index)
            for index, label in enumerate(node_labels.values(), start=1)
        }
        node_cis = {
            label: (value - 0.2, value + 0.2)
            for label, value in node_estimates.items()
        }
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        ci_collections = []
        scatter_point_counts = []
        original_add_collection = matplotlib.axes.Axes.add_collection
        original_plot = matplotlib.axes.Axes.plot
        original_scatter = matplotlib.axes.Axes.scatter

        def fail_ci_plot(self, *args, **kwargs):
            if kwargs.get("color") == "black" and kwargs.get("zorder") == 8:
                raise AssertionError("contMap CI bars should use LineCollection")
            return original_plot(self, *args, **kwargs)

        def capture_collection(self, collection, *args, **kwargs):
            if (
                isinstance(collection, LineCollection)
                and collection.get_zorder() == 8
            ):
                ci_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        def capture_scatter(self, x, y, *args, **kwargs):
            if kwargs.get("c") == "black" and kwargs.get("zorder") == 9:
                try:
                    scatter_point_counts.append(len(x))
                except TypeError:
                    scatter_point_counts.append(1)
            return original_scatter(self, x, y, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_ci_plot)
        monkeypatch.setattr(matplotlib.axes.Axes, "add_collection", capture_collection)
        monkeypatch.setattr(matplotlib.axes.Axes, "scatter", capture_scatter)

        output_path = tmp_path / f"asr_contmap_ci_batched_{circular}.png"
        try:
            svc._plot_contmap(
                tree, node_estimates, node_labels, trait_values,
                "trait", str(output_path), node_cis=node_cis,
            )
        finally:
            plt.close("all")

        ci_segment_counts = sorted(len(collection.get_segments()) for collection in ci_collections)
        if circular:
            assert ci_segment_counts == [len(node_cis)]
        else:
            assert ci_segment_counts == [len(node_cis), 2 * len(node_cis)]
        assert scatter_point_counts == [len(node_cis)]
        assert output_path.exists()

    def test_rectangular_discrete_plot_batches_base_branches(
        self, discrete_args, monkeypatch, tmp_path
    ):
        import matplotlib
        import matplotlib.axes

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection

        discrete_args.circular = False
        discrete_args.no_title = True
        discrete_args.ylabel_fontsize = 0
        svc = AncestralReconstruction(discrete_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clades = list(svc._iter_preorder(tree.root))
        states = ["x", "y"]
        node_posteriors = {
            id(clade): np.array([0.4, 0.6])
            for clade in clades
            if clade.clades
        }
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "y"}
        line_collections = []
        original_add_collection = matplotlib.axes.Axes.add_collection

        def fail_plot(*args, **kwargs):
            raise AssertionError("discrete ASR branches should use LineCollection")

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(matplotlib.axes.Axes, "add_collection", capture_collection)

        output_path = tmp_path / "asr_discrete_base_batched.png"
        try:
            svc._plot_discrete_asr(
                tree, node_posteriors, {}, states, tip_states, str(output_path)
            )
        finally:
            plt.close("all")

        assert len(line_collections) >= 2
        assert output_path.exists()

    @pytest.mark.parametrize("circular", [False, True])
    def test_discrete_plot_batches_node_pies(
        self, discrete_args, monkeypatch, tmp_path, circular
    ):
        import matplotlib
        import matplotlib.axes

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import PatchCollection

        discrete_args.circular = circular
        discrete_args.no_title = True
        discrete_args.ylabel_fontsize = 0
        svc = AncestralReconstruction(discrete_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clades = list(svc._iter_preorder(tree.root))
        states = ["x", "y", "z"]
        posterior = np.array([0.2, 0.3, 0.5])
        node_posteriors = {
            id(clade): posterior
            for clade in clades
            if clade.clades
        }
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "z"}
        pie_collections = []
        original_add_collection = matplotlib.axes.Axes.add_collection

        def fail_add_patch(*args, **kwargs):
            raise AssertionError("discrete ASR node pies should use PatchCollection")

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, PatchCollection):
                pie_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "add_patch", fail_add_patch)
        monkeypatch.setattr(matplotlib.axes.Axes, "add_collection", capture_collection)

        output_path = tmp_path / f"asr_discrete_pies_{circular}.png"
        try:
            svc._plot_discrete_asr(
                tree, node_posteriors, {}, states, tip_states, str(output_path)
            )
        finally:
            plt.close("all")

        assert len(pie_collections) == 1
        assert len(pie_collections[0].get_paths()) == len(node_posteriors) * len(states)
        assert output_path.exists()

    @pytest.mark.parametrize("circular", [False, True])
    def test_discrete_plot_batches_color_clade_overlay(
        self, discrete_args, monkeypatch, tmp_path, circular
    ):
        import matplotlib
        import matplotlib.axes

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection

        color_file = tmp_path / "colors.tsv"
        color_file.write_text("A,B\tclade\t#ff0000\tAB\n")
        discrete_args.circular = circular
        discrete_args.color_file = str(color_file)
        discrete_args.no_title = True
        discrete_args.ylabel_fontsize = 0
        svc = AncestralReconstruction(discrete_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clades = list(svc._iter_preorder(tree.root))
        states = ["x", "y"]
        node_posteriors = {
            id(clade): np.array([0.4, 0.6])
            for clade in clades
            if clade.clades
        }
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "y"}
        overlay_collections = []
        original_add_collection = matplotlib.axes.Axes.add_collection

        def fail_plot(*args, **kwargs):
            raise AssertionError("clade overlay branches should use LineCollection")

        def capture_collection(self, collection, *args, **kwargs):
            if (
                isinstance(collection, LineCollection)
                and collection.get_zorder() == 2
            ):
                overlay_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(matplotlib.axes.Axes, "add_collection", capture_collection)

        output_path = tmp_path / f"asr_discrete_color_{circular}.png"
        try:
            svc._plot_discrete_asr(
                tree, node_posteriors, {}, states, tip_states, str(output_path)
            )
        finally:
            plt.close("all")

        assert overlay_collections
        assert output_path.exists()

    @patch("builtins.print")
    def test_plot_created(self, mocked_print):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            plot_path = f.name

        try:
            args = Namespace(
                tree=TREE_SIMPLE,
                trait_data=TRAITS_FILE,
                trait=None,
                method="fast",
                ci=False,
                plot=plot_path,
                json=False,
            )
            svc = AncestralReconstruction(args)
            svc.run()
            assert os.path.exists(plot_path)
            assert os.path.getsize(plot_path) > 0
        finally:
            if os.path.exists(plot_path):
                os.unlink(plot_path)


# ------------------------------------------------------------------
# Discrete ASR tests
# ------------------------------------------------------------------

@pytest.fixture
def discrete_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=DISCRETE_TRAITS_FILE,
        trait="diet",
        method="fast",
        ci=False,
        plot=None,
        json=False,
        type="discrete",
        model="ER",
    )


@pytest.fixture
def discrete_ard_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=DISCRETE_TRAITS_FILE,
        trait="diet",
        method="fast",
        ci=False,
        plot=None,
        json=False,
        type="discrete",
        model="ARD",
    )


class TestProcessArgsDiscrete:
    def test_discrete_type_and_model(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            trait=None,
            method="fast",
            ci=False,
            plot=None,
            json=False,
            type="discrete",
            model="SYM",
        )
        svc = AncestralReconstruction(args)
        assert svc.trait_type == "discrete"
        assert svc.model == "SYM"

    def test_continuous_defaults_unchanged(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            trait=None,
            method="fast",
            ci=False,
            plot=None,
            json=False,
        )
        svc = AncestralReconstruction(args)
        assert svc.trait_type == "continuous"
        assert svc.model == "ER"


class TestDiscreteTraitParsing:
    def test_multi_col_parsing(self, discrete_args):
        svc = AncestralReconstruction(discrete_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        traits = svc._parse_discrete_trait_data_multi(
            DISCRETE_TRAITS_FILE, tree_tips, "diet"
        )
        assert len(traits) == 8
        assert traits["raccoon"] == "carnivore"
        assert traits["cat"] == "omnivore"

    def test_single_col_parsing(self, discrete_args):
        svc = AncestralReconstruction(discrete_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        # Create a temp 2-col file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write("raccoon\tcarnivore\n")
            f.write("bear\tcarnivore\n")
            f.write("sea_lion\tcarnivore\n")
            f.write("seal\therbivore\n")
            f.write("monkey\therbivore\n")
            f.write("cat\tomnivore\n")
            f.write("weasel\tomnivore\n")
            f.write("dog\tomnivore\n")
            tmp_path = f.name

        try:
            traits = svc._parse_discrete_trait_data_single(tmp_path, tree_tips)
            assert len(traits) == 8
            assert traits["bear"] == "carnivore"
        finally:
            os.unlink(tmp_path)

    def test_single_col_skips_comments_and_blanks(self, discrete_args, tmp_path):
        svc = AncestralReconstruction(discrete_args)
        trait_file = tmp_path / "states.tsv"
        trait_file.write_text(
            "   # ignored\n"
            "\n"
            "raccoon\tcarnivore\n"
            "bear\tomnivore\n"
            "weasel\tcarnivore\n"
        )

        traits = svc._parse_discrete_trait_data_single(
            str(trait_file), ["raccoon", "bear", "weasel"]
        )

        assert traits == {
            "raccoon": "carnivore",
            "bear": "omnivore",
            "weasel": "carnivore",
        }

    def test_single_col_all_shared_emits_no_warnings(
        self, discrete_args, tmp_path, capsys
    ):
        svc = AncestralReconstruction(discrete_args)
        trait_file = tmp_path / "states.tsv"
        trait_file.write_text(
            "raccoon\tomnivore\nbear\tomnivore\nweasel\tcarnivore\n"
        )

        traits = svc._parse_discrete_trait_data_single(
            str(trait_file), ["raccoon", "bear", "weasel"]
        )

        assert traits == {
            "raccoon": "omnivore",
            "bear": "omnivore",
            "weasel": "carnivore",
        }
        assert capsys.readouterr().err == ""

    def test_single_col_wrong_column_count_error(self, discrete_args, tmp_path):
        svc = AncestralReconstruction(discrete_args)
        trait_file = tmp_path / "states.tsv"
        trait_file.write_text(
            "raccoon\tcarnivore\n"
            "bear\tomnivore\textra\n"
            "weasel\tcarnivore\n"
        )

        with pytest.raises(PhykitUserError):
            svc._parse_discrete_trait_data_single(
                str(trait_file), ["raccoon", "bear", "weasel"]
            )

    def test_missing_column_error(self, discrete_args):
        svc = AncestralReconstruction(discrete_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        with pytest.raises(PhykitUserError):
            svc._parse_discrete_trait_data_multi(
                DISCRETE_TRAITS_FILE, tree_tips, "nonexistent"
            )

    def test_multi_col_skips_comments_and_blanks(self, discrete_args, tmp_path):
        svc = AncestralReconstruction(discrete_args)
        trait_file = tmp_path / "states.tsv"
        trait_file.write_text(
            "# ignored before header\n"
            "\n"
            "taxon\tdiet\tactivity\n"
            "\t# ignored before data\n"
            "\n"
            "raccoon\tomnivore\tnocturnal\n"
            "bear\tomnivore\tdiurnal\n"
            "weasel\tcarnivore\tnocturnal\n"
        )

        traits = svc._parse_discrete_trait_data_multi(
            str(trait_file),
            ["raccoon", "bear", "weasel"],
            "activity",
        )

        assert traits == {
            "raccoon": "nocturnal",
            "bear": "diurnal",
            "weasel": "nocturnal",
        }

    def test_multi_col_all_shared_emits_no_warnings(
        self, discrete_args, tmp_path, capsys
    ):
        svc = AncestralReconstruction(discrete_args)
        trait_file = tmp_path / "states.tsv"
        trait_file.write_text(
            "taxon\tdiet\tactivity\n"
            "raccoon\tomnivore\tnocturnal\n"
            "bear\tomnivore\tdiurnal\n"
            "weasel\tcarnivore\tnocturnal\n"
        )

        traits = svc._parse_discrete_trait_data_multi(
            str(trait_file),
            ["raccoon", "bear", "weasel"],
            "activity",
        )

        assert traits == {
            "raccoon": "nocturnal",
            "bear": "diurnal",
            "weasel": "nocturnal",
        }
        assert capsys.readouterr().err == ""

    def test_multi_col_wrong_column_count_error(self, discrete_args, tmp_path):
        svc = AncestralReconstruction(discrete_args)
        trait_file = tmp_path / "states.tsv"
        trait_file.write_text(
            "taxon\tdiet\tactivity\n"
            "raccoon\tomnivore\tnocturnal\n"
            "bear\tomnivore\n"
            "weasel\tcarnivore\tnocturnal\n"
        )

        with pytest.raises(PhykitUserError):
            svc._parse_discrete_trait_data_multi(
                str(trait_file),
                ["raccoon", "bear", "weasel"],
                "activity",
            )

    def test_file_not_found(self, discrete_args):
        svc = AncestralReconstruction(discrete_args)
        with pytest.raises(PhykitUserError):
            svc._parse_discrete_trait_data_single("no_such_file.tsv", [])


class TestMkPrimitives:
    def test_q_matrix_er(self, discrete_args):
        svc = AncestralReconstruction(discrete_args)
        params = np.array([0.5])
        Q = svc._build_q_matrix(params, 3, "ER")
        assert Q.shape == (3, 3)
        # Rows should sum to zero
        for i in range(3):
            assert pytest.approx(np.sum(Q[i, :]), abs=1e-10) == 0.0
        # Off-diagonals should all be 0.5
        for i in range(3):
            for j in range(3):
                if i != j:
                    assert Q[i, j] == pytest.approx(0.5)

    def test_q_matrix_ard(self, discrete_args):
        svc = AncestralReconstruction(discrete_args)
        k = 3
        n_params = k * (k - 1)  # 6
        params = np.arange(1, n_params + 1, dtype=float) * 0.1
        Q = svc._build_q_matrix(params, k, "ARD")
        assert Q.shape == (k, k)
        for i in range(k):
            assert pytest.approx(np.sum(Q[i, :]), abs=1e-10) == 0.0

    def test_felsenstein_pruning_loglik_finite(self, discrete_args):
        svc = AncestralReconstruction(discrete_args)
        tree = svc.read_tree_file()
        import copy
        tree_copy = copy.deepcopy(tree)
        tree_tips = svc.get_tip_names_from_tree(tree_copy)
        tip_states = svc._parse_discrete_trait_data_multi(
            DISCRETE_TRAITS_FILE, tree_tips, "diet"
        )
        states = sorted(set(tip_states.values()))
        k = len(states)
        pi = np.ones(k) / k
        Q = svc._build_q_matrix(np.array([0.01]), k, "ER")
        _, loglik = svc._felsenstein_pruning(tree_copy, tip_states, Q, pi, states)
        assert np.isfinite(loglik)
        assert loglik < 0

    def test_fit_er_loglik(self, discrete_args):
        svc = AncestralReconstruction(discrete_args)
        tree = svc.read_tree_file()
        import copy
        tree_copy = copy.deepcopy(tree)
        tree_tips = svc.get_tip_names_from_tree(tree_copy)
        tip_states = svc._parse_discrete_trait_data_multi(
            DISCRETE_TRAITS_FILE, tree_tips, "diet"
        )
        states = sorted(set(tip_states.values()))
        Q, loglik = svc._fit_q_matrix(tree_copy, tip_states, states, "ER")
        assert np.isfinite(loglik)
        assert loglik == pytest.approx(-8.789, abs=0.5)


class TestDiscreteMarginalPosteriors:
    def _get_posteriors(self, args_fixture):
        import copy
        svc = AncestralReconstruction(args_fixture)
        tree = svc.read_tree_file()
        tree_copy = copy.deepcopy(tree)
        tree_tips = svc.get_tip_names_from_tree(tree_copy)
        tip_states = svc._parse_discrete_trait_data_multi(
            DISCRETE_TRAITS_FILE, tree_tips, "diet"
        )
        states = sorted(set(tip_states.values()))
        Q, loglik = svc._fit_q_matrix(tree_copy, tip_states, states, args_fixture.model)
        posteriors = svc._discrete_marginal_posteriors(
            tree_copy, tip_states, Q, states
        )
        node_labels = svc._label_internal_nodes(tree_copy)
        return posteriors, states, tree_copy, node_labels

    def test_posteriors_sum_to_one(self, discrete_args):
        posteriors, states, tree, _ = self._get_posteriors(discrete_args)
        for node_id, post in posteriors.items():
            assert pytest.approx(np.sum(post), abs=1e-6) == 1.0

    def test_posteriors_non_negative(self, discrete_args):
        posteriors, _, _, _ = self._get_posteriors(discrete_args)
        for node_id, post in posteriors.items():
            assert np.all(post >= 0.0)

    def test_posteriors_correct_shape(self, discrete_args):
        posteriors, states, _, _ = self._get_posteriors(discrete_args)
        k = len(states)
        for node_id, post in posteriors.items():
            assert post.shape == (k,)

    def test_all_internal_nodes_covered(self, discrete_args):
        posteriors, states, tree, node_labels = self._get_posteriors(discrete_args)
        for clade in tree.find_clades(order="preorder"):
            if not clade.is_terminal():
                assert id(clade) in posteriors

    def test_ard_model_valid(self, discrete_ard_args):
        posteriors, states, _, _ = self._get_posteriors(discrete_ard_args)
        for node_id, post in posteriors.items():
            assert pytest.approx(np.sum(post), abs=1e-6) == 1.0
            assert np.all(post >= 0.0)

    def test_map_state_valid(self, discrete_args):
        posteriors, states, _, _ = self._get_posteriors(discrete_args)
        for node_id, post in posteriors.items():
            map_idx = int(np.argmax(post))
            assert 0 <= map_idx < len(states)

    def test_root_posterior_reasonable(self, discrete_args):
        posteriors, states, tree, _ = self._get_posteriors(discrete_args)
        root_post = posteriors[id(tree.root)]
        # Root posterior should have some spread (not all in one state)
        assert np.max(root_post) < 1.0
        assert np.min(root_post) >= 0.0

    def test_discrete_marginal_posteriors_reuse_transition_cache(
        self, discrete_args, monkeypatch
    ):
        svc = AncestralReconstruction(discrete_args)
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "y"}
        states = ["x", "y"]
        Q = np.array([[-1.0, 1.0], [1.0, -1.0]], dtype=float)
        expected = svc._discrete_marginal_posteriors(
            tree,
            tip_states,
            Q,
            states,
        )
        monkeypatch.setattr(
            ancestral_module.np,
            "sum",
            lambda *args, **kwargs: pytest.fail(
                "discrete posterior normalization should use ndarray.sum"
            ),
        )
        calls = []
        original_matrix_exp = svc._matrix_exp

        def wrapped_matrix_exp(Q_arg, t_arg):
            calls.append(t_arg)
            return original_matrix_exp(Q_arg, t_arg)

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("discrete posteriors should use direct traversal")

        monkeypatch.setattr(svc, "_matrix_exp", wrapped_matrix_exp)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        observed = svc._discrete_marginal_posteriors(
            tree,
            tip_states,
            Q,
            states,
        )

        assert calls == [1.0]
        assert observed.keys() == expected.keys()
        for node_id in expected:
            np.testing.assert_allclose(observed[node_id], expected[node_id])


class TestDiscreteRun:
    @patch("builtins.print")
    def test_text_output(self, mocked_print, discrete_args):
        svc = AncestralReconstruction(discrete_args)
        svc.run()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Discrete" in all_output
        assert "Mk" in all_output
        assert "Rate matrix" in all_output
        assert "N1 (root)" in all_output

    @patch("builtins.print")
    def test_json_output_schema(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=DISCRETE_TRAITS_FILE,
            trait="diet",
            method="fast",
            ci=False,
            plot=None,
            json=True,
            type="discrete",
            model="ER",
        )
        svc = AncestralReconstruction(args)
        svc.run()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["method"] == "discrete"
        assert payload["model"] == "ER"
        assert "states" in payload
        assert "q_matrix" in payload
        assert "ancestral_states" in payload
        assert "tip_states" in payload
        assert "N1" in payload["ancestral_states"]
        n1 = payload["ancestral_states"]["N1"]
        assert n1["is_root"] is True
        assert "map_state" in n1
        assert "posteriors" in n1

    def test_run_uses_fast_tip_name_helper_for_discrete_prune_setup(
        self, discrete_args, mocker
    ):
        svc = AncestralReconstruction(discrete_args)
        spy = mocker.spy(svc, "get_tip_names_from_tree")
        mocker.patch("builtins.print")

        svc.run()

        assert spy.call_count == 1

    def test_run_uses_unmodified_tree_read_for_discrete(self, discrete_args, mocker):
        svc = AncestralReconstruction(discrete_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tip_names = ["A", "B", "C", "D"]

        read_unmodified = mocker.patch.object(
            svc, "read_tree_file_unmodified", return_value=tree
        )
        mocker.patch.object(
            svc,
            "read_tree_file",
            side_effect=AssertionError("run should use read_tree_file_unmodified"),
        )
        validate = mocker.patch.object(svc, "validate_tree")
        run_discrete = mocker.patch.object(svc, "_run_discrete")

        svc.run()

        read_unmodified.assert_called_once_with()
        validate.assert_called_once()
        run_discrete.assert_called_once_with(tree, tip_names)

    def test_discrete_run_skips_copy_for_all_shared_tree(
        self, discrete_args, mocker
    ):
        svc = AncestralReconstruction(discrete_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tip_names = ["A", "B", "C", "D"]
        states_by_tip = {"A": "x", "B": "x", "C": "y", "D": "y"}

        mocker.patch.object(
            svc, "_parse_discrete_trait_data_multi", return_value=states_by_tip
        )
        fast_copy = mocker.patch.object(
            svc,
            "_fast_copy",
            side_effect=AssertionError("all-shared analysis should not copy tree"),
        )
        mocker.patch.object(svc, "_label_internal_nodes", return_value={"N1": "N1"})
        fit_q = mocker.patch.object(
            svc, "_fit_q_matrix", return_value=(np.eye(2), -1.0)
        )
        mocker.patch.object(
            svc,
            "_discrete_marginal_posteriors",
            return_value={"N1": np.array([0.5, 0.5])},
        )
        mocker.patch.object(
            svc, "_format_discrete_result", return_value={"ok": True}
        )
        mocker.patch.object(svc, "_print_discrete_text_output")

        svc._run_discrete(tree, tip_names)

        fast_copy.assert_not_called()
        assert fit_q.call_args.args[0] is tree

    def test_discrete_run_copies_before_pruning_missing_tree_tips(
        self, discrete_args, monkeypatch, mocker
    ):
        class OrderedTipStates(dict):
            def __contains__(self, key):
                raise AssertionError("ordered prune path should not scan membership")

        svc = AncestralReconstruction(discrete_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree_copy = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        pruned_tree = Phylo.read(StringIO("((A:1,B:1):1,C:1);"), "newick")
        tip_names = ["A", "B", "C", "D"]
        states_by_tip = OrderedTipStates({"A": "x", "B": "x", "C": "y"})

        monkeypatch.setattr(ancestral_module.Tree, "_ORDERED_MAPPING_PRUNE_MIN_SIZE", 0)
        mocker.patch.object(
            svc, "_parse_discrete_trait_data_multi", return_value=states_by_tip
        )
        fast_copy = mocker.patch.object(svc, "_fast_copy", return_value=tree_copy)
        prune = mocker.patch.object(
            svc, "prune_tree_using_taxa_list", return_value=pruned_tree
        )
        mocker.patch.object(svc, "_label_internal_nodes", return_value={"N1": "N1"})
        fit_q = mocker.patch.object(
            svc, "_fit_q_matrix", return_value=(np.eye(2), -1.0)
        )
        mocker.patch.object(
            svc,
            "_discrete_marginal_posteriors",
            return_value={"N1": np.array([0.5, 0.5])},
        )
        mocker.patch.object(
            svc, "_format_discrete_result", return_value={"ok": True}
        )
        mocker.patch.object(svc, "_print_discrete_text_output")

        svc._run_discrete(tree, tip_names)

        fast_copy.assert_called_once_with(tree)
        prune.assert_called_once_with(tree_copy, ["D"])
        assert fit_q.call_args.args[0] is pruned_tree

    def test_discrete_text_output_uses_cached_descendant_counts(
        self, discrete_args, monkeypatch, mocker
    ):
        svc = AncestralReconstruction(discrete_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        node_labels = svc._label_internal_nodes(tree)
        states = ["x", "y"]
        node_posteriors = {
            node_id: np.array([0.25, 0.75])
            for node_id in node_labels
        }
        Q = np.array([[-1.0, 1.0], [1.0, -1.0]], dtype=float)

        def fail_traversal(*_args, **_kwargs):
            raise AssertionError(
                "discrete text should use cached descendant counts"
            )

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        mocked_print = mocker.patch("builtins.print")

        svc._print_discrete_text_output(
            model="ER",
            trait_name="trait",
            n_tips=4,
            log_likelihood=-1.0,
            states=states,
            Q=Q,
            node_posteriors=node_posteriors,
            node_labels=node_labels,
            tree=tree,
        )

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        mocked_print.assert_called_once()
        assert "Ancestral State Reconstruction (Discrete)" in all_output
        assert "Rate matrix (Q):" in all_output
        assert "Ancestral state posteriors:" in all_output
        assert "N1 (root)" in all_output
        lines = mocked_print.call_args.args[0].splitlines()
        assert lines[11:14] == [
            "                         x           y",
            "             x   -1.000000    1.000000",
            "             y    1.000000   -1.000000",
        ]
        assert lines[16:20] == [
            "  Node          Desc       MAP         x         y",
            "  N1 (root)        4         y    0.2500    0.7500",
            "  N2               2         y    0.2500    0.7500",
            "  N3               2         y    0.2500    0.7500",
        ]

    def test_format_discrete_result_uses_cached_descendant_names(
        self, discrete_args, monkeypatch
    ):
        svc = AncestralReconstruction(discrete_args)
        tree = Phylo.read(StringIO("((B:1,A:1):1,(D:1,C:1):1);"), "newick")
        node_labels = svc._label_internal_nodes(tree)
        states = ["x", "y"]
        node_posteriors = {
            node_id: np.array([0.25, 0.75])
            for node_id in node_labels
        }
        Q = np.array([[-1.0, 1.0], [1.0, -1.0]], dtype=float)

        def fail_traversal(*_args, **_kwargs):
            raise AssertionError(
                "discrete result formatting should use cached descendants"
            )

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        result = svc._format_discrete_result(
            model="ER",
            trait_name="trait",
            n_tips=4,
            log_likelihood=-1.0,
            states=states,
            Q=Q,
            node_posteriors=node_posteriors,
            node_labels=node_labels,
            tree=tree,
            tip_states={"A": "x", "B": "x", "C": "y", "D": "y"},
        )

        assert result["q_matrix"] == {
            "x": {"x": -1.0, "y": 1.0},
            "y": {"x": 1.0, "y": -1.0},
        }
        assert result["ancestral_states"]["N1"]["descendants"] == [
            "A",
            "B",
            "C",
            "D",
        ]

    @patch("builtins.print")
    def test_continuous_regression(self, mocked_print):
        """Ensure continuous mode still works after refactor."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            trait=None,
            method="fast",
            ci=False,
            plot=None,
            json=False,
        )
        svc = AncestralReconstruction(args)
        svc.run()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Ancestral State Reconstruction" in all_output
        assert "Sigma-squared" in all_output

    @patch("builtins.print")
    def test_ard_smoke(self, mocked_print, discrete_ard_args):
        svc = AncestralReconstruction(discrete_ard_args)
        svc.run()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "ARD" in all_output


class TestDiscretePlot:
    @pytest.mark.parametrize("circular", [False, True])
    def test_discrete_plot_uses_direct_tree_traversal(
        self, discrete_args, monkeypatch, tmp_path, circular
    ):
        discrete_args.circular = circular
        svc = AncestralReconstruction(discrete_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clades = list(svc._iter_preorder(tree.root))
        states = ["x", "y"]
        node_posteriors = {
            id(clade): np.array([0.4, 0.6])
            for clade in clades
            if clade.clades
        }
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "y"}

        def fail_traversal(*args, **kwargs):
            raise AssertionError("plot setup should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        output_path = tmp_path / f"asr_discrete_direct_{circular}.png"
        svc._plot_discrete_asr(
            tree, node_posteriors, {}, states, tip_states, str(output_path)
        )

        assert output_path.exists()

    def test_discrete_plot_reuses_clade_lists_for_layout_helpers(
        self, discrete_args, monkeypatch, tmp_path
    ):
        discrete_args.circular = True
        svc = AncestralReconstruction(discrete_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clades = list(svc._iter_preorder(tree.root))
        states = ["x", "y"]
        node_posteriors = {
            id(clade): np.array([0.4, 0.6])
            for clade in clades
            if clade.clades
        }
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "y"}
        layout_calls = {}
        original_node_positions = ancestral_module.compute_node_positions
        original_circular_coords = ancestral_module.compute_circular_coords

        def spy_node_positions(*args, **kwargs):
            layout_calls["node_preorder"] = kwargs.get("preorder_clades")
            return original_node_positions(*args, **kwargs)

        def spy_circular_coords(*args, **kwargs):
            layout_calls["circular_preorder"] = kwargs.get("preorder_clades")
            layout_calls["circular_tips"] = kwargs.get("terminal_clades")
            return original_circular_coords(*args, **kwargs)

        monkeypatch.setattr(
            ancestral_module, "compute_node_positions", spy_node_positions
        )
        monkeypatch.setattr(
            ancestral_module, "compute_circular_coords", spy_circular_coords
        )

        output_path = tmp_path / "asr_discrete_layout_lists.png"
        svc._plot_discrete_asr(
            tree, node_posteriors, {}, states, tip_states, str(output_path)
        )

        expected_preorder = list(svc._iter_preorder(tree.root))
        expected_tips = [clade for clade in expected_preorder if not clade.clades]
        assert layout_calls["circular_preorder"] is layout_calls["node_preorder"]
        assert len(layout_calls["node_preorder"]) == len(expected_preorder)
        assert len(layout_calls["circular_tips"]) == len(expected_tips)
        assert all(
            actual is expected
            for actual, expected in zip(
                layout_calls["node_preorder"], expected_preorder
            )
        )
        assert all(
            actual is expected
            for actual, expected in zip(layout_calls["circular_tips"], expected_tips)
        )
        assert output_path.exists()

    def test_discrete_plot_reuses_cached_rendered_bytes(
        self, discrete_args, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        discrete_args.circular = False
        discrete_args.no_title = True
        discrete_args.ylabel_fontsize = 0
        states = ["x", "y"]
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "y"}
        tree_text = "((A:1,B:1):1,(C:1,D:1):1);"

        def make_plot_inputs():
            svc = AncestralReconstruction(discrete_args)
            tree = Phylo.read(StringIO(tree_text), "newick")
            clades = list(svc._iter_preorder(tree.root))
            posteriors = {
                id(clade): np.array([0.4, 0.6])
                for clade in clades
                if clade.clades
            }
            return svc, tree, posteriors

        svc, tree, node_posteriors = make_plot_inputs()
        first_path = tmp_path / "asr_discrete_first.png"
        svc._plot_discrete_asr(
            tree, node_posteriors, {}, states, tip_states, str(first_path)
        )
        first_bytes = first_path.read_bytes()

        original_import = builtins.__import__

        def fail_matplotlib_import(name, *args, **kwargs):
            if name == "matplotlib" or name.startswith("matplotlib."):
                raise AssertionError("cached discrete ASR plot should skip matplotlib")
            return original_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", fail_matplotlib_import)
        svc, tree, node_posteriors = make_plot_inputs()
        second_path = tmp_path / "asr_discrete_second.png"
        svc._plot_discrete_asr(
            tree, node_posteriors, {}, states, tip_states, str(second_path)
        )

        assert second_path.read_bytes() == first_bytes

    @patch("builtins.print")
    def test_plot_created(self, mocked_print):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            plot_path = f.name

        try:
            args = Namespace(
                tree=TREE_SIMPLE,
                trait_data=DISCRETE_TRAITS_FILE,
                trait="diet",
                method="fast",
                ci=False,
                plot=plot_path,
                json=False,
                type="discrete",
                model="ER",
            )
            svc = AncestralReconstruction(args)
            svc.run()
            assert os.path.exists(plot_path)
            assert os.path.getsize(plot_path) > 0
        finally:
            if os.path.exists(plot_path):
                os.unlink(plot_path)
