import json
import sys
import os
import tempfile

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from phykit.services.tree.ancestral_reconstruction import AncestralReconstruction
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")
DISCRETE_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_discrete_traits.tsv")


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

    def test_missing_column_error(self, discrete_args):
        svc = AncestralReconstruction(discrete_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        with pytest.raises(PhykitUserError):
            svc._parse_discrete_trait_data_multi(
                DISCRETE_TRAITS_FILE, tree_tips, "nonexistent"
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
