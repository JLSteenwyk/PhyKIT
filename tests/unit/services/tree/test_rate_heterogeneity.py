import json
import sys
import os
import tempfile

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from phykit.services.tree.rate_heterogeneity import RateHeterogeneity
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")
REGIMES_FILE = str(SAMPLE_FILES / "tree_simple_regimes.tsv")


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        regime_data=REGIMES_FILE,
        nsim=0,
        seed=None,
        plot=None,
        json=False,
    )


@pytest.fixture
def json_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        regime_data=REGIMES_FILE,
        nsim=0,
        seed=None,
        plot=None,
        json=True,
    )


@pytest.fixture
def sim_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        regime_data=REGIMES_FILE,
        nsim=10,
        seed=42,
        plot=None,
        json=False,
    )


@pytest.fixture
def service(default_args):
    return RateHeterogeneity(default_args)


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            regime_data="r.tsv",
            nsim=0,
            seed=None,
            plot=None,
            json=False,
        )
        svc = RateHeterogeneity.__new__(RateHeterogeneity)
        parsed = svc.process_args(args)
        assert parsed["tree_file_path"] == "t.tre"
        assert parsed["trait_data_path"] == "d.tsv"
        assert parsed["regime_data_path"] == "r.tsv"
        assert parsed["n_sim"] == 0
        assert parsed["seed"] is None
        assert parsed["json_output"] is False

    def test_overrides(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            regime_data="r.tsv",
            nsim=100,
            seed=42,
            plot="out.png",
            json=True,
        )
        svc = RateHeterogeneity.__new__(RateHeterogeneity)
        parsed = svc.process_args(args)
        assert parsed["n_sim"] == 100
        assert parsed["seed"] == 42
        assert parsed["plot_output"] == "out.png"
        assert parsed["json_output"] is True


class TestRegimeParsing:
    def test_valid_file(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        regimes = service._parse_regime_file(REGIMES_FILE, tree_tips)
        assert len(regimes) == 8
        assert regimes["sea_lion"] == "aquatic"
        assert regimes["raccoon"] == "terrestrial"

    def test_comments_and_blanks(self, service, tmp_path):
        regime_file = tmp_path / "regimes.tsv"
        regime_file.write_text(
            "# comment\n\nraccoon\tA\nbear\tA\nsea_lion\tB\n"
            "seal\tB\nmonkey\tA\ncat\tA\nweasel\tA\ndog\tA\n"
        )
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        regimes = service._parse_regime_file(str(regime_file), tree_tips)
        assert len(regimes) == 8
        assert regimes["raccoon"] == "A"

    def test_missing_file(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        with pytest.raises(PhykitUserError):
            service._parse_regime_file("/nonexistent/file.tsv", tree_tips)


class TestFitchParsimony:
    def test_tip_states_preserved(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        parent_map = service._build_parent_map(tree)
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        # Tips should have their correct regime
        for tip in tree.get_terminals():
            if tip.name in regime_assignments:
                assert branch_regimes[id(tip)] == regime_assignments[tip.name]

    def test_correct_branch_assignments(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        parent_map = service._build_parent_map(tree)
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        # All non-root nodes should have a regime
        for clade in tree.find_clades(order="preorder"):
            if clade != tree.root:
                assert id(clade) in branch_regimes
                assert branch_regimes[id(clade)] in {"aquatic", "terrestrial"}


class TestPerRegimeVCV:
    def test_sum_equals_total(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )

        # Sum of per-regime VCVs should approximate total VCV
        C_sum = sum(per_regime_vcv.values())

        # Build total VCV using root-to-tip distances
        n = len(ordered_names)
        C_total = np.zeros((n, n))
        root_to_tip = {}
        for name in ordered_names:
            root_to_tip[name] = tree.distance(tree.root, name)
        for i in range(n):
            for j in range(i, n):
                if i == j:
                    C_total[i, j] = root_to_tip[ordered_names[i]]
                else:
                    d_ij = tree.distance(ordered_names[i], ordered_names[j])
                    shared_path = (
                        root_to_tip[ordered_names[i]]
                        + root_to_tip[ordered_names[j]]
                        - d_ij
                    ) / 2.0
                    C_total[i, j] = shared_path
                    C_total[j, i] = shared_path

        np.testing.assert_allclose(C_sum, C_total, rtol=1e-6)

    def test_individual_values_nonnegative(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )

        for r in regimes:
            assert np.all(per_regime_vcv[r] >= -1e-10)


class TestSingleRate:
    def test_sigma2_positive(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        C_total = sum(per_regime_vcv.values())

        sigma2, anc, ll = service._fit_single_rate(y, C_total)
        assert sigma2 > 0

    def test_log_likelihood_finite(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        C_total = sum(per_regime_vcv.values())

        sigma2, anc, ll = service._fit_single_rate(y, C_total)
        assert np.isfinite(ll)


class TestMultiRate:
    def test_log_likelihood_geq_single(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        C_total = sum(per_regime_vcv.values())

        _, _, ll_single = service._fit_single_rate(y, C_total)
        _, _, ll_multi = service._fit_multi_rate(y, per_regime_vcv, regimes)
        # Multi-rate model should fit at least as well
        assert ll_multi >= ll_single - 1e-6

    def test_per_regime_sigma2_positive(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )

        sigma2s, _, _ = service._fit_multi_rate(y, per_regime_vcv, regimes)
        assert np.all(sigma2s > 0)


class TestLRT:
    def test_chi2_p_in_range(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        C_total = sum(per_regime_vcv.values())

        _, _, ll_single = service._fit_single_rate(y, C_total)
        _, _, ll_multi = service._fit_multi_rate(y, per_regime_vcv, regimes)

        lrt_stat = max(2.0 * (ll_multi - ll_single), 0.0)
        from scipy.stats import chi2
        p = float(chi2.sf(lrt_stat, df=len(regimes) - 1))
        assert 0.0 <= p <= 1.0

    def test_lrt_stat_nonnegative(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        C_total = sum(per_regime_vcv.values())

        _, _, ll_single = service._fit_single_rate(y, C_total)
        _, _, ll_multi = service._fit_multi_rate(y, per_regime_vcv, regimes)

        lrt_stat = max(2.0 * (ll_multi - ll_single), 0.0)
        assert lrt_stat >= 0


class TestRValidation:
    """Validate against R 4.4.0 phytools::brownie.lite with stem=TRUE.

    R code:
        library(phytools)
        tree <- read.tree("tree_simple.tre")
        x <- setNames(c(1.04,2.39,2.30,1.88,0.60,0.56,-0.30,1.18),
                       c("raccoon","bear","sea_lion","seal","monkey","cat","weasel","dog"))
        aquatic_mrca <- findMRCA(tree, c("sea_lion", "seal"))
        tree2 <- paintSubTree(tree, node=aquatic_mrca, state="aquatic",
                              anc.state="terrestrial", stem=TRUE)
        brownie.lite(tree2, x)

    PhyKIT uses Fitch parsimony for regime assignment, which matches
    R's paintSubTree(stem=TRUE) for this tree topology.
    """

    def _get_results(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        C_total = sum(per_regime_vcv.values())

        sigma2_s, anc_s, ll_s = service._fit_single_rate(y, C_total)
        sigma2_m, anc_m, ll_m = service._fit_multi_rate(y, per_regime_vcv, regimes)

        return dict(
            sigma2_single=sigma2_s,
            anc_single=anc_s,
            ll_single=ll_s,
            sigma2_multi=dict(zip(regimes, sigma2_m)),
            anc_multi=anc_m,
            ll_multi=ll_m,
        )

    def test_single_rate_sigma2(self, service):
        r = self._get_results(service)
        # R: 0.0384065703382948
        np.testing.assert_allclose(r["sigma2_single"], 0.0384065703, rtol=1e-6)

    def test_single_rate_log_likelihood(self, service):
        r = self._get_results(service)
        # R: -11.5696774876938
        np.testing.assert_allclose(r["ll_single"], -11.5696774877, rtol=1e-6)

    def test_multi_rate_log_likelihood(self, service):
        r = self._get_results(service)
        # R: -11.204608463063
        np.testing.assert_allclose(r["ll_multi"], -11.2046, rtol=1e-3)


class TestRun:
    @patch("builtins.print")
    def test_text_output(self, mocked_print, default_args):
        svc = RateHeterogeneity(default_args)
        svc.run()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Rate Heterogeneity Test" in all_output
        assert "Single-rate model" in all_output
        assert "Multi-rate model" in all_output
        assert "Likelihood ratio test" in all_output

    @patch("builtins.print")
    def test_json_output(self, mocked_print, json_args):
        svc = RateHeterogeneity(json_args)
        svc.run()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "n_tips" in payload
        assert payload["n_tips"] == 8
        assert "regimes" in payload
        assert "single_rate" in payload
        assert "multi_rate" in payload
        assert "lrt" in payload
        assert "sigma2" in payload["single_rate"]
        assert "sigma2" in payload["multi_rate"]
        assert "chi2_p_value" in payload["lrt"]

    @patch("builtins.print")
    def test_with_simulations(self, mocked_print, sim_args):
        sim_args.json = True
        svc = RateHeterogeneity(sim_args)
        svc.run()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["lrt"]["n_sim"] == 10
        assert payload["lrt"]["sim_p_value"] is not None
        assert 0 <= payload["lrt"]["sim_p_value"] <= 1


class TestEffectSize:
    @patch("builtins.print")
    def test_r2_regime_in_json(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE,
            nsim=0,
            seed=None,
            plot=None,
            json=True,
        )
        svc = RateHeterogeneity(args)
        svc.run()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "r_squared_regime" in payload
        assert isinstance(payload["r_squared_regime"], float)
        assert np.isfinite(payload["r_squared_regime"])

    @patch("builtins.print")
    def test_r2_regime_in_text_output(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE,
            nsim=0,
            seed=None,
            plot=None,
            json=False,
        )
        svc = RateHeterogeneity(args)
        svc.run()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "R2_regime" in all_output

    @patch("builtins.print")
    def test_r2_regime_bounded(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE,
            nsim=0,
            seed=None,
            plot=None,
            json=True,
        )
        svc = RateHeterogeneity(args)
        svc.run()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["r_squared_regime"] <= 1.0

    @patch("builtins.print")
    def test_r2_regime_matches_reference(self, mocked_print):
        """R²_regime must match reference values from Fitch-parsimony regime assignment.

        Reference (PhyKIT, tip-count-weighted):
          sigma2_single = 0.0384065703  (matches R brownie.lite)
          sigma2_aquatic = 0.0088113062
          sigma2_terrestrial = 0.0500190412
          weighted avg = (2/8)*0.00881 + (6/8)*0.05002 = 0.0397170574
          R²_regime = 1 - 0.0397170574 / 0.0384065703 = -0.0341227314

        Note: R brownie.lite uses stochastic mapping (make.simmap) for
        regime assignment, so its per-regime σ² values differ from PhyKIT's
        Fitch parsimony. The single-rate σ² matches exactly.
        """
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE,
            nsim=0,
            seed=None,
            plot=None,
            json=True,
        )
        svc = RateHeterogeneity(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        # Single-rate σ² matches R brownie.lite exactly
        assert payload["single_rate"]["sigma2"] == pytest.approx(
            0.0384065703, abs=1e-6
        )
        # R²_regime is negative (multi-rate weighted avg > single rate)
        assert payload["r_squared_regime"] == pytest.approx(
            -0.0341227314, abs=1e-4
        )


class TestPlot:
    def test_plot_file_created(self, default_args):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            default_args.plot = tmppath
            svc = RateHeterogeneity(default_args)
            svc.run()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_plot_circular(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = Namespace(
                tree=TREE_SIMPLE,
                trait_data=TRAITS_FILE,
                regime_data=REGIMES_FILE,
                nsim=0,
                seed=None,
                plot=tmppath,
                json=False,
                circular=True,
            )
            svc = RateHeterogeneity(args)
            svc.run()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)
