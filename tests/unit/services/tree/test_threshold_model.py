import json
import os
from argparse import Namespace
from io import StringIO
from pathlib import Path

import numpy as np
import pytest
from Bio import Phylo

from phykit.services.tree.threshold_model import ThresholdModel
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_threshold_traits.tsv")


def _make_tree():
    return Phylo.read(
        StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick"
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="traits.tsv",
            traits="habitat,body_mass",
            types="discrete,continuous",
        )
        svc = ThresholdModel(args)
        assert svc.tree_file_path == "t.tre"
        assert svc.traits == ["habitat", "body_mass"]
        assert svc.types == ["discrete", "continuous"]
        assert svc.ngen == 100000
        assert svc.sample == 100
        assert svc.burnin == 0.2
        assert svc.seed is None
        assert svc.plot_output is None
        assert svc.json_output is False

    def test_all_options(self):
        args = Namespace(
            tree="t.tre",
            trait_data="traits.tsv",
            traits="a,b",
            types="continuous,continuous",
            ngen=5000,
            sample=10,
            burnin=0.1,
            seed=42,
            plot_output="trace.png",
            json=True,
        )
        svc = ThresholdModel(args)
        assert svc.ngen == 5000
        assert svc.sample == 10
        assert svc.burnin == 0.1
        assert svc.seed == 42
        assert svc.plot_output == "trace.png"
        assert svc.json_output is True

    def test_json_flag(self):
        args = Namespace(
            tree="t.tre",
            trait_data="traits.tsv",
            traits="a,b",
            types="discrete,discrete",
            json=True,
        )
        svc = ThresholdModel(args)
        assert svc.json_output is True

    def test_missing_traits_raises(self):
        args = Namespace(
            tree="t.tre",
            trait_data="traits.tsv",
            traits=None,
            types="discrete,continuous",
        )
        with pytest.raises(SystemExit):
            ThresholdModel(args)

    def test_invalid_type_raises(self):
        args = Namespace(
            tree="t.tre",
            trait_data="traits.tsv",
            traits="a,b",
            types="discrete,ordinal",
        )
        with pytest.raises(SystemExit):
            ThresholdModel(args)


class TestParseMultiTraitFile:
    def test_valid_parse(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "taxon\tt1\tt2\n"
            "A\t0\t1.5\n"
            "B\t1\t2.3\n"
            "C\t0\t0.8\n"
            "D\t1\t3.1\n"
        )
        tree_tips = ["A", "B", "C", "D"]
        t1, t2, names = ThresholdModel._parse_multi_trait_file(
            str(trait_file), "t1", "t2", "discrete", "continuous", tree_tips
        )
        assert set(names) == {"A", "B", "C", "D"}
        assert t1["A"] == 0.0
        assert t1["B"] == 1.0
        assert abs(t2["A"] - 1.5) < 1e-10

    def test_missing_file_raises(self):
        with pytest.raises(SystemExit):
            ThresholdModel._parse_multi_trait_file(
                "/nonexistent/file.tsv", "t1", "t2",
                "discrete", "continuous", ["A"]
            )

    def test_column_not_found_raises(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("taxon\tt1\tt2\nA\t0\t1.0\n")
        with pytest.raises(SystemExit):
            ThresholdModel._parse_multi_trait_file(
                str(trait_file), "t1", "missing_col",
                "discrete", "continuous", ["A"]
            )

    def test_non_binary_raises(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "taxon\tt1\tt2\n"
            "A\t0\t1.0\n"
            "B\t2\t2.0\n"
            "C\t1\t3.0\n"
        )
        with pytest.raises(SystemExit):
            ThresholdModel._parse_multi_trait_file(
                str(trait_file), "t1", "t2",
                "discrete", "continuous", ["A", "B", "C"]
            )

    def test_taxon_mismatch_warning(self, tmp_path, capsys):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "taxon\tt1\tt2\n"
            "A\t0\t1.0\n"
            "B\t1\t2.0\n"
            "C\t0\t3.0\n"
            "D\t1\t4.0\n"
            "E\t0\t5.0\n"
        )
        tree_tips = ["A", "B", "C", "F"]
        t1, t2, names = ThresholdModel._parse_multi_trait_file(
            str(trait_file), "t1", "t2",
            "discrete", "continuous", tree_tips
        )
        err = capsys.readouterr().err
        assert "Warning" in err
        assert set(names) == {"A", "B", "C"}

    def test_too_few_shared_taxa_raises(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("taxon\tt1\tt2\nA\t0\t1.0\nB\t1\t2.0\n")
        with pytest.raises(SystemExit):
            ThresholdModel._parse_multi_trait_file(
                str(trait_file), "t1", "t2",
                "discrete", "continuous", ["A", "B"]
            )


class TestBuildVCVMatrix:
    def test_symmetric(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        vcv = ThresholdModel._build_vcv_matrix(tree, names)
        np.testing.assert_array_almost_equal(vcv, vcv.T)

    def test_correct_shape(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        vcv = ThresholdModel._build_vcv_matrix(tree, names)
        assert vcv.shape == (4, 4)

    def test_diagonal_root_to_tip(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        vcv = ThresholdModel._build_vcv_matrix(tree, names)
        for i, name in enumerate(names):
            expected = tree.distance(tree.root, name)
            assert abs(vcv[i, i] - expected) < 1e-10


class TestInitializeLiabilities:
    def test_discrete_signs_correct(self):
        rng = np.random.default_rng(42)
        trait_values = {"A": 0.0, "B": 1.0, "C": 0.0, "D": 1.0}
        names = ["A", "B", "C", "D"]
        liabs = ThresholdModel._initialize_liabilities(
            trait_values, "discrete", names, rng
        )
        assert liabs[0] < 0  # A: state 0
        assert liabs[1] > 0  # B: state 1
        assert liabs[2] < 0  # C: state 0
        assert liabs[3] > 0  # D: state 1

    def test_continuous_passthrough(self):
        rng = np.random.default_rng(42)
        trait_values = {"A": 1.5, "B": -0.3, "C": 2.1, "D": 0.7}
        names = ["A", "B", "C", "D"]
        liabs = ThresholdModel._initialize_liabilities(
            trait_values, "continuous", names, rng
        )
        assert abs(liabs[0] - 1.5) < 1e-10
        assert abs(liabs[1] - (-0.3)) < 1e-10


class TestLogLikelihoodBivariateBM:
    def test_r_zero_equals_sum_univariate(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        C_inv = np.linalg.inv(C)
        _, logdet_C = np.linalg.slogdet(C)
        n = len(names)

        x1 = np.array([1.0, 1.5, 0.5, 2.0])
        x2 = np.array([0.3, -0.5, 1.0, 0.8])

        ll_biv = ThresholdModel._log_likelihood_bivariate_bm(
            x1, x2, C_inv, 1.0, 1.0, 0.0, 1.0, 0.5, logdet_C, n
        )

        # At r=0, bivariate likelihood should decompose to two univariates
        r1 = x1 - 1.0
        r2 = x2 - 0.5
        ll1 = -0.5 * (n * np.log(2 * np.pi) + logdet_C + n * np.log(1.0) + float(r1 @ C_inv @ r1))
        ll2 = -0.5 * (n * np.log(2 * np.pi) + logdet_C + n * np.log(1.0) + float(r2 @ C_inv @ r2))

        assert abs(ll_biv - (ll1 + ll2)) < 1e-6

    def test_positive_r_data_higher_ll(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        C_inv = np.linalg.inv(C)
        _, logdet_C = np.linalg.slogdet(C)
        n = len(names)

        # Positively correlated data
        x1 = np.array([1.0, 2.0, 3.0, 4.0])
        x2 = np.array([1.1, 2.1, 2.9, 3.8])

        ll_pos = ThresholdModel._log_likelihood_bivariate_bm(
            x1, x2, C_inv, 1.0, 1.0, 0.8, 2.5, 2.5, logdet_C, n
        )
        ll_neg = ThresholdModel._log_likelihood_bivariate_bm(
            x1, x2, C_inv, 1.0, 1.0, -0.8, 2.5, 2.5, logdet_C, n
        )

        assert ll_pos > ll_neg


class TestSampleLiabilitiesGibbs:
    def test_liabilities_respect_threshold(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        C_inv = np.linalg.inv(C)
        rng = np.random.default_rng(42)

        states = np.array([0, 1, 0, 1])
        liabilities = np.array([-0.5, 0.5, -0.3, 0.8])
        other_liabs = np.array([1.0, 2.0, 0.5, 1.5])

        new_liabs = ThresholdModel._sample_liabilities_gibbs(
            liabilities, states, C, C_inv, 1.0, 0.0,
            other_liabs, 1.0, 0.0, 1.0, rng
        )

        # State 0 must be negative, state 1 must be positive
        assert new_liabs[0] < 0
        assert new_liabs[1] > 0
        assert new_liabs[2] < 0
        assert new_liabs[3] > 0

    def test_deterministic_with_seed(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        C_inv = np.linalg.inv(C)

        states = np.array([0, 1, 0, 1])
        liabilities = np.array([-0.5, 0.5, -0.3, 0.8])
        other_liabs = np.array([1.0, 2.0, 0.5, 1.5])

        rng1 = np.random.default_rng(123)
        result1 = ThresholdModel._sample_liabilities_gibbs(
            liabilities, states, C, C_inv, 1.0, 0.0,
            other_liabs, 1.0, 0.0, 1.0, rng1
        )

        rng2 = np.random.default_rng(123)
        result2 = ThresholdModel._sample_liabilities_gibbs(
            liabilities, states, C, C_inv, 1.0, 0.0,
            other_liabs, 1.0, 0.0, 1.0, rng2
        )

        np.testing.assert_array_almost_equal(result1, result2)


class TestComputeHPD:
    def test_known_normal_samples(self):
        rng = np.random.default_rng(42)
        samples = rng.normal(0, 1, 10000)
        lo, hi = ThresholdModel._compute_hpd(samples, 0.95)
        # 95% HPD of standard normal should be approximately [-1.96, 1.96]
        assert lo < -1.5
        assert hi > 1.5
        assert (hi - lo) < 5.0  # Not absurdly wide

    def test_95_coverage(self):
        rng = np.random.default_rng(42)
        samples = rng.normal(5, 2, 10000)
        lo, hi = ThresholdModel._compute_hpd(samples, 0.95)
        in_interval = np.sum((samples >= lo) & (samples <= hi))
        coverage = in_interval / len(samples)
        assert coverage >= 0.94


class TestRunMCMC:
    def test_expected_keys(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        rng = np.random.default_rng(42)

        t1 = {"A": 0.0, "B": 1.0, "C": 0.0, "D": 1.0}
        t2 = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.8}

        result = ThresholdModel._run_mcmc(
            t1, t2, "discrete", "continuous",
            names, C, 1000, 10, 0.2, rng
        )

        assert "r" in result
        assert "sigma2_1" in result
        assert "sigma2_2" in result
        assert "a1" in result
        assert "a2" in result
        assert "acceptance_rates" in result

    def test_sample_count(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        rng = np.random.default_rng(42)

        t1 = {"A": 0.0, "B": 1.0, "C": 0.0, "D": 1.0}
        t2 = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.8}

        ngen = 2000
        sample = 10
        burnin = 0.5
        result = ThresholdModel._run_mcmc(
            t1, t2, "discrete", "continuous",
            names, C, ngen, sample, burnin, rng
        )

        expected = (ngen - int(ngen * burnin)) // sample
        assert len(result["r"]) == expected

    def test_r_in_bounds(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        rng = np.random.default_rng(42)

        t1 = {"A": 0.0, "B": 1.0, "C": 0.0, "D": 1.0}
        t2 = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.8}

        result = ThresholdModel._run_mcmc(
            t1, t2, "discrete", "continuous",
            names, C, 2000, 10, 0.2, rng
        )

        assert np.all(result["r"] >= -1.0)
        assert np.all(result["r"] <= 1.0)

    def test_sigma2_positive(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        rng = np.random.default_rng(42)

        t1 = {"A": 0.0, "B": 1.0, "C": 0.0, "D": 1.0}
        t2 = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.8}

        result = ThresholdModel._run_mcmc(
            t1, t2, "discrete", "continuous",
            names, C, 2000, 10, 0.2, rng
        )

        assert np.all(result["sigma2_1"] > 0)
        assert np.all(result["sigma2_2"] > 0)

    def test_reproducible_with_seed(self):
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)

        t1 = {"A": 0.0, "B": 1.0, "C": 0.0, "D": 1.0}
        t2 = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.8}

        rng1 = np.random.default_rng(42)
        result1 = ThresholdModel._run_mcmc(
            t1, t2, "discrete", "continuous",
            names, C, 1000, 10, 0.2, rng1
        )

        rng2 = np.random.default_rng(42)
        result2 = ThresholdModel._run_mcmc(
            t1, t2, "discrete", "continuous",
            names, C, 1000, 10, 0.2, rng2
        )

        np.testing.assert_array_almost_equal(result1["r"], result2["r"])

    def test_continuous_continuous(self):
        """Both continuous: should skip Gibbs, just run MH."""
        tree = _make_tree()
        names = sorted([t.name for t in tree.get_terminals()])
        C = ThresholdModel._build_vcv_matrix(tree, names)
        rng = np.random.default_rng(42)

        t1 = {"A": 1.0, "B": 2.0, "C": 0.5, "D": 1.8}
        t2 = {"A": 1.1, "B": 2.1, "C": 0.6, "D": 1.9}

        result = ThresholdModel._run_mcmc(
            t1, t2, "continuous", "continuous",
            names, C, 1000, 10, 0.2, rng
        )

        assert len(result["r"]) > 0
        assert np.all(result["sigma2_1"] > 0)


class TestRun:
    def test_text_output(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            traits="habitat,body_mass",
            types="discrete,continuous",
            ngen=2000,
            sample=10,
            burnin=0.2,
            seed=42,
            plot_output=None,
            json=False,
        )
        svc = ThresholdModel(args)
        svc.run()
        out, _ = capsys.readouterr()
        assert "Trait 1: habitat" in out
        assert "Trait 2: body_mass" in out
        assert "Posterior correlation (r):" in out
        assert "Posterior sigma2_1:" in out
        assert "Acceptance rates:" in out

    def test_json_output(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            traits="habitat,body_mass",
            types="discrete,continuous",
            ngen=2000,
            sample=10,
            burnin=0.2,
            seed=42,
            plot_output=None,
            json=True,
        )
        svc = ThresholdModel(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "metadata" in data
        assert "summary" in data
        assert "posterior_samples" in data
        assert data["metadata"]["trait1"] == "habitat"
        assert "r" in data["summary"]
