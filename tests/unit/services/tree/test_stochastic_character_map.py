import json
import os
import sys
import tempfile

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.stochastic_character_map import StochasticCharacterMap
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
DISCRETE_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_discrete_traits.tsv")


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=DISCRETE_TRAITS_FILE,
        trait="diet",
        model="ER",
        nsim=10,
        seed=42,
        plot=None,
        json=False,
    )


@pytest.fixture
def ard_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=DISCRETE_TRAITS_FILE,
        trait="diet",
        model="ARD",
        nsim=10,
        seed=42,
        plot=None,
        json=False,
    )


@pytest.fixture
def json_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=DISCRETE_TRAITS_FILE,
        trait="diet",
        model="ER",
        nsim=10,
        seed=42,
        plot=None,
        json=True,
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            trait="diet",
        )
        svc = StochasticCharacterMap.__new__(StochasticCharacterMap)
        parsed = svc.process_args(args)
        assert parsed["tree_file_path"] == "t.tre"
        assert parsed["trait_data_path"] == "d.tsv"
        assert parsed["trait_column"] == "diet"
        assert parsed["model"] == "ER"
        assert parsed["nsim"] == 100
        assert parsed["seed"] is None
        assert parsed["plot_output"] is None
        assert parsed["json_output"] is False

    def test_overrides(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            trait="diet",
            model="ARD",
            nsim=500,
            seed=123,
            plot="out.png",
            json=True,
        )
        svc = StochasticCharacterMap.__new__(StochasticCharacterMap)
        parsed = svc.process_args(args)
        assert parsed["model"] == "ARD"
        assert parsed["nsim"] == 500
        assert parsed["seed"] == 123
        assert parsed["plot_output"] == "out.png"
        assert parsed["json_output"] is True


class TestDiscreteTraitParsing:
    def test_basic_parsing(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        assert len(tip_states) == 8
        assert tip_states["raccoon"] == "carnivore"
        assert tip_states["monkey"] == "herbivore"
        assert tip_states["dog"] == "omnivore"

    def test_missing_column(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        with pytest.raises(PhykitUserError):
            svc._parse_discrete_trait_file(
                DISCRETE_TRAITS_FILE, "nonexistent", tree_tips
            )

    def test_comments_and_blanks(self, default_args):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write("# This is a comment\n")
            f.write("\n")
            f.write("taxon\ttrait\n")
            f.write("raccoon\tA\n")
            f.write("bear\tA\n")
            f.write("# Another comment\n")
            f.write("sea_lion\tB\n")
            f.write("seal\tB\n")
            f.write("monkey\tA\n")
            f.write("cat\tB\n")
            f.write("weasel\tA\n")
            f.write("dog\tB\n")
            tmppath = f.name

        try:
            svc = StochasticCharacterMap(default_args)
            tree = svc.read_tree_file()
            tree_tips = svc.get_tip_names_from_tree(tree)
            tip_states = svc._parse_discrete_trait_file(tmppath, "trait", tree_tips)
            assert len(tip_states) == 8
        finally:
            os.unlink(tmppath)


class TestQMatrixFitting:
    """Validated against R 4.4.0 with phytools::fitMk on same tree and traits."""

    def test_er_model_loglik(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))
        Q, loglik = svc._fit_q_matrix(tree, tip_states, states, "ER")

        # R reference: fitMk(tree, trait, model="ER") -> loglik = -8.788898
        # Python finds a marginally better point on a flat likelihood surface
        # (-8.787360 vs -8.788898), both valid ML estimates
        assert loglik == pytest.approx(-8.7889, abs=0.01)

        # Q should have equal off-diagonal rates for ER
        k = len(states)
        off_diag = []
        for i in range(k):
            for j in range(k):
                if i != j:
                    off_diag.append(Q[i, j])
        # All off-diagonal rates should be equal for ER
        assert all(abs(r - off_diag[0]) < 1e-6 for r in off_diag)
        # Rate should be positive
        assert off_diag[0] > 0

    def test_sym_model(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))

        Q_er, loglik_er = svc._fit_q_matrix(tree, tip_states, states, "ER")
        Q_sym, loglik_sym = svc._fit_q_matrix(tree, tip_states, states, "SYM")

        # SYM should give same or better loglik than ER (superset model)
        assert loglik_sym >= loglik_er - 0.01

        # Q should be symmetric
        k = len(states)
        for i in range(k):
            for j in range(k):
                if i != j:
                    assert Q_sym[i, j] == pytest.approx(Q_sym[j, i], abs=1e-4)

    def test_ard_model_loglik(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))

        _, loglik_er = svc._fit_q_matrix(tree, tip_states, states, "ER")
        Q_ard, loglik_ard = svc._fit_q_matrix(tree, tip_states, states, "ARD")

        # R reference: fitMk(tree, trait, model="ARD") -> loglik = -8.430467
        # Python's multi-start optimizer finds a better local optimum (~-8.385)
        # on this 6-parameter surface. Both are valid ML estimates.
        assert loglik_ard >= -8.45  # at least as good as R
        assert np.isfinite(loglik_ard)

        # ARD should give same or better loglik than ER (superset model)
        assert loglik_ard >= loglik_er - 0.01

    def test_q_matrix_structure(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))
        Q, _ = svc._fit_q_matrix(tree, tip_states, states, "ER")

        k = len(states)
        # Rows should sum to 0
        for i in range(k):
            assert abs(Q[i, :].sum()) < 1e-10

        # Off-diagonal should be positive
        for i in range(k):
            for j in range(k):
                if i != j:
                    assert Q[i, j] > 0

        # Diagonal should be negative
        for i in range(k):
            assert Q[i, i] < 0


class TestFelsensteinPruning:
    def test_conditional_likelihoods(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))
        k = len(states)
        pi = np.ones(k) / k

        # Use a simple Q
        Q = np.array([[-2, 1, 1], [1, -2, 1], [1, 1, -2]], dtype=float)
        cond_liks, loglik = svc._felsenstein_pruning(
            tree, tip_states, Q, pi, states
        )

        # Should have conditional likelihoods for all nodes
        node_count = sum(1 for _ in tree.find_clades())
        assert len(cond_liks) == node_count

        # Log-likelihood should be finite and negative
        assert np.isfinite(loglik)
        assert loglik < 0

        # Tip conditional likelihoods should be indicator vectors
        for tip in tree.get_terminals():
            lik = cond_liks[id(tip)]
            assert lik.sum() == pytest.approx(1.0)
            assert all(v in (0.0, 1.0) for v in lik)


class TestStochasticMapping:
    def test_mapping_produces_valid_histories(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))
        state_idx = {s: i for i, s in enumerate(states)}
        k = len(states)
        pi = np.ones(k) / k

        Q, _ = svc._fit_q_matrix(tree, tip_states, states, "ER")
        rng = np.random.default_rng(42)
        mapping = svc._run_single_simulation(
            tree, tip_states, Q, pi, states, rng
        )

        assert "node_states" in mapping
        assert "branch_histories" in mapping

        # Check that tip states match
        for tip in tree.get_terminals():
            if tip.name in tip_states:
                expected = state_idx[tip_states[tip.name]]
                assert mapping["node_states"][id(tip)] == expected

        # Check that histories have valid states
        for clade_id, history in mapping["branch_histories"].items():
            assert len(history) >= 1
            for time_val, state in history:
                assert 0 <= state < k
                assert time_val >= 0

    def test_branch_times_sum_correctly(self, default_args):
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))
        k = len(states)
        pi = np.ones(k) / k

        Q, _ = svc._fit_q_matrix(tree, tip_states, states, "ER")
        rng = np.random.default_rng(42)
        mapping = svc._run_single_simulation(
            tree, tip_states, Q, pi, states, rng
        )

        # Sum dwelling times across all branches
        total_time = 0.0
        for clade in tree.find_clades(order="preorder"):
            if clade == tree.root:
                continue
            clade_id = id(clade)
            if clade_id not in mapping["branch_histories"]:
                continue
            t = clade.branch_length if clade.branch_length else 0.0
            history = mapping["branch_histories"][clade_id]
            branch_total = 0.0
            for h_idx in range(len(history)):
                start_t = history[h_idx][0]
                if h_idx + 1 < len(history):
                    end_t = history[h_idx + 1][0]
                else:
                    end_t = t
                branch_total += end_t - start_t
            assert branch_total == pytest.approx(t, abs=1e-6)
            total_time += branch_total

        # Total tree length = 277.2772 (from R)
        assert total_time == pytest.approx(277.2772, abs=0.01)


class TestSummary:
    def test_dwelling_times_sum_to_tree_length(self, default_args):
        default_args.nsim = 50
        svc = StochasticCharacterMap(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        tip_states = svc._parse_discrete_trait_file(
            DISCRETE_TRAITS_FILE, "diet", tree_tips
        )
        states = sorted(set(tip_states.values()))
        k = len(states)
        pi = np.ones(k) / k

        Q, _ = svc._fit_q_matrix(tree, tip_states, states, "ER")
        rng = np.random.default_rng(42)
        mappings = []
        for _ in range(50):
            m = svc._run_single_simulation(
                tree, tip_states, Q, pi, states, rng
            )
            mappings.append(m)

        summary = svc._summarize_simulations(mappings, states, tree)
        mean_dwelling = summary["mean_dwelling_times"]

        # Total should equal tree length
        assert mean_dwelling.sum() == pytest.approx(277.2772, abs=0.5)


class TestRun:
    def test_text_output(self, default_args, capsys):
        svc = StochasticCharacterMap(default_args)
        svc.run()
        captured = capsys.readouterr()
        assert "Stochastic Character Mapping (SIMMAP)" in captured.out
        assert "Model: ER" in captured.out
        assert "Log-likelihood:" in captured.out
        assert "Mean dwelling times:" in captured.out
        assert "Mean transitions:" in captured.out

    def test_json_output(self, json_args, capsys):
        svc = StochasticCharacterMap(json_args)
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert payload["model"] == "ER"
        assert payload["nsim"] == 10
        assert "q_matrix" in payload
        assert "log_likelihood" in payload
        assert "states" in payload
        assert "mean_dwelling_times" in payload
        assert "mean_transitions" in payload


class TestPlot:
    def test_plot_file_created(self, default_args):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(
            suffix=".png", delete=False
        ) as f:
            tmppath = f.name

        try:
            default_args.plot = tmppath
            svc = StochasticCharacterMap(default_args)
            svc.run()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)
