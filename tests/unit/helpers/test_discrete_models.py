"""
Unit tests for the shared discrete model utilities module.

Tests Q-matrix construction, Felsenstein pruning, parameter counting,
and trait parsing. These functions are the core machinery used by
stochastic_character_map, ancestral_reconstruction, and fit_discrete.
"""
import numpy as np
import pytest
from pathlib import Path

from phykit.helpers.discrete_models import (
    build_q_matrix,
    count_params,
    felsenstein_pruning,
    matrix_exp,
    parse_discrete_traits,
    VALID_DISCRETE_MODELS,
)


class TestBuildQMatrix:
    def test_er_rows_sum_to_zero(self):
        Q = build_q_matrix(np.array([0.5]), 3, "ER")
        for i in range(3):
            assert np.sum(Q[i, :]) == pytest.approx(0.0)

    def test_sym_is_symmetric(self):
        Q = build_q_matrix(np.array([0.1, 0.2, 0.3]), 3, "SYM")
        assert np.allclose(Q, Q.T)

    def test_ard_is_asymmetric(self):
        params = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
        Q = build_q_matrix(params, 3, "ARD")
        assert not np.allclose(Q, Q.T)

    def test_diagonal_negative(self):
        Q = build_q_matrix(np.array([1.0]), 2, "ER")
        assert Q[0, 0] < 0
        assert Q[1, 1] < 0


class TestCountParams:
    def test_er(self):
        assert count_params(2, "ER") == 1
        assert count_params(5, "ER") == 1

    def test_sym(self):
        assert count_params(2, "SYM") == 1
        assert count_params(3, "SYM") == 3
        assert count_params(4, "SYM") == 6

    def test_ard(self):
        assert count_params(2, "ARD") == 2
        assert count_params(3, "ARD") == 6
        assert count_params(4, "ARD") == 12

    def test_invalid_model(self):
        with pytest.raises(SystemExit):
            count_params(3, "FOO")


class TestMatrixExp:
    def test_identity_at_zero_time(self):
        Q = build_q_matrix(np.array([1.0]), 2, "ER")
        P = matrix_exp(Q, 0.0)
        assert np.allclose(P, np.eye(2))

    def test_rows_sum_to_one(self):
        Q = build_q_matrix(np.array([0.5]), 3, "ER")
        P = matrix_exp(Q, 1.0)
        for i in range(3):
            assert np.sum(P[i, :]) == pytest.approx(1.0)


class TestFelsensteinPruning:
    def test_returns_negative_loglik(self, tree_simple):
        tree_simple.root_with_outgroup("dog")
        tip_states = {"dog": "A", "cat": "A", "monkey": "B", "bear": "A",
                       "raccoon": "A", "seal": "B", "sea_lion": "B", "weasel": "A"}
        states = ["A", "B"]
        Q = build_q_matrix(np.array([0.1]), 2, "ER")
        pi = np.ones(2) / 2.0
        _, loglik = felsenstein_pruning(tree_simple, tip_states, Q, pi, states)
        assert loglik < 0

    def test_perfect_data_high_loglik(self, tree_simple):
        """All tips same state should give high likelihood."""
        tree_simple.root_with_outgroup("dog")
        tip_states = {t.name: "A" for t in tree_simple.get_terminals()}
        states = ["A", "B"]
        Q = build_q_matrix(np.array([0.01]), 2, "ER")
        pi = np.ones(2) / 2.0
        _, loglik = felsenstein_pruning(tree_simple, tip_states, Q, pi, states)
        # With low rate and all-same tips, likelihood should be reasonable
        assert loglik > -5.0


class TestParseDiscreteTraits:
    def test_parse_multi_column(self, tmp_path):
        f = tmp_path / "traits.tsv"
        f.write_text("taxon\tdiet\nA\tcarnivore\nB\therbivore\nC\tomnivore\n")
        result = parse_discrete_traits(str(f), ["A", "B", "C", "D"], trait_column="diet")
        assert result == {"A": "carnivore", "B": "herbivore", "C": "omnivore"}

    def test_parse_two_column(self, tmp_path):
        f = tmp_path / "traits.tsv"
        f.write_text("A\tcarnivore\nB\therbivore\nC\tomnivore\n")
        result = parse_discrete_traits(str(f), ["A", "B", "C", "D"])
        assert result == {"A": "carnivore", "B": "herbivore", "C": "omnivore"}

    def test_missing_file_exits(self):
        with pytest.raises(SystemExit):
            parse_discrete_traits("/nonexistent.tsv", ["A", "B", "C"])

    def test_too_few_shared_taxa_exits(self, tmp_path):
        f = tmp_path / "traits.tsv"
        f.write_text("A\tcarnivore\nB\therbivore\n")
        with pytest.raises(SystemExit):
            parse_discrete_traits(str(f), ["X", "Y", "Z"])
