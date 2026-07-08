"""
Unit tests for fit_discrete (discrete trait model comparison).

Fits ER, SYM, and ARD Mk models via maximum likelihood and compares
them using AIC/BIC. Expected log-likelihood values cross-validated
against R's geiger::fitDiscrete() (with multi2di() for multifurcations).
See tests/r_validation/validate_fit_discrete.R for the R script.

Note: PhyKIT and geiger may produce slightly different log-likelihoods
(within ~1.0) due to differences in handling multifurcating trees and
multi-start optimization. The model ranking should agree.
"""
import math
import subprocess
import sys
from io import StringIO
import pytest
import numpy as np
from argparse import Namespace
from Bio import Phylo

from phykit.services.tree.fit_discrete import FitDiscrete
import phykit.services.tree.fit_discrete as fit_discrete_module
from phykit.helpers.discrete_models import (
    build_q_matrix,
    count_params,
    felsenstein_pruning,
    fit_q_matrix,
    parse_discrete_traits,
    VALID_DISCRETE_MODELS,
    _felsenstein_loglik_prepared,
    _prepare_felsenstein_context,
)


def test_module_import_does_not_import_numpy_scipy_or_discrete_helper():
    code = """
import sys
import phykit.services.tree.fit_discrete as module
assert module.VALID_DISCRETE_MODELS == frozenset(["ER", "SYM", "ARD"])
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy.linalg" not in sys.modules
assert "scipy.optimize" not in sys.modules
assert "phykit.helpers.discrete_models" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_discrete_helper_wrappers_cache_after_first_resolution(monkeypatch):
    previous_count_params = fit_discrete_module._COUNT_PARAMS
    fit_discrete_module._COUNT_PARAMS = None

    try:
        assert fit_discrete_module.count_params(3, "SYM") == 3
        cached_count_params = fit_discrete_module._COUNT_PARAMS
        assert cached_count_params is not None

        original_import = __import__

        def fail_discrete_models_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "phykit.helpers.discrete_models":
                raise AssertionError(
                    "discrete helper should be reused after first resolution"
                )
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr("builtins.__import__", fail_discrete_models_import)

        assert fit_discrete_module.count_params(3, "SYM") == 3
        assert fit_discrete_module._COUNT_PARAMS is cached_count_params
    finally:
        fit_discrete_module._COUNT_PARAMS = previous_count_params


@pytest.fixture
def args():
    return Namespace(
        tree="tests/sample_files/tree_simple.tre",
        trait_data="tests/sample_files/tree_simple_discrete_traits.tsv",
        trait="diet",
    )


class TestFitDiscreteInit:
    def test_init_sets_fields(self, args):
        fd = FitDiscrete(args)
        assert fd.tree_file_path == args.tree
        assert fd.trait_data_path == args.trait_data
        assert fd.trait_column == "diet"
        assert fd.selected_models == ["ER", "SYM", "ARD"]
        assert fd.json_output is False

    def test_process_args_custom_models(self):
        args = Namespace(
            tree="t.tre", trait_data="d.tsv", trait="col", models="ER,ARD"
        )
        fd = FitDiscrete(args)
        assert fd.selected_models == ["ER", "ARD"]

    def test_process_args_invalid_model_exits(self):
        args = Namespace(
            tree="t.tre", trait_data="d.tsv", trait="col", models="ER,FOO"
        )
        with pytest.raises(SystemExit):
            FitDiscrete(args)


class TestDiscreteModelsShared:
    def test_default_branch_length_scan_handles_mixed_child_counts(self, monkeypatch):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree scan should not use find_clades")

        monkeypatch.setattr(tree, "find_clades", fail_find_clades)

        assert FitDiscrete._needs_default_branch_lengths(tree) is False
        tree.root.clades[2].clades[1].branch_length = None
        assert FitDiscrete._needs_default_branch_lengths(tree) is True

    def test_build_q_matrix_er(self):
        Q = build_q_matrix([0.5], 3, "ER")
        assert Q.shape == (3, 3)
        # Off-diagonals should all be 0.5
        for i in range(3):
            for j in range(3):
                if i != j:
                    assert Q[i, j] == pytest.approx(0.5)
        # Rows sum to 0
        for i in range(3):
            assert sum(Q[i, :]) == pytest.approx(0.0)

    def test_build_q_matrix_sym(self):
        # k=3 -> 3 params
        Q = build_q_matrix([0.1, 0.2, 0.3], 3, "SYM")
        assert Q[0, 1] == Q[1, 0]  # Symmetric
        assert Q[0, 2] == Q[2, 0]
        assert Q[1, 2] == Q[2, 1]

    def test_build_q_matrix_ard(self):
        # k=3 -> 6 params
        Q = build_q_matrix([0.1, 0.2, 0.3, 0.4, 0.5, 0.6], 3, "ARD")
        assert Q[0, 1] != Q[1, 0]  # Asymmetric

    def test_count_params(self):
        assert count_params(2, "ER") == 1
        assert count_params(3, "ER") == 1
        assert count_params(3, "SYM") == 3
        assert count_params(3, "ARD") == 6
        assert count_params(4, "SYM") == 6
        assert count_params(4, "ARD") == 12

    def test_prepared_pruning_matches_public_pruning(self, args):
        fd = FitDiscrete(args)
        tree = fd.read_tree_file()
        tree_tips = fd.get_tip_names_from_tree(tree)
        tip_states = parse_discrete_traits(
            args.trait_data, tree_tips, trait_column="diet"
        )
        states = sorted(set(tip_states.values()))
        pi = [1.0 / len(states)] * len(states)
        Q = build_q_matrix([0.2], len(states), "ER")

        _, expected = felsenstein_pruning(tree, tip_states, Q, pi, states)
        context = _prepare_felsenstein_context(tree, tip_states, states)
        observed = _felsenstein_loglik_prepared(context, Q, pi)

        assert observed == pytest.approx(expected)

    def test_fit_q_matrix_prepares_pruning_without_generic_traversal(
        self, args, monkeypatch
    ):
        fd = FitDiscrete(args)
        tree = fd.read_tree_file()
        tree_tips = fd.get_tip_names_from_tree(tree)
        tip_states = parse_discrete_traits(
            args.trait_data, tree_tips, trait_column="diet"
        )
        states = sorted(set(tip_states.values()))

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("fit_q_matrix should reuse direct prepared pruning")

        monkeypatch.setattr(tree, "find_clades", fail_find_clades)

        Q, lnL = fit_q_matrix(tree, tip_states, states, "ER")

        assert Q.shape == (len(states), len(states))
        assert math.isfinite(lnL)


class TestFitModel:
    def test_fit_er_returns_valid_result(self, args):
        fd = FitDiscrete(args)
        tree = fd.read_tree_file()
        tree_tips = fd.get_tip_names_from_tree(tree)
        tip_states = parse_discrete_traits(
            args.trait_data, tree_tips, trait_column="diet"
        )
        states = sorted(set(tip_states.values()))
        result = fd._fit_model(tree, tip_states, states, "ER", len(tip_states))

        assert result["model"] == "ER"
        assert result["lnL"] < 0  # Log-likelihood should be negative
        assert result["aic"] > 0
        assert result["bic"] > 0
        assert result["n_params"] == 1
        assert result["n_states"] == 3

    def test_model_comparison_produces_valid_weights(self, args):
        fd = FitDiscrete(args)
        fake_results = [
            {"model": "ER", "lnL": -10.0, "aic": 22.0, "bic": 23.0, "n_params": 1},
            {"model": "SYM", "lnL": -9.5, "aic": 25.0, "bic": 27.0, "n_params": 3},
            {"model": "ARD", "lnL": -8.0, "aic": 28.0, "bic": 32.0, "n_params": 6},
        ]
        results = fd._compute_model_comparison(fake_results)

        # Should be sorted by AIC
        assert results[0]["model"] == "ER"
        assert results[0]["delta_aic"] == 0.0

        # Weights should sum to 1
        total_w = sum(r["aic_weight"] for r in results)
        assert total_w == pytest.approx(1.0)

        # Best model should have highest weight
        assert results[0]["aic_weight"] > results[1]["aic_weight"]


class TestFitDiscreteRun:
    def test_tips_to_prune_for_states_ordered_all_shared_skips_set(self, monkeypatch):
        tree_tips = ["A", "B", "C"]
        tip_states = {"A": "x", "B": "x", "C": "y"}

        def fail_set(*_args, **_kwargs):
            raise AssertionError(
                "ordered all-shared states should not build a set"
            )

        monkeypatch.setattr("builtins.set", fail_set)

        assert FitDiscrete._tips_to_prune_for_states(tree_tips, tip_states) == []

    def test_tips_to_prune_for_states_preserves_partial_pruning(self):
        tree_tips = ["A", "B", "C", "D"]
        tip_states = {"A": "x", "C": "y"}

        assert FitDiscrete._tips_to_prune_for_states(tree_tips, tip_states) == [
            "B",
            "D",
        ]

    def _stub_run_tail(self, monkeypatch, fd, tip_states, captured):
        monkeypatch.setattr(
            "phykit.services.tree.fit_discrete.parse_discrete_traits",
            lambda *_args, **_kwargs: tip_states,
        )
        monkeypatch.setattr(
            "phykit.services.tree.fit_discrete._prepare_felsenstein_context",
            lambda *_args, **_kwargs: object(),
        )

        def fit_model(tree, tip_states_arg, states, model_name, n_obs, pruning_context=None):
            captured.setdefault("trees", []).append(tree)
            return {
                "model": model_name,
                "lnL": -1.0,
                "aic": {"ER": 2.0, "SYM": 4.0, "ARD": 6.0}[model_name],
                "bic": {"ER": 3.0, "SYM": 5.0, "ARD": 7.0}[model_name],
                "n_params": 1,
                "n_states": len(states),
                "q_matrix": [[0.0 for _ in states] for _ in states],
            }

        monkeypatch.setattr(fd, "_fit_model", fit_model)
        monkeypatch.setattr(fd, "_print_text", lambda *_args, **_kwargs: None)

    def test_run_reuses_prepared_pruning_context(self, args, monkeypatch, capsys):
        fd = FitDiscrete(args)
        sentinel_context = object()
        fit_contexts = []

        def fake_prepare(tree, tip_states, states):
            return sentinel_context

        def fake_fit_q_matrix(tree, tip_states, states, model, pruning_context=None):
            fit_contexts.append(pruning_context)
            return np.zeros((len(states), len(states))), -1.0

        monkeypatch.setattr(
            "phykit.services.tree.fit_discrete._prepare_felsenstein_context",
            fake_prepare,
        )
        monkeypatch.setattr(
            "phykit.services.tree.fit_discrete.fit_q_matrix",
            fake_fit_q_matrix,
        )

        fd.run()
        capsys.readouterr()

        assert fit_contexts == [sentinel_context] * 3

    def test_run_all_tips_present_uses_read_only_tree_without_copy(
        self, args, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "y"}
        fd = FitDiscrete(args)
        captured = {}

        monkeypatch.setattr(fd, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(
            fd,
            "read_tree_file",
            lambda: (_ for _ in ()).throw(
                AssertionError("mutable tree read should not be used")
            ),
        )
        monkeypatch.setattr(
            fd,
            "_fast_copy",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("tree should not be copied")
            ),
        )
        self._stub_run_tail(monkeypatch, fd, tip_states, captured)

        fd.run()

        assert captured["trees"] == [tree, tree, tree]
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}

    def test_run_missing_trait_states_copies_before_pruning(self, args, monkeypatch):
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        tip_states = {"A": "x", "B": "x", "C": "y"}
        fd = FitDiscrete(args)
        original_fast_copy = fd._fast_copy
        copied_trees = []
        captured = {}

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        monkeypatch.setattr(fd, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(fd, "_fast_copy", copy_spy)
        self._stub_run_tail(monkeypatch, fd, tip_states, captured)

        fd.run()

        assert len(copied_trees) == 1
        assert captured["trees"] == [copied_trees[0], copied_trees[0], copied_trees[0]]
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}
        assert {tip.name for tip in copied_trees[0].get_terminals()} == {
            "A",
            "B",
            "C",
        }

    def test_run_missing_branch_lengths_copies_before_validation(
        self, args, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("((A,B):1,(C:1,D:1):1);"),
            "newick",
        )
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "y"}
        fd = FitDiscrete(args)
        original_fast_copy = fd._fast_copy
        copied_trees = []
        captured = {}

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        monkeypatch.setattr(fd, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(fd, "_fast_copy", copy_spy)
        self._stub_run_tail(monkeypatch, fd, tip_states, captured)

        fd.run()

        assert len(copied_trees) == 1
        assert captured["trees"] == [copied_trees[0], copied_trees[0], copied_trees[0]]
        assert tree.root.clades[0].clades[0].branch_length is None
        assert copied_trees[0].root.clades[0].clades[0].branch_length == 1e-8

    def test_run_text_output(self, args, capsys):
        fd = FitDiscrete(args)
        fd.run()
        captured = capsys.readouterr()
        assert "ER" in captured.out
        assert "SYM" in captured.out
        assert "ARD" in captured.out
        assert "lnL" in captured.out

    def test_print_text_batches_model_table(self, mocker):
        fd = FitDiscrete.__new__(FitDiscrete)
        printed = mocker.patch("builtins.print")
        results = [
            {
                "model": "ER",
                "lnL": -12.34567,
                "aic": 30.12345,
                "delta_aic": 0.0,
                "aic_weight": 0.75,
                "bic": 33.45678,
                "delta_bic": 1.25,
                "n_params": 2,
            },
            {
                "model": "ARD",
                "lnL": -13.45678,
                "aic": 35.23456,
                "delta_aic": 5.11111,
                "aic_weight": 0.25,
                "bic": 32.12345,
                "delta_bic": 0.0,
                "n_params": 6,
            },
        ]
        header = (
            f"{'Model':<8}{'lnL':>12}{'AIC':>12}{'dAIC':>10}"
            f"{'AIC_w':>10}{'BIC':>12}{'dBIC':>10}{'n_params':>10}"
        )
        expected = "\n".join(
            [
                "Number of taxa: 8",
                "Number of states: 3",
                "States: A, B, C",
                "",
                header,
                "-" * len(header),
                (
                    f"{'ER':<8}{-12.34567:>12.4f}{30.12345:>12.4f}"
                    f"{0.0:>10.4f}{0.75:>10.4f}{33.45678:>12.4f}"
                    f"{1.25:>10.4f}{2:>10d}"
                ),
                (
                    f"{'ARD':<8}{-13.45678:>12.4f}{35.23456:>12.4f}"
                    f"{5.11111:>10.4f}{0.25:>10.4f}{32.12345:>12.4f}"
                    f"{0.0:>10.4f}{6:>10d}"
                ),
            ]
        )

        fd._print_text(results, 8, ["A", "B", "C"])

        printed.assert_called_once_with(expected)

    def test_run_json_output(self, mocker, args):
        args.json = True
        fd = FitDiscrete(args)
        mocked_json = mocker.patch("phykit.services.tree.fit_discrete.print_json")
        fd.run()
        payload = mocked_json.call_args.args[0]
        assert payload["n_taxa"] == 8
        assert payload["n_states"] == 3
        assert len(payload["models"]) == 3

    def test_run_subset_models(self, capsys):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            trait_data="tests/sample_files/tree_simple_discrete_traits.tsv",
            trait="diet",
            models="ER,ARD",
        )
        fd = FitDiscrete(args)
        fd.run()
        captured = capsys.readouterr()
        assert "ER" in captured.out
        assert "ARD" in captured.out
        # SYM should not appear (only 2 models fitted)
        lines = [l for l in captured.out.strip().split("\n") if l.startswith("ER") or l.startswith("SYM") or l.startswith("ARD")]
        assert len(lines) == 2
