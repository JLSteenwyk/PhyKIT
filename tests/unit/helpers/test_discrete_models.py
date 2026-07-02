"""
Unit tests for the shared discrete model utilities module.

Tests Q-matrix construction, Felsenstein pruning, parameter counting,
and trait parsing. These functions are the core machinery used by
stochastic_character_map, ancestral_reconstruction, and fit_discrete.
"""
import builtins
import importlib
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from phykit.errors import PhykitUserError
from phykit.helpers.discrete_models import (
    _felsenstein_loglik_er_rate,
    _felsenstein_loglik_prepared,
    _postorder_clades_direct,
    _prepare_felsenstein_context,
    build_q_matrix,
    count_params,
    felsenstein_pruning,
    fit_q_matrix,
    matrix_exp,
    parse_discrete_traits,
    VALID_DISCRETE_MODELS,
)


def test_module_import_does_not_import_numpy_or_scipy():
    code = """
import sys
import phykit.helpers.discrete_models as module
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy.linalg" not in sys.modules
assert "scipy.optimize" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_module_import_does_not_import_scipy_linalg_or_optimize(monkeypatch):
    module_name = "phykit.helpers.discrete_models"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if (
            name == "scipy.linalg"
            or name.startswith("scipy.linalg.")
            or name == "scipy.optimize"
            or name.startswith("scipy.optimize.")
        ):
            raise AssertionError(
                "discrete_models module import should not import SciPy linalg/optimize"
            )
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", guarded_import)
    try:
        importlib.import_module(module_name)
    finally:
        imported = sys.modules.pop(module_name, None)
        if previous is not None:
            sys.modules[module_name] = previous
        parent_name, _, child_name = module_name.rpartition(".")
        parent = sys.modules.get(parent_name)
        if parent is not None:
            if previous is not None:
                setattr(parent, child_name, previous)
            elif getattr(parent, child_name, None) is imported:
                delattr(parent, child_name)


def test_lazy_numpy_caches_resolved_attributes():
    import phykit.helpers.discrete_models as discrete_models

    lazy_np = discrete_models._LazyNumpy()

    array_attr = lazy_np.array

    assert lazy_np.__dict__["array"] is array_attr
    assert lazy_np.array is array_attr
    assert lazy_np._module is not None


class TestBuildQMatrix:
    def test_er_rows_sum_to_zero(self):
        Q = build_q_matrix(np.array([0.5]), 3, "ER")
        for i in range(3):
            assert np.sum(Q[i, :]) == pytest.approx(0.0)

    def test_er_uses_equal_off_diagonal_rate_and_direct_diagonal(self):
        Q = build_q_matrix(np.array([0.5]), 4, "ER")
        expected = np.array(
            [
                [-1.5, 0.5, 0.5, 0.5],
                [0.5, -1.5, 0.5, 0.5],
                [0.5, 0.5, -1.5, 0.5],
                [0.5, 0.5, 0.5, -1.5],
            ]
        )
        np.testing.assert_allclose(Q, expected)

    def test_er_matrix_construction_avoids_zero_fill(self, monkeypatch):
        import phykit.helpers.discrete_models as discrete_models

        monkeypatch.setattr(
            discrete_models.np,
            "zeros",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("ER construction should not zero-fill first")
            ),
        )
        monkeypatch.setattr(
            discrete_models.np,
            "fill_diagonal",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("ER construction should assign diagonal directly")
            ),
        )

        Q = build_q_matrix(np.array([0.5]), 4, "ER")

        np.testing.assert_allclose(Q.sum(axis=1), np.zeros(4))

    def test_sym_is_symmetric(self):
        Q = build_q_matrix(np.array([0.1, 0.2, 0.3]), 3, "SYM")
        assert np.allclose(Q, Q.T)

    def test_small_sym_matrix_construction_uses_direct_layout(self, monkeypatch):
        import phykit.helpers.discrete_models as discrete_models

        monkeypatch.setattr(
            discrete_models.np,
            "zeros",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("small SYM construction should not zero-fill first")
            ),
        )
        monkeypatch.setattr(
            discrete_models.np,
            "fill_diagonal",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("small SYM construction should assign diagonal directly")
            ),
        )

        sym2 = build_q_matrix(np.array([0.4]), 2, "SYM")
        sym3 = build_q_matrix(np.array([0.1, 0.2, 0.3]), 3, "SYM")

        np.testing.assert_allclose(sym2, np.array([[-0.4, 0.4], [0.4, -0.4]]))
        np.testing.assert_allclose(
            sym3,
            np.array(
                [
                    [-0.3, 0.1, 0.2],
                    [0.1, -0.4, 0.3],
                    [0.2, 0.3, -0.5],
                ],
            ),
        )

    def test_sym_preserves_parameter_order_for_four_states(self):
        params = np.arange(1, 7, dtype=float)
        Q = build_q_matrix(params, 4, "SYM")
        expected = np.array(
            [
                [-6.0, 1.0, 2.0, 3.0],
                [1.0, -10.0, 4.0, 5.0],
                [2.0, 4.0, -12.0, 6.0],
                [3.0, 5.0, 6.0, -14.0],
            ]
        )
        np.testing.assert_allclose(Q, expected)

    def test_ard_is_asymmetric(self):
        params = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
        Q = build_q_matrix(params, 3, "ARD")
        assert not np.allclose(Q, Q.T)

    def test_small_ard_matrix_construction_uses_direct_layout(self, monkeypatch):
        import phykit.helpers.discrete_models as discrete_models

        monkeypatch.setattr(
            discrete_models.np,
            "zeros",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("small ARD construction should not zero-fill first")
            ),
        )
        monkeypatch.setattr(
            discrete_models.np,
            "fill_diagonal",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("small ARD construction should assign diagonal directly")
            ),
        )

        ard2 = build_q_matrix(np.array([0.4, 0.7]), 2, "ARD")
        ard3 = build_q_matrix(
            np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6]), 3, "ARD"
        )

        np.testing.assert_allclose(ard2, np.array([[-0.4, 0.4], [0.7, -0.7]]))
        np.testing.assert_allclose(
            ard3,
            np.array(
                [
                    [-0.3, 0.1, 0.2],
                    [0.3, -0.7, 0.4],
                    [0.5, 0.6, -1.1],
                ],
            ),
        )

    def test_ard_preserves_row_major_off_diagonal_order(self):
        params = np.arange(1, 13, dtype=float)
        Q = build_q_matrix(params, 4, "ARD")
        expected = np.array(
            [
                [-6.0, 1.0, 2.0, 3.0],
                [4.0, -15.0, 5.0, 6.0],
                [7.0, 8.0, -24.0, 9.0],
                [10.0, 11.0, 12.0, -33.0],
            ]
        )
        np.testing.assert_allclose(Q, expected)

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

    def test_two_state_er_uses_analytic_transition_matrix(self, monkeypatch):
        import phykit.helpers.discrete_models as discrete_models

        def fail_expm(*_args, **_kwargs):
            raise AssertionError("two-state ER should not call scipy expm")

        monkeypatch.setattr(discrete_models, "expm", fail_expm)
        monkeypatch.setattr(
            discrete_models.np,
            "exp",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("two-state matrix exponential should use math.exp")
            ),
            raising=False,
        )

        rate = 0.5
        Q = build_q_matrix(np.array([rate]), 2, "ER")
        observed = discrete_models.matrix_exp(Q, 1.0)
        decay = math.exp(-2.0 * rate)
        expected = np.array(
            [
                [0.5 + 0.5 * decay, 0.5 - 0.5 * decay],
                [0.5 - 0.5 * decay, 0.5 + 0.5 * decay],
            ]
        )

        np.testing.assert_allclose(observed, expected)

    def test_two_state_ard_uses_analytic_transition_matrix(self, monkeypatch):
        import phykit.helpers.discrete_models as discrete_models
        from scipy.linalg import expm

        def fail_expm(*_args, **_kwargs):
            raise AssertionError("two-state ARD should not call scipy expm")

        Q = build_q_matrix(np.array([0.35, 1.25]), 2, "ARD")
        expected = expm(Q * 0.7)
        monkeypatch.setattr(discrete_models, "expm", fail_expm)

        observed = discrete_models.matrix_exp(Q, 0.7)

        np.testing.assert_allclose(observed, expected)

    def test_repeated_calls_cache_scipy_expm_import(self, monkeypatch):
        import phykit.helpers.discrete_models as discrete_models

        previous_expm = discrete_models._EXPM
        discrete_models._EXPM = None
        original_import = builtins.__import__
        scipy_linalg_imports = 0

        def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
            nonlocal scipy_linalg_imports
            if name == "scipy.linalg":
                scipy_linalg_imports += 1
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", counting_import)
        try:
            Q = build_q_matrix(np.array([0.5]), 3, "ER")
            discrete_models.matrix_exp(Q, 1.0)
            discrete_models.matrix_exp(Q, 2.0)
        finally:
            discrete_models._EXPM = previous_expm

        assert scipy_linalg_imports == 1

    def test_repeated_calls_cache_scipy_minimize_import(self, monkeypatch):
        import phykit.helpers.discrete_models as discrete_models

        previous_minimize = discrete_models._MINIMIZE
        discrete_models._MINIMIZE = None
        original_import = builtins.__import__
        scipy_optimize_imports = 0

        def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
            nonlocal scipy_optimize_imports
            if name == "scipy.optimize":
                scipy_optimize_imports += 1
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", counting_import)
        try:
            objective = lambda x: (x[0] - 1.0) ** 2
            discrete_models.minimize(
                objective, [0.0], method="Nelder-Mead", options={"maxiter": 1}
            )
            first_call_imports = scipy_optimize_imports
            discrete_models.minimize(
                objective, [0.0], method="Nelder-Mead", options={"maxiter": 1}
            )
        finally:
            discrete_models._MINIMIZE = previous_minimize

        assert first_call_imports > 0
        assert scipy_optimize_imports == first_call_imports

    def test_repeated_calls_cache_scipy_minimize_scalar_import(self, monkeypatch):
        import phykit.helpers.discrete_models as discrete_models

        previous_minimize_scalar = discrete_models._MINIMIZE_SCALAR
        discrete_models._MINIMIZE_SCALAR = None
        original_import = builtins.__import__
        scipy_optimize_imports = 0

        def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
            nonlocal scipy_optimize_imports
            if name == "scipy.optimize":
                scipy_optimize_imports += 1
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", counting_import)
        try:
            objective = lambda x: (x - 1.0) ** 2
            discrete_models.minimize_scalar(
                objective, bounds=(0.0, 2.0), method="bounded"
            )
            first_call_imports = scipy_optimize_imports
            discrete_models.minimize_scalar(
                objective, bounds=(0.0, 2.0), method="bounded"
            )
        finally:
            discrete_models._MINIMIZE_SCALAR = previous_minimize_scalar

        assert first_call_imports > 0
        assert scipy_optimize_imports == first_call_imports


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

    def test_repeated_branch_lengths_reuse_transition_matrices(self, mocker):
        from Bio.Phylo.Newick import Clade, Tree
        import phykit.helpers.discrete_models as discrete_models

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, clades=[
                        Clade(branch_length=1.0, name="A"),
                        Clade(branch_length=1.0, name="B"),
                    ]),
                    Clade(branch_length=1.0, clades=[
                        Clade(branch_length=1.0, name="C"),
                        Clade(branch_length=1.0, name="D"),
                    ]),
                ],
            )
        )
        tip_states = {"A": "0", "B": "1", "C": "0", "D": "1"}
        states = ["0", "1"]
        Q = build_q_matrix(np.array([0.1]), 2, "ER")
        pi = np.ones(2) / 2.0
        spy = mocker.spy(discrete_models, "matrix_exp")

        _, loglik = felsenstein_pruning(tree, tip_states, Q, pi, states)

        assert loglik < 0
        assert spy.call_count == 1

    def test_pruning_root_likelihood_uses_dot_product(self, monkeypatch):
        from Bio.Phylo.Newick import Clade, Tree
        import phykit.helpers.discrete_models as discrete_models

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="A"),
                    Clade(branch_length=1.0, name="B"),
                    Clade(branch_length=1.0, name="C"),
                ],
            )
        )
        tip_states = {"A": "0", "B": "1", "C": "2"}
        states = ["0", "1", "2"]
        Q = build_q_matrix(np.array([0.1]), 3, "ER")
        pi = np.ones(3) / 3.0

        def fail_sum(*_args, **_kwargs):
            raise AssertionError("root likelihood should use np.dot")

        monkeypatch.setattr(discrete_models.np, "sum", fail_sum)

        _, loglik = felsenstein_pruning(tree, tip_states, Q, pi, states)

        assert loglik < 0

    def test_pruning_uses_direct_tree_traversal(self, monkeypatch):
        from Bio.Phylo.BaseTree import TreeMixin
        from Bio.Phylo.Newick import Clade, Tree

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="A"),
                    Clade(branch_length=2.0, name="B"),
                ],
            )
        )

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        _, loglik = felsenstein_pruning(
            tree,
            {"A": "0", "B": "1"},
            build_q_matrix(np.array([0.1]), 2, "ER"),
            np.ones(2) / 2.0,
            ["0", "1"],
        )

        assert loglik < 0

    def test_prepare_context_uses_direct_tree_traversal(self, monkeypatch):
        from Bio.Phylo.Newick import Clade, Tree
        from Bio.Phylo.BaseTree import TreeMixin

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="A"),
                    Clade(branch_length=2.0, name="B"),
                ],
            )
        )

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        context = _prepare_felsenstein_context(
            tree, {"A": "0", "B": "1"}, ["0", "1"]
        )

        assert context["root_index"] == 2
        assert context["child_indices"][2] == (0, 1)
        assert context["branch_lengths"][2] == (1.0, 2.0)
        assert context["unique_branch_lengths"] == (1.0, 2.0)
        assert context["internal_entries_by_length_index"] == [
            (2, (0, 1), (0, 1))
        ]
        assert context["tip_state_indices"].tolist() == [0, 1, -1]
        assert len(context["internal_entries"]) == 1
        assert context["tip_liks"].tolist() == [[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]

    def test_direct_postorder_preserves_left_to_right_order(self):
        from Bio.Phylo.Newick import Clade, Tree

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(clades=[Clade(name="A"), Clade(name="B")]),
                    Clade(clades=[Clade(name="C"), Clade(name="D")]),
                ],
            )
        )

        observed = [
            clade.name if clade.name is not None else "internal"
            for clade in _postorder_clades_direct(tree)
        ]

        assert observed == [
            "A", "B", "internal", "C", "D", "internal", "internal",
        ]

    def test_prepared_loglik_uses_cached_tip_likelihoods(self):
        from Bio.Phylo.Newick import Clade, Tree

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="A"),
                    Clade(branch_length=2.0, name="B"),
                ],
            )
        )
        context = _prepare_felsenstein_context(
            tree, {"A": "0", "B": "1"}, ["0", "1"]
        )
        Q = build_q_matrix(np.array([0.1]), 2, "ER")
        pi = np.ones(2) / 2.0
        expected = _felsenstein_loglik_prepared(context, Q, pi)

        context["tip_state_indices"] = None

        assert _felsenstein_loglik_prepared(context, Q, pi) == expected

    def test_prepared_two_state_er_loglik_uses_scalar_transition_path(
        self, monkeypatch
    ):
        from Bio.Phylo.Newick import Clade, Tree
        import phykit.helpers.discrete_models as discrete_models

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(
                        branch_length=1.0,
                        clades=[
                            Clade(branch_length=0.5, name="A"),
                            Clade(branch_length=0.5, name="B"),
                        ],
                    ),
                    Clade(branch_length=2.0, name="C"),
                ],
            )
        )
        context = _prepare_felsenstein_context(
            tree, {"A": "0", "B": "1", "C": "0"}, ["0", "1"]
        )
        Q = build_q_matrix(np.array([0.1]), 2, "ER")
        pi = np.ones(2) / 2.0
        expected = _felsenstein_loglik_prepared(context, Q, pi)

        def fail_matrix_exp(*_args, **_kwargs):
            raise AssertionError("two-state ER prepared pruning should be scalar")

        monkeypatch.setattr(discrete_models, "matrix_exp", fail_matrix_exp)
        monkeypatch.setattr(
            discrete_models.np,
            "exp",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("two-state ER prepared pruning should use math.exp")
            ),
            raising=False,
        )
        monkeypatch.setattr(
            discrete_models.np,
            "log",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("two-state ER prepared pruning should use math.log")
            ),
            raising=False,
        )

        assert _felsenstein_loglik_prepared(context, Q, pi) == pytest.approx(
            expected
        )

    def test_prepared_two_state_ard_loglik_uses_scalar_transition_path(
        self, monkeypatch
    ):
        from Bio.Phylo.Newick import Clade, Tree
        import phykit.helpers.discrete_models as discrete_models

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(
                        branch_length=1.0,
                        clades=[
                            Clade(branch_length=0.5, name="A"),
                            Clade(branch_length=0.5, name="B"),
                        ],
                    ),
                    Clade(branch_length=2.0, name="C"),
                ],
            )
        )
        context = _prepare_felsenstein_context(
            tree, {"A": "0", "B": "1", "C": "0"}, ["0", "1"]
        )
        Q = build_q_matrix(np.array([0.1, 0.6]), 2, "ARD")
        pi = np.ones(2) / 2.0
        expected = _felsenstein_loglik_prepared(context, Q, pi)

        def fail_matrix_exp(*_args, **_kwargs):
            raise AssertionError("two-state ARD prepared pruning should be scalar")

        monkeypatch.setattr(discrete_models, "matrix_exp", fail_matrix_exp)
        monkeypatch.setattr(
            discrete_models.np,
            "exp",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("two-state ARD prepared pruning should use math.exp")
            ),
            raising=False,
        )
        monkeypatch.setattr(
            discrete_models.np,
            "log",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("two-state ARD prepared pruning should use math.log")
            ),
            raising=False,
        )

        assert _felsenstein_loglik_prepared(context, Q, pi) == pytest.approx(
            expected
        )

    def test_prepared_generic_root_likelihood_uses_dot_product(self, monkeypatch):
        from Bio.Phylo.Newick import Clade, Tree
        import phykit.helpers.discrete_models as discrete_models

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=0.5, name="A"),
                    Clade(branch_length=1.0, name="B"),
                    Clade(branch_length=1.5, name="C"),
                ],
            )
        )
        context = _prepare_felsenstein_context(
            tree, {"A": "0", "B": "1", "C": "2"}, ["0", "1", "2"]
        )
        Q = build_q_matrix(np.array([0.1]), 3, "ER")
        pi = np.ones(3) / 3.0

        def fail_sum(*_args, **_kwargs):
            raise AssertionError("prepared root likelihood should use np.dot")

        monkeypatch.setattr(discrete_models.np, "sum", fail_sum)

        assert _felsenstein_loglik_prepared(context, Q, pi) < 0

    def test_prepared_three_state_sym_loglik_uses_eigendecomp_transition_path(
        self, monkeypatch
    ):
        from Bio.Phylo.Newick import Clade, Tree
        import phykit.helpers.discrete_models as discrete_models

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(
                        branch_length=0.25,
                        clades=[
                            Clade(branch_length=0.5, name="A"),
                            Clade(branch_length=0.75, name="B"),
                        ],
                    ),
                    Clade(branch_length=1.25, name="C"),
                ],
            )
        )
        context = _prepare_felsenstein_context(
            tree, {"A": "0", "B": "1", "C": "2"}, ["0", "1", "2"]
        )
        Q = build_q_matrix(np.array([0.2, 0.3, 0.4]), 3, "SYM")
        pi = np.ones(3) / 3.0
        expected = _felsenstein_loglik_prepared(context, Q, pi)

        def fail_matrix_exp(*_args, **_kwargs):
            raise AssertionError(
                "small prepared pruning should use eigendecomposition"
            )

        monkeypatch.setattr(discrete_models, "matrix_exp", fail_matrix_exp)

        assert _felsenstein_loglik_prepared(context, Q, pi) == pytest.approx(
            expected,
            rel=1e-12,
            abs=1e-12,
        )

    def test_symmetric_eigendecomposition_context_avoids_inverse(
        self, monkeypatch
    ):
        import phykit.helpers.discrete_models as discrete_models

        Q = build_q_matrix(np.array([0.2, 0.3, 0.4]), 3, "SYM")

        def fail_inverse(*_args, **_kwargs):
            raise AssertionError("symmetric matrices should use eigh")

        monkeypatch.setattr(discrete_models.np.linalg, "inv", fail_inverse)
        monkeypatch.setattr(discrete_models.np.linalg, "cond", fail_inverse)

        context = discrete_models._matrix_exp_eigendecomp_context(Q)

        assert context[0] == "symmetric"
        transition = discrete_models._matrix_exp_from_eigendecomp(context, 0.75)
        np.testing.assert_allclose(transition, matrix_exp(Q, 0.75))

    def test_complex_eigendecomposition_context_falls_back(self):
        import phykit.helpers.discrete_models as discrete_models

        Q = build_q_matrix(
            np.array([0.2, 0.3, 0.4, 0.15, 0.05, 0.25]),
            3,
            "ARD",
        )

        assert discrete_models._matrix_exp_eigendecomp_context(Q) is None

    def test_two_state_er_rate_loglik_matches_prepared_q_likelihood(self, monkeypatch):
        from Bio.Phylo.Newick import Clade, Tree
        import phykit.helpers.discrete_models as discrete_models
        from phykit.helpers.discrete_models import (
            _felsenstein_loglik_two_state_er_rate,
        )

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(
                        branch_length=1.0,
                        clades=[
                            Clade(branch_length=0.5, name="A"),
                            Clade(branch_length=0.5, name="B"),
                        ],
                    ),
                    Clade(branch_length=2.0, name="C"),
                ],
            )
        )
        context = _prepare_felsenstein_context(
            tree, {"A": "0", "B": "1", "C": "0"}, ["0", "1"]
        )
        rate = 0.17
        pi = np.ones(2) / 2.0
        Q = build_q_matrix(np.array([rate]), 2, "ER")
        expected = _felsenstein_loglik_prepared(context, Q, pi)

        monkeypatch.setattr(
            discrete_models.np,
            "exp",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("two-state ER rate objective should use math.exp")
            ),
            raising=False,
        )
        monkeypatch.setattr(
            discrete_models.np,
            "log",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("two-state ER rate objective should use math.log")
            ),
            raising=False,
        )

        assert _felsenstein_loglik_two_state_er_rate(
            context,
            rate,
            pi,
        ) == pytest.approx(expected)

    def test_multi_state_er_rate_loglik_matches_prepared_q_likelihood(
        self, monkeypatch
    ):
        from Bio.Phylo.Newick import Clade, Tree
        import phykit.helpers.discrete_models as discrete_models

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(
                        branch_length=1.0,
                        clades=[
                            Clade(branch_length=0.5, name="A"),
                            Clade(branch_length=0.75, name="B"),
                        ],
                    ),
                    Clade(branch_length=2.0, name="C"),
                ],
            )
        )
        context = _prepare_felsenstein_context(
            tree, {"A": "0", "B": "1", "C": "2"}, ["0", "1", "2"]
        )
        rate = 0.17
        pi = np.ones(3) / 3.0
        Q = build_q_matrix(np.array([rate]), 3, "ER")
        expected = _felsenstein_loglik_prepared(context, Q, pi)

        def fail_matrix_exp(*_args, **_kwargs):
            raise AssertionError(
                "multi-state ER rate likelihood should use scalar transitions"
            )

        monkeypatch.setattr(discrete_models, "matrix_exp", fail_matrix_exp)
        monkeypatch.setattr(
            discrete_models.np,
            "exp",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("multi-state ER rate likelihood should use math.exp")
            ),
            raising=False,
        )

        assert _felsenstein_loglik_er_rate(
            context,
            rate,
            pi,
        ) == pytest.approx(expected)

    def test_two_state_er_fit_avoids_objective_q_rebuilds(self, monkeypatch):
        from Bio.Phylo.Newick import Clade, Tree
        import phykit.helpers.discrete_models as discrete_models

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(
                        branch_length=1.0,
                        clades=[
                            Clade(branch_length=0.5, name="A"),
                            Clade(branch_length=0.5, name="B"),
                        ],
                    ),
                    Clade(branch_length=2.0, name="C"),
                ],
            )
        )
        calls = 0
        original = discrete_models.build_q_matrix

        def counting_build_q_matrix(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original(*args, **kwargs)

        monkeypatch.setattr(
            discrete_models,
            "build_q_matrix",
            counting_build_q_matrix,
        )

        Q, loglik = fit_q_matrix(
            tree,
            {"A": "0", "B": "1", "C": "0"},
            ["0", "1"],
            "ER",
        )

        assert calls == 1
        assert Q.shape == (2, 2)
        assert np.isfinite(loglik)

    def test_multi_state_er_fit_avoids_objective_q_rebuilds(self, monkeypatch):
        from Bio.Phylo.Newick import Clade, Tree
        import phykit.helpers.discrete_models as discrete_models

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(
                        branch_length=1.0,
                        clades=[
                            Clade(branch_length=0.5, name="A"),
                            Clade(branch_length=0.5, name="B"),
                        ],
                    ),
                    Clade(branch_length=2.0, name="C"),
                ],
            )
        )
        calls = 0
        original = discrete_models.build_q_matrix

        def counting_build_q_matrix(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original(*args, **kwargs)

        monkeypatch.setattr(
            discrete_models,
            "build_q_matrix",
            counting_build_q_matrix,
        )

        Q, loglik = fit_q_matrix(
            tree,
            {"A": "0", "B": "1", "C": "2"},
            ["0", "1", "2"],
            "ER",
        )

        assert calls == 1
        assert Q.shape == (3, 3)
        assert np.isfinite(loglik)

    def test_two_state_er_fit_uses_scalar_optimizer(self, monkeypatch):
        from Bio.Phylo.Newick import Clade, Tree
        import phykit.helpers.discrete_models as discrete_models

        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="A"),
                    Clade(branch_length=1.0, name="B"),
                    Clade(branch_length=1.0, name="C"),
                ],
            )
        )
        scalar_calls = 0
        original_minimize_scalar = discrete_models.minimize_scalar

        def fail_minimize(*_args, **_kwargs):
            raise AssertionError("two-state ER should use scalar optimization")

        def counting_minimize_scalar(*args, **kwargs):
            nonlocal scalar_calls
            scalar_calls += 1
            return original_minimize_scalar(*args, **kwargs)

        monkeypatch.setattr(discrete_models, "minimize", fail_minimize)
        monkeypatch.setattr(
            discrete_models,
            "minimize_scalar",
            counting_minimize_scalar,
        )

        Q, loglik = fit_q_matrix(
            tree,
            {"A": "0", "B": "1", "C": "0"},
            ["0", "1"],
            "ER",
        )

        assert scalar_calls == 1
        assert Q.shape == (2, 2)
        assert np.isfinite(loglik)


class TestParseDiscreteTraits:
    def test_parse_multi_column(self, tmp_path):
        f = tmp_path / "traits.tsv"
        f.write_text("taxon\tdiet\nA\tcarnivore\nB\therbivore\nC\tomnivore\n")
        result = parse_discrete_traits(str(f), ["A", "B", "C", "D"], trait_column="diet")
        assert result == {"A": "carnivore", "B": "herbivore", "C": "omnivore"}

    def test_parse_multi_column_all_shared_emits_no_warnings(self, tmp_path, capsys):
        f = tmp_path / "traits.tsv"
        f.write_text("taxon\tdiet\nA\tcarnivore\nB\therbivore\nC\tomnivore\n")

        result = parse_discrete_traits(
            str(f), ["A", "B", "C"], trait_column="diet"
        )

        assert result == {"A": "carnivore", "B": "herbivore", "C": "omnivore"}
        assert capsys.readouterr().err == ""

    def test_parse_multi_column_ordered_all_shared_skips_set_validation(
        self, tmp_path, monkeypatch
    ):
        f = tmp_path / "traits.tsv"
        f.write_text("taxon\tdiet\nA\tcarnivore\nB\therbivore\nC\tomnivore\n")

        def fail_set(*_args, **_kwargs):
            raise AssertionError(
                "ordered all-shared validation should not build a set"
            )

        monkeypatch.setattr("builtins.set", fail_set)

        assert parse_discrete_traits(
            str(f), ["A", "B", "C"], trait_column="diet"
        ) == {"A": "carnivore", "B": "herbivore", "C": "omnivore"}

    def test_parse_multi_column_streams_without_readlines(self, monkeypatch):
        class StreamingOnlyFile:
            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, traceback):
                return False

            def __iter__(self):
                return iter(
                    [
                        "   # ignored before header\n",
                        "\n",
                        "taxon\tdiet\thabitat\n",
                        "A\tcarnivore\tforest\n",
                        "\t# ignored between rows\n",
                        "B\therbivore\tplain\n",
                        "C\tomnivore\twetland\n",
                    ]
                )

            def readlines(self):
                raise AssertionError("multi-column trait parsing should stream rows")

        monkeypatch.setattr("builtins.open", lambda path: StreamingOnlyFile())

        result = parse_discrete_traits(
            "traits.tsv", ["A", "B", "C", "D"], trait_column="diet"
        )

        assert result == {"A": "carnivore", "B": "herbivore", "C": "omnivore"}

    def test_parse_multi_column_preserves_logical_line_numbers(self, tmp_path):
        f = tmp_path / "traits.tsv"
        f.write_text(
            "\n".join(
                [
                    "# ignored before header",
                    "taxon\tdiet\thabitat",
                    "A\tcarnivore\tforest",
                    "# ignored between rows",
                    "B\therbivore",
                ]
            )
            + "\n"
        )

        with pytest.raises(PhykitUserError) as excinfo:
            parse_discrete_traits(str(f), ["A", "B", "C"], trait_column="diet")

        assert "Line 3 has 2 columns; expected 3." in excinfo.value.messages

    def test_parse_multi_column_selected_early_validates_trailing_columns(self, tmp_path):
        f = tmp_path / "traits.tsv"
        f.write_text(
            "\n".join(
                [
                    "taxon\tdiet\thabitat\tregion",
                    "A\tcarnivore\tforest\tnorth",
                    "B\therbivore\tplain\tsouth\textra",
                    "C\tomnivore\twetland\teast",
                ]
            )
            + "\n"
        )

        with pytest.raises(PhykitUserError) as excinfo:
            parse_discrete_traits(str(f), ["A", "B", "C"], trait_column="diet")

        assert "Line 3 has 5 columns; expected 4." in excinfo.value.messages

    def test_parse_two_column(self, tmp_path):
        f = tmp_path / "traits.tsv"
        f.write_text("A\tcarnivore\nB\therbivore\nC\tomnivore\n")
        result = parse_discrete_traits(str(f), ["A", "B", "C", "D"])
        assert result == {"A": "carnivore", "B": "herbivore", "C": "omnivore"}

    def test_parse_two_column_all_shared_emits_no_warnings(self, tmp_path, capsys):
        f = tmp_path / "traits.tsv"
        f.write_text("A\tcarnivore\nB\therbivore\nC\tomnivore\n")

        result = parse_discrete_traits(str(f), ["A", "B", "C"])

        assert result == {"A": "carnivore", "B": "herbivore", "C": "omnivore"}
        assert capsys.readouterr().err == ""

    def test_parse_two_column_ordered_all_shared_skips_set_validation(
        self, tmp_path, monkeypatch
    ):
        f = tmp_path / "traits.tsv"
        f.write_text("A\tcarnivore\nB\therbivore\nC\tomnivore\n")

        def fail_set(*_args, **_kwargs):
            raise AssertionError(
                "ordered all-shared validation should not build a set"
            )

        monkeypatch.setattr("builtins.set", fail_set)

        assert parse_discrete_traits(str(f), ["A", "B", "C"]) == {
            "A": "carnivore",
            "B": "herbivore",
            "C": "omnivore",
        }

    def test_parse_two_column_skips_comments_and_preserves_logical_line_numbers(
        self, tmp_path
    ):
        f = tmp_path / "traits.tsv"
        f.write_text(
            "\n".join(
                [
                    "   # ignored before data",
                    "A\tcarnivore",
                    "\t# ignored between rows",
                    "B\therbivore",
                    "C",
                ]
            )
            + "\n"
        )

        with pytest.raises(PhykitUserError) as excinfo:
            parse_discrete_traits(str(f), ["A", "B", "C"])

        assert (
            "Line 3 has 1 columns; expected 2 (taxon, state)."
            in excinfo.value.messages
        )

    def test_parse_two_column_rejects_extra_columns(self, tmp_path):
        f = tmp_path / "traits.tsv"
        f.write_text("A\tcarnivore\nB\therbivore\textra\nC\tomnivore\n")

        with pytest.raises(PhykitUserError) as excinfo:
            parse_discrete_traits(str(f), ["A", "B", "C"])

        assert (
            "Line 2 has 3 columns; expected 2 (taxon, state)."
            in excinfo.value.messages
        )

    def test_missing_file_exits(self):
        with pytest.raises(SystemExit):
            parse_discrete_traits("/nonexistent.tsv", ["A", "B", "C"])

    def test_too_few_shared_taxa_exits(self, tmp_path):
        f = tmp_path / "traits.tsv"
        f.write_text("A\tcarnivore\nB\therbivore\n")
        with pytest.raises(SystemExit):
            parse_discrete_traits(str(f), ["X", "Y", "Z"])
