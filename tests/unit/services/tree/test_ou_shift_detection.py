import builtins
import copy
import importlib
import pickle as stdlib_pickle
import subprocess
import sys
from io import StringIO

import numpy as np
import pytest
from argparse import Namespace
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, TreeMixin

import phykit.services.tree.ou_shift_detection as ou_shift_detection_module
from phykit.services.tree.ou_shift_detection import OUShiftDetection
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")


def test_module_import_does_not_import_numpy_scipy_or_sklearn():
    code = """
import sys
import phykit.services.tree.ou_shift_detection as module

assert callable(module.print_json)
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "pickle" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy.linalg" not in sys.modules
assert "scipy.optimize" not in sys.modules
assert "sklearn" not in sys.modules
assert "sklearn.linear_model" not in sys.modules
assert module._CHO_FACTOR is None
assert module._CHO_SOLVE is None
assert module._SOLVE_TRIANGULAR is None
assert module._MINIMIZE_SCALAR is None
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_module_and_resolved_attributes():
    lazy_np = ou_shift_detection_module._LazyNumpy()

    sqrt = lazy_np.sqrt

    assert lazy_np._module is not None
    assert lazy_np.sqrt is sqrt


def test_lazy_pickle_caches_resolved_copy_helpers(monkeypatch):
    lazy_pickle = ou_shift_detection_module._LazyPickle()

    def cached_dumps(value, **_kwargs):
        return f"cached:{value}".encode("ascii")

    def cached_loads(value):
        return value.decode("ascii").removeprefix("cached:")

    def uncached_dumps(*_args, **_kwargs):
        raise AssertionError("cached dumps should be reused")

    def uncached_loads(*_args, **_kwargs):
        raise AssertionError("cached loads should be reused")

    monkeypatch.setattr(stdlib_pickle, "dumps", cached_dumps)
    monkeypatch.setattr(stdlib_pickle, "loads", cached_loads)

    protocol = lazy_pickle.HIGHEST_PROTOCOL
    assert lazy_pickle.loads(lazy_pickle.dumps("tree", protocol=protocol)) == "tree"

    monkeypatch.setattr(stdlib_pickle, "dumps", uncached_dumps)
    monkeypatch.setattr(stdlib_pickle, "loads", uncached_loads)

    assert lazy_pickle.loads(lazy_pickle.dumps("tree2", protocol=protocol)) == "tree2"
    assert lazy_pickle.__dict__["dumps"] is cached_dumps
    assert lazy_pickle.__dict__["loads"] is cached_loads
    assert lazy_pickle.__dict__["HIGHEST_PROTOCOL"] == stdlib_pickle.HIGHEST_PROTOCOL


def test_lazy_pickle_loads_initializes_matching_dumps_helper():
    lazy_pickle = ou_shift_detection_module._LazyPickle()
    payload = stdlib_pickle.dumps({"tree": "value"})

    assert lazy_pickle.loads(payload) == {"tree": "value"}
    assert lazy_pickle.dumps is stdlib_pickle.dumps


def test_module_import_does_not_import_heavy_optional_packages(monkeypatch):
    module_name = "phykit.services.tree.ou_shift_detection"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if (
            name == "sklearn.linear_model"
            or name.startswith("sklearn.linear_model.")
            or name == "scipy.linalg"
            or name.startswith("scipy.linalg.")
            or name == "scipy.optimize"
            or name.startswith("scipy.optimize.")
        ):
            raise AssertionError(
                "ou_shift_detection module import should not import sklearn or SciPy linalg/optimize"
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


def test_repeated_cholesky_calls_cache_scipy_linalg_imports(monkeypatch):
    previous_cho_factor = ou_shift_detection_module._CHO_FACTOR
    previous_cho_solve = ou_shift_detection_module._CHO_SOLVE
    ou_shift_detection_module._CHO_FACTOR = None
    ou_shift_detection_module._CHO_SOLVE = None
    original_import = builtins.__import__
    scipy_linalg_imports = 0

    def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
        nonlocal scipy_linalg_imports
        if name == "scipy.linalg":
            scipy_linalg_imports += 1
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", counting_import)
    try:
        rng = np.random.default_rng(20260628)
        svc = OUShiftDetection.__new__(OUShiftDetection)
        n = 12
        A = rng.normal(size=(n, n))
        V = A @ A.T + np.eye(n)
        W = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
        x = rng.normal(size=n)

        svc._gls_profile_likelihood_cholesky(x, W, V)
        first_call_imports = scipy_linalg_imports
        svc._gls_profile_likelihood_cholesky(x, W, V)
    finally:
        ou_shift_detection_module._CHO_FACTOR = previous_cho_factor
        ou_shift_detection_module._CHO_SOLVE = previous_cho_solve

    assert first_call_imports > 0
    assert scipy_linalg_imports == first_call_imports


def test_repeated_solve_triangular_calls_cache_scipy_linalg_imports(monkeypatch):
    previous_solve_triangular = ou_shift_detection_module._SOLVE_TRIANGULAR
    ou_shift_detection_module._SOLVE_TRIANGULAR = None
    original_import = builtins.__import__
    scipy_linalg_imports = 0

    def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
        nonlocal scipy_linalg_imports
        if name == "scipy.linalg":
            scipy_linalg_imports += 1
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", counting_import)
    try:
        matrix = np.eye(3)
        values = np.array([1.0, 2.0, 3.0])

        first = ou_shift_detection_module.solve_triangular(
            matrix, values, lower=True, check_finite=False
        )
        first_call_imports = scipy_linalg_imports
        second = ou_shift_detection_module.solve_triangular(
            matrix, values, lower=True, check_finite=False
        )
    finally:
        ou_shift_detection_module._SOLVE_TRIANGULAR = previous_solve_triangular

    np.testing.assert_allclose(first, values)
    np.testing.assert_allclose(second, values)
    assert first_call_imports > 0
    assert scipy_linalg_imports == first_call_imports


def test_repeated_minimize_scalar_calls_cache_scipy_optimize_imports(monkeypatch):
    previous_minimize_scalar = ou_shift_detection_module._MINIMIZE_SCALAR
    ou_shift_detection_module._MINIMIZE_SCALAR = None
    original_import = builtins.__import__
    scipy_optimize_imports = 0

    def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
        nonlocal scipy_optimize_imports
        if name == "scipy.optimize":
            scipy_optimize_imports += 1
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", counting_import)
    try:
        first = ou_shift_detection_module.minimize_scalar(
            lambda value: (value - 0.25) ** 2,
            bounds=(0.0, 1.0),
            method="bounded",
        )
        first_call_imports = scipy_optimize_imports
        second = ou_shift_detection_module.minimize_scalar(
            lambda value: (value - 0.75) ** 2,
            bounds=(0.0, 1.0),
            method="bounded",
        )
    finally:
        ou_shift_detection_module._MINIMIZE_SCALAR = previous_minimize_scalar

    assert first.x == pytest.approx(0.25, abs=1e-4)
    assert second.x == pytest.approx(0.75, abs=1e-4)
    assert first_call_imports > 0
    assert scipy_optimize_imports == first_call_imports


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        criterion="pBIC",
        max_shifts=None,
        json=False,
    )


@pytest.fixture
def svc(default_args):
    return OUShiftDetection(default_args)


@pytest.fixture
def precomputed(svc):
    """Pre-compute common data structures for testing."""
    tree = svc.read_tree_file()
    tree_tips = svc.get_tip_names_from_tree(tree)
    traits = svc._parse_trait_file(TRAITS_FILE, tree_tips)

    shared = set(traits.keys())
    tree_copy = copy.deepcopy(tree)
    tip_names_in_tree = [t.name for t in tree_copy.get_terminals()]
    to_prune = [t for t in tip_names_in_tree if t not in shared]
    if to_prune:
        for t in to_prune:
            tree_copy.prune(t)

    ordered_names = sorted(shared)
    x = np.array([traits[name] for name in ordered_names])

    parent_map = svc._build_parent_map(tree_copy)
    lineage_info = svc._build_lineage_info_no_regime(
        tree_copy, ordered_names, parent_map
    )
    tree_height = max(
        sum(bl for _, bl, _, _ in lineage_info[name])
        for name in ordered_names
    )

    edges = svc._enumerate_edges(tree_copy, parent_map)
    descendant_counts = svc._count_descendants(tree_copy, edges)

    # Set up instance attributes needed by new methods
    svc._ordered_names = ordered_names
    svc._x = x
    svc._n = len(ordered_names)
    svc._lineage_info = lineage_info
    svc._tree_height = tree_height
    svc._edges = edges
    svc._n_edges = len(edges)
    svc._parent_map = parent_map
    svc._S, svc._tip_heights_arr = svc._precompute_shared_path_lengths()

    return dict(
        svc=svc, tree=tree_copy, x=x, ordered_names=ordered_names,
        lineage_info=lineage_info, tree_height=tree_height,
        parent_map=parent_map, edges=edges,
        descendant_counts=descendant_counts,
    )


class TestEnumerateEdges:
    def test_iter_postorder_matches_biopython_order(self, svc):
        tree = Phylo.read(
            StringIO("((A:1,B:2,C:3):4,(D:5,E:6):7,F:8):0;"),
            "newick",
        )

        direct = list(svc._iter_postorder(tree.root))
        reference = list(tree.find_clades(order="postorder"))

        assert [id(clade) for clade in direct] == [
            id(clade) for clade in reference
        ]

    def test_setup_helpers_use_direct_tree_traversal(self, svc, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]

        def fail_traversal(*args, **kwargs):
            raise AssertionError("optimized path should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        parent_map = svc._build_parent_map(tree)
        lineage = svc._build_lineage_info_no_regime(tree, ordered_names, parent_map)
        edges = svc._enumerate_edges(tree, parent_map)
        counts = svc._count_descendants(tree, edges)

        assert len(parent_map) == 6
        assert len(lineage) == 4
        assert len(edges) == 6
        assert sum(counts.values()) == 8
        assert sum(bl for _, bl, _, _ in lineage["A"]) == pytest.approx(2.0)

    def test_parent_map_handles_mixed_child_counts_without_preorder(
        self, svc, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_iter_preorder(*args, **kwargs):
            raise AssertionError("parent map should build directly")

        monkeypatch.setattr(svc, "_iter_preorder", fail_iter_preorder)

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

    def test_iter_preorder_preserves_order_without_reversed(self, svc):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("_iter_preorder should push children directly")

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

    def test_correct_count_excludes_root(self, precomputed):
        edges = precomputed["edges"]
        tree = precomputed["tree"]
        all_clades = list(tree.find_clades())
        expected = len(all_clades) - 1
        assert len(edges) == expected

    def test_root_not_in_edges(self, precomputed):
        edges = precomputed["edges"]
        tree = precomputed["tree"]
        root_id = id(tree.root)
        edge_ids = [eid for eid, _ in edges]
        assert root_id not in edge_ids

    def test_all_edges_have_valid_clades(self, precomputed):
        edges = precomputed["edges"]
        for clade_id, clade in edges:
            assert clade_id == id(clade)
            assert clade is not None


class TestCountDescendants:
    def test_mixed_child_counts(self, svc):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )
        parent_map = svc._build_parent_map(tree)
        edges = svc._enumerate_edges(tree, parent_map)

        counts = svc._count_descendants(tree, edges)

        terminal, binary, trifurcating = tree.root.clades
        assert counts[id(terminal)] == 1
        assert counts[id(binary)] == 2
        assert counts[id(trifurcating)] == 3
        assert all(counts[id(child)] == 1 for child in binary.clades)
        assert all(counts[id(child)] == 1 for child in trifurcating.clades)

    def test_terminal_edges_have_one_descendant(self, precomputed):
        edges = precomputed["edges"]
        descendant_counts = precomputed["descendant_counts"]
        for clade_id, clade in edges:
            if clade.is_terminal():
                assert descendant_counts[clade_id] == 1

    def test_internal_edges_have_multiple_descendants(self, precomputed):
        edges = precomputed["edges"]
        descendant_counts = precomputed["descendant_counts"]
        for clade_id, clade in edges:
            if not clade.is_terminal():
                assert descendant_counts[clade_id] >= 2


class TestBuildShiftWeightMatrix:
    def test_shape(self, precomputed):
        svc = precomputed["svc"]
        ordered_names = precomputed["ordered_names"]
        lineage_info = precomputed["lineage_info"]
        edges = precomputed["edges"]
        tree_height = precomputed["tree_height"]
        n = len(ordered_names)
        p = len(edges)

        W = svc._build_shift_weight_matrix(
            ordered_names, lineage_info, edges, 1.0, tree_height
        )
        assert W.shape == (n, p + 1)

    def test_rows_sum_to_one(self, precomputed):
        svc = precomputed["svc"]
        ordered_names = precomputed["ordered_names"]
        lineage_info = precomputed["lineage_info"]
        edges = precomputed["edges"]
        tree_height = precomputed["tree_height"]

        W = svc._build_shift_weight_matrix(
            ordered_names, lineage_info, edges, 1.0, tree_height
        )
        row_sums = W.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)

    def test_baseline_weight_uses_array_sum(self, precomputed, monkeypatch):
        svc = precomputed["svc"]
        ordered_names = precomputed["ordered_names"]
        lineage_info = precomputed["lineage_info"]
        edges = precomputed["edges"]
        tree_height = precomputed["tree_height"]

        def fail_sum(*_args, **_kwargs):
            raise AssertionError("baseline shift weights should use ndarray.sum")

        monkeypatch.setattr(ou_shift_detection_module.np, "sum", fail_sum)

        W = svc._build_shift_weight_matrix(
            ordered_names, lineage_info, edges, 1.0, tree_height
        )
        np.testing.assert_allclose(W.sum(axis=1), 1.0, atol=1e-10)

    def test_all_weights_nonnegative(self, precomputed):
        svc = precomputed["svc"]
        ordered_names = precomputed["ordered_names"]
        lineage_info = precomputed["lineage_info"]
        edges = precomputed["edges"]
        tree_height = precomputed["tree_height"]

        W = svc._build_shift_weight_matrix(
            ordered_names, lineage_info, edges, 1.0, tree_height
        )
        assert np.all(W >= -1e-15)


class TestBuildIndicatorDesignMatrix:
    def test_shape(self, precomputed):
        svc = precomputed["svc"]
        n = svc._n
        # No shifts: should be n x 1
        W = svc._build_indicator_design_matrix([])
        assert W.shape == (n, 1)
        np.testing.assert_allclose(W[:, 0], 1.0)

    def test_shape_with_shifts(self, precomputed):
        svc = precomputed["svc"]
        n = svc._n
        W = svc._build_indicator_design_matrix([0, 1])
        assert W.shape == (n, 3)

    def test_intercept_is_ones(self, precomputed):
        svc = precomputed["svc"]
        W = svc._build_indicator_design_matrix([0, 1, 2])
        np.testing.assert_allclose(W[:, 0], 1.0)

    def test_values_are_binary(self, precomputed):
        svc = precomputed["svc"]
        edges = precomputed["edges"]
        all_indices = list(range(len(edges)))
        W = svc._build_indicator_design_matrix(all_indices)
        # All values should be 0 or 1
        for col in range(1, W.shape[1]):
            unique_vals = np.unique(W[:, col])
            assert all(v in (0.0, 1.0) for v in unique_vals)

    def test_terminal_edge_marks_one_tip(self, precomputed):
        svc = precomputed["svc"]
        edges = precomputed["edges"]
        # Find a terminal edge
        for j, (clade_id, clade) in enumerate(edges):
            if clade.is_terminal():
                W = svc._build_indicator_design_matrix([j])
                # Only one tip should have a 1 in the shift column
                assert W[:, 1].sum() == 1.0
                break

    def test_uses_cached_lineage_rows(self, precomputed):
        svc = precomputed["svc"]
        shift_indices = list(range(min(3, len(svc._edges))))

        expected = np.zeros((svc._n, len(shift_indices) + 1))
        expected[:, 0] = 1.0
        shift_clade_ids = {
            svc._edges[edge_idx][0]: col_idx + 1
            for col_idx, edge_idx in enumerate(shift_indices)
        }
        for row_idx, name in enumerate(svc._ordered_names):
            for clade_id, *_ in svc._lineage_info[name]:
                if clade_id in shift_clade_ids:
                    expected[row_idx, shift_clade_ids[clade_id]] = 1.0

        W = svc._build_indicator_design_matrix(shift_indices)

        np.testing.assert_array_equal(W, expected)
        assert svc._lineage_rows_by_clade_id is not None


class TestExtractLassoConfigs:
    def test_column_l2_norms_matches_linalg_norm(self):
        matrix = np.array(
            [
                [3.0, 1.0, 0.0],
                [4.0, 2.0, 0.0],
                [0.0, 2.0, 0.0],
            ]
        )

        observed = ou_shift_detection_module._column_l2_norms(matrix)
        expected = np.linalg.norm(matrix, axis=0)

        np.testing.assert_allclose(observed, expected)

    def test_uses_flat_nonzero_for_coefficient_path(self, monkeypatch):
        svc = OUShiftDetection.__new__(OUShiftDetection)
        svc._n = 5
        X_star = np.column_stack(
            [
                np.ones(5),
                np.array([0.0, 1.0, 0.0, 1.0, 0.0]),
                np.array([1.0, 0.0, 1.0, 0.0, 1.0]),
                np.array([0.0, 0.0, 1.0, 1.0, 0.0]),
            ]
        )
        y_star = np.array([0.2, 1.0, 0.4, 1.2, 0.6])
        coefs_path = np.array(
            [
                [0.0, 0.5, 0.5, 0.0],
                [0.0, 0.0, 0.25, 0.25],
                [0.0, 0.0, 0.0, 0.5],
            ]
        )

        def fake_lars_path(_X, _y, method, max_iter):
            assert method == "lasso"
            assert max_iter == 3
            return None, None, coefs_path

        def fail_where(*_args, **_kwargs):
            raise AssertionError("lasso config extraction should use flat indices")

        monkeypatch.setattr(ou_shift_detection_module, "_lars_path", fake_lars_path)
        monkeypatch.setattr(ou_shift_detection_module.np, "where", fail_where)

        configs = svc._extract_lasso_configs(X_star, y_star, max_shifts=3)

        assert configs == [[0], [0, 1], [1, 2]]

    def test_extract_lasso_configs_avoids_linalg_norm(self, monkeypatch):
        svc = OUShiftDetection.__new__(OUShiftDetection)
        svc._n = 5
        X_star = np.column_stack(
            [
                np.ones(5),
                np.array([0.0, 1.0, 0.0, 1.0, 0.0]),
                np.array([1.0, 0.0, 1.0, 0.0, 1.0]),
            ]
        )
        y_star = np.array([0.2, 1.0, 0.4, 1.2, 0.6])
        coefs_path = np.array(
            [
                [0.0, 0.5, 0.5],
                [0.0, 0.0, 0.25],
            ]
        )

        def fake_lars_path(_X, _y, method, max_iter):
            assert method == "lasso"
            assert max_iter == 2
            return None, None, coefs_path

        def fail_norm(*_args, **_kwargs):
            raise AssertionError("LASSO column scaling should use column L2 helper")

        monkeypatch.setattr(ou_shift_detection_module, "_lars_path", fake_lars_path)
        monkeypatch.setattr(ou_shift_detection_module.np.linalg, "norm", fail_norm)

        configs = svc._extract_lasso_configs(X_star, y_star, max_shifts=2)

        assert configs == [[0], [0, 1]]


class TestComputePBIC:
    def test_returns_finite_for_valid_model(self, precomputed):
        """pBIC should return a finite value for a valid model fit."""
        svc = precomputed["svc"]
        n = svc._n
        x = svc._x

        # Fit null model and check pBIC is finite
        alpha = 1.0
        V = svc._build_ou_vcv_fast(alpha)
        V_inv = np.linalg.inv(V)
        W = np.ones((n, 1))
        theta = svc._gls_theta_hat(x, W, V_inv)
        e = x - W @ theta
        sig2 = float(e @ V_inv @ e) / n
        sign, logdet = np.linalg.slogdet(V)
        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)

        pbic = svc._compute_pbic(n, ll, 0, W, V_inv, sig2, x)
        assert np.isfinite(pbic)

    def test_vcv_pbic_matches_inverse_pbic(self, precomputed):
        svc = precomputed["svc"]
        n = svc._n
        x = svc._x

        alpha = 1.0
        V = svc._build_ou_vcv_fast(alpha)
        V_inv = np.linalg.inv(V)
        W = svc._build_indicator_design_matrix([0, 1])
        theta = svc._gls_theta_hat(x, W, V_inv)
        e = x - W @ theta
        sig2 = float(e @ V_inv @ e) / n
        sign, logdet = np.linalg.slogdet(V)
        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)

        expected = svc._compute_pbic(n, ll, 2, W, V_inv, sig2, x)
        observed = svc._compute_pbic_from_vcv(n, ll, 2, W, V, sig2, x)

        np.testing.assert_allclose(observed, expected, rtol=1e-10, atol=1e-10)

    def test_pbic_from_info_matches_scaled_inverse_formula(self):
        rng = np.random.default_rng(20260628)
        svc = OUShiftDetection.__new__(OUShiftDetection)
        svc._n_edges = 80
        n = 60
        k = 4
        d = k + 1
        x = rng.normal(size=n)
        A = rng.normal(size=(d, d))
        info = A @ A.T + np.eye(d)
        sigma2 = 1.25
        ll = -12.5

        var_y = np.var(x, ddof=1)
        scaled = sigma2 * np.linalg.inv(info) * (n - d) / (var_y * n)
        _, ld = np.linalg.slogdet(scaled)
        expected = (
            2.0 * k * np.log(svc._n_edges - 1)
            - 2.0 * ll
            + 2.0 * np.log(n)
            - ld
        )

        observed = svc._compute_pbic_from_info(n, ll, k, info, sigma2, x)

        assert observed == pytest.approx(expected)

    def test_pbic_from_info_avoids_inverse(self, monkeypatch):
        rng = np.random.default_rng(20260628)
        svc = OUShiftDetection.__new__(OUShiftDetection)
        svc._n_edges = 80
        n = 60
        k = 3
        d = k + 1
        x = rng.normal(size=n)
        A = rng.normal(size=(d, d))
        info = A @ A.T + np.eye(d)

        def fail_inverse(_matrix):
            raise AssertionError("pBIC should not invert the information matrix")

        monkeypatch.setattr(np.linalg, "inv", fail_inverse)

        observed = svc._compute_pbic_from_info(n, -12.5, k, info, 1.25, x)

        assert np.isfinite(observed)

    def test_more_shifts_higher_penalty(self, precomputed):
        """pBIC penalty should increase with number of shifts."""
        svc = precomputed["svc"]
        n = svc._n
        x = svc._x

        # Create two models: 0 shifts and 2 shifts (with same LL)
        alpha = 1.0
        V = svc._build_ou_vcv_fast(alpha)
        V_inv = np.linalg.inv(V)

        # Null model
        W0 = np.ones((n, 1))
        theta0 = svc._gls_theta_hat(x, W0, V_inv)
        e0 = x - W0 @ theta0
        sig2_0 = float(e0 @ V_inv @ e0) / n
        sign, logdet = np.linalg.slogdet(V)
        ll_0 = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2_0) + logdet + n)

        # 2-shift model (using first 2 edges)
        edges2 = svc._edges[:2]
        W2 = svc._build_shift_weight_matrix(
            svc._ordered_names, svc._lineage_info, edges2, alpha, svc._tree_height
        )
        theta2 = svc._gls_theta_hat(x, W2, V_inv)
        e2 = x - W2 @ theta2
        sig2_2 = float(e2 @ V_inv @ e2) / n
        ll_2 = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2_2) + logdet + n)

        pbic_0 = svc._compute_pbic(n, ll_0, 0, W0, V_inv, sig2_0, x)
        pbic_2 = svc._compute_pbic(n, ll_2, 2, W2, V_inv, sig2_2, x)

        # Both should be finite
        assert np.isfinite(pbic_0)
        assert np.isfinite(pbic_2)


class TestComputeBIC:
    def test_value(self, svc):
        n = 10
        ll = -20.0
        k = 2
        bic = svc._compute_bic(n, ll, k)
        # R-style: -2*LL + (2k+3)*log(n)
        expected = 40.0 + 7 * np.log(10)
        np.testing.assert_allclose(bic, expected, atol=1e-10)


class TestComputeAICc:
    def test_value(self, svc):
        n = 20
        ll = -30.0
        k = 2
        aicc = svc._compute_aicc(n, ll, k)
        # R-style: p = 2k+3 = 7
        p = 7
        expected = 60.0 + 2.0 * p + 2.0 * p * 8 / (20 - p - 1)
        np.testing.assert_allclose(aicc, expected, atol=1e-10)

    def test_small_n_gives_inf(self, svc):
        # p = 2*2+3 = 7, need n > p+1 = 8
        aicc = svc._compute_aicc(3, -10.0, 1)
        assert aicc == float("inf")


class TestTransformToIndependent:
    def test_identity_covariance_after_transform(self, precomputed):
        svc = precomputed["svc"]
        ordered_names = precomputed["ordered_names"]
        lineage_info = precomputed["lineage_info"]
        edges = precomputed["edges"]
        tree_height = precomputed["tree_height"]
        x = precomputed["x"]

        alpha = 1.0
        V = svc._build_ou_vcv_single_alpha_no_regime(
            ordered_names, lineage_info, alpha, 1.0
        )
        W = svc._build_shift_weight_matrix(
            ordered_names, lineage_info, edges, alpha, tree_height
        )
        X_star, y_star = svc._transform_to_independent(x, W, V)

        L = np.linalg.cholesky(V)
        L_inv = np.linalg.inv(L)
        transformed_V = L_inv @ V @ L_inv.T
        np.testing.assert_allclose(transformed_V, np.eye(len(x)), atol=1e-10)

    def test_regularizes_positive_semidefinite_covariance(self, svc):
        x = np.array([1.0, 2.0])
        W = np.ones((2, 1))
        V = np.ones((2, 2))

        X_star, y_star = svc._transform_to_independent(x, W, V)

        assert np.all(np.isfinite(X_star))
        assert np.all(np.isfinite(y_star))


class TestGLSProfileLikelihood:
    def test_cholesky_matches_inverse(self, precomputed):
        svc = precomputed["svc"]
        x = precomputed["x"]
        V = svc._build_ou_vcv_fast(alpha=1.0)
        W = svc._build_indicator_design_matrix([0, 1])

        theta_fast, sig2_fast, ll_fast = svc._gls_profile_likelihood_cholesky(
            x, W, V
        )
        theta_inv, sig2_inv, ll_inv = svc._gls_profile_likelihood_inverse(
            x, W, V
        )

        np.testing.assert_allclose(theta_fast, theta_inv)
        np.testing.assert_allclose(sig2_fast, sig2_inv)
        np.testing.assert_allclose(ll_fast, ll_inv)

    def test_cholesky_uses_single_solve(self, precomputed, monkeypatch):
        svc = precomputed["svc"]
        x = precomputed["x"]
        V = svc._build_ou_vcv_fast(alpha=1.0)
        W = svc._build_indicator_design_matrix([0, 1])
        original = ou_shift_detection_module.cho_solve
        calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original(*args, **kwargs)

        monkeypatch.setattr(
            ou_shift_detection_module, "cho_solve", counting_cho_solve
        )

        theta_fast, sig2_fast, ll_fast = svc._gls_profile_likelihood_cholesky(
            x, W, V
        )
        theta_inv, sig2_inv, ll_inv = svc._gls_profile_likelihood_inverse(
            x, W, V
        )

        assert calls == 1
        np.testing.assert_allclose(theta_fast, theta_inv)
        np.testing.assert_allclose(sig2_fast, sig2_inv)
        np.testing.assert_allclose(ll_fast, ll_inv)

    def test_profile_likelihood_falls_back_to_inverse(self, svc, monkeypatch):
        expected = (np.array([1.0]), 0.5, -2.0)
        monkeypatch.setattr(
            svc,
            "_gls_profile_likelihood_cholesky",
            lambda *_args: (_ for _ in ()).throw(ValueError("not positive definite")),
        )
        monkeypatch.setattr(
            svc, "_gls_profile_likelihood_inverse", lambda *_args: expected
        )

        observed = svc._gls_profile_likelihood(
            np.array([1.0]), np.ones((1, 1)), np.eye(1)
        )

        assert observed is expected

    def test_cholesky_profile_uses_lstsq_for_collinear_predictors(self, svc):
        x = np.array([0.0, 1.0, 4.0])
        W = np.ones((3, 2))

        theta, sigma2, log_likelihood = svc._gls_profile_likelihood_cholesky(
            x, W, np.eye(3)
        )

        assert theta.shape == (2,)
        assert sigma2 > 0
        assert np.isfinite(log_likelihood)

    def test_profile_likelihood_rejects_nonpositive_residual_variance(self, svc):
        x = np.ones(3)
        W = np.ones((3, 1))

        _, sigma2_cholesky, ll_cholesky = svc._gls_profile_likelihood_cholesky(
            x, W, np.eye(3)
        )
        _, sigma2_inverse, ll_inverse = svc._gls_profile_likelihood_inverse(
            x, W, np.eye(3)
        )

        assert sigma2_cholesky <= 0
        assert ll_cholesky == float("-inf")
        assert sigma2_inverse <= 0
        assert ll_inverse == float("-inf")

    def test_inverse_profile_rejects_indefinite_covariance(self, svc):
        x = np.array([0.0, 1.0])
        W = np.zeros((2, 1))
        V = np.diag([-1.0, 1.0])

        _, sigma2, log_likelihood = svc._gls_profile_likelihood_inverse(x, W, V)

        assert sigma2 > 0
        assert log_likelihood == float("-inf")

    def test_covariance_inverse_uses_fallbacks(self, svc, monkeypatch):
        monkeypatch.setattr(
            ou_shift_detection_module,
            "cho_factor",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(ValueError("factor")),
        )

        np.testing.assert_allclose(svc._invert_vcv(np.eye(2)), np.eye(2))

        monkeypatch.setattr(
            np.linalg,
            "inv",
            lambda *_args: (_ for _ in ()).throw(np.linalg.LinAlgError("inverse")),
        )
        assert svc._invert_vcv(np.eye(2)) is None

    def test_theta_hat_uses_lstsq_for_singular_information(self, svc):
        x = np.array([1.0, 2.0, 3.0])
        W = np.ones((3, 2))

        theta = svc._gls_theta_hat(x, W, np.eye(3))

        np.testing.assert_allclose(W @ theta, np.full(3, 2.0))


class TestFitOU1ForAlpha:
    def test_returns_positive_alpha(self, precomputed):
        svc = precomputed["svc"]
        x = precomputed["x"]
        ordered_names = precomputed["ordered_names"]
        lineage_info = precomputed["lineage_info"]
        tree_height = precomputed["tree_height"]

        alpha = svc._fit_ou1_for_alpha(x, ordered_names, lineage_info, tree_height)
        assert alpha > 0


class TestPrecomputeSharedPathLengths:
    def test_diagonal_equals_tip_heights(self, precomputed):
        svc = precomputed["svc"]
        S = svc._S
        tip_heights = svc._tip_heights_arr
        n = svc._n
        for i in range(n):
            np.testing.assert_allclose(S[i, i], tip_heights[i], atol=1e-12)

    def test_shared_path_symmetric(self, precomputed):
        svc = precomputed["svc"]
        S = svc._S
        np.testing.assert_allclose(S, S.T, atol=1e-12)

    def test_precompute_populates_lineage_row_cache(self, precomputed):
        svc = precomputed["svc"]
        rows_by_clade_id = svc._get_lineage_rows_by_clade_id()

        assert rows_by_clade_id is svc._lineage_rows_by_clade_id
        for row_idx, name in enumerate(svc._ordered_names):
            for clade_id, *_ in svc._lineage_info[name]:
                assert row_idx in rows_by_clade_id[clade_id]


class TestBuildOUVCVFast:
    def test_matches_slow_path(self, precomputed):
        """Fast VCV should match the loop-based VCV."""
        svc = precomputed["svc"]
        ordered_names = precomputed["ordered_names"]
        lineage_info = precomputed["lineage_info"]

        alpha = 1.5
        V_slow = svc._build_ou_vcv_single_alpha_no_regime(
            ordered_names, lineage_info, alpha, 1.0
        )
        V_fast = svc._build_ou_vcv_fast(alpha, 1.0)
        np.testing.assert_allclose(V_fast, V_slow, atol=1e-10)

    def test_no_regime_vcv_uses_block_shared_paths(self):
        class IterOnlyPath:
            def __init__(self, items):
                self.items = items

            def __iter__(self):
                return iter(self.items)

            def __getitem__(self, _index):
                raise AssertionError(
                    "no-regime OU VCV should not use pairwise path indexing"
                )

        svc = OUShiftDetection.__new__(OUShiftDetection)
        ordered_names = ["a", "b"]
        lineage_info = {
            "a": IterOnlyPath([(0, 1.0, 0.0, 1.0), (1, 2.0, 1.0, 3.0)]),
            "b": IterOnlyPath([(0, 1.0, 0.0, 1.0), (2, 4.0, 1.0, 5.0)]),
        }
        alpha = 0.5
        sigma2 = 2.0
        shared = np.array([[3.0, 1.0], [1.0, 5.0]])
        heights = np.array([3.0, 5.0])
        distance = heights[:, None] + heights[None, :] - 2.0 * shared
        expected = (sigma2 / (2.0 * alpha)) * np.exp(-alpha * distance) * (
            1.0 - np.exp(-2.0 * alpha * shared)
        )

        observed = svc._build_ou_vcv_single_alpha_no_regime(
            ordered_names,
            lineage_info,
            alpha,
            sigma2,
        )

        np.testing.assert_allclose(observed, expected, atol=1e-12)

    def test_bm_limit(self, precomputed):
        """Alpha near 0 should give BM covariance."""
        svc = precomputed["svc"]
        V_bm = svc._build_ou_vcv_fast(1e-12)
        np.testing.assert_allclose(V_bm, svc._S, atol=1e-6)

    def test_slow_path_bm_limit(self, precomputed):
        svc = precomputed["svc"]

        observed = svc._build_ou_vcv_single_alpha_no_regime(
            precomputed["ordered_names"],
            precomputed["lineage_info"],
            alpha=1e-12,
            sigma2=2.0,
        )

        np.testing.assert_allclose(observed, 2.0 * svc._S)


class TestNumericalFallbacks:
    def test_lineage_row_cache_rebuilds_from_paths(self, precomputed):
        svc = precomputed["svc"]
        svc._lineage_rows_by_clade_id = None

        rows_by_clade = svc._get_lineage_rows_by_clade_id()

        assert rows_by_clade is svc._lineage_rows_by_clade_id
        assert rows_by_clade
        assert all(rows.dtype == np.intp for rows in rows_by_clade.values())

    def test_zero_length_bm_shift_keeps_background_weight(self, svc):
        lineage = {"A": [(101, 0.0, 0.0, 0.0)]}

        weights = svc._build_shift_weight_matrix(
            ["A"], lineage, [(101, object())], alpha=0.0, tree_height=0.0
        )

        np.testing.assert_allclose(weights, [[1.0, 0.0]])

    def test_alpha_fit_slow_path_penalizes_invalid_likelihood(self, svc, monkeypatch):
        svc._x = np.array([1.0, 2.0, 3.0])
        svc._tree_height = 1.0
        monkeypatch.setattr(
            svc,
            "_build_ou_vcv_single_alpha_no_regime",
            lambda *_args: np.eye(3),
        )
        monkeypatch.setattr(
            svc,
            "_gls_profile_likelihood",
            lambda *_args: (np.array([0.0]), 0.0, float("-inf")),
        )

        def exercise_objective(objective, _bounds):
            assert objective(0.5) == 1e20
            return 0.5, -1e20

        monkeypatch.setattr(svc, "_optimize_parameter", exercise_objective)

        assert svc._fit_ou1_for_alpha(
            ordered_names=["A", "B", "C"],
            lineage_info={"A": [], "B": [], "C": []},
            tree_height=1.0,
        ) == pytest.approx(0.5)

    @pytest.mark.parametrize("exception", [ValueError("bad"), RuntimeError("bad")])
    def test_alpha_optimizer_returns_midpoint_on_failure(
        self, svc, monkeypatch, exception
    ):
        monkeypatch.setattr(
            ou_shift_detection_module,
            "minimize_scalar",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(exception),
        )

        assert svc._optimize_alpha(lambda _alpha: 0.0, (2.0, 6.0)) == 4.0

    def test_multibracket_optimizer_tolerates_failed_brackets(self, svc, monkeypatch):
        monkeypatch.setattr(
            ou_shift_detection_module,
            "minimize_scalar",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("bad")),
        )

        parameter, score = svc._optimize_parameter(
            lambda _alpha: 0.0, (2.0, 6.0), niter=2
        )

        assert parameter == 4.0
        assert score == float("-inf")

    def test_multibracket_optimizer_ignores_empty_brackets(self, svc):
        parameter, score = svc._optimize_parameter(
            lambda _alpha: 0.0, (2.0, 2.0), niter=2
        )

        assert parameter == 2.0
        assert score == float("-inf")

    def test_lasso_config_extraction_handles_degenerate_designs(
        self, svc, monkeypatch
    ):
        svc._n = 3
        y = np.arange(3.0)

        assert svc._extract_lasso_configs(np.ones((3, 1)), y, max_shifts=2) == []
        assert svc._extract_lasso_configs(np.ones((3, 2)), y, max_shifts=0) == []

        monkeypatch.setattr(
            ou_shift_detection_module,
            "_lars_path",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(ValueError("lasso")),
        )
        assert svc._extract_lasso_configs(np.eye(3), y, max_shifts=2) == []

    def test_pbic_vcv_uses_inverse_fallback(self, svc, monkeypatch):
        monkeypatch.setattr(
            ou_shift_detection_module,
            "cho_factor",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(ValueError("factor")),
        )
        svc._n_edges = 5
        x = np.array([1.0, 2.0, 4.0])
        W = np.ones((3, 1))

        observed = svc._compute_pbic_from_vcv(
            3, -2.0, 0, W, np.eye(3), 1.0, x
        )

        assert np.isfinite(observed)

        monkeypatch.setattr(
            np.linalg,
            "inv",
            lambda *_args: (_ for _ in ()).throw(np.linalg.LinAlgError("inverse")),
        )
        assert svc._compute_pbic_from_vcv(
            3, -2.0, 0, W, np.eye(3), 1.0, x
        ) is None

    def test_pbic_rejects_degenerate_information(self, svc):
        svc._n_edges = 5

        assert svc._compute_pbic_from_info(
            3, -2.0, 0, np.eye(1), 1.0, np.ones(3)
        ) == float("inf")
        assert svc._compute_pbic_from_info(
            2, -2.0, 1, np.eye(2), 1.0, np.array([1.0, 2.0])
        ) == float("inf")
        assert svc._compute_pbic_from_info(
            3, -2.0, 0, np.eye(1), 0.0, np.array([1.0, 2.0, 4.0])
        ) == float("inf")
        assert svc._compute_pbic_from_info(
            3, -2.0, 0, -np.eye(1), 1.0, np.array([1.0, 2.0, 4.0])
        ) == float("inf")


class TestEndToEnd:
    @staticmethod
    def _stub_run_setup(svc, mocker, edges):
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:1);"), "newick")
        traits = {"A": 1.0, "B": 2.0, "C": 3.0}
        lineage = {name: [(0, 1.0, 0.0, 1.0)] for name in traits}
        mocker.patch.object(svc, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(
            svc, "get_tip_names_from_tree", return_value=list(traits)
        )
        mocker.patch.object(svc, "_parse_trait_file", return_value=traits)
        mocker.patch.object(svc, "_build_parent_map", return_value={})
        mocker.patch.object(
            svc, "_build_lineage_info_no_regime", return_value=lineage
        )
        mocker.patch.object(svc, "_enumerate_edges", return_value=edges)
        mocker.patch.object(
            svc,
            "_precompute_shared_path_lengths",
            return_value=(np.eye(3), np.ones(3)),
        )
        build_result = mocker.patch.object(
            svc,
            "_build_final_result",
            side_effect=lambda best, _tree: {"best": best},
        )
        mocker.patch.object(svc, "_output")
        return build_result

    def test_pipeline_runs_with_setup(self, precomputed):
        """Pipeline should complete with pre-computed data."""
        svc = precomputed["svc"]
        # Use _fit_and_score_config for null model
        null = svc._fit_and_score_config([])
        assert null is not None
        assert null["n_shifts"] == 0
        assert "log_likelihood" in null
        assert "pBIC" in null

    def test_fit_with_shifts(self, precomputed):
        """Fitting with shift edges should complete."""
        svc = precomputed["svc"]
        # Fit with first 2 edges as shifts
        result = svc._fit_and_score_config([0, 1])
        assert result is not None
        assert result["n_shifts"] == 2
        assert result["alpha"] > 0
        assert np.isfinite(result["pBIC"])

    def test_run_completes(self, default_args):
        """Full run completes without error."""
        from unittest.mock import patch
        svc = OUShiftDetection(default_args)
        with patch("builtins.print"):
            svc.run()

    def test_run_uses_fast_tip_name_helper_for_prune_setup(self, default_args, mocker):
        svc = OUShiftDetection(default_args)
        spy = mocker.spy(svc, "get_tip_names_from_tree")
        mocker.patch("builtins.print")

        svc.run()

        assert spy.call_count == 1

    def test_prepare_shared_trait_data_skips_prune_scan_when_all_shared(self):
        tree_tips = ["A", "B", "C"]
        traits = {"A": 1.0, "B": 2.0, "C": 3.0}

        to_prune, ordered_names = OUShiftDetection._prepare_shared_trait_data(
            tree_tips,
            traits,
        )

        assert to_prune == []
        assert ordered_names == ["A", "B", "C"]

    def test_prepare_shared_trait_data_preserves_partial_prune_order(self):
        tree_tips = ["A", "B", "C", "D"]
        traits = {"A": 1.0, "C": 3.0, "B": 2.0}

        to_prune, ordered_names = OUShiftDetection._prepare_shared_trait_data(
            tree_tips,
            traits,
        )

        assert to_prune == ["D"]
        assert ordered_names == ["A", "B", "C"]

    def test_run_uses_unmodified_tree_read(self, default_args, mocker):
        svc = OUShiftDetection(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        traits = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        lineage = {name: [(0, 1.0, 0, 1)] for name in traits}
        best = {"shifts": [], "alpha": 0.0, "pBIC": 1.0, "BIC": 1.0, "AICc": 1.0}

        read_unmodified = mocker.patch.object(
            svc, "read_tree_file_unmodified", return_value=tree
        )
        mocker.patch.object(
            svc,
            "read_tree_file",
            side_effect=AssertionError("run should use read_tree_file_unmodified"),
        )
        validate = mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(svc, "get_tip_names_from_tree", return_value=list(traits))
        mocker.patch.object(svc, "_parse_trait_file", return_value=traits)
        mocker.patch.object(svc, "_build_parent_map", return_value={})
        mocker.patch.object(
            svc, "_build_lineage_info_no_regime", return_value=lineage
        )
        mocker.patch.object(svc, "_enumerate_edges", return_value=[])
        mocker.patch.object(
            svc, "_precompute_shared_path_lengths", return_value=(None, None)
        )
        mocker.patch.object(svc, "_fit_and_score_config", return_value=best)
        build_result = mocker.patch.object(
            svc, "_build_final_result", return_value={"ok": True}
        )
        mocker.patch.object(svc, "_output")
        fast_copy = mocker.patch.object(
            svc,
            "_fast_copy",
            side_effect=AssertionError("all-shared analysis should not copy tree"),
        )

        svc.run()

        read_unmodified.assert_called_once_with()
        validate.assert_called_once()
        fast_copy.assert_not_called()
        assert build_result.call_args.args[1] is tree

    def test_run_copies_before_pruning_missing_tree_tips(self, default_args, mocker):
        svc = OUShiftDetection(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree_copy = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        traits = {"A": 1.0, "B": 2.0, "C": 3.0}
        lineage = {name: [(0, 1.0, 0, 1)] for name in traits}
        best = {"shifts": [], "alpha": 0.0, "pBIC": 1.0, "BIC": 1.0, "AICc": 1.0}

        mocker.patch.object(svc, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(
            svc, "get_tip_names_from_tree", return_value=["A", "B", "C", "D"]
        )
        mocker.patch.object(svc, "_parse_trait_file", return_value=traits)
        fast_copy = mocker.patch.object(svc, "_fast_copy", return_value=tree_copy)
        prune = mocker.spy(tree_copy, "prune")
        mocker.patch.object(svc, "_build_parent_map", return_value={})
        mocker.patch.object(
            svc, "_build_lineage_info_no_regime", return_value=lineage
        )
        mocker.patch.object(svc, "_enumerate_edges", return_value=[])
        mocker.patch.object(
            svc, "_precompute_shared_path_lengths", return_value=(None, None)
        )
        mocker.patch.object(svc, "_fit_and_score_config", return_value=best)
        build_result = mocker.patch.object(
            svc, "_build_final_result", return_value={"ok": True}
        )
        mocker.patch.object(svc, "_output")

        svc.run()

        fast_copy.assert_called_once_with(tree)
        prune.assert_called_once_with("D")
        assert build_result.call_args.args[1] is tree_copy

    def test_run_uses_empty_result_when_zero_edge_model_fails(
        self, default_args, mocker
    ):
        svc = OUShiftDetection(default_args)
        build_result = self._stub_run_setup(svc, mocker, edges=[])
        mocker.patch.object(svc, "_fit_and_score_config", return_value=None)

        svc.run()

        fallback = build_result.call_args.args[0]
        assert fallback["n_shifts"] == 0
        assert fallback["log_likelihood"] == float("-inf")

    def test_run_uses_empty_result_when_lasso_and_null_models_fail(
        self, default_args, mocker
    ):
        svc = OUShiftDetection(default_args)
        build_result = self._stub_run_setup(svc, mocker, edges=[(1, object())])
        mocker.patch.object(svc, "_lasso_pass_with_selection", return_value=None)
        mocker.patch.object(svc, "_fit_and_score_config", return_value=None)

        svc.run()

        assert build_result.call_args.args[0]["n_shifts"] == 0

    def test_run_selects_improved_second_pass_model(self, default_args, mocker):
        svc = OUShiftDetection(default_args)
        build_result = self._stub_run_setup(svc, mocker, edges=[(1, object())])
        first = {"alpha": 0.5, "pBIC": 10.0}
        second = {"alpha": 0.75, "pBIC": 5.0}
        lasso = mocker.patch.object(
            svc, "_lasso_pass_with_selection", side_effect=[first, second]
        )

        svc.run()

        assert lasso.call_count == 2
        assert build_result.call_args.args[0] is second


class TestOutput:
    def test_print_text_output_batches_detected_shifts(self, svc, mocker):
        printed = mocker.patch("builtins.print")
        result = {
            "n_tips": 8,
            "n_shifts": 2,
            "criterion": "pBIC",
            "alpha": 1.234567,
            "sigma2": 0.987654,
            "theta_root": 2.345678,
            "log_likelihood": -123.456789,
            "pBIC": 456.789,
            "BIC": 567.891,
            "AICc": 345.678,
            "shifts": [
                {"description": "terminal branch to A", "optimum": 1.25},
                {"description": "stem of clade B,C", "optimum": -0.5},
            ],
        }

        svc._print_text_output(result)

        printed.assert_called_once_with(
            "============================================================\n"
            "OU Shift Detection (l1ou)\n"
            "============================================================\n"
            "Number of tips:       8\n"
            "Number of shifts:     2\n"
            "Selection criterion:  pBIC\n"
            "Alpha (OU strength):  1.234567\n"
            "Sigma² (BM rate):     0.987654\n"
            "Root optimum (θ₀):    2.345678\n"
            "Log-likelihood:       -123.4568\n"
            "pBIC:                 456.7890\n"
            "BIC:                  567.8910\n"
            "AICc:                 345.6780\n"
            "\n"
            "Detected shifts:\n"
            "------------------------------------------------------------\n"
            "  Shift 1: terminal branch to A\n"
            "           New optimum: 1.250000\n"
            "  Shift 2: stem of clade B,C\n"
            "           New optimum: -0.500000\n"
            "============================================================"
        )

    def test_print_text_output_reports_no_shifts(self, svc, capsys):
        result = {
            "n_tips": 3,
            "n_shifts": 0,
            "criterion": "pBIC",
            "alpha": 0.5,
            "sigma2": 1.0,
            "theta_root": 2.0,
            "log_likelihood": -3.0,
            "pBIC": 10.0,
            "BIC": 11.0,
            "AICc": 12.0,
            "shifts": [],
        }

        svc._print_text_output(result)

        assert "No shifts detected" in capsys.readouterr().out

    def test_print_text_output_allows_closed_pipe(self, svc, monkeypatch):
        result = {
            "n_tips": 3,
            "n_shifts": 0,
            "criterion": "pBIC",
            "alpha": 0.5,
            "sigma2": 1.0,
            "theta_root": 2.0,
            "log_likelihood": -3.0,
            "pBIC": 10.0,
            "BIC": 11.0,
            "AICc": 12.0,
            "shifts": [],
        }
        monkeypatch.setattr(
            "builtins.print",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(BrokenPipeError),
        )

        svc._print_text_output(result)


class TestDescribeEdge:
    def test_terminal(self, precomputed):
        svc = precomputed["svc"]
        tree = precomputed["tree"]
        edges = precomputed["edges"]
        for clade_id, clade in edges:
            if clade.is_terminal():
                desc = svc._describe_edge(clade, tree)
                assert "terminal branch to" in desc
                break

    def test_internal(self, precomputed):
        svc = precomputed["svc"]
        tree = precomputed["tree"]
        edges = precomputed["edges"]
        for clade_id, clade in edges:
            if not clade.is_terminal():
                desc = svc._describe_edge(clade, tree)
                assert "stem of" in desc
                break

    def test_uses_direct_terminal_name_traversal(self, precomputed, monkeypatch):
        from Bio.Phylo.BaseTree import Clade

        svc = precomputed["svc"]
        tree = precomputed["tree"]

        def fail_get_terminals(self):
            raise AssertionError("edge descriptions should use direct traversal")

        monkeypatch.setattr(Clade, "get_terminals", fail_get_terminals)

        desc = svc._describe_edge(tree.root, tree)
        assert desc.startswith("stem of")

    def test_falls_back_to_biopython_terminal_traversal(self, svc, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:1);"), "newick")
        monkeypatch.setattr(
            ou_shift_detection_module.Tree,
            "calculate_terminal_names_fast",
            staticmethod(lambda _clade: None),
        )

        description = svc._describe_edge(tree.root.clades[0], tree)

        assert description == "stem of (A, B)"


class TestProcessArgs:
    def test_invalid_criterion(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            criterion="INVALID",
            max_shifts=None,
            json=False,
        )
        with pytest.raises(PhykitUserError):
            OUShiftDetection(args)

    def test_valid_criteria(self):
        for crit in ("pBIC", "BIC", "AICc"):
            args = Namespace(
                tree=TREE_SIMPLE,
                trait_data=TRAITS_FILE,
                criterion=crit,
                max_shifts=None,
                json=False,
            )
            svc = OUShiftDetection(args)
            assert svc.criterion == crit


class TestTraitParsing:
    def test_comments_and_blanks(self, tmp_path, default_args):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "# comment\n\nraccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\n"
            "seal\t4.0\nmonkey\t5.0\ncat\t6.0\nweasel\t7.0\ndog\t8.0\n"
        )
        svc = OUShiftDetection(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        traits = svc._parse_trait_file(str(trait_file), tree_tips)

        assert len(traits) == 8
        assert traits["raccoon"] == pytest.approx(1.0)

    def test_all_shared_trait_file_emits_no_warnings(
        self, tmp_path, default_args, capsys
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "raccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\nseal\t4.0\n"
            "monkey\t5.0\ncat\t6.0\nweasel\t7.0\ndog\t8.0\n"
        )
        svc = OUShiftDetection(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        traits = svc._parse_trait_file(str(trait_file), tree_tips)

        stderr = capsys.readouterr().err
        assert set(traits) == set(tree_tips)
        assert stderr == ""

    def test_ordered_all_shared_trait_file_skips_sets(
        self, tmp_path, default_args, monkeypatch
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("A\t1.0\nB\t2.0\nC\t3.0\n")
        svc = OUShiftDetection(default_args)
        original_set = builtins.set

        def fail_set(*args, **kwargs):
            raise AssertionError("ordered exact path should not build sets")

        monkeypatch.setattr(builtins, "set", fail_set)
        traits = svc._parse_trait_file(str(trait_file), ["A", "B", "C"])

        assert traits == {"A": 1.0, "B": 2.0, "C": 3.0}
        assert builtins.set is fail_set
        monkeypatch.setattr(builtins, "set", original_set)

    def test_extra_columns_error(self, tmp_path, default_args):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "raccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\ndog\t4.0\nextra\t1.0\t2.0\n"
        )
        svc = OUShiftDetection(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_trait_file(str(trait_file), tree_tips)

        assert (
            "Line 5 in trait file has 3 columns; expected 2."
            in exc_info.value.messages
        )

    def test_non_numeric_trait_error(self, tmp_path, default_args):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("raccoon\t1.0\nbear\tbad\nsea_lion\t3.0\ndog\t4.0\n")
        svc = OUShiftDetection(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_trait_file(str(trait_file), tree_tips)

        assert (
            "Non-numeric trait value 'bad' for taxon 'bear' on line 2."
            in exc_info.value.messages
        )

    def test_missing_trait_file_error(self, tmp_path, default_args):
        missing = tmp_path / "missing_traits.tsv"
        svc = OUShiftDetection(default_args)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_trait_file(str(missing), ["A", "B", "C"])

        assert exc_info.value.code == 2
        assert str(missing) in exc_info.value.messages[0]

    def test_unordered_complete_trait_file_is_accepted(self, tmp_path, default_args):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("C\t3.0\nA\t1.0\nB\t2.0\n")
        svc = OUShiftDetection(default_args)

        traits = svc._parse_trait_file(str(trait_file), ["A", "B", "C"])

        assert traits == {"C": 3.0, "A": 1.0, "B": 2.0}

    def test_partial_overlap_warns_and_filters_traits(
        self, tmp_path, default_args, capsys
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("A\t1.0\nB\t2.0\nC\t3.0\nE\t5.0\n")
        svc = OUShiftDetection(default_args)

        traits = svc._parse_trait_file(
            str(trait_file), ["A", "B", "C", "D"]
        )

        stderr = capsys.readouterr().err
        assert traits == {"A": 1.0, "B": 2.0, "C": 3.0}
        assert "1 taxa in tree but not in trait file: D" in stderr
        assert "1 taxa in trait file but not in tree: E" in stderr

    def test_fewer_than_three_shared_taxa_is_rejected(
        self, tmp_path, default_args
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("A\t1.0\nB\t2.0\nD\t4.0\n")
        svc = OUShiftDetection(default_args)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_trait_file(str(trait_file), ["A", "B", "C"])

        assert "Only 2 shared taxa" in exc_info.value.messages[0]


class TestModelFailureHandling:
    def test_lasso_pass_skips_unfittable_candidate(self, precomputed, monkeypatch):
        svc = precomputed["svc"]
        null = {"pBIC": 10.0, "n_shifts": 0}
        monkeypatch.setattr(
            svc,
            "_build_shift_weight_matrix",
            lambda *_args: np.ones((svc._n, len(svc._edges) + 1)),
        )
        monkeypatch.setattr(svc, "_build_ou_vcv_fast", lambda *_args: np.eye(svc._n))
        monkeypatch.setattr(
            svc,
            "_transform_to_independent",
            lambda x, W, _V: (W, x),
        )
        monkeypatch.setattr(svc, "_extract_lasso_configs", lambda *_args: [[0]])
        monkeypatch.setattr(
            svc, "_fit_and_score_config", lambda config: null if not config else None
        )

        assert svc._lasso_pass_with_selection(0.5, 1) is null

    def test_shift_model_rejects_invalid_final_fit(self, precomputed, monkeypatch):
        svc = precomputed["svc"]
        monkeypatch.setattr(svc, "_optimize_alpha", lambda *_args: 0.5)
        monkeypatch.setattr(
            svc,
            "_gls_profile_likelihood",
            lambda *_args: (np.zeros(2), 0.0, float("-inf")),
        )

        assert svc._fit_and_score_config([0]) is None

    def test_shift_model_rejects_unavailable_pbic(self, precomputed, monkeypatch):
        svc = precomputed["svc"]
        monkeypatch.setattr(svc, "_optimize_alpha", lambda *_args: 0.5)
        monkeypatch.setattr(
            svc,
            "_gls_profile_likelihood",
            lambda *_args: (np.array([1.0, 0.5]), 1.0, -2.0),
        )
        monkeypatch.setattr(svc, "_compute_pbic_from_vcv", lambda *_args: None)

        assert svc._fit_and_score_config([0]) is None

    @pytest.mark.parametrize("pbic_available", [False, True])
    def test_null_model_rejects_invalid_fit_or_pbic(
        self, precomputed, monkeypatch, pbic_available
    ):
        svc = precomputed["svc"]
        monkeypatch.setattr(svc, "_fit_ou1_for_alpha", lambda: 0.5)
        if pbic_available:
            monkeypatch.setattr(
                svc,
                "_gls_profile_likelihood",
                lambda *_args: (np.array([1.0]), 0.0, float("-inf")),
            )
        else:
            monkeypatch.setattr(
                svc,
                "_gls_profile_likelihood",
                lambda *_args: (np.array([1.0]), 1.0, -2.0),
            )
            monkeypatch.setattr(svc, "_compute_pbic_from_vcv", lambda *_args: None)

        assert svc._fit_null_model() is None

    def test_backward_elimination_keeps_model_when_reduced_fit_fails(
        self, svc, monkeypatch
    ):
        fitted = {"shift_edge_indices": [0], "pBIC": 10.0}
        monkeypatch.setattr(svc, "_fit_and_score_config", lambda _config: None)

        assert svc._backward_eliminate_single_pass(fitted) is fitted

    def test_cholesky_inverse_path_returns_matrix(self, svc):
        np.testing.assert_allclose(svc._invert_vcv(np.eye(2)), np.eye(2))
