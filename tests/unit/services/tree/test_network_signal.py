import builtins
import importlib
import json
import subprocess
import sys
import pytest
import numpy as np
from argparse import Namespace
from io import StringIO

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

from phykit.services.tree.network_signal import NetworkSignal
import phykit.services.tree.network_signal as network_signal_module
from phykit.errors import PhykitUserError


def test_module_import_does_not_import_numpy_or_scipy():
    code = """
import sys
import phykit.services.tree.network_signal as module
assert callable(module.print_json)
assert hasattr(module.np, "__getattr__")
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy.linalg" not in sys.modules
assert "scipy.optimize" not in sys.modules
assert module._CHO_FACTOR is None
assert module._CHO_SOLVE is None
assert module._MINIMIZE_SCALAR is None
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_module_import_does_not_import_scipy_linalg_or_optimize(monkeypatch):
    module_name = "phykit.services.tree.network_signal"
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
                "network_signal module import should not import SciPy linalg/optimize"
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
    previous_cho_factor = network_signal_module._CHO_FACTOR
    previous_cho_solve = network_signal_module._CHO_SOLVE
    network_signal_module._CHO_FACTOR = None
    network_signal_module._CHO_SOLVE = None
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
        n = 12
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        x = rng.normal(size=n)

        NetworkSignal._log_likelihood_cholesky(x, C)
        first_call_imports = scipy_linalg_imports
        NetworkSignal._log_likelihood_cholesky(x, C)
    finally:
        network_signal_module._CHO_FACTOR = previous_cho_factor
        network_signal_module._CHO_SOLVE = previous_cho_solve

    assert first_call_imports > 0
    assert scipy_linalg_imports == first_call_imports


def test_repeated_minimize_scalar_calls_cache_scipy_optimize_imports(monkeypatch):
    previous_minimize_scalar = network_signal_module._MINIMIZE_SCALAR
    network_signal_module._MINIMIZE_SCALAR = None
    original_import = builtins.__import__
    scipy_optimize_imports = 0

    def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
        nonlocal scipy_optimize_imports
        if name == "scipy.optimize":
            scipy_optimize_imports += 1
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", counting_import)
    try:
        first = network_signal_module.minimize_scalar(
            lambda value: (value - 0.25) ** 2,
            bounds=(0.0, 1.0),
            method="bounded",
        )
        first_call_imports = scipy_optimize_imports
        second = network_signal_module.minimize_scalar(
            lambda value: (value - 0.75) ** 2,
            bounds=(0.0, 1.0),
            method="bounded",
        )
    finally:
        network_signal_module._MINIMIZE_SCALAR = previous_minimize_scalar

    assert first.x == pytest.approx(0.25, abs=1e-5)
    assert second.x == pytest.approx(0.75, abs=1e-5)
    assert first_call_imports > 0
    assert scipy_optimize_imports == first_call_imports


def test_inverse_from_spd_matrix_cholesky_matches_inverse():
    rng = np.random.default_rng(20260628)
    A = rng.normal(size=(8, 8))
    matrix = A @ A.T + np.eye(8) * 0.25

    observed = network_signal_module._inverse_from_spd_matrix(matrix)

    np.testing.assert_allclose(observed, np.linalg.inv(matrix))


def test_inverse_from_spd_matrix_avoids_inverse_for_spd_matrix(monkeypatch):
    rng = np.random.default_rng(20260628)
    A = rng.normal(size=(6, 6))
    matrix = A @ A.T + np.eye(6)

    def fail_inverse(_matrix):
        raise AssertionError("SPD covariance inversion should use Cholesky")

    monkeypatch.setattr(network_signal_module.np.linalg, "inv", fail_inverse)

    observed = network_signal_module._inverse_from_spd_matrix(matrix)

    assert observed.shape == matrix.shape
    assert np.all(np.isfinite(observed))


def test_inverse_from_spd_matrix_keeps_inverse_fallback():
    matrix = np.array([[1.0, 2.0], [2.0, 1.0]])

    observed = network_signal_module._inverse_from_spd_matrix(matrix)

    np.testing.assert_allclose(observed, np.linalg.inv(matrix))


def _make_tree(newick):
    return Phylo.read(StringIO(newick), "newick")


class TestValidateTree:
    def test_validate_tree_uses_direct_traversal(self, monkeypatch):
        svc = NetworkSignal.__new__(NetworkSignal)
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")

        def fail_traversal(*args, **kwargs):
            raise AssertionError("tree validation should use direct traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)

        svc._validate_tree(tree)

    def test_validate_tree_rejects_short_tree_before_branch_lengths(self):
        svc = NetworkSignal.__new__(NetworkSignal)
        tree = _make_tree("(A,B);")

        with pytest.raises(PhykitUserError) as exc_info:
            svc._validate_tree(tree)

        assert "at least 3 tips" in exc_info.value.messages[0]

    def test_validate_tree_rejects_missing_branch_lengths(self):
        svc = NetworkSignal.__new__(NetworkSignal)
        tree = _make_tree("(A:1,B,C:1);")

        with pytest.raises(PhykitUserError) as exc_info:
            svc._validate_tree(tree)

        assert "must have lengths" in exc_info.value.messages[0]

    def test_validate_tree_warns_for_non_root_trifurcation_only(self, capsys):
        svc = NetworkSignal.__new__(NetworkSignal)

        svc._validate_tree(_make_tree("(A:1,B:1,C:1);"))
        assert capsys.readouterr().err == ""

        svc._validate_tree(_make_tree("((A:1,B:1,C:1):1,D:1);"))
        assert "1 polytomous node" in capsys.readouterr().err

    def test_validate_tree_handles_mixed_child_counts(self, monkeypatch, capsys):
        svc = NetworkSignal.__new__(NetworkSignal)
        tree = _make_tree("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);")

        def fail_traversal(*args, **kwargs):
            raise AssertionError("tree validation should use direct traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)

        svc._validate_tree(tree)

        assert "1 polytomous node" in capsys.readouterr().err


class TestParseTraitFile:
    def test_parse_trait_file_skips_comments_and_blanks(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "# comment\n"
            "\n"
            "A\t1.5\n"
            "B\t2.5\n"
            "C\t3.5\n"
            "D\t4.5\n"
        )
        svc = NetworkSignal.__new__(NetworkSignal)

        traits = svc._parse_trait_file(str(trait_file), ["A", "B", "C"])

        assert traits == {"A": 1.5, "B": 2.5, "C": 3.5}

    def test_parse_trait_file_all_shared_emits_no_warnings(self, tmp_path, capsys):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("A\t1.5\nB\t2.5\nC\t3.5\n")
        svc = NetworkSignal.__new__(NetworkSignal)

        traits = svc._parse_trait_file(str(trait_file), ["A", "B", "C"])

        stderr = capsys.readouterr().err
        assert traits == {"A": 1.5, "B": 2.5, "C": 3.5}
        assert stderr == ""

    def test_parse_trait_file_rejects_wrong_column_count(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("A\t1.0\nB\t2.0\textra\nC\t3.0\n")
        svc = NetworkSignal.__new__(NetworkSignal)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_trait_file(str(trait_file), ["A", "B", "C"])

        assert "Line 2 in trait file has 3 columns; expected 2." in (
            exc_info.value.messages[0]
        )

    def test_parse_trait_file_rejects_non_numeric_value(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("A\t1.0\nB\tmissing\nC\t3.0\n")
        svc = NetworkSignal.__new__(NetworkSignal)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_trait_file(str(trait_file), ["A", "B", "C"])

        assert (
            "Non-numeric trait value 'missing' for taxon 'B' on line 2."
            in exc_info.value.messages[0]
        )


# ---------------------------------------------------------------------------
# Helpers to build DAGs without going through the full service
# ---------------------------------------------------------------------------

def _build_simple_dag():
    """Build DAG from ((A:1,B:1):0.5,C:1.5)"""
    tree = _make_tree("((A:1,B:1):0.5,C:1.5);")
    return NetworkSignal._tree_to_dag(tree)


def _build_five_taxon_dag():
    """Build DAG from (((A:1,B:1):0.5,C:1.5):0.5,(D:1,E:1):1)"""
    tree = _make_tree("(((A:1,B:1):0.5,C:1.5):0.5,(D:1,E:1):1);")
    return NetworkSignal._tree_to_dag(tree)


class TestRun:
    def test_run_uses_unmodified_tree_read(self, mocker):
        args = Namespace(
            tree="tree.tre",
            trait_data="traits.tsv",
            method="blombergs_k",
            permutations=0,
            hybrid=None,
            quartet_json=None,
            json=False,
        )
        svc = NetworkSignal(args)
        tree = _make_tree("((A:1,B:1):1,C:1);")
        traits = {"A": 1.0, "B": 2.0, "C": 3.0}
        nodes = [0, 1, 2, 3]
        parents = {
            0: [],
            1: [(0, 1.0, 1.0)],
            2: [(0, 1.0, 1.0)],
            3: [(0, 1.0, 1.0)],
        }
        tip_map = {1: "A", 2: "B", 3: "C"}

        read_unmodified = mocker.patch.object(
            svc, "read_tree_file_unmodified", return_value=tree
        )
        mocker.patch.object(
            svc,
            "read_tree_file",
            side_effect=AssertionError("run should use read_tree_file_unmodified"),
        )
        validate = mocker.patch.object(svc, "_validate_tree")
        mocker.patch.object(svc, "get_tip_names_from_tree", return_value=list(traits))
        mocker.patch.object(svc, "_parse_trait_file", return_value=traits)
        tree_to_dag = mocker.patch.object(
            svc, "_tree_to_dag", return_value=(nodes, parents, tip_map)
        )
        mocker.patch.object(svc, "_compute_network_vcv", return_value=np.eye(3))
        mocker.patch.object(svc, "_blombergs_k", return_value={"K": 1.0})
        mocker.patch.object(svc, "_output")

        svc.run()

        read_unmodified.assert_called_once_with()
        validate.assert_called_once_with(tree)
        assert tree_to_dag.call_args.args[0] is tree

    def test_output_batches_hybrid_edges_with_indexed_donor_lookup(self, mocker):
        svc = NetworkSignal.__new__(NetworkSignal)
        svc.json_output = False
        svc.method = "both"
        tip_map = {
            1: "Recipient",
            2: "Donor",
            3: "Other",
        }
        parents = {
            1: [(10, 1.0, 0.7), (20, 1.0, 0.3)],
            2: [(20, 1.0, 1.0)],
            3: [(10, 1.0, 1.0)],
        }
        results = {
            "blombergs_k": {"K": 1.25, "p_value": 0.0123},
            "pagels_lambda": {
                "lambda": 0.5,
                "log_likelihood": -10.0,
                "p_value": 0.0456,
            },
        }
        mocked_print = mocker.patch("builtins.print")

        svc._output(results, ["Recipient", "Donor", "Other"], None, parents, tip_map)

        mocked_print.assert_called_once()
        output = mocked_print.call_args.args[0]
        assert "Hybrid edge: Donor -> Recipient (gamma=0.3000)" in output
        assert "Network taxa: 3" in output
        assert "Blomberg's K: 1.2500    p-value: 0.0123" in output
        assert (
            "Pagel's lambda: 0.5000    log-likelihood: -10.0000    p-value: 0.0456"
            in output
        )


# ===========================================================================
# TestTreeToNetwork
# ===========================================================================

class TestTreeToNetwork:
    def test_simple_three_taxon(self):
        """((A:1,B:1):0.5,C:1.5) -> 5 nodes, 3 tips, root has no parents.

        Biopython parses this as: root (2 children: internal node + C).
        Total = 1 root + 1 internal + 3 tips = 5 nodes.
        """
        nodes, parents, tip_map = _build_simple_dag()
        assert len(nodes) == 5
        assert len(tip_map) == 3
        assert set(tip_map.values()) == {"A", "B", "C"}
        # Root is first in topological order and has no parents
        root = nodes[0]
        assert root not in parents or len(parents[root]) == 0

    def test_five_taxon(self):
        """(((A:1,B:1):0.5,C:1.5):0.5,(D:1,E:1):1) -> 9 nodes, 5 tips."""
        nodes, parents, tip_map = _build_five_taxon_dag()
        assert len(nodes) == 9
        assert len(tip_map) == 5
        assert set(tip_map.values()) == {"A", "B", "C", "D", "E"}


# ===========================================================================
# TestAddHybridEdge
# ===========================================================================

class TestAddHybridEdge:
    def test_add_single_hybrid(self):
        """B->C gamma=0.3: C gets 2 parents with gammas 0.7 and 0.3."""
        nodes, parents, tip_map = _build_simple_dag()
        # Find node IDs for B and C by name
        name_to_id = {v: k for k, v in tip_map.items()}
        donor_id = name_to_id["B"]
        recipient_id = name_to_id["C"]

        NetworkSignal._add_hybrid_edge(
            nodes, parents, tip_map, "B", "C", gamma=0.3
        )

        # Recipient C should now have exactly 2 parents
        assert len(parents[recipient_id]) == 2
        gammas = sorted([entry[2] for entry in parents[recipient_id]])
        assert gammas == pytest.approx([0.3, 0.7])

    def test_invalid_donor_raises(self):
        """Nonexistent donor raises PhykitUserError (SystemExit)."""
        nodes, parents, tip_map = _build_simple_dag()
        with pytest.raises(SystemExit):
            NetworkSignal._add_hybrid_edge(
                nodes, parents, tip_map, "NONEXISTENT", "C", gamma=0.3
            )

    def test_gamma_out_of_range_raises(self):
        """gamma=0.6 raises PhykitUserError."""
        nodes, parents, tip_map = _build_simple_dag()
        with pytest.raises(SystemExit):
            NetworkSignal._add_hybrid_edge(
                nodes, parents, tip_map, "B", "C", gamma=0.6
            )


# ===========================================================================
# TestNetworkVCV
# ===========================================================================

class TestNetworkVCV:
    @staticmethod
    def _compute_network_vcv_scalar_reference(nodes, parents, tip_map, ordered_tips):
        n_nodes = len(nodes)
        node_to_idx = {nid: i for i, nid in enumerate(nodes)}
        V = np.zeros((n_nodes, n_nodes))

        for nid in nodes:
            idx = node_to_idx[nid]
            parent_list = parents.get(nid, [])

            if len(parent_list) == 0:
                V[idx, idx] = 0.0
            elif len(parent_list) == 1:
                p_id, edge_len, _ = parent_list[0]
                p_idx = node_to_idx[p_id]
                V[idx, idx] = V[p_idx, p_idx] + edge_len
                for j in range(idx):
                    V[idx, j] = V[p_idx, j]
                    V[j, idx] = V[p_idx, j]
            elif len(parent_list) == 2:
                p1_id, l1, g1 = parent_list[0]
                p2_id, l2, g2 = parent_list[1]
                p1_idx = node_to_idx[p1_id]
                p2_idx = node_to_idx[p2_id]
                V[idx, idx] = (
                    g1 * g1 * (V[p1_idx, p1_idx] + l1)
                    + g2 * g2 * (V[p2_idx, p2_idx] + l2)
                    + 2 * g1 * g2 * V[p1_idx, p2_idx]
                )
                for j in range(idx):
                    cov = g1 * V[p1_idx, j] + g2 * V[p2_idx, j]
                    V[idx, j] = cov
                    V[j, idx] = cov

        name_to_id = {v: k for k, v in tip_map.items()}
        tip_ids = [name_to_id[tip_name] for tip_name in ordered_tips]
        vcv = np.zeros((len(ordered_tips), len(ordered_tips)))
        for i, tid_i in enumerate(tip_ids):
            for j, tid_j in enumerate(tip_ids):
                vcv[i, j] = V[node_to_idx[tid_i], node_to_idx[tid_j]]
        return vcv

    def test_network_vcv_matches_scalar_reference_with_custom_tip_order(self):
        nodes, parents, tip_map = _build_five_taxon_dag()
        NetworkSignal._add_hybrid_edge(nodes, parents, tip_map, "B", "D", gamma=0.3)
        ordered_tips = ["E", "C", "A", "D", "B"]

        fast = NetworkSignal._compute_network_vcv(
            nodes, parents, tip_map, ordered_tips
        )
        reference = self._compute_network_vcv_scalar_reference(
            nodes, parents, tip_map, ordered_tips
        )

        np.testing.assert_array_almost_equal(fast, reference)

    def test_tree_only_matches_standard(self):
        """No hybrid edge -> VCV matches expected tree VCV.

        For ((A:1,B:1):0.5,C:1.5):
        V[A,A]=V[B,B]=V[C,C]=1.5, V[A,B]=0.5, V[A,C]=V[B,C]=0.0
        """
        nodes, parents, tip_map = _build_simple_dag()
        ordered_tips = sorted(tip_map.values())  # ['A', 'B', 'C']
        vcv = NetworkSignal._compute_network_vcv(nodes, parents, tip_map, ordered_tips)

        assert vcv.shape == (3, 3)
        # Diagonal: root-to-tip distances
        np.testing.assert_almost_equal(vcv[0, 0], 1.5)  # A
        np.testing.assert_almost_equal(vcv[1, 1], 1.5)  # B
        np.testing.assert_almost_equal(vcv[2, 2], 1.5)  # C
        # Off-diagonal: shared path lengths
        np.testing.assert_almost_equal(vcv[0, 1], 0.5)  # A-B share 0.5
        np.testing.assert_almost_equal(vcv[0, 2], 0.0)  # A-C share 0
        np.testing.assert_almost_equal(vcv[1, 2], 0.0)  # B-C share 0

    def test_hybrid_changes_vcv(self):
        """Adding a hybrid edge changes VCV; C's covariance with B increases."""
        nodes, parents, tip_map = _build_simple_dag()
        ordered_tips = sorted(tip_map.values())
        vcv_tree = NetworkSignal._compute_network_vcv(nodes, parents, tip_map, ordered_tips)

        # Re-build and add hybrid
        nodes2, parents2, tip_map2 = _build_simple_dag()
        NetworkSignal._add_hybrid_edge(nodes2, parents2, tip_map2, "B", "C", gamma=0.3)
        vcv_net = NetworkSignal._compute_network_vcv(nodes2, parents2, tip_map2, ordered_tips)

        # C's covariance with B should increase (B is donor for C)
        idx_b = ordered_tips.index("B")
        idx_c = ordered_tips.index("C")
        assert vcv_net[idx_b, idx_c] > vcv_tree[idx_b, idx_c]

    def test_hybrid_vcv_symmetry(self):
        """Network VCV is symmetric with 5-taxon tree and cross-clade hybrid B->D."""
        nodes, parents, tip_map = _build_five_taxon_dag()
        NetworkSignal._add_hybrid_edge(nodes, parents, tip_map, "B", "D", gamma=0.3)
        ordered_tips = sorted(tip_map.values())
        vcv = NetworkSignal._compute_network_vcv(nodes, parents, tip_map, ordered_tips)
        np.testing.assert_array_almost_equal(vcv, vcv.T)

    def test_hybrid_vcv_positive_definite(self):
        """All eigenvalues of the network VCV are positive with 5-taxon cross-clade hybrid."""
        nodes, parents, tip_map = _build_five_taxon_dag()
        NetworkSignal._add_hybrid_edge(nodes, parents, tip_map, "B", "D", gamma=0.3)
        ordered_tips = sorted(tip_map.values())
        vcv = NetworkSignal._compute_network_vcv(nodes, parents, tip_map, ordered_tips)
        eigenvalues = np.linalg.eigvalsh(vcv)
        assert all(ev > 0 for ev in eigenvalues)

    def test_gamma_zero_approaches_tree(self):
        """gamma=0.001 -> VCV approximately equals tree VCV (within decimal=2)."""
        nodes_tree, parents_tree, tip_map_tree = _build_simple_dag()
        ordered_tips = sorted(tip_map_tree.values())
        vcv_tree = NetworkSignal._compute_network_vcv(
            nodes_tree, parents_tree, tip_map_tree, ordered_tips
        )

        nodes_net, parents_net, tip_map_net = _build_simple_dag()
        NetworkSignal._add_hybrid_edge(
            nodes_net, parents_net, tip_map_net, "B", "C", gamma=0.001
        )
        vcv_net = NetworkSignal._compute_network_vcv(
            nodes_net, parents_net, tip_map_net, ordered_tips
        )

        np.testing.assert_almost_equal(vcv_net, vcv_tree, decimal=2)


# ===========================================================================
# TestQuartetJsonInference
# ===========================================================================

class TestQuartetJsonInference:
    def test_identify_swap_pair(self):
        """{A,B}|{D,E} dom + {A,D}|{B,E} minor -> swap = {B,D}."""
        dominant = (frozenset({"A", "B"}), frozenset({"D", "E"}))
        minor = (frozenset({"A", "D"}), frozenset({"B", "E"}))
        swap = NetworkSignal._identify_swap_pair(dominant, minor)
        assert swap == frozenset({"B", "D"})

    def test_top_two_topology_indices_preserves_stable_ties(self):
        assert NetworkSignal._top_two_topology_indices([5, 5, 4]) == (0, 1)
        assert NetworkSignal._top_two_topology_indices([5, 4, 5]) == (0, 2)
        assert NetworkSignal._top_two_topology_indices([4, 5, 5]) == (1, 2)
        assert NetworkSignal._top_two_topology_indices([3, 5, 4]) == (1, 2)

    def test_infer_from_quartet_json_caches_repeated_topology_strings(self, mocker):
        quartet_data = {
            "quartets": [
                {
                    "taxa": ["A", "B", "C", "D"],
                    "classification": "hybrid",
                    "counts": [40, 40, 10],
                    "dominant_topology": "{A, B} | {C, D}",
                },
                {
                    "taxa": ["A", "B", "C", "D"],
                    "classification": "hybrid",
                    "counts": [40, 40, 10],
                    "dominant_topology": "{A, B} | {C, D}",
                },
            ]
        }
        parser_spy = mocker.spy(NetworkSignal, "_parse_topology_string")

        edges = NetworkSignal._infer_hybrid_edges(quartet_data)

        assert parser_spy.call_count == 1
        assert edges == [
            {"donor": "B", "recipient": "C", "gamma": 0.4444, "n_quartets": 2}
        ]

    def test_infer_from_quartet_json_uses_plain_dict_pair_counts(self, monkeypatch):
        quartet_data = {
            "quartets": [
                {
                    "taxa": ["A", "B", "C", "D"],
                    "classification": "hybrid",
                    "counts": [40, 40, 10],
                    "dominant_topology": "{A, B} | {C, D}",
                },
            ]
        }

        def fail_counter(*_args, **_kwargs):
            raise AssertionError("hybrid edge inference should not instantiate Counter")

        monkeypatch.setattr(network_signal_module, "Counter", fail_counter, raising=False)

        edges = NetworkSignal._infer_hybrid_edges(quartet_data)

        assert edges == [
            {"donor": "B", "recipient": "C", "gamma": 0.4444, "n_quartets": 1}
        ]

    def test_infer_from_quartet_json_still_validates_topology_strings(self):
        quartet_data = {
            "quartets": [
                {
                    "taxa": ["A", "B", "C", "D"],
                    "classification": "hybrid",
                    "counts": [40, 40, 10],
                    "dominant_topology": "not a topology",
                },
            ]
        }

        with pytest.raises(PhykitUserError):
            NetworkSignal._infer_hybrid_edges(quartet_data)

    def test_infer_from_quartet_json(self):
        """Load quartet data fixture with a hybrid quartet, get >=1 edge."""
        # Create fake quartet_data with at least one hybrid classification
        quartet_data = {
            "quartets": [
                {
                    "taxa": ["A", "B", "D", "E"],
                    "classification": "hybrid",
                    "counts": [45, 35, 20],
                    "dominant_topology": "{A, B} | {D, E}",
                },
                {
                    "taxa": ["A", "B", "C", "D"],
                    "classification": "tree",
                    "counts": [8, 1, 1],
                    "dominant_topology": "{A, B} | {C, D}",
                },
            ]
        }
        edges = NetworkSignal._infer_hybrid_edges(quartet_data)
        assert len(edges) >= 1
        # The swap pair for the hybrid quartet above: dominant is {A,B}|{D,E}
        # The two elevated topologies are topo0 (dominant, counts[0]=45) and topo1 (counts[1]=35)
        # Topo 0: {A,B}|{D,E}, Topo 1: {A,D}|{B,E}
        # Swap pair from dominant {A,B}|{D,E} to minor {A,D}|{B,E}: B and D switch
        any_has_b_or_d = any(
            "B" in (e["donor"], e["recipient"]) or "D" in (e["donor"], e["recipient"])
            for e in edges
        )
        assert any_has_b_or_d
        # Each edge should include n_quartets key from aggregation
        for e in edges:
            assert "n_quartets" in e
            assert e["n_quartets"] >= 1

    def test_infer_from_quartet_json_uses_stable_tie_order(self):
        quartet_data = {
            "quartets": [
                {
                    "taxa": ["A", "B", "C", "D"],
                    "classification": "hybrid",
                    "counts": [40, 40, 10],
                    "dominant_topology": "{A, B} | {C, D}",
                },
            ]
        }

        edges = NetworkSignal._infer_hybrid_edges(quartet_data)

        assert edges == [
            {"donor": "B", "recipient": "C", "gamma": 0.4444, "n_quartets": 1}
        ]

    def test_infer_from_quartet_json_handles_all_topology_index_pairs(self):
        quartet_data = {
            "quartets": [
                {
                    "taxa": ["A", "B", "C", "D"],
                    "classification": "hybrid",
                    "counts": counts,
                    "dominant_topology": topology,
                }
                for counts, topology in [
                    ([5, 4, 1], "{A, B} | {C, D}"),
                    ([5, 1, 4], "{A, B} | {C, D}"),
                    ([4, 5, 1], "{A, C} | {B, D}"),
                    ([1, 5, 4], "{A, C} | {B, D}"),
                    ([4, 1, 5], "{A, D} | {B, C}"),
                    ([1, 4, 5], "{A, D} | {B, C}"),
                ]
            ]
        }

        edges = NetworkSignal._infer_hybrid_edges(quartet_data)

        assert edges == [
            {"donor": "B", "recipient": "C", "gamma": 0.4, "n_quartets": 2},
            {"donor": "B", "recipient": "D", "gamma": 0.4, "n_quartets": 2},
            {"donor": "C", "recipient": "D", "gamma": 0.4, "n_quartets": 2},
        ]

    def test_no_hybrid_quartets_returns_empty(self):
        """All 'tree' quartets -> empty list."""
        quartet_data = {
            "quartets": [
                {
                    "taxa": ["A", "B", "C", "D"],
                    "classification": "tree",
                    "counts": [8, 1, 1],
                    "dominant_topology": "{A, B} | {C, D}",
                },
                {
                    "taxa": ["A", "B", "C", "E"],
                    "classification": "tree",
                    "counts": [9, 0, 1],
                    "dominant_topology": "{A, B} | {C, E}",
                },
            ]
        }
        edges = NetworkSignal._infer_hybrid_edges(quartet_data)
        assert edges == []


# ===========================================================================
# TestBlombergsKNetwork
# ===========================================================================

class TestBlombergsKNetwork:
    def test_k_on_tree_vcv(self):
        """K works on a tree VCV, returns K > 0 and a p_value."""
        nodes, parents, tip_map = _build_simple_dag()
        ordered_tips = sorted(tip_map.values())
        vcv = NetworkSignal._compute_network_vcv(nodes, parents, tip_map, ordered_tips)
        x = np.array([1.0, 1.5, 3.0])  # trait values for A, B, C
        result = NetworkSignal._blombergs_k(x, vcv, n_perm=100)
        assert result["K"] > 0
        assert 0 <= result["p_value"] <= 1
        assert result["permutations"] == 100

    def test_k_changes_with_hybrid(self):
        """K differs between tree and network VCV."""
        # Tree VCV
        nodes_t, parents_t, tip_map_t = _build_simple_dag()
        ordered_tips = sorted(tip_map_t.values())
        vcv_tree = NetworkSignal._compute_network_vcv(
            nodes_t, parents_t, tip_map_t, ordered_tips
        )

        # Network VCV
        nodes_n, parents_n, tip_map_n = _build_simple_dag()
        NetworkSignal._add_hybrid_edge(
            nodes_n, parents_n, tip_map_n, "B", "C", gamma=0.3
        )
        vcv_net = NetworkSignal._compute_network_vcv(
            nodes_n, parents_n, tip_map_n, ordered_tips
        )

        x = np.array([1.0, 1.5, 3.0])
        k_tree = NetworkSignal._blombergs_k(x, vcv_tree, n_perm=100)
        k_net = NetworkSignal._blombergs_k(x, vcv_net, n_perm=100)
        assert k_tree["K"] != pytest.approx(k_net["K"], abs=1e-6)

    def test_batched_permutations_match_scalar_loop(self):
        """Batched K permutations preserve the previous scalar calculations."""
        nodes, parents, tip_map = _build_five_taxon_dag()
        ordered_tips = sorted(tip_map.values())
        vcv = NetworkSignal._compute_network_vcv(nodes, parents, tip_map, ordered_tips)
        x = np.array([1.0, 1.5, 3.0, 2.0, 2.5])

        C_inv = np.linalg.inv(vcv)
        ones = np.ones(len(x))
        sum_C_inv = float(ones @ C_inv @ ones)
        expected_ratio = (np.trace(vcv) - len(x) / sum_C_inv) / (len(x) - 1)

        rng = np.random.default_rng(seed=42)
        scalar = np.empty(25)
        for perm_i in range(25):
            x_perm = rng.permutation(x)
            a_p = float((ones @ C_inv @ x_perm) / sum_C_inv)
            e_p = x_perm - a_p
            obs_p = float(e_p @ e_p) / float(e_p @ C_inv @ e_p)
            scalar[perm_i] = obs_p / expected_ratio

        batched = NetworkSignal._blombergs_k_permutations(
            x,
            C_inv,
            sum_C_inv,
            expected_ratio,
            n_perm=25,
            batch_size=7,
        )

        assert batched == pytest.approx(scalar)

    def test_batched_permutations_use_one_quadratic_einsum_per_batch(self, monkeypatch):
        nodes, parents, tip_map = _build_five_taxon_dag()
        ordered_tips = sorted(tip_map.values())
        vcv = NetworkSignal._compute_network_vcv(nodes, parents, tip_map, ordered_tips)
        x = np.array([1.0, 1.5, 3.0, 2.0, 2.5])
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(len(x))
        sum_C_inv = float(ones @ C_inv @ ones)
        expected_ratio = (np.trace(vcv) - len(x) / sum_C_inv) / (len(x) - 1)
        einsum_calls = []
        real_einsum = np.einsum

        def tracking_einsum(*args, **kwargs):
            einsum_calls.append(args[0])
            return real_einsum(*args, **kwargs)

        monkeypatch.setattr(network_signal_module.np, "einsum", tracking_einsum)

        NetworkSignal._blombergs_k_permutations(
            x,
            C_inv,
            sum_C_inv,
            expected_ratio,
            n_perm=8,
            batch_size=3,
        )

        assert einsum_calls == ["ij,ij->i", "ij,ij->i", "ij,ij->i"]


# ===========================================================================
# TestPagelsLambdaNetwork
# ===========================================================================

class TestPagelsLambdaNetwork:
    def test_lambda_on_tree_vcv(self):
        """Lambda works on a tree VCV, returns 0 <= lambda <= 1."""
        nodes, parents, tip_map = _build_simple_dag()
        ordered_tips = sorted(tip_map.values())
        vcv = NetworkSignal._compute_network_vcv(nodes, parents, tip_map, ordered_tips)
        x = np.array([1.0, 1.5, 3.0])
        result = NetworkSignal._pagels_lambda(x, vcv, max_lambda=1.0)
        assert 0 <= result["lambda"] <= 1
        assert "log_likelihood" in result
        assert "p_value" in result

    def test_lambda_p_value_does_not_import_scipy_stats(self, monkeypatch):
        nodes, parents, tip_map = _build_simple_dag()
        ordered_tips = sorted(tip_map.values())
        vcv = NetworkSignal._compute_network_vcv(nodes, parents, tip_map, ordered_tips)
        x = np.array([1.0, 1.5, 3.0])
        real_import = __import__

        def fail_scipy_stats_import(name, *args, **kwargs):
            if name == "scipy.stats" or name.startswith("scipy.stats."):
                raise AssertionError("lambda p-value should not import scipy.stats")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr("builtins.__import__", fail_scipy_stats_import)

        result = NetworkSignal._pagels_lambda(x, vcv, max_lambda=1.0)

        assert 0 <= result["p_value"] <= 1

    def test_log_likelihood_cholesky_matches_inverse(self):
        """Cholesky likelihood preserves the inverse-based likelihood."""
        nodes, parents, tip_map = _build_five_taxon_dag()
        ordered_tips = sorted(tip_map.values())
        vcv = NetworkSignal._compute_network_vcv(nodes, parents, tip_map, ordered_tips)
        x = np.array([1.0, 1.5, 3.0, 2.0, 2.5])

        ll_chol, a_chol = NetworkSignal._log_likelihood_cholesky(x, vcv)
        ll_inv, a_inv = NetworkSignal._log_likelihood_inverse(x, vcv)

        assert ll_chol == pytest.approx(ll_inv)
        assert a_chol == pytest.approx(a_inv)

    def test_log_likelihood_cholesky_uses_single_solve(self, monkeypatch):
        rng = np.random.default_rng(20260628)
        n = 12
        a_matrix = rng.normal(size=(n, n))
        vcv = a_matrix @ a_matrix.T + np.eye(n)
        x = rng.normal(size=n)
        original_cho_solve = network_signal_module.cho_solve
        solve_calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal solve_calls
            solve_calls += 1
            return original_cho_solve(*args, **kwargs)

        monkeypatch.setattr(network_signal_module, "cho_solve", counting_cho_solve)

        ll_chol, a_chol = NetworkSignal._log_likelihood_cholesky(x, vcv)
        ll_inv, a_inv = NetworkSignal._log_likelihood_inverse(x, vcv)

        assert solve_calls == 1
        assert ll_chol == pytest.approx(ll_inv)
        assert a_chol == pytest.approx(a_inv)

    def test_lambda_changes_with_hybrid(self):
        """Lambda differs between tree and network VCV."""
        # Tree VCV
        nodes_t, parents_t, tip_map_t = _build_simple_dag()
        ordered_tips = sorted(tip_map_t.values())
        vcv_tree = NetworkSignal._compute_network_vcv(
            nodes_t, parents_t, tip_map_t, ordered_tips
        )

        # Network VCV
        nodes_n, parents_n, tip_map_n = _build_simple_dag()
        NetworkSignal._add_hybrid_edge(
            nodes_n, parents_n, tip_map_n, "B", "C", gamma=0.3
        )
        vcv_net = NetworkSignal._compute_network_vcv(
            nodes_n, parents_n, tip_map_n, ordered_tips
        )

        x = np.array([1.0, 1.5, 3.0])
        lam_tree = NetworkSignal._pagels_lambda(x, vcv_tree, max_lambda=1.0)
        lam_net = NetworkSignal._pagels_lambda(x, vcv_net, max_lambda=1.0)
        # The lambda values or log-likelihoods should differ
        assert (
            lam_tree["lambda"] != pytest.approx(lam_net["lambda"], abs=1e-6)
            or lam_tree["log_likelihood"]
            != pytest.approx(lam_net["log_likelihood"], abs=1e-6)
        )


class TestPolytomyHandling:
    """Polytomous nodes are represented as star topologies in the DAG."""

    def test_star_tree_creates_valid_dag(self):
        """A 4-way star tree produces a valid DAG."""
        tree = _make_tree("(A:1,B:1,C:1,D:1);")
        nodes, parents, tip_map = NetworkSignal._tree_to_dag(tree)
        assert len(tip_map) == 4
        assert set(tip_map.values()) == {"A", "B", "C", "D"}
        # All tips should have exactly one parent (the root)
        for tip_id in tip_map:
            assert len(parents[tip_id]) == 1

    def test_polytomy_vcv_equal_covariance(self):
        """Tips in a star polytomy should have equal off-diagonal covariance."""
        tree = _make_tree("(A:1,B:1,C:1);")
        nodes, parents, tip_map = NetworkSignal._tree_to_dag(tree)
        ordered_tips = sorted(tip_map.values())
        vcv = NetworkSignal._compute_network_vcv(
            nodes, parents, tip_map, ordered_tips
        )
        # All diagonal elements should be equal (same branch length)
        assert vcv[0, 0] == pytest.approx(vcv[1, 1])
        assert vcv[1, 1] == pytest.approx(vcv[2, 2])
        # All off-diagonal elements should be equal (same shared ancestry)
        assert vcv[0, 1] == pytest.approx(vcv[0, 2])
        assert vcv[0, 1] == pytest.approx(vcv[1, 2])
