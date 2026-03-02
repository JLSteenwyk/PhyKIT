import json
import pytest
import numpy as np
from io import StringIO

from Bio import Phylo

from phykit.services.tree.network_signal import NetworkSignal
from phykit.errors import PhykitUserError


def _make_tree(newick):
    return Phylo.read(StringIO(newick), "newick")


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
