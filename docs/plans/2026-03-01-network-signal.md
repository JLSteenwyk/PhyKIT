# Network Phylogenetic Signal Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Compute Blomberg's K and Pagel's lambda on phylogenetic networks using the Bastide et al. 2018 network VCV algorithm.

**Architecture:** Convert a species tree + hybrid edge specifications into an internal DAG, compute the network VCV matrix recursively, then feed it into the existing K/lambda formulas from `phylogenetic_signal.py`. Two input modes: explicit hybrid edges (`--hybrid`) and auto-inference from quartet_network JSON (`--quartet-json`).

**Tech Stack:** Python 3, numpy, scipy (optimize, stats), Bio.Phylo (tree parsing only)

---

### Task 1: Create test fixture files

**Files:**
- Create: `tests/sample_files/network_signal_tree.tre`
- Create: `tests/sample_files/network_signal_traits.tsv`
- Create: `tests/sample_files/network_signal_quartet.json`

**Step 1: Create species tree**

Create a rooted 5-taxon tree file:
```
(((A:1.0,B:1.0):0.5,C:1.5):0.5,(D:1.0,E:1.0):1.0);
```

**Step 2: Create trait data**

Create `network_signal_traits.tsv`:
```
A	2.5
B	2.8
C	1.2
D	5.1
E	4.9
```

**Step 3: Create quartet_network JSON fixture**

Create a minimal `network_signal_quartet.json` with 5 quartets (C(5,4)=5), where one is classified as "hybrid" with A and D as the swap pair. The hybrid quartet should have taxa [A, B, D, E] with:
- counts: [6, 3, 1] (dominant ab|de, elevated minor ad|be, lower ae|bd)
- cfs: [0.6, 0.3, 0.1]
- classification: "hybrid"
- dominant_topology: "{A, B} | {D, E}"

All other quartets should be "tree" classification.

Full JSON structure matching quartet_network output:
```json
{
  "input_tree_count": 10,
  "taxa": ["A", "B", "C", "D", "E"],
  "alpha": 0.05,
  "beta": 0.95,
  "total_quartets": 5,
  "tree_count": 4,
  "hybrid_count": 1,
  "unresolved_count": 0,
  "pruned_to_shared_taxa": false,
  "quartets": [
    {
      "taxa": ["A", "B", "C", "D"],
      "counts": [8, 1, 1],
      "cfs": [0.8, 0.1, 0.1],
      "classification": "tree",
      "p_star": 0.001,
      "p_tree": 0.5,
      "dominant_topology": "{A, B} | {C, D}"
    },
    {
      "taxa": ["A", "B", "C", "E"],
      "counts": [8, 1, 1],
      "cfs": [0.8, 0.1, 0.1],
      "classification": "tree",
      "p_star": 0.001,
      "p_tree": 0.5,
      "dominant_topology": "{A, B} | {C, E}"
    },
    {
      "taxa": ["A", "B", "D", "E"],
      "counts": [6, 3, 1],
      "cfs": [0.6, 0.3, 0.1],
      "classification": "hybrid",
      "p_star": 0.01,
      "p_tree": 0.02,
      "dominant_topology": "{A, B} | {D, E}"
    },
    {
      "taxa": ["A", "C", "D", "E"],
      "counts": [7, 2, 1],
      "cfs": [0.7, 0.2, 0.1],
      "classification": "tree",
      "p_star": 0.001,
      "p_tree": 0.3,
      "dominant_topology": "{A, C} | {D, E}"
    },
    {
      "taxa": ["B", "C", "D", "E"],
      "counts": [7, 2, 1],
      "cfs": [0.7, 0.2, 0.1],
      "classification": "tree",
      "p_star": 0.001,
      "p_tree": 0.3,
      "dominant_topology": "{B, C} | {D, E}"
    }
  ]
}
```

**Step 4: Commit**

```bash
git add tests/sample_files/network_signal_*
git commit -m "test: add fixture files for network_signal command"
```

---

### Task 2: Unit tests + Network VCV computation (TDD)

**Files:**
- Create: `tests/unit/services/tree/test_network_signal.py`
- Create: `phykit/services/tree/network_signal.py`

This task implements the core network VCV algorithm and the supporting network construction code, driven by unit tests.

**Unit test classes to create:**

#### TestTreeToNetwork

Test converting a Biopython tree to the internal DAG representation.

```python
import pytest
import numpy as np
from io import StringIO
from Bio import Phylo

from phykit.services.tree.network_signal import NetworkSignal


class TestTreeToNetwork:
    def test_simple_three_taxon(self):
        """Convert ((A:1,B:1):0.5,C:1.5); to DAG."""
        tree = Phylo.read(StringIO("((A:1,B:1):0.5,C:1.5);"), "newick")
        nodes, parents, tip_map = NetworkSignal._tree_to_dag(tree)
        # Root has no parents
        root_id = nodes[0]
        assert parents[root_id] == []
        # Tips exist
        assert "A" in tip_map.values()
        assert "B" in tip_map.values()
        assert "C" in tip_map.values()
        # All non-root nodes have exactly one parent
        for nid in nodes[1:]:
            assert len(parents[nid]) == 1

    def test_five_taxon(self):
        """Convert (((A:1,B:1):0.5,C:1.5):0.5,(D:1,E:1):1); to DAG."""
        tree = Phylo.read(
            StringIO("(((A:1,B:1):0.5,C:1.5):0.5,(D:1,E:1):1);"), "newick"
        )
        nodes, parents, tip_map = NetworkSignal._tree_to_dag(tree)
        assert len(tip_map) == 5
        # 5 tips + 3 internal + 1 root = 9 nodes
        assert len(nodes) == 9
```

#### TestAddHybridEdge

Test adding hybrid edges to the DAG.

```python
class TestAddHybridEdge:
    def test_add_single_hybrid(self):
        """Add hybrid edge B->C with gamma=0.3."""
        tree = Phylo.read(StringIO("((A:1,B:1):0.5,C:1.5);"), "newick")
        nodes, parents, tip_map = NetworkSignal._tree_to_dag(tree)
        name_to_id = {v: k for k, v in tip_map.items()}

        NetworkSignal._add_hybrid_edge(
            nodes, parents, tip_map, "B", "C", 0.3
        )

        # C should now have 2 parents (hybrid node)
        c_id = name_to_id["C"]
        assert len(parents[c_id]) == 2
        gammas = [p[2] for p in parents[c_id]]
        assert sorted(gammas) == pytest.approx([0.3, 0.7])

    def test_invalid_donor_raises(self):
        """Donor taxon not in tree should raise error."""
        tree = Phylo.read(StringIO("((A:1,B:1):0.5,C:1.5);"), "newick")
        nodes, parents, tip_map = NetworkSignal._tree_to_dag(tree)
        with pytest.raises(SystemExit):
            NetworkSignal._add_hybrid_edge(
                nodes, parents, tip_map, "Z", "C", 0.3
            )

    def test_gamma_out_of_range_raises(self):
        """Gamma >= 0.5 or <= 0 should raise error."""
        tree = Phylo.read(StringIO("((A:1,B:1):0.5,C:1.5);"), "newick")
        nodes, parents, tip_map = NetworkSignal._tree_to_dag(tree)
        with pytest.raises(SystemExit):
            NetworkSignal._add_hybrid_edge(
                nodes, parents, tip_map, "B", "C", 0.6
            )
```

#### TestNetworkVCV

Test the network VCV computation.

```python
class TestNetworkVCV:
    def test_tree_only_matches_standard(self):
        """With no hybrid edges, network VCV should match tree VCV."""
        tree = Phylo.read(StringIO("((A:1,B:1):0.5,C:1.5);"), "newick")
        nodes, parents, tip_map = NetworkSignal._tree_to_dag(tree)
        ordered_tips = sorted(tip_map.values())
        vcv = NetworkSignal._compute_network_vcv(nodes, parents, tip_map, ordered_tips)

        # Expected tree VCV for ((A:1,B:1):0.5,C:1.5):
        # V[A,A] = 1.5, V[B,B] = 1.5, V[C,C] = 1.5
        # V[A,B] = 0.5 (shared root-to-AB-ancestor = 0.5)
        # V[A,C] = V[B,C] = 0.0 (diverge at root)
        expected = np.array([
            [1.5, 0.5, 0.0],
            [0.5, 1.5, 0.0],
            [0.0, 0.0, 1.5],
        ])
        np.testing.assert_array_almost_equal(vcv, expected)

    def test_hybrid_changes_vcv(self):
        """Adding a hybrid edge should change the VCV matrix."""
        tree = Phylo.read(StringIO("((A:1,B:1):0.5,C:1.5);"), "newick")
        nodes, parents, tip_map = NetworkSignal._tree_to_dag(tree)
        ordered_tips = sorted(tip_map.values())

        vcv_tree = NetworkSignal._compute_network_vcv(
            nodes, parents, tip_map, ordered_tips
        )

        # Add hybrid: B -> C with gamma=0.3
        NetworkSignal._add_hybrid_edge(nodes, parents, tip_map, "B", "C", 0.3)
        vcv_net = NetworkSignal._compute_network_vcv(
            nodes, parents, tip_map, ordered_tips
        )

        # VCV should differ
        assert not np.allclose(vcv_tree, vcv_net)

        # C's variance should change (hybrid lineage adds variance)
        # C's covariance with A and B should increase (shared ancestry via B's parent)
        c_idx = ordered_tips.index("C")
        a_idx = ordered_tips.index("A")
        b_idx = ordered_tips.index("B")
        assert vcv_net[c_idx, b_idx] > vcv_tree[c_idx, b_idx]

    def test_hybrid_vcv_symmetry(self):
        """Network VCV should always be symmetric."""
        tree = Phylo.read(
            StringIO("(((A:1,B:1):0.5,C:1.5):0.5,(D:1,E:1):1);"), "newick"
        )
        nodes, parents, tip_map = NetworkSignal._tree_to_dag(tree)
        NetworkSignal._add_hybrid_edge(nodes, parents, tip_map, "B", "D", 0.2)
        ordered_tips = sorted(tip_map.values())
        vcv = NetworkSignal._compute_network_vcv(nodes, parents, tip_map, ordered_tips)
        np.testing.assert_array_almost_equal(vcv, vcv.T)

    def test_hybrid_vcv_positive_definite(self):
        """Network VCV should be positive definite."""
        tree = Phylo.read(
            StringIO("(((A:1,B:1):0.5,C:1.5):0.5,(D:1,E:1):1);"), "newick"
        )
        nodes, parents, tip_map = NetworkSignal._tree_to_dag(tree)
        NetworkSignal._add_hybrid_edge(nodes, parents, tip_map, "B", "D", 0.2)
        ordered_tips = sorted(tip_map.values())
        vcv = NetworkSignal._compute_network_vcv(nodes, parents, tip_map, ordered_tips)
        eigenvalues = np.linalg.eigvalsh(vcv)
        assert all(ev > 0 for ev in eigenvalues)

    def test_gamma_zero_equals_tree(self):
        """With gamma approaching 0, network VCV should approach tree VCV."""
        tree = Phylo.read(StringIO("((A:1,B:1):0.5,C:1.5);"), "newick")

        # Tree VCV
        nodes1, parents1, tip_map1 = NetworkSignal._tree_to_dag(tree)
        ordered_tips = sorted(tip_map1.values())
        vcv_tree = NetworkSignal._compute_network_vcv(
            nodes1, parents1, tip_map1, ordered_tips
        )

        # Network with tiny gamma
        tree2 = Phylo.read(StringIO("((A:1,B:1):0.5,C:1.5);"), "newick")
        nodes2, parents2, tip_map2 = NetworkSignal._tree_to_dag(tree2)
        NetworkSignal._add_hybrid_edge(nodes2, parents2, tip_map2, "B", "C", 0.001)
        vcv_net = NetworkSignal._compute_network_vcv(
            nodes2, parents2, tip_map2, ordered_tips
        )

        np.testing.assert_array_almost_equal(vcv_tree, vcv_net, decimal=2)
```

#### TestQuartetJsonInference

Test inferring hybrid edges from quartet_network JSON.

```python
class TestQuartetJsonInference:
    def test_identify_swap_pair(self):
        """Identify swap pair from dominant and elevated minor topologies."""
        # Dominant: {A, B} | {D, E}, Minor: {A, D} | {B, E}
        # A stays left, B moves right->left for D, D moves right->left for B
        # Swap pair = (B, D) — B and D switch sides
        dominant = (frozenset({"A", "B"}), frozenset({"D", "E"}))
        minor = (frozenset({"A", "D"}), frozenset({"B", "E"}))
        swap = NetworkSignal._identify_swap_pair(dominant, minor)
        assert swap == frozenset({"B", "D"})

    def test_infer_from_quartet_json(self):
        """Infer hybrid edges from quartet_network JSON with one hybrid quartet."""
        import json
        with open("tests/sample_files/network_signal_quartet.json") as f:
            data = json.load(f)
        edges = NetworkSignal._infer_hybrid_edges(data)
        # Should find one hybrid edge
        assert len(edges) >= 1
        # The swap pair for {A,B}|{D,E} dominant, {A,D}|{B,E} minor is (B,D)
        donors_recipients = {(e["donor"], e["recipient"]) for e in edges}
        swap_taxa = set()
        for d, r in donors_recipients:
            swap_taxa.add(d)
            swap_taxa.add(r)
        # B and D should be involved
        assert "B" in swap_taxa or "D" in swap_taxa

    def test_no_hybrid_quartets_returns_empty(self):
        """If no hybrid quartets, return empty list."""
        data = {
            "quartets": [
                {
                    "taxa": ["A", "B", "C", "D"],
                    "counts": [8, 1, 1],
                    "cfs": [0.8, 0.1, 0.1],
                    "classification": "tree",
                    "dominant_topology": "{A, B} | {C, D}",
                }
            ]
        }
        edges = NetworkSignal._infer_hybrid_edges(data)
        assert edges == []
```

#### TestBlombergsKNetwork and TestPagelsLambdaNetwork

Test that K and lambda work with network VCV (reusing existing formulas).

```python
class TestBlombergsKNetwork:
    def test_k_on_tree_vcv(self):
        """K should work on a standard tree VCV."""
        tree = Phylo.read(StringIO("((A:1,B:1):0.5,C:1.5);"), "newick")
        nodes, parents, tip_map = NetworkSignal._tree_to_dag(tree)
        ordered_tips = sorted(tip_map.values())
        vcv = NetworkSignal._compute_network_vcv(nodes, parents, tip_map, ordered_tips)
        x = np.array([1.0, 1.2, 5.0])  # A, B, C
        result = NetworkSignal._blombergs_k(x, vcv, 100)
        assert "K" in result
        assert "p_value" in result
        assert result["K"] > 0

    def test_k_changes_with_hybrid(self):
        """K should differ between tree and network."""
        x = np.array([1.0, 1.2, 5.0])

        tree1 = Phylo.read(StringIO("((A:1,B:1):0.5,C:1.5);"), "newick")
        nodes1, parents1, tip_map1 = NetworkSignal._tree_to_dag(tree1)
        ordered = sorted(tip_map1.values())
        vcv1 = NetworkSignal._compute_network_vcv(nodes1, parents1, tip_map1, ordered)
        k1 = NetworkSignal._blombergs_k(x, vcv1, 100)

        tree2 = Phylo.read(StringIO("((A:1,B:1):0.5,C:1.5);"), "newick")
        nodes2, parents2, tip_map2 = NetworkSignal._tree_to_dag(tree2)
        NetworkSignal._add_hybrid_edge(nodes2, parents2, tip_map2, "B", "C", 0.3)
        vcv2 = NetworkSignal._compute_network_vcv(nodes2, parents2, tip_map2, ordered)
        k2 = NetworkSignal._blombergs_k(x, vcv2, 100)

        assert k1["K"] != pytest.approx(k2["K"], abs=0.01)


class TestPagelsLambdaNetwork:
    def test_lambda_on_tree_vcv(self):
        """Lambda should work on a standard tree VCV."""
        tree = Phylo.read(StringIO("((A:1,B:1):0.5,C:1.5);"), "newick")
        nodes, parents, tip_map = NetworkSignal._tree_to_dag(tree)
        ordered_tips = sorted(tip_map.values())
        vcv = NetworkSignal._compute_network_vcv(nodes, parents, tip_map, ordered_tips)
        x = np.array([1.0, 1.2, 5.0])
        result = NetworkSignal._pagels_lambda(x, vcv)
        assert "lambda" in result
        assert "p_value" in result
        assert 0 <= result["lambda"] <= 1.0

    def test_lambda_changes_with_hybrid(self):
        """Lambda should differ between tree and network."""
        x = np.array([1.0, 1.2, 5.0])

        tree1 = Phylo.read(StringIO("((A:1,B:1):0.5,C:1.5);"), "newick")
        nodes1, parents1, tip_map1 = NetworkSignal._tree_to_dag(tree1)
        ordered = sorted(tip_map1.values())
        vcv1 = NetworkSignal._compute_network_vcv(nodes1, parents1, tip_map1, ordered)
        lam1 = NetworkSignal._pagels_lambda(x, vcv1)

        tree2 = Phylo.read(StringIO("((A:1,B:1):0.5,C:1.5);"), "newick")
        nodes2, parents2, tip_map2 = NetworkSignal._tree_to_dag(tree2)
        NetworkSignal._add_hybrid_edge(nodes2, parents2, tip_map2, "B", "C", 0.3)
        vcv2 = NetworkSignal._compute_network_vcv(nodes2, parents2, tip_map2, ordered)
        lam2 = NetworkSignal._pagels_lambda(x, vcv2)

        # Lambda values should differ (network accounts for hybrid ancestry)
        assert lam1["lambda"] != pytest.approx(lam2["lambda"], abs=0.01)
```

**Step 2: Implement `network_signal.py`**

Create `phykit/services/tree/network_signal.py` with these methods:

```python
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class NetworkSignal(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.hybrid_specs = parsed["hybrid_specs"]
        self.quartet_json_path = parsed["quartet_json_path"]
        self.method = parsed["method"]
        self.permutations = parsed["permutations"]
        self.json_output = parsed["json_output"]
        self.verbose = parsed["verbose"]

    def process_args(self, args) -> Dict:
        hybrid_specs = []
        if getattr(args, "hybrid", None):
            for spec in args.hybrid:
                parts = spec.split(":")
                if len(parts) != 3:
                    raise PhykitUserError(
                        [
                            f"Invalid hybrid spec '{spec}'.",
                            "Format: donor:recipient:gamma (e.g., B:C:0.3)",
                        ]
                    )
                donor, recipient, gamma_str = parts
                try:
                    gamma = float(gamma_str)
                except ValueError:
                    raise PhykitUserError(
                        [f"Invalid gamma value '{gamma_str}' in hybrid spec '{spec}'."]
                    )
                hybrid_specs.append((donor, recipient, gamma))

        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            hybrid_specs=hybrid_specs,
            quartet_json_path=getattr(args, "quartet_json", None),
            method=getattr(args, "method", "both"),
            permutations=getattr(args, "permutations", 1000),
            json_output=getattr(args, "json", False),
            verbose=getattr(args, "verbose", False),
        )

    # ---------------------------------------------------------------
    # Tree to DAG conversion
    # ---------------------------------------------------------------

    @staticmethod
    def _tree_to_dag(tree):
        """Convert a Biopython tree to a DAG representation.

        Returns:
            nodes: list of node IDs in topological order (root first)
            parents: dict mapping node_id -> list of (parent_id, edge_length, gamma)
            tip_map: dict mapping node_id -> taxon_name (tips only)
        """
        all_clades = list(tree.find_clades(order="level"))
        clade_to_id = {}
        nodes = []
        for i, clade in enumerate(all_clades):
            clade_to_id[id(clade)] = i
            nodes.append(i)

        parents = {}
        for clade in all_clades:
            cid = clade_to_id[id(clade)]
            if clade == tree.root:
                parents[cid] = []
            # Build parent links from children
            for child in clade.clades:
                child_id = clade_to_id[id(child)]
                length = child.branch_length if child.branch_length is not None else 0.0
                parents[child_id] = [(cid, length, 1.0)]

        tip_map = {}
        for clade in tree.get_terminals():
            tip_map[clade_to_id[id(clade)]] = clade.name

        return nodes, parents, tip_map

    # ---------------------------------------------------------------
    # Add hybrid edges
    # ---------------------------------------------------------------

    @staticmethod
    def _add_hybrid_edge(nodes, parents, tip_map, donor, recipient, gamma):
        """Add a hybrid edge: recipient inherits gamma from donor's parent.

        The recipient node becomes a hybrid node with two parents:
        - Original tree parent (weight 1-gamma)
        - Donor's tree parent (weight gamma)
        """
        name_to_id = {v: k for k, v in tip_map.items()}

        if donor not in name_to_id:
            raise PhykitUserError(
                [f"Donor taxon '{donor}' not found in tree."]
            )
        if recipient not in name_to_id:
            raise PhykitUserError(
                [f"Recipient taxon '{recipient}' not found in tree."]
            )
        if gamma <= 0 or gamma >= 0.5:
            raise PhykitUserError(
                [
                    f"Gamma must be between 0 and 0.5 (exclusive), got {gamma}.",
                    "Gamma represents the minority inheritance from the donor lineage.",
                ]
            )

        donor_id = name_to_id[donor]
        recipient_id = name_to_id[recipient]

        # Get donor's tree parent
        if not parents[donor_id]:
            raise PhykitUserError([f"Donor '{donor}' has no parent (is it the root?)."])
        donor_parent_id, donor_edge_len, _ = parents[donor_id][0]

        # Get recipient's tree parent
        if not parents[recipient_id]:
            raise PhykitUserError(
                [f"Recipient '{recipient}' has no parent (is it the root?)."]
            )
        recip_parent_id, recip_edge_len, _ = parents[recipient_id][0]

        # Recipient becomes hybrid: two parents
        parents[recipient_id] = [
            (recip_parent_id, recip_edge_len, 1.0 - gamma),  # tree parent
            (donor_parent_id, donor_edge_len, gamma),          # hybrid parent
        ]

    # ---------------------------------------------------------------
    # Network VCV computation (Bastide et al. 2018)
    # ---------------------------------------------------------------

    @staticmethod
    def _compute_network_vcv(nodes, parents, tip_map, ordered_tips):
        """Compute the VCV matrix for tips under BM on a network.

        Uses the Bastide et al. 2018 recursive algorithm.
        Processes nodes in topological order, maintaining a full
        node-by-node covariance matrix.
        """
        n_nodes = len(nodes)
        node_idx = {nid: i for i, nid in enumerate(nodes)}

        V = np.zeros((n_nodes, n_nodes))

        for nid in nodes:
            i = node_idx[nid]
            p_list = parents.get(nid, [])

            if len(p_list) == 0:
                # Root
                V[i, i] = 0.0

            elif len(p_list) == 1:
                # Tree node
                parent_id, length, _ = p_list[0]
                pi = node_idx[parent_id]
                V[i, i] = V[pi, pi] + length
                for j in range(i):
                    V[i, j] = V[pi, j]
                    V[j, i] = V[pi, j]

            elif len(p_list) == 2:
                # Hybrid node
                p1_id, l1, gamma = p_list[0]
                p2_id, l2, gamma2 = p_list[1]
                pi1 = node_idx[p1_id]
                pi2 = node_idx[p2_id]

                V[i, i] = (
                    gamma ** 2 * (V[pi1, pi1] + l1)
                    + gamma2 ** 2 * (V[pi2, pi2] + l2)
                    + 2 * gamma * gamma2 * V[pi1, pi2]
                )

                for j in range(i):
                    v = gamma * V[pi1, j] + gamma2 * V[pi2, j]
                    V[i, j] = v
                    V[j, i] = v

        # Extract tip submatrix
        name_to_id = {v: k for k, v in tip_map.items()}
        tip_indices = [node_idx[name_to_id[name]] for name in ordered_tips]
        return V[np.ix_(tip_indices, tip_indices)]

    # ---------------------------------------------------------------
    # Quartet JSON inference
    # ---------------------------------------------------------------

    @staticmethod
    def _parse_topology_string(topo_str):
        """Parse '{A, B} | {C, D}' into (frozenset, frozenset)."""
        parts = topo_str.split("|")
        left = frozenset(
            t.strip() for t in parts[0].strip().strip("{}").split(",")
        )
        right = frozenset(
            t.strip() for t in parts[1].strip().strip("{}").split(",")
        )
        return (left, right)

    @staticmethod
    def _identify_swap_pair(dominant, minor):
        """Identify the swap pair between dominant and minor topologies.

        The swap pair consists of the two taxa that change sides.
        """
        dom_left, dom_right = dominant
        min_left, min_right = minor

        # Taxa that moved from dominant-left to minor-right (or vice versa)
        moved_from_left = dom_left - min_left
        moved_from_right = dom_right - min_right

        return frozenset(moved_from_left | moved_from_right)

    @staticmethod
    def _infer_hybrid_edges(quartet_data):
        """Infer hybrid edges from quartet_network JSON output.

        For each hybrid quartet:
        1. Find the dominant and elevated minor topologies
        2. Identify the swap pair
        3. Aggregate evidence across quartets
        4. Estimate gamma from CF ratios
        """
        swap_counts = Counter()
        swap_gammas = {}  # swap_pair -> list of gamma estimates

        for q in quartet_data.get("quartets", []):
            if q["classification"] != "hybrid":
                continue

            taxa = q["taxa"]
            cfs = q["cfs"]
            dominant_str = q["dominant_topology"]

            # Parse dominant topology
            dominant = NetworkSignal._parse_topology_string(dominant_str)

            # Find indices sorted by CF (descending)
            sorted_indices = sorted(range(3), key=lambda i: cfs[i], reverse=True)
            dominant_idx = sorted_indices[0]
            elevated_idx = sorted_indices[1]
            lower_idx = sorted_indices[2]

            # Build the three topologies for this quartet
            a, b, c, d = sorted(taxa)
            topologies = [
                (frozenset({a, b}), frozenset({c, d})),  # ab|cd
                (frozenset({a, c}), frozenset({b, d})),  # ac|bd
                (frozenset({a, d}), frozenset({b, c})),  # ad|bc
            ]

            minor_topo = topologies[elevated_idx]
            swap = NetworkSignal._identify_swap_pair(
                topologies[dominant_idx], minor_topo
            )

            if len(swap) == 2:
                swap_counts[swap] += 1
                # Estimate gamma from CF ratio
                cf_elevated = cfs[elevated_idx]
                cf_lower = cfs[lower_idx]
                if cf_elevated + cf_lower > 0:
                    gamma_est = cf_elevated / (cf_elevated + cf_lower)
                    # Clamp to (0, 0.5)
                    gamma_est = min(gamma_est, 0.499)
                    gamma_est = max(gamma_est, 0.001)
                    swap_gammas.setdefault(swap, []).append(gamma_est)

        if not swap_counts:
            return []

        # Select pairs with strongest evidence
        edges = []
        for swap_pair, count in swap_counts.most_common():
            gammas = swap_gammas.get(swap_pair, [0.1])
            avg_gamma = sum(gammas) / len(gammas)
            # Ensure gamma < 0.5
            avg_gamma = min(avg_gamma, 0.499)
            taxa_list = sorted(swap_pair)
            edges.append({
                "donor": taxa_list[0],
                "recipient": taxa_list[1],
                "gamma": round(avg_gamma, 4),
                "n_quartets": count,
            })

        return edges

    # ---------------------------------------------------------------
    # Trait file parsing (same pattern as phylogenetic_signal.py)
    # ---------------------------------------------------------------

    def _parse_trait_file(self, path, tree_tips):
        """Parse tab-separated trait file: taxon<TAB>value."""
        try:
            with open(path) as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ]
            )

        traits = {}
        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise PhykitUserError(
                    [
                        f"Line {line_num} in trait file has {len(parts)} columns; expected 2.",
                        "Each line should be: taxon_name<tab>trait_value",
                    ]
                )
            taxon, value_str = parts
            try:
                traits[taxon] = float(value_str)
            except ValueError:
                raise PhykitUserError(
                    [f"Non-numeric trait value '{value_str}' for taxon '{taxon}' on line {line_num}."]
                )

        tree_tip_set = set(tree_tips)
        trait_taxa_set = set(traits.keys())
        shared = tree_tip_set & trait_taxa_set

        tree_only = tree_tip_set - trait_taxa_set
        trait_only = trait_taxa_set - tree_tip_set
        if tree_only:
            print(
                f"Warning: {len(tree_only)} taxa in tree but not in trait file: "
                f"{', '.join(sorted(tree_only))}",
                file=sys.stderr,
            )
        if trait_only:
            print(
                f"Warning: {len(trait_only)} taxa in trait file but not in tree: "
                f"{', '.join(sorted(trait_only))}",
                file=sys.stderr,
            )

        if len(shared) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa between tree and trait file.",
                    "At least 3 shared taxa are required.",
                ]
            )

        return {taxon: traits[taxon] for taxon in shared}

    # ---------------------------------------------------------------
    # K and lambda (same math as phylogenetic_signal.py)
    # ---------------------------------------------------------------

    @staticmethod
    def _log_likelihood(x, C):
        """Concentrated log-likelihood with sigma^2 profiled out."""
        n = len(x)
        ones = np.ones(n)
        C_inv = np.linalg.inv(C)
        a_hat = float((ones @ C_inv @ x) / (ones @ C_inv @ ones))
        e = x - a_hat
        sig2 = float(e @ C_inv @ e / n)
        sign, logdet_C = np.linalg.slogdet(C)
        logdet_sig2C = n * np.log(sig2) + logdet_C
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet_sig2C + n)
        return ll, a_hat

    @staticmethod
    def _blombergs_k(x, vcv, n_perm):
        """Blomberg's K following phytools::phylosig(method='K')."""
        n = len(x)
        ones = np.ones(n)
        C = vcv.copy()
        C_inv = np.linalg.inv(C)
        sum_C_inv = float(ones @ C_inv @ ones)
        a_hat = float((ones @ C_inv @ x) / sum_C_inv)
        e = x - a_hat
        observed_ratio = float(e @ e) / float(e @ C_inv @ e)
        expected_ratio = (np.trace(C) - n / sum_C_inv) / (n - 1)
        K = observed_ratio / expected_ratio

        rng = np.random.default_rng(seed=42)
        k_perm = np.empty(n_perm)
        for perm_i in range(n_perm):
            x_perm = rng.permutation(x)
            a_p = float((ones @ C_inv @ x_perm) / sum_C_inv)
            e_p = x_perm - a_p
            obs_p = float(e_p @ e_p) / float(e_p @ C_inv @ e_p)
            k_perm[perm_i] = obs_p / expected_ratio

        p_value = float(np.mean(k_perm >= K))
        return dict(K=float(K), p_value=p_value, permutations=n_perm)

    @staticmethod
    def _pagels_lambda(x, vcv, max_lambda=1.0):
        """Pagel's lambda with multi-interval optimization."""
        n = len(x)
        diag_vals = np.diag(vcv).copy()
        niter = 10

        def neg_ll(lam):
            C_lam = vcv * lam
            np.fill_diagonal(C_lam, diag_vals)
            try:
                ll, _ = NetworkSignal._log_likelihood(x, C_lam)
                return -ll
            except (np.linalg.LinAlgError, FloatingPointError, ValueError):
                return 1e10

        bounds_lo = np.linspace(0, max_lambda - max_lambda / niter, niter)
        bounds_hi = np.linspace(max_lambda / niter, max_lambda, niter)

        best_ll = -np.inf
        lambda_hat = 0.0
        for lo, hi in zip(bounds_lo, bounds_hi):
            res = minimize_scalar(neg_ll, bounds=(lo, hi), method="bounded")
            ll_val = -res.fun
            if ll_val > best_ll:
                best_ll = ll_val
                lambda_hat = res.x

        C_fitted = vcv * lambda_hat
        np.fill_diagonal(C_fitted, diag_vals)
        ll_fitted, _ = NetworkSignal._log_likelihood(x, C_fitted)

        C_zero = vcv * 0.0
        np.fill_diagonal(C_zero, diag_vals)
        ll_zero, _ = NetworkSignal._log_likelihood(x, C_zero)

        lr_stat = max(2.0 * (ll_fitted - ll_zero), 0.0)
        p_value = float(chi2.sf(lr_stat, df=1))

        return dict(
            **{"lambda": float(lambda_hat)},
            log_likelihood=float(ll_fitted),
            p_value=p_value,
        )

    # ---------------------------------------------------------------
    # Main entry point
    # ---------------------------------------------------------------

    def run(self):
        tree = self.read_tree_file()

        # Validate tree
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for network signal analysis."]
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                raise PhykitUserError(
                    ["All branches in the tree must have lengths."]
                )

        # Parse traits
        tree_tips = self.get_tip_names_from_tree(tree)
        traits = self._parse_trait_file(self.trait_data_path, tree_tips)
        ordered_names = sorted(traits.keys())
        x = np.array([traits[name] for name in ordered_names])

        # Build network DAG
        nodes, parents, tip_map = self._tree_to_dag(tree)

        # Determine hybrid edges
        hybrid_edges = []
        if self.hybrid_specs:
            for donor, recipient, gamma in self.hybrid_specs:
                self._add_hybrid_edge(nodes, parents, tip_map, donor, recipient, gamma)
                hybrid_edges.append(
                    {"donor": donor, "recipient": recipient, "gamma": gamma}
                )
        elif self.quartet_json_path:
            qpath = Path(self.quartet_json_path)
            if not qpath.exists():
                raise PhykitUserError(
                    [f"{self.quartet_json_path} not found."]
                )
            with open(qpath) as f:
                quartet_data = json.load(f)
            inferred = self._infer_hybrid_edges(quartet_data)
            if not inferred:
                raise PhykitUserError(
                    [
                        "No hybrid quartets found in quartet JSON.",
                        "Cannot infer hybrid edges. Use --hybrid to specify manually.",
                    ]
                )
            for edge in inferred:
                self._add_hybrid_edge(
                    nodes, parents, tip_map,
                    edge["donor"], edge["recipient"], edge["gamma"],
                )
            hybrid_edges = inferred
        else:
            raise PhykitUserError(
                ["Provide --hybrid or --quartet-json to specify network reticulations."]
            )

        # Compute network VCV
        vcv = self._compute_network_vcv(nodes, parents, tip_map, ordered_names)

        # Compute signal metrics
        results = {"hybrid_edges": hybrid_edges, "n_taxa": len(ordered_names)}

        if self.method in ("both", "blombergs_k"):
            results["blombergs_k"] = self._blombergs_k(x, vcv, self.permutations)

        if self.method in ("both", "lambda"):
            results["pagels_lambda"] = self._pagels_lambda(x, vcv)

        # Output
        self._output(results)

    def _output(self, results):
        if self.json_output:
            print_json(results)
            return

        try:
            for edge in results["hybrid_edges"]:
                print(
                    f"Hybrid edge: {edge['donor']} -> {edge['recipient']} "
                    f"(gamma={edge['gamma']:.4f})"
                )
            print(f"Network taxa: {results['n_taxa']}")
            print("---")

            if "blombergs_k" in results:
                k = results["blombergs_k"]
                print(
                    f"Blomberg's K: {k['K']:.4f}\t"
                    f"p-value: {k['p_value']:.4f}"
                )

            if "pagels_lambda" in results:
                lam = results["pagels_lambda"]
                print(
                    f"Pagel's lambda: {lam['lambda']:.4f}\t"
                    f"log-likelihood: {lam['log_likelihood']:.4f}\t"
                    f"p-value: {lam['p_value']:.4f}"
                )
        except BrokenPipeError:
            pass
```

**Step 3: Run tests**

```bash
python -m pytest tests/unit/services/tree/test_network_signal.py -v
```

All tests should pass.

**Step 4: Commit**

```bash
git add phykit/services/tree/network_signal.py tests/unit/services/tree/test_network_signal.py
git commit -m "feat: add network_signal service with VCV computation and K/lambda"
```

---

### Task 3: Register command in PhyKIT CLI

**Files:**
- Modify: `phykit/service_factories.py` (line ~79)
- Modify: `phykit/cli_registry.py` (line ~116)
- Modify: `phykit/services/tree/__init__.py`
- Modify: `phykit/phykit.py` (handler + wrapper + help text)
- Modify: `setup.py` (entry points)

**Step 1: Add service factory**

In `phykit/service_factories.py`, after the `RelativeRateTest` line, add:
```python
NetworkSignal = _LazyServiceFactory("phykit.services.tree.network_signal", "NetworkSignal")
```

**Step 2: Add CLI aliases**

In `phykit/cli_registry.py`, add in the Tree aliases section:
```python
"netsig": "network_signal",
"net_signal": "network_signal",
```

**Step 3: Add to tree __init__.py**

In `phykit/services/tree/__init__.py`, add to `_EXPORTS`:
```python
"NetworkSignal": "network_signal",
```

**Step 4: Add handler in phykit.py**

Add a static method `network_signal` to the `Phykit` class (near the other tree handlers):

```python
@staticmethod
def network_signal(argv=None):
    """Phylogenetic signal on a network (Blomberg's K and Pagel's lambda)."""
    from .service_factories import NetworkSignal

    description = (
        "Compute phylogenetic signal (Blomberg's K and/or Pagel's lambda) "
        "on a phylogenetic network using the Bastide et al. (2018) VCV algorithm."
    )
    parser = _new_parser(description=description)
    parser.add_argument("-t", "--tree", type=str, required=True, help=SUPPRESS)
    parser.add_argument("-d", "--trait-data", type=str, required=True, help=SUPPRESS)

    hybrid_group = parser.add_mutually_exclusive_group(required=True)
    hybrid_group.add_argument(
        "--hybrid", nargs="+", type=str, help=SUPPRESS
    )
    hybrid_group.add_argument(
        "--quartet-json", type=str, help=SUPPRESS
    )

    parser.add_argument(
        "--method",
        type=str,
        choices=["both", "blombergs_k", "lambda"],
        default="both",
        required=False,
        help=SUPPRESS,
    )
    parser.add_argument(
        "--permutations", type=int, default=1000, required=False, help=SUPPRESS
    )
    parser.add_argument("-v", "--verbose", action="store_true", help=SUPPRESS)
    _add_json_argument(parser)
    _run_service(parser, argv, NetworkSignal)
```

Add module-level wrapper function:
```python
def network_signal(argv=None):
    Phykit.network_signal(sys.argv[1:])
```

Add help text in `__init__` banner:
```
network_signal (alias: netsig; net_signal)
    - Phylogenetic signal on a network
```

**Step 5: Add entry points in setup.py**

```python
"pk_network_signal = phykit.phykit:network_signal",
"pk_netsig = phykit.phykit:network_signal",
"pk_net_signal = phykit.phykit:network_signal",
```

**Step 6: Run existing tests to verify no regressions**

```bash
python -m pytest tests/unit/services/tree/test_network_signal.py -v
```

**Step 7: Commit**

```bash
git add phykit/service_factories.py phykit/cli_registry.py phykit/services/tree/__init__.py phykit/phykit.py setup.py
git commit -m "feat: register network_signal command with aliases and CLI handler"
```

---

### Task 4: Integration tests

**Files:**
- Create: `tests/integration/tree/test_network_signal_integration.py`

**Integration tests:**

```python
import json
import sys

import pytest
from mock import patch

from phykit.phykit import Phykit


@pytest.mark.integration
class TestNetworkSignalIntegration:
    @patch("builtins.print")
    def test_explicit_hybrid(self, mocked_print, tmp_path):
        """Run with explicit --hybrid specification."""
        tree = tmp_path / "tree.tre"
        tree.write_text("(((A:1,B:1):0.5,C:1.5):0.5,(D:1,E:1):1);")
        traits = tmp_path / "traits.tsv"
        traits.write_text("A\t2.5\nB\t2.8\nC\t1.2\nD\t5.1\nE\t4.9\n")

        testargs = [
            "phykit", "network_signal",
            "-t", str(tree), "-d", str(traits),
            "--hybrid", "B:D:0.3",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Hybrid edge: B -> D" in output
        assert "Blomberg's K:" in output
        assert "Pagel's lambda:" in output

    @patch("builtins.print")
    def test_netsig_alias(self, mocked_print, tmp_path):
        """Test netsig alias."""
        tree = tmp_path / "tree.tre"
        tree.write_text("((A:1,B:1):0.5,C:1.5);")
        traits = tmp_path / "traits.tsv"
        traits.write_text("A\t1.0\nB\t1.2\nC\t5.0\n")

        testargs = [
            "phykit", "netsig",
            "-t", str(tree), "-d", str(traits),
            "--hybrid", "A:C:0.2",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Blomberg's K:" in output

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        """Test JSON output mode."""
        tree = tmp_path / "tree.tre"
        tree.write_text("((A:1,B:1):0.5,C:1.5);")
        traits = tmp_path / "traits.tsv"
        traits.write_text("A\t1.0\nB\t1.2\nC\t5.0\n")

        testargs = [
            "phykit", "netsig",
            "-t", str(tree), "-d", str(traits),
            "--hybrid", "A:C:0.2",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "blombergs_k" in payload
        assert "pagels_lambda" in payload
        assert "hybrid_edges" in payload

    @patch("builtins.print")
    def test_quartet_json_mode(self, mocked_print, tmp_path):
        """Test --quartet-json mode."""
        tree = tmp_path / "tree.tre"
        tree.write_text("(((A:1,B:1):0.5,C:1.5):0.5,(D:1,E:1):1);")
        traits = tmp_path / "traits.tsv"
        traits.write_text("A\t2.5\nB\t2.8\nC\t1.2\nD\t5.1\nE\t4.9\n")

        # Use sample fixture
        import shutil
        qjson = tmp_path / "quartets.json"
        shutil.copy("tests/sample_files/network_signal_quartet.json", str(qjson))

        testargs = [
            "phykit", "network_signal",
            "-t", str(tree), "-d", str(traits),
            "--quartet-json", str(qjson),
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Hybrid edge:" in output
        assert "Blomberg's K:" in output

    @patch("builtins.print")
    def test_method_k_only(self, mocked_print, tmp_path):
        """Test --method blombergs_k."""
        tree = tmp_path / "tree.tre"
        tree.write_text("((A:1,B:1):0.5,C:1.5);")
        traits = tmp_path / "traits.tsv"
        traits.write_text("A\t1.0\nB\t1.2\nC\t5.0\n")

        testargs = [
            "phykit", "netsig",
            "-t", str(tree), "-d", str(traits),
            "--hybrid", "A:C:0.2",
            "--method", "blombergs_k",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Blomberg's K:" in output
        assert "Pagel's lambda:" not in output
```

**Step 2: Run all tests**

```bash
python -m pytest tests/unit/services/tree/test_network_signal.py tests/integration/tree/test_network_signal_integration.py -v
```

**Step 3: Commit**

```bash
git add tests/integration/tree/test_network_signal_integration.py
git commit -m "test: add integration tests for network_signal command"
```

---

### Task 5: Documentation and version bump

**Files:**
- Modify: `docs/usage/index.rst`
- Modify: `docs/change_log/index.rst`
- Modify: `phykit/version.py`

**Step 1: Add usage docs**

Add a section for `network_signal` in `docs/usage/index.rst` following the pattern of other commands. Include:
- Algorithm description (Bastide et al. 2018)
- Both input modes (--hybrid and --quartet-json)
- Example usage
- Output interpretation
- References

**Step 2: Add changelog entry**

Add to `docs/change_log/index.rst`:
```
2.1.26
======
* Added ``network_signal`` (alias: ``netsig``; ``net_signal``): computes Blomberg's K and
  Pagel's lambda on phylogenetic networks using the Bastide et al. (2018) variance-covariance
  algorithm. Accepts explicit hybrid edge specifications or auto-infers from ``quartet_network``
  JSON output. Blomberg's K on networks is a novel implementation — first available in any tool.
```

**Step 3: Bump version**

Update `phykit/version.py` to `2.1.26`.

**Step 4: Commit**

```bash
git add docs/usage/index.rst docs/change_log/index.rst phykit/version.py
git commit -m "docs: add network_signal usage docs and changelog for v2.1.26"
```

---

## Verification

```bash
# Run all network_signal tests
python -m pytest tests/unit/services/tree/test_network_signal.py tests/integration/tree/test_network_signal_integration.py -v

# Manual test with explicit hybrid
python -c "
from argparse import Namespace
from phykit.services.tree.network_signal import NetworkSignal
svc = NetworkSignal(Namespace(
    tree='tests/sample_files/network_signal_tree.tre',
    trait_data='tests/sample_files/network_signal_traits.tsv',
    hybrid=['B:D:0.3'],
    quartet_json=None,
    method='both',
    permutations=1000,
    json=False,
    verbose=False,
))
svc.run()
"

# Manual test with quartet JSON
python -c "
from argparse import Namespace
from phykit.services.tree.network_signal import NetworkSignal
svc = NetworkSignal(Namespace(
    tree='tests/sample_files/network_signal_tree.tre',
    trait_data='tests/sample_files/network_signal_traits.tsv',
    hybrid=None,
    quartet_json='tests/sample_files/network_signal_quartet.json',
    method='both',
    permutations=1000,
    json=False,
    verbose=False,
))
svc.run()
"

# Verify no-hybrid-edge case matches phylogenetic_signal
# (tree-only VCV should produce same K as phylogenetic_signal command)
```
