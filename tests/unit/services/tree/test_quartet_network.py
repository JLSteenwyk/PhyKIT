import itertools
from argparse import Namespace
from io import StringIO
import math
import subprocess
import sys

import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

from phykit.errors import PhykitUserError
import phykit.services.tree.quartet_network as quartet_network_module
from phykit.services.tree.quartet_network import QuartetNetwork, _position_extent


class SinglePassPositions:
    def __init__(self, positions):
        self.positions = positions
        self.iterations = 0

    def __iter__(self):
        self.iterations += 1
        if self.iterations > 1:
            raise AssertionError("positions should be scanned once")
        return iter(self.positions)


def test_position_extent_scans_positions_once():
    positions = SinglePassPositions([(1.0, -2.0), (4.0, 3.0), (-1.0, 0.0)])

    assert _position_extent(positions) == 5.0
    assert positions.iterations == 1
    assert _position_extent([]) == 0


def _write(path, content):
    path.write_text(content)
    return str(path)


def _make_tree(newick):
    return Phylo.read(StringIO(newick), "newick")


def _legacy_compute_split_directions(ordering, circular_splits):
    n = len(ordering)
    angles = {taxon: 2 * math.pi * i / n for i, taxon in enumerate(ordering)}
    directions = {}
    for split, count, freq in circular_splits:
        gap_positions = []
        for i in range(n):
            curr = ordering[i]
            nxt = ordering[(i + 1) % n]
            if (curr in split) != (nxt in split):
                gap_positions.append(i)
        if len(gap_positions) != 2:
            continue
        g1 = math.pi * (2 * gap_positions[0] + 1) / n
        g2 = math.pi * (2 * gap_positions[1] + 1) / n
        cx = math.cos(g2) - math.cos(g1)
        cy = math.sin(g2) - math.sin(g1)
        dx = -cy
        dy = cx
        length = math.sqrt(dx * dx + dy * dy)
        if length > 1e-10:
            dx /= length
            dy /= length
        else:
            dx, dy = 1.0, 0.0
        cx_split = sum(math.cos(angles[t]) for t in split) / len(split)
        cy_split = sum(math.sin(angles[t]) for t in split) / len(split)
        if dx * cx_split + dy * cy_split < 0:
            dx = -dx
            dy = -dy
        directions[split] = (dx, dy)
    return directions


def test_module_import_does_not_import_biophylo_or_numpy():
    code = """
import sys
import phykit.services.tree.quartet_network as module

assert hasattr(module.Phylo, "read")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "json" not in sys.modules
assert "numpy" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


class TestExtractBipartitions:
    def test_balanced_tree(self):
        tree = _make_tree("((A,B),(C,D));")
        all_taxa = frozenset({"A", "B", "C", "D"})
        bips = QuartetNetwork._extract_bipartitions(tree, all_taxa)
        # Both internal clades produce a bipartition for {A,B}|{C,D}
        assert len(bips) >= 1
        all_sides = set()
        for side_a, side_b in bips:
            all_sides.add(frozenset(side_a))
        assert frozenset({"A", "B"}) in all_sides or frozenset({"C", "D"}) in all_sides

    def test_caterpillar_tree(self):
        tree = _make_tree("(A,(B,(C,D)));")
        all_taxa = frozenset({"A", "B", "C", "D"})
        bips = QuartetNetwork._extract_bipartitions(tree, all_taxa)
        # Non-trivial bipartition: {C,D} | {A,B}
        assert len(bips) == 1
        sides = {frozenset(bips[0][0]), frozenset(bips[0][1])}
        assert frozenset({"C", "D"}) in sides

    def test_six_taxa_tree(self):
        tree = _make_tree("((A,B),((C,D),(E,F)));")
        all_taxa = frozenset({"A", "B", "C", "D", "E", "F"})
        bips = QuartetNetwork._extract_bipartitions(tree, all_taxa)
        # Non-trivial bipartitions include {A,B}, {C,D,E,F}, {C,D}, {E,F}
        assert len(bips) >= 3
        all_sides = {frozenset(b[0]) for b in bips}
        assert frozenset({"C", "D"}) in all_sides
        assert frozenset({"E", "F"}) in all_sides

    def test_extract_bipartitions_uses_direct_standard_tree_traversal(
        self, monkeypatch
    ):
        tree = _make_tree("((A,B),((C,D),(E,F)));")
        all_taxa = frozenset({"A", "B", "C", "D", "E", "F"})

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard helper should not use find_clades")

        def fail_get_nonterminals(*_args, **_kwargs):
            raise AssertionError("standard helper should not use get_nonterminals")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)
        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_get_nonterminals)

        bips = QuartetNetwork._extract_bipartitions(tree, all_taxa)

        all_sides = {frozenset(bip[0]) for bip in bips}
        assert frozenset({"A", "B"}) in all_sides
        assert frozenset({"C", "D"}) in all_sides
        assert frozenset({"E", "F"}) in all_sides


class TestDetermineQuartetTopology:
    def test_topology_0_ab_cd(self):
        tree = _make_tree("((A,B),(C,D));")
        all_taxa = frozenset({"A", "B", "C", "D"})
        bips = QuartetNetwork._extract_bipartitions(tree, all_taxa)
        result = QuartetNetwork._determine_quartet_topology(("A", "B", "C", "D"), bips)
        assert result == 0  # ab|cd

    def test_topology_1_ac_bd(self):
        tree = _make_tree("((A,C),(B,D));")
        all_taxa = frozenset({"A", "B", "C", "D"})
        bips = QuartetNetwork._extract_bipartitions(tree, all_taxa)
        result = QuartetNetwork._determine_quartet_topology(("A", "B", "C", "D"), bips)
        assert result == 1  # ac|bd

    def test_topology_2_ad_bc(self):
        tree = _make_tree("((A,D),(B,C));")
        all_taxa = frozenset({"A", "B", "C", "D"})
        bips = QuartetNetwork._extract_bipartitions(tree, all_taxa)
        result = QuartetNetwork._determine_quartet_topology(("A", "B", "C", "D"), bips)
        assert result == 2  # ad|bc

    def test_star_topology_returns_none(self):
        tree = _make_tree("(A,B,C,D);")
        all_taxa = frozenset({"A", "B", "C", "D"})
        bips = QuartetNetwork._extract_bipartitions(tree, all_taxa)
        result = QuartetNetwork._determine_quartet_topology(("A", "B", "C", "D"), bips)
        assert result is None


class TestStarTest:
    def test_uniform_counts_high_p(self):
        # Near-uniform counts => high p-value (fail to reject star)
        p = QuartetNetwork._star_test([35, 33, 32])
        assert p > 0.05

    def test_skewed_counts_low_p(self):
        # Highly skewed => low p-value (reject star)
        p = QuartetNetwork._star_test([70, 15, 15])
        assert p < 0.05

    def test_empty_counts(self):
        p = QuartetNetwork._star_test([0, 0, 0])
        assert p == 1.0

    def test_matches_r_reference(self):
        # R: chisq.test(c(8,0,2))$p.value = 0.005517
        p = QuartetNetwork._star_test([8, 0, 2])
        assert abs(p - 0.005517) < 0.001


class TestTreeTest:
    def test_perfect_tree_p_one(self):
        # All counts in one topology => perfect tree, p=1
        p = QuartetNetwork._tree_test([10, 0, 0])
        assert p == 1.0

    def test_equal_minors_high_p(self):
        # Equal minor topologies => consistent with T3 model
        p = QuartetNetwork._tree_test([70, 15, 15])
        assert p > 0.05

    def test_asymmetric_minors_low_p(self):
        # Asymmetric minors => reject T3 model (hybrid signal)
        # R: quartetTreeTest(c(45,35,20), model="T3", lambda=0, method="conservative")
        #    p = 0.0418
        p = QuartetNetwork._tree_test([45, 35, 20])
        assert p < 0.05
        assert abs(p - 0.0418) < 0.005

    def test_empty_counts(self):
        p = QuartetNetwork._tree_test([0, 0, 0])
        assert p == 1.0


class TestClassifyQuartet:
    def test_tree_like(self):
        # R reference: p_star=7.3e-14 (<<0.95), p_T3=1.0 (>0.05) => tree
        counts = [70, 15, 15]
        result = QuartetNetwork._classify_quartet(counts, alpha=0.05, beta=0.95)
        assert result["classification"] == "tree"
        assert result["p_star"] < 0.95  # star rejected
        assert result["p_tree"] > 0.05  # tree not rejected

    def test_hybrid(self):
        # R reference: p_star=0.0087 (<<0.95), p_T3=0.042 (<0.05) => hybrid
        counts = [45, 35, 20]
        result = QuartetNetwork._classify_quartet(counts, alpha=0.05, beta=0.95)
        assert result["classification"] == "hybrid"
        assert result["p_star"] < 0.95
        assert result["p_tree"] < 0.05

    def test_classification_does_not_import_scipy_stats(self, monkeypatch):
        real_import = __import__

        def fail_scipy_stats_import(name, *args, **kwargs):
            if name == "scipy.stats" or name.startswith("scipy.stats."):
                raise AssertionError("quartet classification should not import scipy.stats")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr("builtins.__import__", fail_scipy_stats_import)

        result = QuartetNetwork._classify_quartet([45, 35, 20], alpha=0.05, beta=0.95)

        assert result["classification"] == "hybrid"
        assert result["p_star"] == pytest.approx(0.008651695203120634)
        assert result["p_tree"] == pytest.approx(0.04180220628706992)

    def test_unresolved_with_high_beta(self):
        # R reference: p_star=0.932. With beta=0.90 => p_star > beta => star
        counts = [35, 33, 32]
        result = QuartetNetwork._classify_quartet(counts, alpha=0.05, beta=0.90)
        assert result["classification"] == "unresolved"
        assert result["p_star"] > 0.90

    def test_tree_with_default_beta(self):
        # R reference: p_star=0.932 < beta=0.95 => star rejected => check T3
        # p_T3 > 0.05 => tree
        counts = [35, 33, 32]
        result = QuartetNetwork._classify_quartet(counts, alpha=0.05, beta=0.95)
        assert result["classification"] == "tree"

    def test_empty_counts(self):
        counts = [0, 0, 0]
        result = QuartetNetwork._classify_quartet(counts, alpha=0.05, beta=0.95)
        assert result["classification"] == "unresolved"
        assert result["p_star"] == 1.0
        assert result["p_tree"] == 1.0

    def test_matches_nanuq_sample_data(self):
        # From R NANUQ on gene_trees_for_network.nwk:
        # counts=(8,0,2): p_star=0.0055, p_T3=0.214 => tree
        result = QuartetNetwork._classify_quartet([8, 0, 2], alpha=0.05, beta=0.95)
        assert result["classification"] == "tree"

        # counts=(10,0,0): p_star=4.5e-5, p_T3=1.0 => tree
        result = QuartetNetwork._classify_quartet([10, 0, 0], alpha=0.05, beta=0.95)
        assert result["classification"] == "tree"

        # counts=(4,1,5): p_star=0.273, p_T3=0.183 => tree
        result = QuartetNetwork._classify_quartet([4, 1, 5], alpha=0.05, beta=0.95)
        assert result["classification"] == "tree"


class TestNanuqDistance:
    def test_topology_index_helpers_preserve_stable_ties(self):
        assert QuartetNetwork._dominant_topology_index([5, 5, 4]) == 0
        assert QuartetNetwork._dominant_topology_index([4, 5, 5]) == 1
        assert QuartetNetwork._dominant_topology_index([4, 3, 5]) == 2

        assert QuartetNetwork._top_two_topology_indices([5, 5, 4]) == (0, 1)
        assert QuartetNetwork._top_two_topology_indices([5, 4, 5]) == (0, 2)
        assert QuartetNetwork._top_two_topology_indices([4, 5, 5]) == (1, 2)
        assert QuartetNetwork._top_two_topology_indices([3, 5, 4]) == (1, 2)

    def test_tree_quartet_adds_cross_pairs(self):
        """Tree-like quartet ab|cd should add 1 to cross-split pairs only."""
        all_taxa = frozenset({"A", "B", "C", "D"})
        quartet_results = {
            ("A", "B", "C", "D"): {
                "classification": "tree",
                "counts": [5, 0, 0],
                "p_star": 0.001,
                "p_tree": 1.0,
            }
        }
        taxa_list, dist, taxa_idx = QuartetNetwork._compute_nanuq_distance(
            all_taxa, quartet_results
        )
        ia, ib = taxa_idx["A"], taxa_idx["B"]
        ic, id_ = taxa_idx["C"], taxa_idx["D"]
        # Cross-split pairs (A,C), (A,D), (B,C), (B,D) get +1
        assert dist[ia][ic] == 1.0
        assert dist[ia][id_] == 1.0
        assert dist[ib][ic] == 1.0
        assert dist[ib][id_] == 1.0
        # Same-side pairs (A,B), (C,D) stay 0
        assert dist[ia][ib] == 0.0
        assert dist[ic][id_] == 0.0

    def test_star_quartet_adds_all_pairs(self):
        """Unresolved quartet should add 1 to all 6 pairwise distances."""
        all_taxa = frozenset({"A", "B", "C", "D"})
        quartet_results = {
            ("A", "B", "C", "D"): {
                "classification": "unresolved",
                "counts": [3, 3, 4],
                "p_star": 0.99,
                "p_tree": 1.0,
            }
        }
        _, dist, idx = QuartetNetwork._compute_nanuq_distance(all_taxa, quartet_results)
        for t1 in ["A", "B", "C", "D"]:
            for t2 in ["A", "B", "C", "D"]:
                if t1 != t2:
                    assert dist[idx[t1]][idx[t2]] == 1.0

    def test_hybrid_quartet_adds_half(self):
        """Hybrid quartet should add 0.5 for each of the two elevated topologies."""
        all_taxa = frozenset({"A", "B", "C", "D"})
        # Topology 0 (ab|cd) dominant, topology 1 (ac|bd) elevated
        quartet_results = {
            ("A", "B", "C", "D"): {
                "classification": "hybrid",
                "counts": [45, 35, 20],
                "p_star": 0.009,
                "p_tree": 0.04,
            }
        }
        _, dist, idx = QuartetNetwork._compute_nanuq_distance(all_taxa, quartet_results)
        # Topo 0 cross pairs: (A,C),(A,D),(B,C),(B,D) get +0.5
        # Topo 1 cross pairs: (A,B),(A,D),(C,B),(C,D) get +0.5
        # (A,D) appears in both → gets 1.0
        assert dist[idx["A"]][idx["D"]] == 1.0
        # (A,C) only in topo 0 → 0.5
        assert dist[idx["A"]][idx["C"]] == 0.5


class TestNeighborJoiningOrder:
    def test_four_taxa_tree_ordering(self):
        """NJ on a tree distance should group cherries adjacently."""
        # Tree ((A,B),(C,D)) with d(A,B)=2, d(C,D)=2, cross=4
        taxa = ["A", "B", "C", "D"]
        dist = [
            [0, 2, 4, 4],
            [2, 0, 4, 4],
            [4, 4, 0, 2],
            [4, 4, 2, 0],
        ]
        ordering = QuartetNetwork._neighbor_joining_order(taxa, dist)
        assert set(ordering) == {"A", "B", "C", "D"}
        # A and B should be adjacent in circular ordering
        ia = ordering.index("A")
        ib = ordering.index("B")
        assert abs(ia - ib) == 1 or abs(ia - ib) == 3  # adjacent circularly

    def test_two_taxa(self):
        ordering = QuartetNetwork._neighbor_joining_order(["X", "Y"], [[0, 1], [1, 0]])
        assert ordering == ["X", "Y"]


class TestCircularSplitWeights:
    def test_tree_distance_one_positive_split(self):
        """A pure tree distance should produce one positive internal split."""
        ordering = ["A", "B", "C", "D"]
        taxa_idx = {"A": 0, "B": 1, "C": 2, "D": 3}
        # Tree ((A,B),(C,D)): d(A,B)=2, d(C,D)=2, cross=4
        dist = [
            [0, 2, 4, 4],
            [2, 0, 4, 4],
            [4, 4, 0, 2],
            [4, 4, 2, 0],
        ]
        splits = QuartetNetwork._compute_circular_split_weights(
            ordering, dist, taxa_idx
        )
        # Should have exactly one positive split: {A,B}|{C,D}
        assert len(splits) == 1
        split, weight = splits[0]
        assert split == frozenset({"A", "B"}) or split == frozenset({"C", "D"})
        assert weight > 0

    def test_equal_distances_no_splits(self):
        """Uniform distances should produce no positive splits."""
        ordering = ["A", "B", "C", "D"]
        taxa_idx = {"A": 0, "B": 1, "C": 2, "D": 3}
        dist = [
            [0, 1, 1, 1],
            [1, 0, 1, 1],
            [1, 1, 0, 1],
            [1, 1, 1, 0],
        ]
        splits = QuartetNetwork._compute_circular_split_weights(
            ordering, dist, taxa_idx
        )
        assert len(splits) == 0

    def test_equal_size_circular_split_tiebreak_avoids_sorting(self, monkeypatch):
        ordering = ["D", "C", "B", "A"]
        taxa_idx = {taxon: idx for idx, taxon in enumerate(ordering)}
        dist = [
            [0, 2, 4, 4],
            [2, 0, 4, 4],
            [4, 4, 0, 2],
            [4, 4, 2, 0],
        ]

        def fail_sorted(*_args, **_kwargs):
            raise AssertionError("equal-size split tiebreak should use min")

        monkeypatch.setattr(
            quartet_network_module,
            "sorted",
            fail_sorted,
            raising=False,
        )

        splits = QuartetNetwork._compute_circular_split_weights(
            ordering,
            dist,
            taxa_idx,
        )

        assert splits == [(frozenset({"A", "B"}), 2.0)]


class TestBuildSplitsGraph:
    class CountingOrdering(list):
        def __init__(self, values):
            super().__init__(values)
            self.index_count = 0

        def __getitem__(self, index):
            self.index_count += 1
            return super().__getitem__(index)

    def test_is_circular_split_handles_wraparound_boundaries(self):
        ordering = ["A", "B", "C", "D", "E", "F"]

        assert QuartetNetwork._is_circular_split(
            frozenset(["E", "F", "A"]),
            ordering,
        )
        assert QuartetNetwork._is_circular_split(
            frozenset(["B", "C", "D"]),
            ordering,
        )
        assert not QuartetNetwork._is_circular_split(
            frozenset(["A", "C", "E"]),
            ordering,
        )
        assert not QuartetNetwork._is_circular_split(frozenset(), [])

    def test_is_circular_split_stops_after_third_boundary(self):
        ordering = self.CountingOrdering(["A", "B", "C", "D", "E", "F"])

        assert not QuartetNetwork._is_circular_split(
            frozenset(["A", "C", "E"]),
            ordering,
        )
        assert ordering.index_count < len(ordering)

    def test_compute_split_directions_matches_legacy_wraparound_boundaries(self):
        ordering = ["A", "B", "C", "D", "E", "F"]
        circular_splits = [
            (frozenset(["E", "F", "A"]), 3, 0.3),
            (frozenset(["B", "C", "D"]), 2, 0.2),
            (frozenset(["A", "C", "E"]), 1, 0.1),
        ]

        observed = QuartetNetwork._compute_split_directions(
            ordering,
            circular_splits,
        )
        expected = _legacy_compute_split_directions(ordering, circular_splits)

        assert observed.keys() == expected.keys()
        for split, direction in expected.items():
            assert observed[split] == pytest.approx(direction)

    def test_compute_split_directions_uses_supplied_gap_positions(self, monkeypatch):
        ordering = ["A", "B", "C", "D", "E", "F"]
        circular_splits = [
            (frozenset(["E", "F", "A"]), 3, 0.3),
            (frozenset(["B", "C", "D"]), 2, 0.2),
        ]
        gap_positions_by_split = {
            split: QuartetNetwork._circular_gap_positions(split, ordering)
            for split, _count, _freq in circular_splits
        }

        def fail_gap_scan(*args, **kwargs):
            raise AssertionError("cached gap positions should be reused")

        monkeypatch.setattr(QuartetNetwork, "_circular_gap_positions", fail_gap_scan)

        observed = QuartetNetwork._compute_split_directions(
            ordering,
            circular_splits,
            gap_positions_by_split,
        )
        expected = _legacy_compute_split_directions(ordering, circular_splits)

        assert observed.keys() == expected.keys()
        for split, direction in expected.items():
            assert observed[split] == pytest.approx(direction)

    @staticmethod
    def _legacy_build_splits_graph(circular_splits, all_taxa):
        splits_list = [s[0] for s in circular_splits]
        weights = [s[2] for s in circular_splits]
        n_splits = len(splits_list)
        if n_splits == 0:
            return {}, set(), [], [], []
        taxon_signs = {}
        for taxon in all_taxa:
            signs = tuple(1 if taxon in sp else -1 for sp in splits_list)
            taxon_signs[taxon] = signs
        forbidden = {}
        for i in range(n_splits):
            for j in range(i + 1, n_splits):
                si_pos = splits_list[i]
                si_neg = all_taxa - si_pos
                sj_pos = splits_list[j]
                sj_neg = all_taxa - sj_pos
                fb = set()
                if not (si_pos & sj_pos):
                    fb.add((1, 1))
                if not (si_pos & sj_neg):
                    fb.add((1, -1))
                if not (si_neg & sj_pos):
                    fb.add((-1, 1))
                if not (si_neg & sj_neg):
                    fb.add((-1, -1))
                if fb:
                    forbidden[(i, j)] = fb
        valid_nodes = set()
        for combo in range(2 ** n_splits):
            signs = tuple(1 if (combo >> j) & 1 else -1 for j in range(n_splits))
            valid = True
            for (i, j), fb_pairs in forbidden.items():
                if (signs[i], signs[j]) in fb_pairs:
                    valid = False
                    break
            if valid:
                valid_nodes.add(signs)
        for signs in taxon_signs.values():
            valid_nodes.add(signs)
        edges = []
        valid_list = list(valid_nodes)
        for i in range(len(valid_list)):
            for j in range(i + 1, len(valid_list)):
                s1 = valid_list[i]
                s2 = valid_list[j]
                diff = [k for k in range(n_splits) if s1[k] != s2[k]]
                if len(diff) == 1:
                    edges.append((s1, s2, diff[0]))
        return taxon_signs, valid_nodes, edges, splits_list, weights

    @staticmethod
    def _canonical_edges(edges):
        return {
            (min(s1, s2), max(s1, s2), split_idx)
            for s1, s2, split_idx in edges
        }

    def test_build_splits_graph_matches_exhaustive_reference(self):
        ordering = [f"T{i}" for i in range(10)]
        all_taxa = frozenset(ordering)
        circular_splits = [
            (frozenset({"T0", "T1"}), 7, 0.7),
            (frozenset({"T2", "T3", "T4"}), 6, 0.6),
            (frozenset({"T5", "T6"}), 5, 0.5),
            (frozenset({"T7", "T8", "T9"}), 4, 0.4),
            (frozenset({"T1", "T2", "T3"}), 3, 0.3),
            (frozenset({"T4", "T5"}), 2, 0.2),
        ]

        observed = QuartetNetwork._build_splits_graph(circular_splits, all_taxa)
        expected = self._legacy_build_splits_graph(circular_splits, all_taxa)

        assert observed[0] == expected[0]
        assert observed[1] == expected[1]
        assert self._canonical_edges(observed[2]) == self._canonical_edges(expected[2])
        assert observed[3] == expected[3]
        assert observed[4] == expected[4]


class TestNetworkPlot:
    def test_draw_quartet_network_batches_edges(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection
        from phykit.helpers.plot_config import PlotConfig

        service = QuartetNetwork.__new__(QuartetNetwork)
        service.plot_config = PlotConfig(show_title=False)
        ordering = ["A", "B", "C", "D"]
        all_taxa = frozenset(ordering)
        taxa_idx = {taxon: idx for idx, taxon in enumerate(ordering)}

        monkeypatch.setattr(
            service,
            "_compute_nanuq_distance",
            lambda all_taxa, quartet_results: (ordering, [], taxa_idx),
        )
        monkeypatch.setattr(
            service,
            "_neighbor_joining_order",
            lambda taxa_list, dist_matrix: ordering,
        )
        monkeypatch.setattr(
            service,
            "_compute_circular_split_weights",
            lambda ordering_arg, dist_matrix, taxa_idx_arg: [
                (frozenset({"A", "B"}), 1.0),
                (frozenset({"B", "C"}), 0.5),
            ],
        )

        line_collections = []
        original_add_collection = matplotlib.axes.Axes.add_collection

        def fail_plot(*args, **kwargs):
            raise AssertionError("quartet network edges should use LineCollection")

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(matplotlib.axes.Axes, "add_collection", capture_collection)

        output_path = tmp_path / "quartet_net.png"
        service._draw_quartet_network(all_taxa, {}, str(output_path))

        assert len(line_collections) >= 2
        assert output_path.exists()


class TestComputeQuartetCFs:
    def test_bitmask_topology_matches_public_bipartition_topology(self):
        tree = _make_tree("((A,B),((C,D),(E,F)));")
        all_taxa = frozenset({"A", "B", "C", "D", "E", "F"})
        taxa_sorted = sorted(all_taxa)
        taxa_index = {taxon: idx for idx, taxon in enumerate(taxa_sorted)}

        bips = QuartetNetwork._extract_bipartitions(tree, all_taxa)
        masks = QuartetNetwork._extract_bipartition_masks(
            tree,
            all_taxa,
            taxa_index,
        )

        for quartet in itertools.combinations(taxa_sorted, 4):
            expected = QuartetNetwork._determine_quartet_topology(quartet, bips)
            quartet_bits = tuple(1 << taxa_index[taxon] for taxon in quartet)
            observed = QuartetNetwork._determine_quartet_topology_from_masks(
                quartet_bits,
                masks,
            )
            assert observed == expected

    def test_direct_postorder_matches_standard_tree_traversal(self, monkeypatch):
        tree = _make_tree("((A,B),((C,D),(E,F)));")
        expected = list(tree.find_clades(order="postorder"))

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("direct postorder helper should not call find_clades")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        assert QuartetNetwork._postorder_clades_direct(tree) == expected

    def test_extract_bipartition_masks_uses_single_postorder_pass(self, monkeypatch):
        tree = _make_tree("((A,B),((C,D),(E,F)));")
        all_taxa = frozenset({"A", "B", "C", "D", "E", "F"})
        taxa_index = {taxon: idx for idx, taxon in enumerate(sorted(all_taxa))}

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("mask extraction should use direct postorder")

        def fail_get_nonterminals(*args, **kwargs):
            raise AssertionError("mask extraction should not call get_nonterminals")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)
        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_get_nonterminals)

        masks = QuartetNetwork._extract_bipartition_masks(
            tree,
            all_taxa,
            taxa_index,
        )

        assert set(masks) == {
            (1 << taxa_index["A"]) | (1 << taxa_index["B"]),
            (1 << taxa_index["C"]) | (1 << taxa_index["D"]),
            (1 << taxa_index["E"]) | (1 << taxa_index["F"]),
            (
                (1 << taxa_index["C"])
                | (1 << taxa_index["D"])
                | (1 << taxa_index["E"])
                | (1 << taxa_index["F"])
            ),
        }

    def test_extract_bipartition_masks_handles_mixed_child_counts(self):
        tree = _make_tree("(((A):1,((B,C):1,D):1):1,E:1,F:1);")
        all_taxa = frozenset({"A", "B", "C", "D", "E", "F"})
        taxa_index = {taxon: idx for idx, taxon in enumerate(sorted(all_taxa))}

        masks = QuartetNetwork._extract_bipartition_masks(
            tree,
            all_taxa,
            taxa_index,
        )

        assert set(masks) == {
            (1 << taxa_index["B"]) | (1 << taxa_index["C"]),
            (
                (1 << taxa_index["B"])
                | (1 << taxa_index["C"])
                | (1 << taxa_index["D"])
            ),
            (
                (1 << taxa_index["A"])
                | (1 << taxa_index["B"])
                | (1 << taxa_index["C"])
                | (1 << taxa_index["D"])
            ),
        }

    def test_four_taxa_three_trees(self):
        trees = [
            _make_tree("((A,B),(C,D));"),
            _make_tree("((A,B),(C,D));"),
            _make_tree("((A,C),(B,D));"),
        ]
        all_taxa = frozenset({"A", "B", "C", "D"})
        results = QuartetNetwork._compute_quartet_cfs(trees, all_taxa)
        # Only one quartet: (A,B,C,D)
        assert len(results) == 1
        quartet = ("A", "B", "C", "D")
        counts = results[quartet]
        assert counts[0] == 2  # ab|cd appears in 2 trees
        assert counts[1] == 1  # ac|bd appears in 1 tree
        assert counts[2] == 0  # ad|bc appears in 0 trees

    def test_precomputed_quartet_masks_match_existing_helper(self):
        quartet_bits = (1, 2, 4, 8)
        quartet_data = (
            1 | 2 | 4 | 8,
            1 | 2,
            1 | 4,
            1 | 8,
            4 | 8,
            2 | 8,
            2 | 4,
        )
        split_masks = [1 | 4, 2 | 8]

        expected = QuartetNetwork._determine_quartet_topology_from_masks(
            quartet_bits,
            split_masks,
        )
        observed = QuartetNetwork._determine_quartet_topology_from_precomputed_masks(
            quartet_data,
            split_masks,
        )

        assert observed == expected == 1

    def test_compute_quartet_cfs_uses_precomputed_masks(self, monkeypatch):
        trees = [
            _make_tree("((A,B),(C,D));"),
            _make_tree("((A,C),(B,D));"),
        ]
        all_taxa = frozenset({"A", "B", "C", "D"})

        def fail_old_helper(*_args, **_kwargs):
            raise AssertionError("precomputed quartet masks should be used")

        monkeypatch.setattr(
            QuartetNetwork,
            "_determine_quartet_topology_from_masks",
            fail_old_helper,
        )

        results = QuartetNetwork._compute_quartet_cfs(trees, all_taxa)

        assert results[("A", "B", "C", "D")] == [1, 1, 0]


class TestPruneToTaxa:
    def test_normalize_taxa_uses_fast_terminal_names_for_identical_parsed_trees(
        self, monkeypatch
    ):
        svc = QuartetNetwork(
            Namespace(
                trees="unused",
                alpha=0.05,
                beta=0.95,
                missing_taxa="error",
                plot_output=None,
                json=False,
                fig_width=None,
                fig_height=None,
                dpi=300,
                no_title=False,
                title=None,
                legend_position=None,
                ylabel_fontsize=None,
                xlabel_fontsize=None,
                title_fontsize=None,
                axis_fontsize=None,
                colors=None,
                ladderize=False,
                cladogram=False,
                circular=False,
                color_file=None,
            )
        )
        trees = [
            _make_tree("((A,B),(C,D));"),
            _make_tree("((A,C),(B,D));"),
        ]

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("parsed trees should use fast terminal names")

        def fail_prune(*_args, **_kwargs):
            raise AssertionError("identical taxon sets should not be pruned")

        for tree in trees:
            monkeypatch.setattr(tree, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(QuartetNetwork, "_prune_to_taxa", fail_prune)

        normalized, pruned, taxa = svc._normalize_taxa(trees)

        assert normalized is trees
        assert pruned is False
        assert taxa == {"A", "B", "C", "D"}

    def test_all_tip_sets_identical_does_not_slice_rows(self):
        class NoSliceList(list):
            def __getitem__(self, key):
                if isinstance(key, slice):
                    raise AssertionError("tip-set scan should not slice")
                return super().__getitem__(key)

        tip_sets = NoSliceList([
            {"A", "B", "C", "D"},
            {"A", "B", "C", "D"},
            {"A", "B", "C", "D"},
        ])

        assert quartet_network_module._all_tip_sets_identical(tip_sets) is True

    def test_uses_direct_terminal_collection(self, monkeypatch):
        tree = _make_tree("((A,B),(C,D));")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("standard tree should use direct terminal traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        QuartetNetwork._prune_to_taxa(tree, {"A", "C"})

        assert QuartetNetwork._tips(tree) == {"A", "C"}

    def test_uses_batch_prune_for_standard_tree(self, monkeypatch):
        tree = _make_tree("((A:1,B:2):3,(C:4,D:5):6);")

        def fail_prune(*_args, **_kwargs):
            raise AssertionError("per-target prune should not be used")

        monkeypatch.setattr(TreeMixin, "prune", fail_prune)

        QuartetNetwork._prune_to_taxa(tree, {"B", "D"})

        assert [tip.name for tip in tree.get_terminals()] == ["B", "D"]
        assert [tip.branch_length for tip in tree.get_terminals()] == [5.0, 11.0]


class TestQuartetNetworkRun:
    def test_parse_trees_skips_comments_and_blank_lines(self, tmp_path):
        tree_file = tmp_path / "trees.nwk"
        _write(
            tree_file,
            "# ignored\n\n  ((A,B),(C,D));  \n# also ignored\n((A,C),(B,D));\n",
        )
        svc = QuartetNetwork(
            Namespace(
                trees=str(tree_file),
                alpha=0.05,
                beta=0.95,
                missing_taxa="error",
                plot_output=None,
                json=False,
            )
        )

        trees = svc._parse_trees_from_source(str(tree_file))

        assert len(trees) == 2
        assert [sorted(QuartetNetwork._tips(tree)) for tree in trees] == [
            ["A", "B", "C", "D"],
            ["A", "B", "C", "D"],
        ]

    def test_parse_tree_path_list_avoids_per_row_path_objects(self, tmp_path, mocker):
        tree_a = tmp_path / "a.nwk"
        tree_b = tmp_path / "b.nwk"
        _write(tree_a, "((A,B),(C,D));\n")
        _write(tree_b, "((A,C),(B,D));\n")
        tree_list = tmp_path / "tree-list.txt"
        _write(tree_list, "a.nwk\nb.nwk\n")
        path_calls = 0
        original_path = __import__(
            "phykit.services.tree.quartet_network",
            fromlist=["Path"],
        ).Path

        def counting_path(*args, **kwargs):
            nonlocal path_calls
            path_calls += 1
            return original_path(*args, **kwargs)

        mocker.patch(
            "phykit.services.tree.quartet_network.Path",
            side_effect=counting_path,
        )
        svc = QuartetNetwork(
            Namespace(
                trees=str(tree_list),
                alpha=0.05,
                beta=0.95,
                missing_taxa="error",
                plot_output=None,
                json=False,
            )
        )

        trees = svc._parse_trees_from_source(str(tree_list))

        assert len(trees) == 2
        assert path_calls == 1

    def test_text_output(self, tmp_path, capsys):
        tree_file = tmp_path / "trees.nwk"
        _write(
            tree_file,
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n"
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,B),(C,D));\n",
        )

        svc = QuartetNetwork(
            Namespace(
                trees=str(tree_file),
                alpha=0.05,
                beta=0.95,
                missing_taxa="error",
                plot_output=None,
                json=False,
            )
        )
        svc.run()
        captured = capsys.readouterr()
        assert "Number of input trees: 6" in captured.out
        assert "Number of taxa: 4" in captured.out
        assert "Total quartets: 1" in captured.out

    def test_text_output_batches_quartet_report(self, tmp_path, mocker):
        tree_file = tmp_path / "trees.nwk"
        _write(tree_file, "((A,B),(C,D));\n")
        svc = QuartetNetwork(
            Namespace(
                trees=str(tree_file),
                alpha=0.05,
                beta=0.95,
                missing_taxa="shared",
                plot_output=None,
                json=False,
            )
        )
        mocker.patch.object(svc, "_parse_trees_from_source", return_value=["t1", "t2", "t3"])
        mocker.patch.object(
            svc,
            "_normalize_taxa",
            return_value=(["t1", "t2", "t3"], True, {"A", "B", "C", "D"}),
        )
        mocker.patch.object(
            svc,
            "_compute_quartet_cfs",
            return_value={
                ("A", "B", "C", "D"): [3, 0, 0],
                ("A", "C", "B", "D"): [1, 2, 0],
            },
        )
        mocker.patch.object(
            svc,
            "_classify_quartet",
            side_effect=[
                {"classification": "tree", "p_star": 0.012345, "p_tree": 0.98765},
                {"classification": "hybrid", "p_star": 0.11111, "p_tree": 0.22222},
            ],
        )
        printed = mocker.patch("builtins.print")

        svc.run()

        printed.assert_called_once_with(
            "Number of input trees: 3\n"
            "Number of taxa: 4\n"
            "Significance level (alpha): 0.05\n"
            "Star tree threshold (beta): 0.95\n"
            "Pruned to shared taxa: yes\n"
            "Total quartets: 2\n"
            "Tree-like: 1 (50.0%)\n"
            "Hybrid: 1 (50.0%)\n"
            "Unresolved: 0 (0.0%)\n"
            "---\n"
            "{A, B} | {C, D}\t1.0000\t0.0000\t0.0000\ttree\t"
            "p_star=0.0123\tp_tree=0.9877\n"
            "{A, B} | {C, D}\t0.3333\t0.6667\t0.0000\thybrid\t"
            "p_star=0.1111\tp_tree=0.2222"
        )

    def test_json_output(self, tmp_path, monkeypatch):
        tree_file = tmp_path / "trees.nwk"
        _write(
            tree_file,
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n"
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,B),(C,D));\n",
        )

        captured = {}
        svc = QuartetNetwork(
            Namespace(
                trees=str(tree_file),
                alpha=0.05,
                beta=0.95,
                missing_taxa="error",
                plot_output=None,
                json=True,
            )
        )
        monkeypatch.setattr(
            "phykit.services.tree.quartet_network.print_json",
            lambda payload: captured.setdefault("payload", payload),
        )
        svc.run()

        payload = captured["payload"]
        assert payload["input_tree_count"] == 6
        assert payload["total_quartets"] == 1
        assert len(payload["quartets"]) == 1
        assert payload["quartets"][0]["classification"] in ("tree", "hybrid", "unresolved")
        assert "p_star" in payload["quartets"][0]
        assert "p_tree" in payload["quartets"][0]
        assert payload["beta"] == 0.95

    def test_fewer_than_four_taxa_raises(self, tmp_path):
        tree_file = tmp_path / "trees.nwk"
        _write(tree_file, "((A,B),C);\n((A,B),C);\n")

        svc = QuartetNetwork(
            Namespace(
                trees=str(tree_file),
                alpha=0.05,
                beta=0.95,
                missing_taxa="error",
                plot_output=None,
                json=False,
            )
        )
        with pytest.raises(PhykitUserError):
            svc.run()

    def test_missing_taxa_error_mode_raises(self, tmp_path):
        tree_file = tmp_path / "trees.nwk"
        _write(tree_file, "((A,B),(C,D));\n((A,B),(C,E));\n")

        svc = QuartetNetwork(
            Namespace(
                trees=str(tree_file),
                alpha=0.05,
                beta=0.95,
                missing_taxa="error",
                plot_output=None,
                json=False,
            )
        )
        with pytest.raises(PhykitUserError):
            svc.run()

    def test_missing_taxa_shared_mode(self, tmp_path, capsys):
        tree_file = tmp_path / "trees.nwk"
        _write(
            tree_file,
            "((A,B),(C,(D,E)));\n((A,B),(C,(D,F)));\n"
            "((A,B),(C,(D,E)));\n((A,B),(C,(D,F)));\n"
            "((A,B),(C,(D,E)));\n((A,B),(C,(D,F)));\n",
        )

        svc = QuartetNetwork(
            Namespace(
                trees=str(tree_file),
                alpha=0.05,
                beta=0.95,
                missing_taxa="shared",
                plot_output=None,
                json=False,
            )
        )
        svc.run()
        captured = capsys.readouterr()
        assert "Pruned to shared taxa: yes" in captured.out

    def test_plot_output_creates_file(self, tmp_path):
        tree_file = tmp_path / "trees.nwk"
        _write(
            tree_file,
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n"
            "((A,B),(C,D));\n((A,B),(C,D));\n((A,B),(C,D));\n",
        )
        plot_path = str(tmp_path / "quartet_net.png")

        svc = QuartetNetwork(
            Namespace(
                trees=str(tree_file),
                alpha=0.05,
                beta=0.95,
                missing_taxa="error",
                plot_output=plot_path,
                json=False,
            )
        )
        svc.run()
        assert (tmp_path / "quartet_net.png").exists()

    def test_sample_data_all_tree_like(self, tmp_path, capsys):
        """Validate against R's NANUQ: all 15 quartets in
        gene_trees_for_network.nwk should be classified as 'tree'."""
        svc = QuartetNetwork(
            Namespace(
                trees="tests/sample_files/gene_trees_for_network.nwk",
                alpha=0.05,
                beta=0.95,
                missing_taxa="error",
                plot_output=None,
                json=False,
            )
        )
        svc.run()
        captured = capsys.readouterr()
        assert "Tree-like: 15 (100.0%)" in captured.out
        assert "Hybrid: 0 (0.0%)" in captured.out
        assert "Unresolved: 0 (0.0%)" in captured.out


class TestPolytomyHandling:
    """Polytomous nodes (collapsed branches) should be skipped in bipartition
    extraction, so quartets spanning a polytomy are treated as unresolved."""

    def test_star_tree_no_bipartitions(self):
        """A 4-way star tree (A,B,C,D) has no informative bipartitions."""
        tree = _make_tree("(A:1,B:1,C:1,D:1);")
        all_taxa = frozenset(["A", "B", "C", "D"])
        bips = QuartetNetwork._extract_bipartitions(tree, all_taxa)
        assert len(bips) == 0

    def test_trifurcating_root_allowed(self):
        """A trifurcating root (standard unrooted) should NOT be skipped."""
        tree = _make_tree("((A:1,B:1):1,C:1,D:1);")
        all_taxa = frozenset(["A", "B", "C", "D"])
        bips = QuartetNetwork._extract_bipartitions(tree, all_taxa)
        # The resolved subclade (A,B) produces one bipartition
        assert len(bips) >= 1

    def test_internal_polytomy_skipped(self):
        """An internal 4-way polytomy should produce no bipartitions."""
        tree = _make_tree("((A:1,B:1,C:1,D:1):1,E:1);")
        all_taxa = frozenset(["A", "B", "C", "D", "E"])
        bips = QuartetNetwork._extract_bipartitions(tree, all_taxa)
        # The polytomous node (A,B,C,D) is skipped; only the root edge remains
        # but it's trivial (4 vs 1), so no bipartitions
        assert len(bips) == 0

    def test_partial_polytomy_preserves_resolved(self):
        """Resolved subclades of a polytomy still produce bipartitions."""
        tree = _make_tree("((A:1,B:1):1,C:1,D:1,E:1);")
        all_taxa = frozenset(["A", "B", "C", "D", "E"])
        bips = QuartetNetwork._extract_bipartitions(tree, all_taxa)
        # (A,B) is resolved and produces {A,B}|{C,D,E}
        ab = frozenset(["A", "B"])
        cde = frozenset(["C", "D", "E"])
        assert (ab, cde) in bips or (cde, ab) in bips

    def test_polytomy_quartet_returns_none(self):
        """A quartet spanning an unresolved polytomy returns None."""
        tree = _make_tree("(A:1,B:1,C:1,D:1);")
        all_taxa = frozenset(["A", "B", "C", "D"])
        bips = QuartetNetwork._extract_bipartitions(tree, all_taxa)
        topo = QuartetNetwork._determine_quartet_topology(
            ("A", "B", "C", "D"), bips
        )
        assert topo is None

    def test_polytomy_trees_not_counted(self):
        """Polytomous gene trees should not contribute to quartet counts."""
        resolved = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        polytomy = _make_tree("(A:1,B:1,C:1,D:1);")
        all_taxa = frozenset(["A", "B", "C", "D"])
        cfs = QuartetNetwork._compute_quartet_cfs(
            [resolved, polytomy], all_taxa
        )
        quartet = ("A", "B", "C", "D")
        counts = cfs[quartet]
        # Only the resolved tree contributes (topology 0: ab|cd)
        assert sum(counts) == 1
        assert counts[0] == 1
