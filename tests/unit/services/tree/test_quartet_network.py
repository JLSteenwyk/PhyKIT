from argparse import Namespace
from io import StringIO

import pytest
from Bio import Phylo

from phykit.errors import PhykitUserError
from phykit.services.tree.quartet_network import QuartetNetwork


def _write(path, content):
    path.write_text(content)
    return str(path)


def _make_tree(newick):
    return Phylo.read(StringIO(newick), "newick")


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


class TestComputeQuartetCFs:
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


class TestQuartetNetworkRun:
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
