from argparse import Namespace
from io import StringIO

import pytest
from Bio import Phylo

from phykit.errors import PhykitUserError
from phykit.services.tree.consensus_network import ConsensusNetwork


def _write(path, content):
    path.write_text(content)
    return str(path)


def _make_tree(newick):
    return Phylo.read(StringIO(newick), "newick")


class TestCanonicalSplit:
    def test_returns_smaller_side(self):
        all_taxa = frozenset({"A", "B", "C", "D"})
        result = ConsensusNetwork._canonical_split(frozenset({"A", "B", "C"}), all_taxa)
        assert result == frozenset({"D"})

    def test_returns_smaller_side_when_given_small(self):
        all_taxa = frozenset({"A", "B", "C", "D"})
        result = ConsensusNetwork._canonical_split(frozenset({"A"}), all_taxa)
        assert result == frozenset({"A"})

    def test_equal_size_returns_lexicographically_smaller(self):
        all_taxa = frozenset({"A", "B", "C", "D"})
        result1 = ConsensusNetwork._canonical_split(frozenset({"A", "B"}), all_taxa)
        result2 = ConsensusNetwork._canonical_split(frozenset({"C", "D"}), all_taxa)
        assert result1 == result2
        assert result1 == frozenset({"A", "B"})


class TestExtractSplitsFromTree:
    def test_simple_balanced_tree(self):
        tree = _make_tree("((A,B),(C,D));")
        all_taxa = frozenset({"A", "B", "C", "D"})
        splits = ConsensusNetwork._extract_splits_from_tree(tree, all_taxa)
        # Only non-trivial split: {A,B} (or equivalently {C,D})
        assert splits == {frozenset({"A", "B"})}

    def test_caterpillar_tree(self):
        tree = _make_tree("(A,(B,(C,D)));")
        all_taxa = frozenset({"A", "B", "C", "D"})
        splits = ConsensusNetwork._extract_splits_from_tree(tree, all_taxa)
        # Non-trivial split: {C,D} canonicalizes to {A,B} (equal size, lexicographic)
        # {B,C,D} is all-but-one so trivial
        assert splits == {frozenset({"A", "B"})}

    def test_six_taxa_tree(self):
        tree = _make_tree("((A,B),((C,D),(E,F)));")
        all_taxa = frozenset({"A", "B", "C", "D", "E", "F"})
        splits = ConsensusNetwork._extract_splits_from_tree(tree, all_taxa)
        # Non-trivial splits: {A,B}, {C,D}, {E,F}, {C,D,E,F}
        # {C,D,E,F} canonicalizes to {A,B} since |{A,B}| < |{C,D,E,F}|
        # So we get: {A,B}, {C,D}, {E,F}
        assert frozenset({"A", "B"}) in splits
        assert frozenset({"C", "D"}) in splits
        assert frozenset({"E", "F"}) in splits


class TestCountSplits:
    def test_counts_across_trees(self):
        trees = [
            _make_tree("((A,B),(C,D));"),
            _make_tree("((A,B),(C,D));"),
            _make_tree("((A,C),(B,D));"),
        ]
        all_taxa = frozenset({"A", "B", "C", "D"})
        counts = ConsensusNetwork._count_splits(trees, all_taxa)
        assert counts[frozenset({"A", "B"})] == 2
        assert counts[frozenset({"A", "C"})] == 1


class TestFilterSplits:
    def test_threshold_filtering(self):
        from collections import Counter
        counts = Counter({
            frozenset({"A", "B"}): 8,
            frozenset({"C", "D"}): 5,
            frozenset({"E", "F"}): 2,
        })
        result = ConsensusNetwork._filter_splits(counts, 10, 0.5)
        splits_in_result = [r[0] for r in result]
        assert frozenset({"A", "B"}) in splits_in_result
        assert frozenset({"C", "D"}) in splits_in_result
        assert frozenset({"E", "F"}) not in splits_in_result

    def test_sorted_by_frequency(self):
        from collections import Counter
        counts = Counter({
            frozenset({"A", "B"}): 5,
            frozenset({"C", "D"}): 8,
        })
        result = ConsensusNetwork._filter_splits(counts, 10, 0.1)
        assert result[0][0] == frozenset({"C", "D"})
        assert result[1][0] == frozenset({"A", "B"})

    def test_zero_threshold_includes_all(self):
        from collections import Counter
        counts = Counter({
            frozenset({"A", "B"}): 1,
        })
        result = ConsensusNetwork._filter_splits(counts, 10, 0.0)
        assert len(result) == 1


class TestComputeCircularOrdering:
    def test_contains_all_taxa(self):
        trees = [_make_tree("((A,B),(C,D));")]
        all_taxa = frozenset({"A", "B", "C", "D"})
        ordering = ConsensusNetwork._compute_circular_ordering(trees, all_taxa)
        assert set(ordering) == set(all_taxa)
        assert len(ordering) == len(all_taxa)


class TestConsensusNetworkRun:
    def test_text_output(self, tmp_path, capsys):
        tree_file = tmp_path / "trees.nwk"
        _write(tree_file, "((A,B),(C,D));\n((A,B),(C,D));\n((A,C),(B,D));\n")

        svc = ConsensusNetwork(
            Namespace(
                trees=str(tree_file),
                threshold=0.1,
                missing_taxa="error",
                plot_output=None,
                json=False,
            )
        )
        svc.run()
        captured = capsys.readouterr()
        assert "Number of input trees: 3" in captured.out
        assert "{A, B}" in captured.out

    def test_json_output(self, tmp_path, monkeypatch):
        tree_file = tmp_path / "trees.nwk"
        _write(tree_file, "((A,B),(C,D));\n((A,B),(C,D));\n")

        captured = {}
        svc = ConsensusNetwork(
            Namespace(
                trees=str(tree_file),
                threshold=0.1,
                missing_taxa="error",
                plot_output=None,
                json=True,
            )
        )
        monkeypatch.setattr(
            "phykit.services.tree.consensus_network.print_json",
            lambda payload: captured.setdefault("payload", payload),
        )
        svc.run()

        payload = captured["payload"]
        assert payload["input_tree_count"] == 2
        assert payload["taxa_count"] == 4
        assert len(payload["splits"]) > 0

    def test_single_tree(self, tmp_path, capsys):
        tree_file = tmp_path / "trees.nwk"
        _write(tree_file, "((A,B),(C,D));\n")

        svc = ConsensusNetwork(
            Namespace(
                trees=str(tree_file),
                threshold=0.1,
                missing_taxa="error",
                plot_output=None,
                json=False,
            )
        )
        svc.run()
        captured = capsys.readouterr()
        assert "Number of input trees: 1" in captured.out

    def test_all_identical_trees(self, tmp_path, capsys):
        tree_file = tmp_path / "trees.nwk"
        _write(tree_file, "((A,B),(C,D));\n" * 5)

        svc = ConsensusNetwork(
            Namespace(
                trees=str(tree_file),
                threshold=0.1,
                missing_taxa="error",
                plot_output=None,
                json=False,
            )
        )
        svc.run()
        captured = capsys.readouterr()
        assert "1.0000" in captured.out

    def test_missing_taxa_error_mode_raises(self, tmp_path):
        tree_file = tmp_path / "trees.nwk"
        _write(tree_file, "((A,B),(C,D));\n((A,B),(C,E));\n")

        svc = ConsensusNetwork(
            Namespace(
                trees=str(tree_file),
                threshold=0.1,
                missing_taxa="error",
                plot_output=None,
                json=False,
            )
        )
        with pytest.raises(PhykitUserError):
            svc.run()

    def test_missing_taxa_shared_mode(self, tmp_path, capsys):
        tree_file = tmp_path / "trees.nwk"
        _write(tree_file, "((A,B),(C,D));\n((A,B),(C,E));\n")

        svc = ConsensusNetwork(
            Namespace(
                trees=str(tree_file),
                threshold=0.1,
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
        _write(tree_file, "((A,B),(C,D));\n((A,B),(C,D));\n")
        plot_path = str(tmp_path / "network.png")

        svc = ConsensusNetwork(
            Namespace(
                trees=str(tree_file),
                threshold=0.1,
                missing_taxa="error",
                plot_output=plot_path,
                json=False,
            )
        )
        svc.run()
        assert (tmp_path / "network.png").exists()


class TestPolytomyHandling:
    """Polytomous nodes should be skipped in split extraction."""

    def test_star_tree_no_splits(self):
        """A 4-way star (A,B,C,D) produces no non-trivial splits."""
        tree = _make_tree("(A:1,B:1,C:1,D:1);")
        all_taxa = frozenset(["A", "B", "C", "D"])
        splits = ConsensusNetwork._extract_splits_from_tree(tree, all_taxa)
        assert len(splits) == 0

    def test_trifurcating_root_allowed(self):
        """A trifurcating root (standard unrooted) is not skipped."""
        tree = _make_tree("((A:1,B:1):1,C:1,D:1);")
        all_taxa = frozenset(["A", "B", "C", "D"])
        splits = ConsensusNetwork._extract_splits_from_tree(tree, all_taxa)
        assert len(splits) >= 1

    def test_internal_polytomy_skipped(self):
        """An internal polytomy produces no splits from that node."""
        tree = _make_tree("((A:1,B:1,C:1,D:1):1,E:1,F:1);")
        all_taxa = frozenset(["A", "B", "C", "D", "E", "F"])
        splits = ConsensusNetwork._extract_splits_from_tree(tree, all_taxa)
        # The 4-way polytomy (A,B,C,D) is skipped; only resolved subclades
        # contribute. No resolved subclades here, so no splits.
        assert len(splits) == 0

    def test_partial_polytomy_preserves_resolved(self):
        """Resolved subclades within a polytomy still produce splits."""
        tree = _make_tree("(((A:1,B:1):1,C:1,D:1):1,E:1,F:1);")
        all_taxa = frozenset(["A", "B", "C", "D", "E", "F"])
        splits = ConsensusNetwork._extract_splits_from_tree(tree, all_taxa)
        # (A,B) is resolved and should produce a split
        assert len(splits) >= 1
