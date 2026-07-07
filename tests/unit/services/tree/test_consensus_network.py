from argparse import Namespace
from io import StringIO
from pathlib import Path
import math
import subprocess
import sys

import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

from phykit.errors import PhykitUserError
from phykit.services.tree.consensus_network import ConsensusNetwork, _position_extent
import phykit.services.tree.consensus_network as consensus_network_module


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
import phykit.services.tree.consensus_network as module

assert hasattr(module.Phylo, "read")
assert hasattr(module.Consensus, "majority_consensus")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


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

    def test_empty_split_still_canonicalizes_to_empty_set(self):
        assert ConsensusNetwork._canonical_split(frozenset(), frozenset()) == frozenset()

    def test_canonical_split_mask_equal_size_uses_lexicographic_tiebreak(self):
        names = ("A", "B", "C", "D")
        all_mask = (1 << len(names)) - 1
        ab_mask = (1 << 0) | (1 << 1)
        cd_mask = (1 << 2) | (1 << 3)

        assert (
            ConsensusNetwork._canonical_split_mask(
                cd_mask,
                all_mask,
                len(names),
                names,
            )
            == ab_mask
        )
        assert (
            ConsensusNetwork._canonical_split_mask(
                ab_mask,
                all_mask,
                len(names),
                names,
            )
            == ab_mask
        )


class TestTips:
    def test_uses_fast_terminal_names_for_parsed_tree(self, monkeypatch):
        tree = _make_tree("((A,B),(C,D));")

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("get_terminals fallback should not be called")

        monkeypatch.setattr(tree, "get_terminals", fail_get_terminals)

        assert ConsensusNetwork._tips(tree) == {"A", "B", "C", "D"}

    def test_prune_to_taxa_uses_direct_terminal_collection(self, monkeypatch):
        tree = _make_tree("((A:1,B:2):3,(C:4,D:5):6);")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("standard tree should use direct terminal traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        ConsensusNetwork._prune_to_taxa(tree, {"B", "D"})

        assert ConsensusNetwork._tips(tree) == {"B", "D"}


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

    def test_extract_splits_matches_legacy_terminal_scan(self):
        tree = _make_tree("((A,B),((C,D),(E,F)));")
        all_taxa = frozenset({"A", "B", "C", "D", "E", "F"})
        expected = set()
        for clade in tree.get_nonterminals():
            n_children = len(clade.clades)
            if n_children > 2:
                is_root = clade == tree.root
                if not (is_root and n_children == 3):
                    continue
            tips = frozenset(tip.name for tip in clade.get_terminals())
            if len(tips) <= 1 or len(tips) >= len(all_taxa) - 1:
                continue
            if tips == all_taxa:
                continue
            expected.add(ConsensusNetwork._canonical_split(tips, all_taxa))

        assert ConsensusNetwork._extract_splits_from_tree(tree, all_taxa) == expected

    def test_extract_splits_avoids_repeated_terminal_scans(self, monkeypatch):
        from Bio.Phylo.BaseTree import Clade

        tree = _make_tree("((A,B),((C,D),(E,F)));")
        all_taxa = frozenset({"A", "B", "C", "D", "E", "F"})

        def fail_get_terminals(self):
            raise AssertionError("get_terminals should not be called")

        monkeypatch.setattr(Clade, "get_terminals", fail_get_terminals)

        splits = ConsensusNetwork._extract_splits_from_tree(tree, all_taxa)

        assert frozenset({"A", "B"}) in splits

    def test_extract_splits_uses_direct_mask_traversal(self, monkeypatch):
        tree = _make_tree("((A,B),((C,D),(E,F)));")
        all_taxa = frozenset({"A", "B", "C", "D", "E", "F"})

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("generic clade traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        splits = ConsensusNetwork._extract_splits_from_tree(tree, all_taxa)

        assert frozenset({"A", "B"}) in splits

    def test_direct_mask_traversal_handles_mixed_child_counts(self, monkeypatch):
        tree = _make_tree("(A,(B,C),(D,(E,F),G),H);")
        all_taxa = frozenset({"A", "B", "C", "D", "E", "F", "G", "H"})
        expected = ConsensusNetwork._extract_splits_from_tree_legacy(
            tree,
            all_taxa,
        )

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("generic clade traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        observed = ConsensusNetwork._extract_splits_from_tree(tree, all_taxa)

        assert observed == expected


class TestCountSplits:
    def test_counts_across_trees(self):
        trees = [
            _make_tree("((A,B),(C,D));"),
            _make_tree("((A,B),(C,D));"),
            _make_tree("((A,C),(B,D));"),
        ]
        all_taxa = frozenset({"A", "B", "C", "D"})
        counts, possible = ConsensusNetwork._count_splits(trees, all_taxa)
        assert counts[frozenset({"A", "B"})] == 2
        assert counts[frozenset({"A", "C"})] == 1

    def test_counts_shared_taxa_with_direct_mask_counter(self, monkeypatch):
        trees = [
            _make_tree("((A,B),(C,D));"),
            _make_tree("((A,B),(C,D));"),
            _make_tree("((A,C),(B,D));"),
        ]
        all_taxa = frozenset({"A", "B", "C", "D"})

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("generic clade traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        counts, possible = ConsensusNetwork._count_splits(trees, all_taxa)

        assert counts[frozenset({"A", "B"})] == 2
        assert counts[frozenset({"A", "C"})] == 1
        assert possible[frozenset({"A", "B"})] == len(trees)

    def test_allow_mode_uses_fast_taxon_sets(self, monkeypatch):
        trees = [
            _make_tree("((A,B),(C,D));"),
            _make_tree("((A,B),(C,E));"),
        ]
        all_taxa = frozenset({"A", "B", "C", "D", "E"})

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("allow-mode taxon setup should use fast names")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        counts, possible = ConsensusNetwork._count_splits(
            trees,
            all_taxa,
            allow_mode=True,
        )

        assert counts[frozenset({"A", "B"})] == 2
        assert possible[frozenset({"A", "B"})] == len(trees)

    def test_allow_mode_batches_split_counts_with_counter_update(
        self, monkeypatch
    ):
        class RecordingCounter(dict):
            instances = []

            def __init__(self, *args, **kwargs):
                super().__init__()
                self.update_calls = []
                self.__class__.instances.append(self)
                if args or kwargs:
                    self.update(*args, **kwargs)

            def __missing__(self, key):
                return 0

            def update(self, iterable=None, **kwargs):
                if iterable is None:
                    items = ()
                elif hasattr(iterable, "items"):
                    items = tuple(iterable.items())
                    for key, value in items:
                        self[key] = value
                    self.update_calls.append(items)
                    return
                else:
                    items = tuple(iterable)
                self.update_calls.append(frozenset(items))
                for key in items:
                    self[key] = self.get(key, 0) + 1
                for key, value in kwargs.items():
                    self[key] = value

        split_ab = frozenset({"A", "B"})
        split_ac = frozenset({"A", "C"})
        trees = [object(), object(), object()]

        monkeypatch.setattr(consensus_network_module, "Counter", RecordingCounter)
        monkeypatch.setattr(
            ConsensusNetwork,
            "_tips",
            staticmethod(lambda _tree: {"A", "B", "C", "D"}),
        )
        monkeypatch.setattr(
            ConsensusNetwork,
            "_extract_splits_from_tree",
            staticmethod(
                lambda tree, _taxa: (
                    {split_ab, split_ac} if tree is trees[0] else {split_ab}
                )
            ),
        )

        counts, possible = ConsensusNetwork._count_splits(
            trees,
            frozenset({"A", "B", "C", "D"}),
            allow_mode=True,
        )

        assert counts == {split_ab: 3, split_ac: 1}
        assert possible == {split_ab: 3, split_ac: 3}
        assert RecordingCounter.instances[0].update_calls == [
            frozenset({split_ab, split_ac}),
            frozenset({split_ab}),
            frozenset({split_ab}),
        ]


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

    def test_uses_direct_terminal_name_traversal(self, monkeypatch):
        consensus = _make_tree("((C,D),(A,B));")
        trees = [consensus]
        all_taxa = frozenset({"A", "B", "C", "D"})

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("consensus ordering should use fast terminal names")

        monkeypatch.setattr(
            "phykit.services.tree.consensus_network.Consensus.majority_consensus",
            lambda *_args, **_kwargs: consensus,
        )
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        ordering = ConsensusNetwork._compute_circular_ordering(trees, all_taxa)
        assert ordering == ["C", "D", "A", "B"]


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

        assert ConsensusNetwork._is_circular_split(
            frozenset(["E", "F", "A"]),
            ordering,
        )
        assert ConsensusNetwork._is_circular_split(
            frozenset(["B", "C", "D"]),
            ordering,
        )
        assert not ConsensusNetwork._is_circular_split(
            frozenset(["A", "C", "E"]),
            ordering,
        )
        assert not ConsensusNetwork._is_circular_split(frozenset(), [])

    def test_is_circular_split_stops_after_third_boundary(self):
        ordering = self.CountingOrdering(["A", "B", "C", "D", "E", "F"])

        assert not ConsensusNetwork._is_circular_split(
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

        observed = ConsensusNetwork._compute_split_directions(
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
            split: ConsensusNetwork._circular_gap_positions(split, ordering)
            for split, _count, _freq in circular_splits
        }

        def fail_gap_scan(*args, **kwargs):
            raise AssertionError("cached gap positions should be reused")

        monkeypatch.setattr(ConsensusNetwork, "_circular_gap_positions", fail_gap_scan)

        observed = ConsensusNetwork._compute_split_directions(
            ordering,
            circular_splits,
            gap_positions_by_split,
        )
        expected = _legacy_compute_split_directions(ordering, circular_splits)

        assert observed.keys() == expected.keys()
        for split, direction in expected.items():
            assert observed[split] == pytest.approx(direction)

    def test_compute_split_directions_reuses_taxon_coordinates(self, monkeypatch):
        ordering = [f"T{i}" for i in range(12)]
        splits = [
            frozenset(ordering[start:start + size])
            for start, size in (
                (0, 2),
                (1, 3),
                (2, 4),
                (4, 3),
                (5, 5),
            )
        ]
        circular_splits = [
            (split, idx + 1, 0.1 * (idx + 1))
            for idx, split in enumerate(splits)
        ]
        gap_positions_by_split = {
            split: ConsensusNetwork._circular_gap_positions(split, ordering)
            for split in splits
        }
        real_cos = math.cos
        cos_calls = 0

        def counted_cos(value):
            nonlocal cos_calls
            cos_calls += 1
            return real_cos(value)

        monkeypatch.setattr(consensus_network_module.math, "cos", counted_cos)

        directions = ConsensusNetwork._compute_split_directions(
            ordering,
            circular_splits,
            gap_positions_by_split,
        )

        assert len(directions) == len(circular_splits)
        assert cos_calls == len(ordering) + 2 * len(circular_splits)

    def test_compute_split_directions_accumulates_split_center_once(self):
        class CountingSplit(frozenset):
            def __new__(cls, values):
                obj = super().__new__(cls, values)
                obj.iterations = 0
                return obj

            def __iter__(self):
                self.iterations += 1
                return super().__iter__()

        ordering = ["A", "B", "C", "D", "E", "F"]
        split = CountingSplit(["A", "B", "C"])
        circular_splits = [(split, 1, 0.5)]
        gap_positions_by_split = {
            split: ConsensusNetwork._circular_gap_positions(split, ordering)
        }
        split.iterations = 0

        directions = ConsensusNetwork._compute_split_directions(
            ordering,
            circular_splits,
            gap_positions_by_split,
        )

        assert split in directions
        assert split.iterations == 1

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

        observed = ConsensusNetwork._build_splits_graph(circular_splits, all_taxa)
        expected = self._legacy_build_splits_graph(circular_splits, all_taxa)

        assert observed[0] == expected[0]
        assert observed[1] == expected[1]
        assert self._canonical_edges(observed[2]) == self._canonical_edges(expected[2])
        assert observed[3] == expected[3]
        assert observed[4] == expected[4]

    def test_build_splits_graph_independent_splits_has_hypercube_edges(self):
        n_splits = 4
        taxa = [f"T{mask:04b}" for mask in range(1 << n_splits)]
        all_taxa = frozenset(taxa)
        circular_splits = [
            (
                frozenset(
                    taxa[mask]
                    for mask in range(1 << n_splits)
                    if mask & (1 << split_idx)
                ),
                1,
                1.0,
            )
            for split_idx in range(n_splits)
        ]

        taxon_signs, valid_nodes, edges, splits_list, weights = (
            ConsensusNetwork._build_splits_graph(circular_splits, all_taxa)
        )

        assert len(taxon_signs) == len(taxa)
        assert len(valid_nodes) == 1 << n_splits
        expected_edges = n_splits * (1 << (n_splits - 1))
        assert len(self._canonical_edges(edges)) == expected_edges
        assert len(edges) == expected_edges
        assert splits_list == [split[0] for split in circular_splits]
        assert weights == [1.0] * n_splits


class TestNetworkPlot:
    def test_draw_network_batches_edges(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection
        from phykit.helpers.plot_config import PlotConfig

        service = ConsensusNetwork.__new__(ConsensusNetwork)
        service.plot_config = PlotConfig(show_title=False)
        ordering = ["A", "B", "C", "D"]
        all_taxa = frozenset(ordering)
        filtered_splits = [
            (frozenset({"A", "B"}), 7, 0.7),
            (frozenset({"B", "C"}), 5, 0.5),
        ]
        line_collections = []
        original_add_collection = matplotlib.axes.Axes.add_collection

        def fail_plot(*args, **kwargs):
            raise AssertionError("network edges should use LineCollection")

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(matplotlib.axes.Axes, "add_collection", capture_collection)

        output_path = str(tmp_path / "consensus_network.png")
        service._draw_network(ordering, filtered_splits, all_taxa, output_path)

        assert len(line_collections) >= 2
        assert Path(output_path).exists()

    def test_draw_network_skips_redundant_tight_layout(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.pyplot as plt
        from phykit.helpers.plot_config import PlotConfig

        service = ConsensusNetwork.__new__(ConsensusNetwork)
        service.plot_config = PlotConfig(show_title=False)
        ordering = ["A", "B", "C", "D"]
        all_taxa = frozenset(ordering)
        filtered_splits = [
            (frozenset({"A", "B"}), 7, 0.7),
            (frozenset({"B", "C"}), 5, 0.5),
        ]

        def fail_tight_layout(*args, **kwargs):
            raise AssertionError("bbox_inches='tight' handles saved bounds")

        monkeypatch.setattr(plt, "tight_layout", fail_tight_layout)

        output_path = str(tmp_path / "consensus_network_no_tight_layout.png")
        service._draw_network(ordering, filtered_splits, all_taxa, output_path)

        assert Path(output_path).exists()
        assert Path(output_path).stat().st_size > 0

    def test_draw_network_batches_unlabeled_fallback_points(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from phykit.helpers.plot_config import PlotConfig

        service = ConsensusNetwork.__new__(ConsensusNetwork)
        service.plot_config = PlotConfig(show_title=False, ylabel_fontsize=0)
        ordering = ["A", "B", "C", "D", "E"]
        all_taxa = frozenset(ordering)

        original_scatter = matplotlib.axes.Axes.scatter
        scatter_calls = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("unlabeled fallback points should use scatter")

        def capture_scatter(self, x, y, *args, **kwargs):
            scatter_calls.append((x, y, kwargs))
            return original_scatter(self, x, y, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(matplotlib.axes.Axes, "scatter", capture_scatter)

        output_path = str(tmp_path / "consensus_network_points.png")
        service._draw_network(ordering, [], all_taxa, output_path)

        assert len(scatter_calls) == 1
        x, y, kwargs = scatter_calls[0]
        assert len(x) == len(y) == len(ordering)
        assert kwargs["color"] == "black"
        assert Path(output_path).exists()


class TestPruneToTaxa:
    def test_avoids_generic_terminal_collection(self, monkeypatch):
        tree = _make_tree("((A,B),(C,D));")
        original_get_terminals = tree.get_terminals
        calls = 0

        def counting_get_terminals():
            nonlocal calls
            calls += 1
            if calls > 1:
                raise AssertionError("get_terminals should only be called once")
            return original_get_terminals()

        monkeypatch.setattr(tree, "get_terminals", counting_get_terminals)

        ConsensusNetwork._prune_to_taxa(tree, {"A", "C"})

        assert calls == 0
        assert {tip.name for tip in original_get_terminals()} == {"A", "C"}

    def test_uses_batch_prune_for_standard_tree(self, monkeypatch):
        tree = _make_tree("((A:1,B:2):3,(C:4,D:5):6);")

        def fail_prune(*_args, **_kwargs):
            raise AssertionError("per-target prune should not be used")

        monkeypatch.setattr(TreeMixin, "prune", fail_prune)

        ConsensusNetwork._prune_to_taxa(tree, {"B", "D"})

        assert [tip.name for tip in tree.get_terminals()] == ["B", "D"]
        assert [tip.branch_length for tip in tree.get_terminals()] == [5.0, 11.0]


class TestConsensusNetworkRun:
    def test_parse_trees_skips_comments_and_blank_lines(self, tmp_path):
        tree_file = tmp_path / "trees.nwk"
        _write(
            tree_file,
            "   # ignored\n\n  ((A,B),(C,D));  \n\t# also ignored\n((A,C),(B,D));\n",
        )
        svc = ConsensusNetwork(
            Namespace(
                trees=str(tree_file),
                threshold=0.1,
                missing_taxa="error",
                plot_output=None,
                json=False,
            )
        )

        trees = svc._parse_trees_from_source(str(tree_file))

        assert len(trees) == 2
        assert [sorted(ConsensusNetwork._tips(tree)) for tree in trees] == [
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
            "phykit.services.tree.consensus_network",
            fromlist=["Path"],
        ).Path

        def counting_path(*args, **kwargs):
            nonlocal path_calls
            path_calls += 1
            return original_path(*args, **kwargs)

        mocker.patch(
            "phykit.services.tree.consensus_network.Path",
            side_effect=counting_path,
        )
        svc = ConsensusNetwork(
            Namespace(
                trees=str(tree_list),
                threshold=0.1,
                missing_taxa="error",
                plot_output=None,
                json=False,
            )
        )

        trees = svc._parse_trees_from_source(str(tree_list))

        assert len(trees) == 2
        assert path_calls == 1

    def test_normalize_taxa_identical_sets_uses_fast_path(self, monkeypatch):
        svc = ConsensusNetwork(
            Namespace(
                trees="unused",
                threshold=0.1,
                missing_taxa="error",
                plot_output=None,
                json=False,
            )
        )
        trees = [
            _make_tree("((A,B),(C,D));"),
            _make_tree("((A,C),(B,D));"),
        ]

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("parsed trees should use fast terminal names")

        def fail_prune(*_args, **_kwargs):
            raise AssertionError("identical taxon sets should not be pruned")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(ConsensusNetwork, "_prune_to_taxa", fail_prune)

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
            {"A", "B", "C"},
            {"A", "B", "C"},
            {"A", "B", "C"},
        ])

        assert consensus_network_module._all_tip_sets_identical(tip_sets) is True

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

    def test_text_output_batches_filtered_splits(self, tmp_path, mocker):
        tree_file = tmp_path / "trees.nwk"
        tree_file.write_text("((A,B),(C,D));\n")
        svc = ConsensusNetwork(
            Namespace(
                trees=str(tree_file),
                threshold=0.25,
                missing_taxa="shared",
                max_splits=30,
                histogram=None,
                plot_output=None,
                json=False,
            )
        )
        split_counts = {
            frozenset({"A", "B"}): 3,
            frozenset({"A", "C"}): 2,
        }
        filtered = [
            (frozenset({"A", "B"}), 3, 1.0),
            (frozenset({"A", "C"}), 2, 2 / 3),
        ]

        mocker.patch.object(
            svc,
            "_parse_trees_from_source",
            return_value=["t1", "t2", "t3"],
        )
        mocker.patch.object(
            svc,
            "_normalize_taxa",
            return_value=(["t1", "t2", "t3"], True, {"A", "B", "C", "D"}),
        )
        mocker.patch.object(svc, "_count_splits", return_value=(split_counts, None))
        mocker.patch.object(svc, "_filter_splits", return_value=filtered)
        printed = mocker.patch("builtins.print")

        svc.run()

        printed.assert_called_once_with(
            "Number of input trees: 3\n"
            "Number of taxa: 4\n"
            "Threshold: 0.25\n"
            "Pruned to shared taxa: yes\n"
            "Total unique splits: 2\n"
            "Splits above threshold: 2\n"
            "---\n"
            "{A, B}\t3/3\t1.0000\n"
            "{A, C}\t2/3\t0.6667"
        )

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
