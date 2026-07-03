import csv
import math
import os
import subprocess
import sys
from argparse import Namespace
from io import StringIO

import numpy as np
import pytest

from phykit.errors import PhykitUserError
from phykit.services.tree.neighbor_net import NeighborNet, _position_extent
import phykit.services.tree.neighbor_net as neighbor_net_module


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


def test_module_import_does_not_import_numpy_or_biopython_fasta_parser():
    code = """
import sys
import phykit.services.tree.neighbor_net as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_proxy_caches_resolved_attributes():
    lazy_np = neighbor_net_module._LazyNumpy()

    first_empty = lazy_np.empty
    second_empty = lazy_np.empty

    assert first_empty is second_empty
    assert lazy_np.__dict__["empty"] is first_empty


def _legacy_compute_split_directions(ordering, circular_splits):
    n = len(ordering)
    angles = {taxon: 2 * math.pi * i / n for i, taxon in enumerate(ordering)}
    directions = {}
    for split, weight in circular_splits:
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


SAMPLE_ALIGNMENT = os.path.join(
    os.path.dirname(__file__),
    os.pardir,
    os.pardir,
    os.pardir,
    "sample_files",
    "12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit",
)


def _write_csv_matrix(path, taxa, D):
    """Helper to write a distance matrix as CSV."""
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([""] + taxa)
        for i, taxon in enumerate(taxa):
            writer.writerow([taxon] + [f"{D[i, j]:.6f}" for j in range(len(taxa))])


def _make_simple_distance_matrix():
    """4-taxon distance matrix with clear structure."""
    taxa = ["A", "B", "C", "D"]
    D = np.array(
        [
            [0.0, 0.1, 0.5, 0.5],
            [0.1, 0.0, 0.5, 0.5],
            [0.5, 0.5, 0.0, 0.1],
            [0.5, 0.5, 0.1, 0.0],
        ]
    )
    return taxa, D


class TestDistanceMatrixFromAlignment:
    def test_read_alignment_preserves_taxa_order_uppercases_and_keeps_last_duplicate(
        self, tmp_path
    ):
        path = tmp_path / "alignment.fa"
        path.write_text(
            ">A description\nac gt\n"
            ">B\nn- x\n?\n"
            ">A duplicate\nTGCA\n"
        )

        taxa, seqs = NeighborNet._read_alignment(str(path))

        assert taxa == ["A", "B", "A"]
        assert seqs == {"A": "TGCA", "B": "N-X?"}

    def test_vectorized_distance_matrix_matches_legacy_loop(self, tmp_path):
        output = str(tmp_path / "network.png")
        args = Namespace(
            alignment=SAMPLE_ALIGNMENT,
            distance_matrix=None,
            output=output,
            metric="p-distance",
            max_splits=30,
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
        service = NeighborNet(args)
        taxa = ["A", "B", "C", "D"]
        seqs = {
            "A": "ACGT-NX?",
            "B": "ACGTACGT",
            "C": "TGCAACGT",
            "D": "NNNN----",
        }

        for metric in ("p-distance", "identity", "jc"):
            service.metric = metric
            fast = service._compute_distance_matrix_from_equal_length_sequences(
                taxa, seqs
            )
            legacy = service._compute_distance_matrix_from_sequences_legacy(
                taxa, seqs
            )
            np.testing.assert_allclose(fast, legacy)

    def test_equal_length_distance_matrix_uses_float_count_products(
        self, mocker, tmp_path
    ):
        output = str(tmp_path / "network.png")
        args = Namespace(
            alignment=SAMPLE_ALIGNMENT,
            distance_matrix=None,
            output=output,
            metric="p-distance",
            max_splits=30,
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
        service = NeighborNet(args)
        zeros_like_spy = mocker.spy(neighbor_net_module.np, "zeros_like")

        service._compute_distance_matrix_from_equal_length_sequences(
            ["A", "B"],
            {"A": "ACGT", "B": "AGGT"},
        )

        assert any(
            call.kwargs.get("dtype") is neighbor_net_module.np.float64
            for call in zeros_like_spy.call_args_list
        )

    def test_clean_ascii_distance_matrix_uses_direct_comparison(
        self, mocker, monkeypatch, tmp_path
    ):
        output = str(tmp_path / "network.png")
        args = Namespace(
            alignment=SAMPLE_ALIGNMENT,
            distance_matrix=None,
            output=output,
            metric="p-distance",
            max_splits=30,
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
        service = NeighborNet(args)
        monkeypatch.setattr(neighbor_net_module, "_NO_SKIP_DIRECT_MIN_LENGTH", 1)
        monkeypatch.setattr(neighbor_net_module, "_NO_SKIP_DIRECT_MAX_TAXA", 10)
        direct_path = mocker.spy(
            service,
            "_compute_distance_matrix_from_no_skip_matrix",
        )
        mocker.patch.object(
            neighbor_net_module.np,
            "isin",
            side_effect=AssertionError(
                "clean direct comparison should skip validity-mask setup"
            ),
        )

        observed = service._compute_distance_matrix_from_equal_length_sequences(
            ["A", "B", "C"],
            {
                "A": "ACGTACGT",
                "B": "ACGTTCGT",
                "C": "TCGTACGA",
            },
        )

        np.testing.assert_allclose(
            observed,
            np.array(
                [
                    [0.0, 0.125, 0.25],
                    [0.125, 0.0, 0.375],
                    [0.25, 0.375, 0.0],
                ]
            ),
        )
        direct_path.assert_called_once()

    def test_identity_no_overlap_stays_zero(self, tmp_path):
        output = str(tmp_path / "network.png")
        args = Namespace(
            alignment=SAMPLE_ALIGNMENT,
            distance_matrix=None,
            output=output,
            metric="identity",
            max_splits=30,
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
        service = NeighborNet(args)
        taxa = ["A", "B"]
        seqs = {"A": "NNNN", "B": "----"}

        D = service._compute_distance_matrix_from_sequences(taxa, seqs)

        assert D[0, 1] == 0.0
        assert D[1, 0] == 0.0

    def test_distance_matrix_symmetric(self, tmp_path):
        output = str(tmp_path / "network.png")
        args = Namespace(
            alignment=SAMPLE_ALIGNMENT,
            distance_matrix=None,
            output=output,
            metric="p-distance",
            max_splits=30,
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
        service = NeighborNet(args)
        taxa, D = service._compute_distances_from_alignment()
        # Symmetry check
        np.testing.assert_array_almost_equal(D, D.T)

    def test_distance_matrix_diagonal_zero(self, tmp_path):
        output = str(tmp_path / "network.png")
        args = Namespace(
            alignment=SAMPLE_ALIGNMENT,
            distance_matrix=None,
            output=output,
            metric="p-distance",
            max_splits=30,
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
        service = NeighborNet(args)
        taxa, D = service._compute_distances_from_alignment()
        for i in range(len(taxa)):
            assert D[i, i] == 0.0

    def test_jc_metric(self, tmp_path):
        output = str(tmp_path / "network.png")
        args = Namespace(
            alignment=SAMPLE_ALIGNMENT,
            distance_matrix=None,
            output=output,
            metric="jc",
            max_splits=30,
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
        service = NeighborNet(args)
        taxa, D_jc = service._compute_distances_from_alignment()
        # JC distances should be >= 0 and symmetric
        assert np.all(D_jc >= 0)
        np.testing.assert_array_almost_equal(D_jc, D_jc.T)
        # JC distances should generally be >= p-distances
        service.metric = "p-distance"
        _, D_p = service._compute_distances_from_alignment()
        for i in range(len(taxa)):
            for j in range(i + 1, len(taxa)):
                if D_jc[i, j] < 10.0:  # not saturated
                    assert D_jc[i, j] >= D_p[i, j] - 1e-10


class TestCircularOrdering:
    def test_direct_terminal_names_preserves_biopython_order(self):
        from Bio import Phylo

        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1,E:1):1);"), "newick")

        assert NeighborNet._terminal_names_direct(tree) == [
            tip.name for tip in tree.get_terminals()
        ]

    def test_nj_circular_ordering_uses_direct_terminal_names(self, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_get_terminals(self):
            raise AssertionError("get_terminals should not be called")

        monkeypatch.setattr(DistanceTreeConstructor, "nj", lambda self, dm: tree)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        taxa, D = _make_simple_distance_matrix()

        assert NeighborNet._nj_circular_ordering(D, taxa) == ["A", "B", "C", "D"]

    def test_nj_circular_ordering_uses_lower_triangle_rows(self, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo import TreeConstruction

        captured = {}
        tree = Phylo.read(StringIO("(A:1,B:1,C:1);"), "newick")

        class FakeDistanceMatrix:
            def __init__(self, names, matrix):
                captured["names"] = names
                captured["matrix"] = matrix

        monkeypatch.setattr(TreeConstruction, "DistanceMatrix", FakeDistanceMatrix)
        monkeypatch.setattr(
            TreeConstruction.DistanceTreeConstructor,
            "nj",
            lambda self, dm: tree,
        )

        taxa = ["A", "B", "C"]
        D = np.array(
            [
                [0.0, 1.5, 2.5],
                [1.5, 0.0, 3.5],
                [2.5, 3.5, 0.0],
            ]
        )

        assert NeighborNet._nj_circular_ordering(D, taxa) == taxa
        assert captured == {
            "names": taxa,
            "matrix": [[0.0], [1.5, 0.0], [2.5, 3.5, 0.0]],
        }

    def test_contains_all_taxa(self):
        taxa, D = _make_simple_distance_matrix()
        ordering = NeighborNet._nj_circular_ordering(D, taxa)
        assert set(ordering) == set(taxa)

    def test_ordering_length(self):
        taxa, D = _make_simple_distance_matrix()
        ordering = NeighborNet._nj_circular_ordering(D, taxa)
        assert len(ordering) == len(taxa)


class TestCircularSplits:
    class CountingOrdering(list):
        def __init__(self, values):
            super().__init__(values)
            self.index_count = 0

        def __getitem__(self, index):
            self.index_count += 1
            return super().__getitem__(index)

    def test_splits_are_contiguous(self):
        ordering = ["A", "B", "C", "D", "E"]
        splits = NeighborNet._enumerate_circular_splits(ordering)
        n = len(ordering)
        for split in splits:
            # Verify the split forms a contiguous arc
            assert NeighborNet._is_circular_split(split, ordering)

    def test_is_circular_split_handles_wraparound_boundaries(self):
        ordering = ["A", "B", "C", "D", "E", "F"]

        assert NeighborNet._is_circular_split(frozenset(["E", "F", "A"]), ordering)
        assert NeighborNet._is_circular_split(frozenset(["B", "C", "D"]), ordering)
        assert not NeighborNet._is_circular_split(frozenset(["A", "C", "E"]), ordering)
        assert not NeighborNet._is_circular_split(frozenset(), [])

    def test_is_circular_split_stops_after_third_boundary(self):
        ordering = self.CountingOrdering(["A", "B", "C", "D", "E", "F"])

        assert not NeighborNet._is_circular_split(
            frozenset(["A", "C", "E"]),
            ordering,
        )
        assert ordering.index_count < len(ordering)

    def test_compute_split_directions_matches_legacy_wraparound_boundaries(self):
        ordering = ["A", "B", "C", "D", "E", "F"]
        circular_splits = [
            (frozenset(["E", "F", "A"]), 1.0),
            (frozenset(["B", "C", "D"]), 0.5),
            (frozenset(["A", "C", "E"]), 0.25),
        ]

        observed = NeighborNet._compute_split_directions(ordering, circular_splits)
        expected = _legacy_compute_split_directions(ordering, circular_splits)

        assert observed.keys() == expected.keys()
        for split, direction in expected.items():
            assert observed[split] == pytest.approx(direction)

    def test_compute_split_directions_uses_supplied_gap_positions(self, monkeypatch):
        ordering = ["A", "B", "C", "D", "E", "F"]
        circular_splits = [
            (frozenset(["E", "F", "A"]), 1.0),
            (frozenset(["B", "C", "D"]), 0.5),
        ]
        gap_positions_by_split = {
            split: NeighborNet._circular_gap_positions(split, ordering)
            for split, _weight in circular_splits
        }

        def fail_gap_scan(*args, **kwargs):
            raise AssertionError("cached gap positions should be reused")

        monkeypatch.setattr(NeighborNet, "_circular_gap_positions", fail_gap_scan)

        observed = NeighborNet._compute_split_directions(
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
        splits = NeighborNet._enumerate_circular_splits(ordering)[:8]
        circular_splits = [(split, 1.0) for split in splits]
        gap_positions_by_split = {
            split: NeighborNet._circular_gap_positions(split, ordering)
            for split in splits
        }
        real_cos = math.cos
        cos_calls = 0

        def counted_cos(value):
            nonlocal cos_calls
            cos_calls += 1
            return real_cos(value)

        monkeypatch.setattr(neighbor_net_module.math, "cos", counted_cos)

        directions = NeighborNet._compute_split_directions(
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
        circular_splits = [(split, 1.0)]
        gap_positions_by_split = {
            split: NeighborNet._circular_gap_positions(split, ordering)
        }
        split.iterations = 0

        directions = NeighborNet._compute_split_directions(
            ordering,
            circular_splits,
            gap_positions_by_split,
        )

        assert split in directions
        assert split.iterations == 1

    def test_split_count_four_taxa(self):
        ordering = ["A", "B", "C", "D"]
        splits = NeighborNet._enumerate_circular_splits(ordering)
        # For 4 taxa in a circle, non-trivial circular splits:
        # size-2 arcs: {A,B}, {B,C}, {C,D}, {D,A} = 4 splits
        # but canonical deduplication: {A,B}=={C,D}, {B,C}=={A,D}
        # So we should get 2 unique canonical splits
        assert len(splits) == 2

    def test_no_trivial_splits(self):
        ordering = ["A", "B", "C", "D", "E"]
        splits = NeighborNet._enumerate_circular_splits(ordering)
        for split in splits:
            assert len(split) >= 2
            assert len(split) <= len(ordering) - 2

    def test_enumerate_circular_splits_avoids_full_order_rescans(self):
        class CountingOrdering(list):
            def __init__(self, values):
                super().__init__(values)
                self.iter_count = 0

            def __iter__(self):
                self.iter_count += 1
                return super().__iter__()

        ordering = CountingOrdering(["A", "B", "C", "D", "E", "F"])

        splits = NeighborNet._enumerate_circular_splits(ordering)

        assert splits
        assert ordering.iter_count == 0

    def test_enumerate_circular_splits_matches_nested_index_reference(self):
        ordering = [f"T{i}" for i in range(12)]
        seen = set()
        expected = []
        for start in range(len(ordering)):
            for size in range(2, len(ordering) - 1):
                complement_size = len(ordering) - size
                side = frozenset(
                    ordering[(start + offset) % len(ordering)]
                    for offset in range(size)
                )
                if size < complement_size:
                    canonical = side
                elif size > complement_size:
                    canonical = frozenset(
                        ordering[(start + size + offset) % len(ordering)]
                        for offset in range(complement_size)
                    )
                else:
                    complement = frozenset(
                        ordering[(start + size + offset) % len(ordering)]
                        for offset in range(complement_size)
                    )
                    if sorted(side) < sorted(complement):
                        canonical = side
                    else:
                        canonical = complement
                if canonical not in seen:
                    seen.add(canonical)
                    expected.append(canonical)

        observed = NeighborNet._enumerate_circular_splits(ordering)

        assert observed == expected

    def test_equal_size_split_uses_sorted_lexicographic_tiebreak(self):
        ordering = ["D", "C", "B", "A"]

        splits = NeighborNet._enumerate_circular_splits(ordering)

        assert splits == [frozenset({"A", "B"}), frozenset({"A", "D"})]

    def test_canonical_split_equal_size_uses_minimum_taxon_tiebreak(self):
        all_taxa = frozenset({"A", "B", "C", "D"})

        assert NeighborNet._canonical_split(
            frozenset({"C", "D"}), all_taxa
        ) == frozenset({"A", "B"})
        assert NeighborNet._canonical_split(
            frozenset({"A", "B"}), all_taxa
        ) == frozenset({"A", "B"})


class TestNNLSWeights:
    def test_weights_nonneg(self):
        taxa, D = _make_simple_distance_matrix()
        ordering = NeighborNet._nj_circular_ordering(D, taxa)
        splits = NeighborNet._enumerate_circular_splits(ordering)
        weights = NeighborNet._estimate_split_weights(D, splits, ordering, taxa)
        assert np.all(weights >= -1e-10)

    def test_weights_match_explicit_design_matrix(self):
        from scipy.optimize import nnls

        taxa, D = _make_simple_distance_matrix()
        ordering = ["A", "B", "C", "D"]
        splits = [
            frozenset({"A", "B"}),
            frozenset({"B", "C"}),
            frozenset({"A", "C"}),
        ]

        observed = NeighborNet._estimate_split_weights(D, splits, ordering, taxa)

        pairs = [(i, j) for i in range(len(taxa)) for j in range(i + 1, len(taxa))]
        A = np.zeros((len(pairs), len(splits)))
        b = np.zeros(len(pairs))
        for pair_idx, (i, j) in enumerate(pairs):
            b[pair_idx] = D[i, j]
            for split_idx, split in enumerate(splits):
                A[pair_idx, split_idx] = (taxa[i] in split) != (taxa[j] in split)
        expected, _ = nnls(A, b)

        np.testing.assert_allclose(observed, expected)


class TestSplitsGraph:
    @staticmethod
    def _legacy_build_splits_graph(circular_splits, all_taxa):
        splits_list = [s[0] for s in circular_splits]
        weights = [s[1] for s in circular_splits]
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
        ordering = ["A", "B", "C", "D", "E", "F", "G"]
        all_taxa = frozenset(ordering)
        splits = NeighborNet._enumerate_circular_splits(ordering)[:8]
        circular_splits = [(split, float(idx + 1)) for idx, split in enumerate(splits)]

        observed = NeighborNet._build_splits_graph(circular_splits, all_taxa)
        expected = self._legacy_build_splits_graph(circular_splits, all_taxa)

        assert observed[0] == expected[0]
        assert observed[1] == expected[1]
        assert self._canonical_edges(observed[2]) == self._canonical_edges(expected[2])
        assert observed[3] == expected[3]
        assert observed[4] == expected[4]


class TestNetworkPlot:
    def test_draw_network_batches_edges(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection
        from phykit.helpers.plot_config import PlotConfig

        service = NeighborNet.__new__(NeighborNet)
        service.output_path = str(tmp_path / "neighbor_net.png")
        service.plot_config = PlotConfig(show_title=False)
        ordering = ["A", "B", "C", "D"]
        all_taxa = frozenset(ordering)
        plot_splits = [
            (frozenset({"A", "B"}), 1.0),
            (frozenset({"B", "C"}), 0.5),
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

        service._draw_network(ordering, plot_splits, all_taxa)

        assert len(line_collections) >= 2
        assert os.path.exists(service.output_path)

    def test_draw_network_batches_unlabeled_fallback_points(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from phykit.helpers.plot_config import PlotConfig

        service = NeighborNet.__new__(NeighborNet)
        service.output_path = str(tmp_path / "neighbor_net_points.png")
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

        service._draw_network(ordering, [], all_taxa)

        assert len(scatter_calls) == 1
        x, y, kwargs = scatter_calls[0]
        assert len(x) == len(y) == len(ordering)
        assert kwargs["color"] == "black"
        assert os.path.exists(service.output_path)


class TestFromDistanceMatrixCSV:
    def test_reads_csv(self, tmp_path):
        taxa, D = _make_simple_distance_matrix()
        csv_path = str(tmp_path / "dist.csv")
        _write_csv_matrix(csv_path, taxa, D)
        output = str(tmp_path / "network.png")

        args = Namespace(
            alignment=None,
            distance_matrix=csv_path,
            output=output,
            metric="p-distance",
            max_splits=30,
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
        service = NeighborNet(args)
        read_taxa, read_D = service._read_distance_matrix()
        assert read_taxa == taxa
        np.testing.assert_array_almost_equal(read_D, D)

    def test_reads_csv_without_label_cell_and_preserves_short_row_zeros(
        self, tmp_path
    ):
        csv_path = str(tmp_path / "dist_no_label.csv")
        with open(csv_path, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["A", "B", "C"])
            writer.writerow(["0.0", "0.1"])
            writer.writerow(["0.1", "0.0", "0.2"])
            writer.writerow(["0.3", "0.2", "0.0"])

        args = Namespace(
            alignment=None,
            distance_matrix=csv_path,
            output=str(tmp_path / "network.png"),
            metric="p-distance",
            max_splits=30,
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
        service = NeighborNet(args)

        taxa, D = service._read_distance_matrix()

        assert taxa == ["A", "B", "C"]
        np.testing.assert_array_equal(
            D,
            np.array(
                [
                    [0.0, 0.1, 0.0],
                    [0.1, 0.0, 0.2],
                    [0.3, 0.2, 0.0],
                ]
            ),
        )

    def test_read_distance_matrix_falls_back_for_quoted_csv(self, tmp_path):
        csv_path = str(tmp_path / "dist_quoted.csv")
        with open(csv_path, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["", "A,one", "B"])
            writer.writerow(["A,one", "0.0", "0.1"])
            writer.writerow(["B", "0.1", "0.0"])

        args = Namespace(
            alignment=None,
            distance_matrix=csv_path,
            output=str(tmp_path / "network.png"),
            metric="p-distance",
            max_splits=30,
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
        service = NeighborNet(args)

        taxa, D = service._read_distance_matrix()

        assert taxa == ["A,one", "B"]
        np.testing.assert_array_equal(D, np.array([[0.0, 0.1], [0.1, 0.0]]))

    def test_read_distance_matrix_uses_row_slice_assignments(
        self, tmp_path, monkeypatch
    ):
        csv_path = str(tmp_path / "dist.csv")
        with open(csv_path, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["", "A", "B"])
            writer.writerow(["A", "0.0", "0.1"])
            writer.writerow(["B", "0.1", "0.0"])

        class FakeMatrix:
            def __init__(self, shape):
                self.shape = shape
                self.assignments = []

            def __setitem__(self, key, value):
                row_key, col_key = key
                if not isinstance(col_key, slice):
                    raise AssertionError("distance rows should be assigned by slice")
                self.assignments.append((row_key, col_key, tuple(value)))

        fake_matrix = FakeMatrix((2, 2))

        def fake_zeros(shape):
            assert shape == (2, 2)
            return fake_matrix

        monkeypatch.setattr(neighbor_net_module.np, "zeros", fake_zeros)
        args = Namespace(
            alignment=None,
            distance_matrix=csv_path,
            output=str(tmp_path / "network.png"),
            metric="p-distance",
            max_splits=30,
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
        service = NeighborNet(args)

        taxa, D = service._read_distance_matrix()

        assert taxa == ["A", "B"]
        assert D is fake_matrix
        assert fake_matrix.assignments == [
            (0, slice(None, 2, None), (0.0, 0.1)),
            (1, slice(None, 2, None), (0.1, 0.0)),
        ]

    def test_text_output_batches_positive_split_rows(self, tmp_path, mocker):
        output = str(tmp_path / "network.png")
        args = Namespace(
            alignment=None,
            distance_matrix="dist.csv",
            output=output,
            metric="p-distance",
            max_splits=30,
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
        service = NeighborNet(args)
        taxa = ["A", "B", "C", "D"]
        splits = [frozenset({"A", "B"}), frozenset({"B", "C"})]

        mocker.patch.object(
            service,
            "_read_distance_matrix",
            return_value=(taxa, np.zeros((4, 4))),
        )
        mocker.patch.object(service, "_nj_circular_ordering", return_value=taxa)
        mocker.patch.object(
            service,
            "_enumerate_circular_splits",
            return_value=splits,
        )
        mocker.patch.object(
            service,
            "_estimate_split_weights",
            return_value=np.array([0.5, 0.25]),
        )
        mocker.patch.object(service, "_draw_network")
        printed = mocker.patch("builtins.print")

        service.run()

        printed.assert_called_once_with(
            "NeighborNet Analysis\n"
            "Taxa: 4\n"
            "Distance metric: p-distance\n"
            "Circular splits with positive weight: 2\n"
            "---\n"
            "{A, B}\t0.500000\n"
            "{B, C}\t0.250000\n"
            f"Output: {output}"
        )


class TestCreatesPNG:
    def test_creates_output_file(self, tmp_path):
        taxa, D = _make_simple_distance_matrix()
        csv_path = str(tmp_path / "dist.csv")
        _write_csv_matrix(csv_path, taxa, D)
        output = str(tmp_path / "network.png")

        args = Namespace(
            alignment=None,
            distance_matrix=csv_path,
            output=output,
            metric="p-distance",
            max_splits=30,
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
        service = NeighborNet(args)
        service.run()
        assert os.path.exists(output)


class TestFromAlignment:
    def test_alignment_run(self, tmp_path):
        output = str(tmp_path / "network.png")
        args = Namespace(
            alignment=SAMPLE_ALIGNMENT,
            distance_matrix=None,
            output=output,
            metric="p-distance",
            max_splits=30,
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
        service = NeighborNet(args)
        service.run()
        assert os.path.exists(output)


class TestJsonOutput:
    def test_json_structure(self, tmp_path, capsys):
        import json

        taxa, D = _make_simple_distance_matrix()
        csv_path = str(tmp_path / "dist.csv")
        _write_csv_matrix(csv_path, taxa, D)
        output = str(tmp_path / "network.png")

        args = Namespace(
            alignment=None,
            distance_matrix=csv_path,
            output=output,
            metric="p-distance",
            max_splits=30,
            json=True,
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
        service = NeighborNet(args)
        service.run()

        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert payload["taxa_count"] == 4
        assert "circular_ordering" in payload
        assert "splits" in payload
        assert "output" in payload
        assert set(payload["taxa"]) == {"A", "B", "C", "D"}


class TestRequiresInput:
    def test_no_input_raises_error(self, tmp_path):
        output = str(tmp_path / "network.png")
        args = Namespace(
            alignment=None,
            distance_matrix=None,
            output=output,
            metric="p-distance",
            max_splits=30,
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
        with pytest.raises(PhykitUserError):
            NeighborNet(args)
