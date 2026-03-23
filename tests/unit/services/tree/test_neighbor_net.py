import csv
import os
from argparse import Namespace

import numpy as np
import pytest

from phykit.errors import PhykitUserError
from phykit.services.tree.neighbor_net import NeighborNet


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
    def test_contains_all_taxa(self):
        taxa, D = _make_simple_distance_matrix()
        ordering = NeighborNet._nj_circular_ordering(D, taxa)
        assert set(ordering) == set(taxa)

    def test_ordering_length(self):
        taxa, D = _make_simple_distance_matrix()
        ordering = NeighborNet._nj_circular_ordering(D, taxa)
        assert len(ordering) == len(taxa)


class TestCircularSplits:
    def test_splits_are_contiguous(self):
        ordering = ["A", "B", "C", "D", "E"]
        splits = NeighborNet._enumerate_circular_splits(ordering)
        n = len(ordering)
        for split in splits:
            # Verify the split forms a contiguous arc
            assert NeighborNet._is_circular_split(split, ordering)

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


class TestNNLSWeights:
    def test_weights_nonneg(self):
        taxa, D = _make_simple_distance_matrix()
        ordering = NeighborNet._nj_circular_ordering(D, taxa)
        splits = NeighborNet._enumerate_circular_splits(ordering)
        weights = NeighborNet._estimate_split_weights(D, splits, ordering, taxa)
        assert np.all(weights >= -1e-10)


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
