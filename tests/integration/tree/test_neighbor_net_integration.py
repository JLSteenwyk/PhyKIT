import csv
import json
import os
import sys

import numpy as np
import pytest
from mock import patch

from phykit.phykit import Phykit


SAMPLE_ALIGNMENT = os.path.join(
    os.path.dirname(__file__),
    os.pardir,
    os.pardir,
    "sample_files",
    "12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit",
)


def _write_csv_matrix(path, taxa, D):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([""] + taxa)
        for i, taxon in enumerate(taxa):
            writer.writerow([taxon] + [f"{D[i, j]:.6f}" for j in range(len(taxa))])


def _simple_matrix():
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


@pytest.mark.integration
class TestNeighborNetIntegration:
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print, tmp_path):
        output = str(tmp_path / "network.png")
        testargs = [
            "phykit", "neighbor_net",
            "-a", SAMPLE_ALIGNMENT,
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert os.path.exists(output)
        full_output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "NeighborNet Analysis" in full_output
        assert "Taxa: 12" in full_output

    @patch("builtins.print")
    def test_alias_nnet(self, mocked_print, tmp_path):
        output = str(tmp_path / "network.png")
        testargs = [
            "phykit", "nnet",
            "-a", SAMPLE_ALIGNMENT,
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert os.path.exists(output)
        full_output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "NeighborNet Analysis" in full_output

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        output = str(tmp_path / "network.png")
        testargs = [
            "phykit", "neighbor_net",
            "-a", SAMPLE_ALIGNMENT,
            "-o", output,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert os.path.exists(output)
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["taxa_count"] == 12
        assert "circular_ordering" in payload
        assert "splits" in payload

    @patch("builtins.print")
    def test_from_csv_matrix(self, mocked_print, tmp_path):
        taxa, D = _simple_matrix()
        csv_path = str(tmp_path / "dist.csv")
        _write_csv_matrix(csv_path, taxa, D)
        output = str(tmp_path / "network.png")

        testargs = [
            "phykit", "neighbor_net",
            "--distance-matrix", csv_path,
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert os.path.exists(output)
        full_output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Taxa: 4" in full_output

    @patch("builtins.print")
    def test_jc_metric(self, mocked_print, tmp_path):
        output = str(tmp_path / "network.png")
        testargs = [
            "phykit", "neighbor_net",
            "-a", SAMPLE_ALIGNMENT,
            "-o", output,
            "--metric", "jc",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        assert os.path.exists(output)
        full_output = "\n".join(str(call) for call in mocked_print.call_args_list)
        assert "Distance metric: jc" in full_output
