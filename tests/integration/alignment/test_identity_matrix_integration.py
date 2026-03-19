import json
import os
import sys
import pytest
from mock import patch
from pathlib import Path

from phykit.phykit import Phykit

here = Path(__file__)
sample_files = here.parent.parent.parent / "sample_files"
ALIGNMENT = str(sample_files / "12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit")


@pytest.mark.integration
class TestIdentityMatrixIntegration:
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print, tmp_path):
        """Basic end-to-end run produces text output and output file."""
        out_file = str(tmp_path / "identity_matrix.png")
        testargs = [
            "phykit",
            "identity_matrix",
            "-a", ALIGNMENT,
            "-o", out_file,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        # Check that expected text was printed
        printed = " ".join(str(c) for c in mocked_print.mock_calls)
        assert "Sequence Identity Matrix" in printed
        assert "Taxa:" in printed
        assert "Mean pairwise identity:" in printed
        assert os.path.exists(out_file)
        assert os.path.getsize(out_file) > 0

    @patch("builtins.print")
    def test_alias_id_matrix(self, mocked_print, tmp_path):
        """The id_matrix alias should work."""
        out_file = str(tmp_path / "id_matrix.png")
        testargs = [
            "phykit",
            "id_matrix",
            "-a", ALIGNMENT,
            "-o", out_file,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        printed = " ".join(str(c) for c in mocked_print.mock_calls)
        assert "Sequence Identity Matrix" in printed
        assert os.path.exists(out_file)

    @patch("builtins.print")
    def test_alias_seqid(self, mocked_print, tmp_path):
        """The seqid alias should work."""
        out_file = str(tmp_path / "seqid.png")
        testargs = [
            "phykit",
            "seqid",
            "-a", ALIGNMENT,
            "-o", out_file,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        printed = " ".join(str(c) for c in mocked_print.mock_calls)
        assert "Sequence Identity Matrix" in printed
        assert os.path.exists(out_file)

    @patch("builtins.print")
    def test_p_distance(self, mocked_print, tmp_path):
        """--metric p-distance should work and report p-distance in output."""
        out_file = str(tmp_path / "pdist.png")
        testargs = [
            "phykit",
            "identity_matrix",
            "-a", ALIGNMENT,
            "-o", out_file,
            "--metric", "p-distance",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        printed = " ".join(str(c) for c in mocked_print.mock_calls)
        assert "Metric: p-distance" in printed
        assert os.path.exists(out_file)

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        """--json should produce valid JSON with correct structure."""
        out_file = str(tmp_path / "json_out.png")
        testargs = [
            "phykit",
            "identity_matrix",
            "-a", ALIGNMENT,
            "-o", out_file,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()

        # The last print call should be JSON
        payload = json.loads(mocked_print.call_args.args[0])
        assert "n_taxa" in payload
        assert "alignment_length" in payload
        assert "metric" in payload
        assert payload["metric"] == "identity"
        assert "mean_identity" in payload
        assert "min_identity" in payload
        assert "max_identity" in payload
        assert "taxa_order" in payload
        assert "output_file" in payload
        assert payload["n_taxa"] == len(payload["taxa_order"])
        assert os.path.exists(out_file)
