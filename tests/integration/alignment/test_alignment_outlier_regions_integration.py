import json
import sys
from pathlib import Path
from unittest.mock import patch

from Bio import SeqIO
import pytest

from phykit.phykit import Phykit


def _write_local_error_alignment(path: Path) -> Path:
    lines = []
    for index in range(9):
        lines.extend([f">normal_{index}", "A" * 120])
    lines.extend(
        [
            ">local_error",
            ("A" * 50) + ("C" * 10) + ("A" * 60),
        ]
    )
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.mark.integration
class TestAlignmentOutlierRegions:
    @patch("builtins.print")
    def test_tsv_output_and_alias(self, mocked_print, tmp_path):
        alignment = _write_local_error_alignment(tmp_path / "local_error.fa")
        testargs = ["phykit", "aor", str(alignment)]

        with patch.object(sys, "argv", testargs):
            Phykit()

        lines = mocked_print.call_args.args[0].splitlines()
        assert lines[0].startswith("taxon\talignment_start")
        assert lines[1].startswith("local_error\t49\t58\t49\t58\t10")

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        alignment = _write_local_error_alignment(tmp_path / "local_error.fa")
        testargs = [
            "phykit",
            "alignment_outlier_regions",
            str(alignment),
            "--json",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["doi"] == "10.1111/2041-210X.13696"
        assert payload["summary"] == {
            "affected_taxa": 1,
            "masked_residues": 10,
            "regions": 1,
        }
        assert payload["regions"][0]["taxon"] == "local_error"

    @patch("builtins.print")
    def test_report_and_masked_alignment_outputs(self, mocked_print, tmp_path):
        alignment = _write_local_error_alignment(tmp_path / "local_error.fa")
        report = tmp_path / "regions.tsv"
        masked = tmp_path / "masked.fa"
        testargs = [
            "phykit",
            "outlier_regions",
            str(alignment),
            "--report",
            str(report),
            "--mask-output",
            str(masked),
            "--mask-character",
            "N",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        assert report.read_text().splitlines()[1].startswith("local_error\t49\t58")
        records = list(SeqIO.parse(masked, "fasta"))
        assert len(records) == 10
        assert str(records[-1].seq).count("N") == 10
        printed = [call.args[0] for call in mocked_print.call_args_list]
        assert f"report\t{report}" in printed
        assert f"masked_alignment\t{masked}" in printed

    def test_help_cites_taper(self, capsys):
        with pytest.raises(SystemExit) as error:
            Phykit.alignment_outlier_regions(["--help"])

        assert error.value.code == 0
        output = capsys.readouterr().out
        assert "TAPER" in output
        assert "doi:10.1111/2041-210X.13696" in output

