import json
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

from phykit.phykit import Phykit


SAMPLE_DIR = Path(__file__).parents[2] / "sample_files"
REFERENCE_ALIGNMENT = SAMPLE_DIR / "codon_dnds_reference.fa"


@pytest.mark.integration
class TestCodonDnDs:
    def test_canonical_command_verbose(self, capsys):
        testargs = [
            "phykit",
            "codon_dnds",
            str(REFERENCE_ALIGNMENT),
            "--verbose",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        lines = capsys.readouterr().out.splitlines()
        assert lines[0].startswith("taxon_a\ttaxon_b\tdN\tdS\tomega")
        assert lines[1].startswith("rna1\trna2\t0.067715135\t0.20119799")
        assert lines[1].endswith("\t23\t0\tok")

    def test_kaks_alias_json(self, capsys):
        testargs = [
            "phykit",
            "kaks",
            str(REFERENCE_ALIGNMENT),
            "--method",
            "LWL85",
            "--json",
        ]

        with patch.object(sys, "argv", testargs):
            Phykit()

        payload = json.loads(capsys.readouterr().out)
        assert payload["method"] == "LWL85"
        assert payload["summary"]["pairs_total"] == 1
        assert payload["summary"]["mean_dN"] == pytest.approx(0.0687284919)

    def test_invalid_frame_reports_user_error(self, tmp_path, capsys):
        alignment = tmp_path / "invalid_frame.fa"
        alignment.write_text(">a\nATGG\n>b\nATGG\n")
        testargs = ["phykit", "dnds", str(alignment)]

        with pytest.raises(SystemExit) as error:
            with patch.object(sys, "argv", testargs):
                Phykit()

        assert error.value.code == 2
        assert "not divisible by 3" in capsys.readouterr().out
