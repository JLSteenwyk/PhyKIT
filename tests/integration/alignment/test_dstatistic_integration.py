import pytest
import sys
import json
from mock import patch
from pathlib import Path

from phykit.phykit import Phykit


def _write_alignment(path, seqs):
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


# Standard test alignment: 4 ABBA, 2 BABA, 4 invariant = 10 sites
# D = (4-2)/(4+2) = 0.3333
STANDARD_SEQS = {
    "P1":       "AAAACCAAAA",
    "P2":       "CCCCAAAAAA",
    "P3":       "CCCCCCAAAA",
    "Outgroup": "AAAAAAAAAA",
}


@pytest.mark.integration
class TestDstatisticIntegration:
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), STANDARD_SEQS)
        testargs = [
            "phykit",
            "dstatistic",
            "-a", str(aln),
            "--p1", "P1",
            "--p2", "P2",
            "--p3", "P3",
            "--outgroup", "Outgroup",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        output = " ".join(str(c) for c in mocked_print.call_args_list)
        assert "ABBA sites: 4" in output
        assert "BABA sites: 2" in output
        assert "D-statistic: 0.3333" in output

    @patch("builtins.print")
    def test_alias_dstat(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), STANDARD_SEQS)
        testargs = [
            "phykit",
            "dstat",
            "-a", str(aln),
            "--p1", "P1",
            "--p2", "P2",
            "--p3", "P3",
            "--outgroup", "Outgroup",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        output = " ".join(str(c) for c in mocked_print.call_args_list)
        assert "D-statistic: 0.3333" in output

    @patch("builtins.print")
    def test_alias_abba_baba(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), STANDARD_SEQS)
        testargs = [
            "phykit",
            "abba_baba",
            "-a", str(aln),
            "--p1", "P1",
            "--p2", "P2",
            "--p3", "P3",
            "--outgroup", "Outgroup",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        output = " ".join(str(c) for c in mocked_print.call_args_list)
        assert "D-statistic: 0.3333" in output

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), STANDARD_SEQS)
        testargs = [
            "phykit",
            "dstatistic",
            "-a", str(aln),
            "--p1", "P1",
            "--p2", "P2",
            "--p3", "P3",
            "--outgroup", "Outgroup",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        # JSON output is a single print call
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["p1"] == "P1"
        assert payload["p2"] == "P2"
        assert payload["p3"] == "P3"
        assert payload["outgroup"] == "Outgroup"
        assert payload["alignment_length"] == 10
        assert payload["informative_sites"] == 6
        assert payload["abba_count"] == 4
        assert payload["baba_count"] == 2
        assert abs(payload["d_statistic"] - 0.3333) < 0.01

    @patch("builtins.print")
    def test_custom_block_size(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), STANDARD_SEQS)
        testargs = [
            "phykit",
            "dstatistic",
            "-a", str(aln),
            "--p1", "P1",
            "--p2", "P2",
            "--p3", "P3",
            "--outgroup", "Outgroup",
            "--block-size", "5",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        output = " ".join(str(c) for c in mocked_print.call_args_list)
        assert "Block jackknife (block size: 5):" in output
        assert "Standard error:" in output

    @patch("builtins.print")
    def test_missing_taxon_error(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), STANDARD_SEQS)
        testargs = [
            "phykit",
            "dstatistic",
            "-a", str(aln),
            "--p1", "nonexistent",
            "--p2", "P2",
            "--p3", "P3",
            "--outgroup", "Outgroup",
        ]
        with pytest.raises(SystemExit):
            with patch.object(sys, "argv", testargs):
                Phykit()
