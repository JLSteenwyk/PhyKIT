import pytest
import sys
import json
from mock import patch

from phykit.phykit import Phykit


def _write_alignment(path, seqs):
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


def _build_alignment_from_patterns(patterns):
    """Build sequences from pattern strings (P1, P2, P3, P4, O positions).
    A -> nucleotide 'A' (ancestral), B -> nucleotide 'C' (derived).
    """
    p1_seq = ''.join('A' if p[0] == 'A' else 'C' for p in patterns)
    p2_seq = ''.join('A' if p[1] == 'A' else 'C' for p in patterns)
    p3_seq = ''.join('A' if p[2] == 'A' else 'C' for p in patterns)
    p4_seq = ''.join('A' if p[3] == 'A' else 'C' for p in patterns)
    o_seq = ''.join('A' if p[4] == 'A' else 'C' for p in patterns)
    return {
        "P1": p1_seq,
        "P2": p2_seq,
        "P3": p3_seq,
        "P4": p4_seq,
        "Outgroup": o_seq,
    }


# Standard test alignment
STANDARD_PATTERNS = ['AAABA'] * 3 + ['AABAA'] * 2 + ['ABABA'] * 1
STANDARD_SEQS = _build_alignment_from_patterns(STANDARD_PATTERNS)


@pytest.mark.integration
class TestDfoilIntegration:
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), STANDARD_SEQS)
        testargs = [
            "phykit",
            "dfoil",
            "-a", str(aln),
            "--p1", "P1",
            "--p2", "P2",
            "--p3", "P3",
            "--p4", "P4",
            "--outgroup", "Outgroup",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        output = " ".join(str(c) for c in mocked_print.call_args_list)
        assert "DFOIL Test" in output
        assert "DFO:" in output
        assert "DIL:" in output
        assert "DFI:" in output
        assert "DOL:" in output
        assert "Sign pattern:" in output

    @patch("builtins.print")
    def test_alias_dfoil_test(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), STANDARD_SEQS)
        testargs = [
            "phykit",
            "dfoil_test",
            "-a", str(aln),
            "--p1", "P1",
            "--p2", "P2",
            "--p3", "P3",
            "--p4", "P4",
            "--outgroup", "Outgroup",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        output = " ".join(str(c) for c in mocked_print.call_args_list)
        assert "DFOIL Test" in output

    @patch("builtins.print")
    def test_json_output(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), STANDARD_SEQS)
        testargs = [
            "phykit",
            "dfoil",
            "-a", str(aln),
            "--p1", "P1",
            "--p2", "P2",
            "--p3", "P3",
            "--p4", "P4",
            "--outgroup", "Outgroup",
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["p1"] == "P1"
        assert payload["p2"] == "P2"
        assert payload["p3"] == "P3"
        assert payload["p4"] == "P4"
        assert payload["outgroup"] == "Outgroup"
        assert "alignment_length" in payload
        assert "informative_sites" in payload
        assert "pattern_counts" in payload
        assert "dfo" in payload
        assert "dil" in payload
        assert "dfi" in payload
        assert "dol" in payload
        assert "sign_pattern" in payload
        assert "interpretation" in payload

    @patch("builtins.print")
    def test_missing_taxon_error(self, mocked_print, tmp_path):
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), STANDARD_SEQS)
        testargs = [
            "phykit",
            "dfoil",
            "-a", str(aln),
            "--p1", "nonexistent",
            "--p2", "P2",
            "--p3", "P3",
            "--p4", "P4",
            "--outgroup", "Outgroup",
        ]
        with pytest.raises(SystemExit):
            with patch.object(sys, "argv", testargs):
                Phykit()
