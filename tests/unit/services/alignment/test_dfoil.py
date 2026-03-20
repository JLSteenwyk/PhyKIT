import io
import json
import sys

import pytest
from argparse import Namespace
from math import isclose

from phykit.services.alignment.dfoil import Dfoil
from phykit.errors import PhykitUserError


def _write_alignment(path, seqs):
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


def _make_args(alignment, p1="P1", p2="P2", p3="P3", p4="P4",
               outgroup="Outgroup", json_output=False):
    return Namespace(
        alignment=alignment,
        p1=p1,
        p2=p2,
        p3=p3,
        p4=p4,
        outgroup=outgroup,
        json=json_output,
    )


def _build_alignment_from_patterns(patterns):
    """Build a FASTA alignment from a list of 5-character pattern strings.

    Each pattern has positions: P1, P2, P3, P4, O.
    A = matches outgroup (use nucleotide 'A'), B = differs (use nucleotide 'C').
    Outgroup is always 'A'.
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


def _run_and_capture(svc):
    captured = io.StringIO()
    old_stdout = sys.stdout
    sys.stdout = captured
    try:
        svc.run()
    finally:
        sys.stdout = old_stdout
    return captured.getvalue()


class TestDfoilPatternCounting:
    def test_pattern_counting(self, tmp_path):
        """Verify correct pattern counts for a known alignment."""
        # Create an alignment with specific patterns
        patterns = [
            'AAABA', 'AAABA',  # 2x AAABA
            'AABAA',           # 1x AABAA
            'ABABA',           # 1x ABABA
            'BABAA',           # 1x BABAA
            'BAAAA',           # 1x BAAAA
        ]
        seqs = _build_alignment_from_patterns(patterns)
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), json_output=True)
        svc = Dfoil(args)
        output = _run_and_capture(svc)
        payload = json.loads(output)

        assert payload["pattern_counts"]["AAABA"] == 2
        assert payload["pattern_counts"]["AABAA"] == 1
        assert payload["pattern_counts"]["ABABA"] == 1
        assert payload["pattern_counts"]["BABAA"] == 1
        assert payload["pattern_counts"]["BAAAA"] == 1
        assert payload["informative_sites"] == 6

    def test_all_invariant_zero(self, tmp_path):
        """All same alleles should yield D = 0 for all statistics."""
        aln = tmp_path / "test.fa"
        seqs = {
            "P1":       "AAAAAAAAAA",
            "P2":       "AAAAAAAAAA",
            "P3":       "AAAAAAAAAA",
            "P4":       "AAAAAAAAAA",
            "Outgroup": "AAAAAAAAAA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), json_output=True)
        svc = Dfoil(args)
        output = _run_and_capture(svc)
        payload = json.loads(output)

        assert payload["informative_sites"] == 0
        assert payload["dfo"]["value"] == 0.0
        assert payload["dil"]["value"] == 0.0
        assert payload["dfi"]["value"] == 0.0
        assert payload["dol"]["value"] == 0.0
        assert payload["sign_pattern"] == "0000"
        assert payload["interpretation"] == "No significant introgression detected"

    def test_gaps_skipped(self, tmp_path):
        """Sites with gaps should be excluded from counting."""
        patterns = ['AAABA', 'AABAA']
        seqs = _build_alignment_from_patterns(patterns)
        # Insert a gap at site 0 for P1
        seqs["P1"] = "-" + seqs["P1"][1:]
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), json_output=True)
        svc = Dfoil(args)
        output = _run_and_capture(svc)
        payload = json.loads(output)

        # Only site 1 (AABAA) should be counted
        assert payload["informative_sites"] == 1
        assert payload["pattern_counts"]["AABAA"] == 1

    def test_non_biallelic_skipped(self, tmp_path):
        """Sites with 3+ alleles should be excluded."""
        aln = tmp_path / "test.fa"
        # Site 0: P1=A, P2=C, P3=G, P4=A, O=A -> 3 alleles -> skip
        # Site 1: AAABA pattern -> count
        seqs = {
            "P1":       "AA",
            "P2":       "CA",
            "P3":       "GA",
            "P4":       "AC",
            "Outgroup": "AA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), json_output=True)
        svc = Dfoil(args)
        output = _run_and_capture(svc)
        payload = json.loads(output)

        assert payload["informative_sites"] == 1
        assert payload["pattern_counts"]["AAABA"] == 1

    def test_missing_taxon_raises(self, tmp_path):
        """Error if any of the 5 taxa not found in alignment."""
        aln = tmp_path / "test.fa"
        seqs = {
            "P1":       "AAAAAAAAAA",
            "P2":       "AAAAAAAAAA",
            "P3":       "AAAAAAAAAA",
            "P4":       "AAAAAAAAAA",
            "Outgroup": "AAAAAAAAAA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), p1="missing_taxon")
        svc = Dfoil(args)
        with pytest.raises(PhykitUserError):
            svc.run()


class TestDfoilComputation:
    def test_dfo_computation(self, tmp_path):
        """Verify DFO formula with known pattern counts.

        DFO_left  = AAABA + ABABA + BABAA + BBBAA
        DFO_right = AABAA + ABBAA + BAABA + BBABA
        DFO = (left - right) / (left + right)
        """
        # Create patterns: 3 on left side, 1 on right side
        patterns = [
            'AAABA',  # DFO left
            'ABABA',  # DFO left
            'BABAA',  # DFO left
            'AABAA',  # DFO right
        ]
        seqs = _build_alignment_from_patterns(patterns)
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), json_output=True)
        svc = Dfoil(args)
        output = _run_and_capture(svc)
        payload = json.loads(output)

        assert payload["dfo"]["left"] == 3
        assert payload["dfo"]["right"] == 1
        # DFO = (3-1)/(3+1) = 0.5
        assert isclose(payload["dfo"]["value"], 0.5, rel_tol=0.01)

    def test_sign_pattern_positive(self, tmp_path):
        """Verify sign classification for a significantly positive D."""
        # Create a large number of DFO-left patterns to ensure significance
        patterns = ['AAABA'] * 50 + ['AABAA'] * 5
        seqs = _build_alignment_from_patterns(patterns)
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), json_output=True)
        svc = Dfoil(args)
        output = _run_and_capture(svc)
        payload = json.loads(output)

        # DFO should be significantly positive
        assert payload["dfo"]["value"] > 0
        assert payload["dfo"]["p_value"] < 0.05
        # First character of sign pattern should be '+'
        assert payload["sign_pattern"][0] == '+'

    def test_sign_pattern_zero_when_equal(self, tmp_path):
        """Equal left/right counts should produce sign '0'."""
        # Equal DFO-left and DFO-right
        patterns = ['AAABA'] * 10 + ['AABAA'] * 10
        seqs = _build_alignment_from_patterns(patterns)
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), json_output=True)
        svc = Dfoil(args)
        output = _run_and_capture(svc)
        payload = json.loads(output)

        assert isclose(payload["dfo"]["value"], 0.0, abs_tol=0.001)
        # First character should be '0'
        assert payload["sign_pattern"][0] == '0'

    def test_interpretation_lookup_no_introgression(self, tmp_path):
        """Verify that '0000' sign pattern produces no introgression message."""
        # All invariant -> all D = 0 -> sign pattern '0000'
        aln = tmp_path / "test.fa"
        seqs = {
            "P1":       "AAAAAAAAAA",
            "P2":       "AAAAAAAAAA",
            "P3":       "AAAAAAAAAA",
            "P4":       "AAAAAAAAAA",
            "Outgroup": "AAAAAAAAAA",
        }
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), json_output=True)
        svc = Dfoil(args)
        output = _run_and_capture(svc)
        payload = json.loads(output)

        assert payload["sign_pattern"] == "0000"
        assert payload["interpretation"] == "No significant introgression detected"

    def test_interpretation_lookup_known_pattern(self, tmp_path):
        """Verify that '++00' sign pattern maps to ancestor(P1,P2) <-> P3."""
        # To get '++00' we need DFO > 0 (sig), DIL > 0 (sig), DFI ~ 0, DOL ~ 0
        # DFO_left  = AAABA + ABABA + BABAA + BBBAA
        # DFO_right = AABAA + ABBAA + BAABA + BBABA
        # DIL_left  = AAABA + ABBAA + BAABA + BBBAA
        # DIL_right = AABAA + ABABA + BABAA + BBABA
        # DFI_left  = ABAAA + ABABA + BABAA + BABBA
        # DFI_right = BAAAA + ABBAA + BAABA + ABBBA
        # DOL_left  = ABAAA + ABBAA + BAABA + BABBA
        # DOL_right = BAAAA + ABABA + BABAA + ABBBA
        #
        # For DFO positive: AAABA heavy
        # For DIL positive: AAABA heavy (also in DIL_left)
        # For DFI ~ 0: need DFI_left ~ DFI_right
        # For DOL ~ 0: need DOL_left ~ DOL_right
        #
        # AAABA contributes to DFO_left and DIL_left only
        # We need enough of it to make DFO and DIL significant
        patterns = ['AAABA'] * 100
        seqs = _build_alignment_from_patterns(patterns)
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), json_output=True)
        svc = Dfoil(args)
        output = _run_and_capture(svc)
        payload = json.loads(output)

        assert payload["sign_pattern"] == "++00"
        assert "ancestor of (P1,P2) <-> P3" in payload["interpretation"]


class TestDfoilOutput:
    def test_json_output_has_all_fields(self, tmp_path):
        """JSON output must contain all expected fields."""
        patterns = ['AAABA', 'AABAA', 'ABABA', 'BAAAA']
        seqs = _build_alignment_from_patterns(patterns)
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), json_output=True)
        svc = Dfoil(args)
        output = _run_and_capture(svc)
        payload = json.loads(output)

        # Top-level fields
        assert payload["p1"] == "P1"
        assert payload["p2"] == "P2"
        assert payload["p3"] == "P3"
        assert payload["p4"] == "P4"
        assert payload["outgroup"] == "Outgroup"
        assert "alignment_length" in payload
        assert "informative_sites" in payload
        assert "pattern_counts" in payload

        # D-statistic sub-objects
        for d_name in ["dfo", "dil", "dfi", "dol"]:
            assert d_name in payload
            d_obj = payload[d_name]
            assert "value" in d_obj
            assert "left" in d_obj
            assert "right" in d_obj
            assert "chi2" in d_obj
            assert "p_value" in d_obj

        assert "sign_pattern" in payload
        assert "interpretation" in payload

    def test_text_output_structure(self, tmp_path):
        """Text output should contain key sections."""
        patterns = ['AAABA', 'AABAA', 'ABABA', 'BAAAA']
        seqs = _build_alignment_from_patterns(patterns)
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), json_output=False)
        svc = Dfoil(args)
        output = _run_and_capture(svc)

        assert "DFOIL Test (Pease & Hahn 2015)" in output
        assert "Topology:" in output
        assert "Alignment length:" in output
        assert "Informative sites:" in output
        assert "Site pattern counts:" in output
        assert "D-statistics:" in output
        assert "DFO:" in output
        assert "DIL:" in output
        assert "DFI:" in output
        assert "DOL:" in output
        assert "Sign pattern:" in output
        assert "Interpretation:" in output
