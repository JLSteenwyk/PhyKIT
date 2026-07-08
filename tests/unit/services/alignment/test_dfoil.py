import io
import json
import subprocess
import sys

import pytest
from argparse import Namespace
from math import isclose

from phykit.services.alignment.dfoil import Dfoil, PATTERNS
import phykit.services.alignment.dfoil as dfoil_module
from phykit.errors import PhykitUserError


def test_module_import_does_not_import_biopython_fasta_parser():
    code = """
import sys
import phykit.services.alignment.dfoil as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_proxy_caches_resolved_attributes():
    lazy_np = dfoil_module._LazyNumpy()

    assert lazy_np.frombuffer is lazy_np.frombuffer
    assert lazy_np._module is not None


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
    def test_read_fasta_uses_first_header_token_uppercases_and_keeps_last_duplicate(
        self, tmp_path
    ):
        aln = tmp_path / "test.fa"
        aln.write_text(
            ">P1 description\nac gt\nn-\n"
            ">P2 second description\ncc gg\nn-\n"
            ">P1 replacement\ntttt\n"
        )
        assert Dfoil._read_fasta(str(aln)) == {
            "P1": "TTTT",
            "P2": "CCGGN-",
        }

    def test_vectorized_counting_matches_scalar_reference(self):
        """Fast pattern counting must preserve the original scalar semantics."""
        columns = [
            ("A", "A", "A", "C", "A"),  # AAABA
            ("C", "A", "A", "A", "A"),  # BAAAA
            ("C", "C", "C", "C", "A"),  # BBBBA, biallelic but uninformative
            ("A", "A", "A", "A", "A"),  # invariant, skipped by scalar code
            ("A", "C", "G", "A", "A"),  # non-biallelic
            ("C", "G", "A", "A", "A"),  # two derived alleles, skipped
            ("T", "T", "G", "T", "T"),  # AABAA with a non-A outgroup
            ("T", "G", "T", "G", "T"),  # ABABA with a non-A outgroup
            ("?", "A", "A", "C", "A"),  # skipped character
        ]
        seqs = ["".join(column[i] for column in columns) for i in range(5)]

        fast_counts = Dfoil._count_site_patterns(*seqs)
        scalar_counts = Dfoil._count_site_patterns_scalar(*seqs)

        assert fast_counts == scalar_counts
        assert fast_counts["AAABA"] == 1
        assert fast_counts["BAAAA"] == 1
        assert fast_counts["BBBBA"] == 1
        assert fast_counts["AABAA"] == 1
        assert fast_counts["ABABA"] == 1
        assert fast_counts["AAAAA"] == 0

    def test_scalar_counting_handles_unicode_and_biallelic_semantics(self):
        columns = [
            ("Ω", "A", "A", "A", "A"),  # BAAAA
            ("Ω", "Ω", "Ω", "Ω", "A"),  # BBBBA
            ("A", "A", "A", "A", "A"),  # invariant, skipped
            ("Ω", "B", "A", "A", "A"),  # non-biallelic, skipped
            ("-", "Ω", "A", "A", "A"),  # skipped character
        ]
        seqs = ["".join(column[i] for column in columns) for i in range(5)]

        counts = Dfoil._count_site_patterns_scalar(*seqs)

        assert counts["BAAAA"] == 1
        assert counts["BBBBA"] == 1
        assert counts["AAAAA"] == 0
        for pattern in PATTERNS:
            if pattern not in {"BAAAA", "BBBBA"}:
                assert counts[pattern] == 0

    def test_scalar_counting_uses_shared_skip_character_constant(self, monkeypatch):
        monkeypatch.setattr(dfoil_module, "_SCALAR_SKIP_CHARS", "")

        counts = Dfoil._count_site_patterns_scalar("-", "A", "A", "A", "A")

        assert counts["BAAAA"] == 1

    def test_all_valid_ascii_counting_skips_validity_mask(self, mocker):
        mocker.patch.object(
            dfoil_module.np,
            "ones",
            side_effect=AssertionError(
                "all-valid ASCII DFOIL counts should not build a valid mask"
            ),
        )
        seqs = _build_alignment_from_patterns(
            ["AAABA", "AABAA", "ABABA", "BABAA", "BBBBA"]
        )

        counts = Dfoil._count_site_patterns(
            seqs["P1"],
            seqs["P2"],
            seqs["P3"],
            seqs["P4"],
            seqs["Outgroup"],
        )

        assert counts["AAABA"] == 1
        assert counts["AABAA"] == 1
        assert counts["ABABA"] == 1
        assert counts["BABAA"] == 1
        assert counts["BBBBA"] == 1
        assert counts["AAAAA"] == 0
        assert list(counts) == PATTERNS

    def test_small_ascii_with_skips_uses_lookup_validity_mask(self, mocker):
        mocker.patch.object(
            dfoil_module.np,
            "ones",
            side_effect=AssertionError(
                "small DFOIL alignments with skips should use lookup validity mask"
            ),
        )
        lookup_spy = mocker.spy(dfoil_module, "_get_skip_lookup")
        seq_dict = _build_alignment_from_patterns(["AAABA", "AABAA"])
        seqs = (
            seq_dict["P1"] + "?",
            seq_dict["P2"] + "A",
            seq_dict["P3"] + "A",
            seq_dict["P4"] + "C",
            seq_dict["Outgroup"] + "A",
        )

        counts = Dfoil._count_site_patterns(*seqs)

        assert counts == Dfoil._count_site_patterns_scalar(*seqs)
        assert counts["AAABA"] == 1
        assert counts["AABAA"] == 1
        assert lookup_spy.call_count == 1

    def test_skip_code_counts_build_result_without_indexing_bincounts(
        self,
        monkeypatch,
    ):
        class IterableCounts:
            def __iter__(self):
                return iter([0, 1, 2, 3] + [0] * (len(PATTERNS) - 4))

            def __getitem__(self, _index):
                raise AssertionError("skip-code result should iterate counts")

        monkeypatch.setattr(
            dfoil_module.np,
            "bincount",
            lambda *_args, **_kwargs: IterableCounts(),
        )

        counts = Dfoil._count_site_patterns("?AAA", "AAAA", "AAAA", "CAAA", "AAAA")

        assert list(counts) == PATTERNS
        assert counts["AAABA"] == 1
        assert counts["AABAA"] == 2
        assert counts["AABBA"] == 3

    def test_large_ascii_with_skips_keeps_loop_validity_mask(self, monkeypatch):
        monkeypatch.setattr(dfoil_module, "_SKIP_LOOKUP_SMALL_ALIGNMENT_MAX", 1)
        monkeypatch.setattr(
            dfoil_module,
            "_get_skip_lookup",
            lambda: (_ for _ in ()).throw(
                AssertionError("large DFOIL alignments should not use lookup mask")
            ),
        )
        seq_dict = _build_alignment_from_patterns(["AAABA", "AABAA"])
        seqs = (
            seq_dict["P1"] + "?",
            seq_dict["P2"] + "A",
            seq_dict["P3"] + "A",
            seq_dict["P4"] + "C",
            seq_dict["Outgroup"] + "A",
        )

        counts = Dfoil._count_site_patterns(*seqs)

        assert counts == Dfoil._count_site_patterns_scalar(*seqs)
        assert counts["AAABA"] == 1
        assert counts["AABAA"] == 1

    def test_all_invariant_sequences_skip_numpy_counting(self, mocker):
        mocker.patch.object(
            dfoil_module.np,
            "frombuffer",
            side_effect=AssertionError(
                "all-invariant DFOIL counts should return before NumPy setup"
            ),
        )

        counts = Dfoil._count_site_patterns(
            "ACGT" * 10,
            "ACGT" * 10,
            "ACGT" * 10,
            "ACGT" * 10,
            "ACGT" * 10,
        )

        assert counts == {pattern: 0 for pattern in PATTERNS}
        counts["AAAAA"] = 1
        fresh_counts = Dfoil._count_site_patterns(
            "ACGT" * 10,
            "ACGT" * 10,
            "ACGT" * 10,
            "ACGT" * 10,
            "ACGT" * 10,
        )
        assert fresh_counts["AAAAA"] == 0

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

    def test_chi_square_p_values_do_not_import_scipy(self, tmp_path, monkeypatch):
        patterns = ['AAABA'] * 50 + ['AABAA'] * 5
        seqs = _build_alignment_from_patterns(patterns)
        aln = tmp_path / "test.fa"
        _write_alignment(str(aln), seqs)
        args = _make_args(str(aln), json_output=True)
        svc = Dfoil(args)
        real_import = __import__

        def fail_scipy_import(name, *args, **kwargs):
            if name == "scipy" or name.startswith("scipy."):
                raise AssertionError("DFOIL chi-square p-values should not import scipy")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr("builtins.__import__", fail_scipy_import)

        output = _run_and_capture(svc)
        payload = json.loads(output)

        assert payload["dfo"]["p_value"] < 0.05

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

    def test_print_text_output_batches_report(self, monkeypatch):
        """Text report output should be emitted in a single print call."""
        svc = Dfoil.__new__(Dfoil)
        svc.p1 = "P1"
        svc.p2 = "P2"
        svc.p3 = "P3"
        svc.p4 = "P4"
        svc.outgroup = "Outgroup"
        counts = {pattern: i for i, pattern in enumerate(PATTERNS)}
        printed = []

        def fake_print(*args, **kwargs):
            printed.append((args, kwargs))

        monkeypatch.setattr("builtins.print", fake_print)

        svc._print_text_output(
            1000,
            120,
            counts,
            0.1234,
            -0.2345,
            0.0,
            0.6789,
            0.0004,
            0.007,
            0.04,
            0.5,
            "+-00",
            "Ambiguous or complex introgression pattern",
        )

        expected = "\n".join([
            "DFOIL Test (Pease & Hahn 2015)",
            "================================",
            "Topology: ((P1, P2), (P3, P4), Outgroup)",
            "P1: P1, P2: P2, P3: P3, P4: P4, Outgroup: Outgroup",
            "",
            "Alignment length: 1000",
            "Informative sites: 120",
            "",
            "Site pattern counts:",
            "  AAABA: 1  AABAA: 2  AABBA: 3  ABAAA: 4",
            "  ABABA: 5  ABBAA: 6  ABBBA: 7  BAAAA: 8",
            "  BAABA: 9  BABAA: 10  BABBA: 11  BBAAA: 12",
            "  BBABA: 13  BBBAA: 14",
            "",
            "D-statistics:",
            "  DFO:  0.1234  (p = 0.000400 ***)",
            "  DIL:  -0.2345  (p = 0.007000 **)",
            "  DFI:  0.0000  (p = 0.040000 *)",
            "  DOL:  0.6789  (p = 0.500000)",
            "",
            "Sign pattern: +-00",
            "Interpretation: Ambiguous or complex introgression pattern",
        ])
        assert printed == [((expected,), {})]
