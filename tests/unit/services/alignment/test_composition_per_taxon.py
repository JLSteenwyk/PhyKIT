import pytest
import subprocess
import sys
from argparse import Namespace
from math import isclose
from types import SimpleNamespace
import numpy as np
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.composition_per_taxon import CompositionPerTaxon
import phykit.services.alignment.composition_per_taxon as composition_per_taxon_module


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestCompositionPerTaxon(object):
    def test_composition_mask_and_sop_modules_defer_heavy_imports(self):
        modules = [
            "phykit.services.alignment.composition_per_taxon",
            "phykit.services.alignment.mask_alignment",
            "phykit.services.alignment.sum_of_pairs_score",
        ]
        code = f"""
import sys
modules = {modules!r}
for module_name in modules:
    module = __import__(module_name, fromlist=["*"])
    assert hasattr(module.np, "__getattr__")
    assert callable(module.print_json)

assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_init_sets_alignment_file_path(self, args):
        composition = CompositionPerTaxon(args)
        assert composition.alignment_file_path == args.alignment

    def test_composition_per_taxon(self, alignment_simple, args):
        composition = CompositionPerTaxon(args)
        symbols, rows = composition.calculate_composition_per_taxon(
            alignment_simple, is_protein=False
        )

        assert symbols == ["A", "C", "G", "T"]
        by_taxon = {taxon: vals for taxon, vals in rows}

        assert isclose(by_taxon["1"][0], 0.4, rel_tol=0.001)  # A
        assert isclose(by_taxon["1"][1], 0.0, rel_tol=0.001)  # C
        assert isclose(by_taxon["1"][2], 0.2, rel_tol=0.001)  # G
        assert isclose(by_taxon["1"][3], 0.4, rel_tol=0.001)  # T

        assert isclose(by_taxon["4"][0], 0.6, rel_tol=0.001)  # A
        assert isclose(by_taxon["4"][1], 0.0, rel_tol=0.001)  # C
        assert isclose(by_taxon["4"][2], 0.2, rel_tol=0.001)  # G
        assert isclose(by_taxon["4"][3], 0.2, rel_tol=0.001)  # T

    def test_process_args_defaults_json_false(self):
        parsed = CompositionPerTaxon(Namespace(alignment="x.fa")).process_args(
            Namespace(alignment="x.fa")
        )
        assert parsed["json_output"] is False

    def test_calculate_composition_per_taxon_all_invalid(self, args):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("----"), id="t1"),
                SeqRecord(Seq("NNNN"), id="t2"),
            ]
        )
        symbols, rows = svc.calculate_composition_per_taxon(alignment, is_protein=False)
        assert symbols == []
        assert rows == []

    def test_calculate_composition_per_taxon_single_valid_symbol_shortcut(
        self, args, mocker
    ):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AAAA"), id="t1"),
                SeqRecord(Seq("A-A?"), id="t2"),
                SeqRecord(Seq("----"), id="t3"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.composition_per_taxon.np.array",
            side_effect=AssertionError("single-symbol composition should skip counts"),
        )

        symbols, rows = svc.calculate_composition_per_taxon(
            alignment,
            is_protein=False,
        )
        by_taxon = {taxon: vals for taxon, vals in rows}

        assert symbols == ["A"]
        assert np.array_equal(by_taxon["t1"], np.array([1.0]))
        assert np.array_equal(by_taxon["t2"], np.array([1.0]))
        assert np.array_equal(by_taxon["t3"], np.array([0.0]))

    def test_calculate_composition_per_taxon_ignores_ambiguous_sites_and_uppercases(self, args):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("aCGN"), id="t1"),
                SeqRecord(Seq("a--t"), id="t2"),
            ]
        )

        symbols, rows = svc.calculate_composition_per_taxon(alignment, is_protein=False)
        by_taxon = {taxon: vals for taxon, vals in rows}

        assert symbols == ["A", "C", "G", "T"]
        assert np.allclose(by_taxon["t1"], np.array([1 / 3, 1 / 3, 1 / 3, 0.0]))
        assert np.allclose(by_taxon["t2"], np.array([0.5, 0.0, 0.0, 0.5]))

    def test_calculate_composition_per_taxon_matrix_path_ascii(self, args):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("aaCC"), id="t1"),
                SeqRecord(Seq("GGtt"), id="t2"),
            ]
        )

        symbols, rows = svc.calculate_composition_per_taxon(alignment, is_protein=False)
        by_taxon = {taxon: vals for taxon, vals in rows}

        assert symbols == ["A", "C", "G", "T"]
        assert np.allclose(by_taxon["t1"], np.array([0.5, 0.5, 0.0, 0.0]))
        assert np.allclose(by_taxon["t2"], np.array([0.0, 0.0, 0.5, 0.5]))

    def test_calculate_composition_per_taxon_ascii_small_alphabet_skips_bincount(self, args, mocker):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AACCGGTT"), id="t1"),
                SeqRecord(Seq("AA--NN??"), id="t2"),
            ]
        )
        bincount_spy = mocker.spy(composition_per_taxon_module.np, "bincount")

        symbols, rows = svc.calculate_composition_per_taxon(alignment, is_protein=False)
        by_taxon = {taxon: vals for taxon, vals in rows}

        bincount_spy.assert_not_called()
        assert symbols == ["A", "C", "G", "T"]
        assert np.allclose(by_taxon["t1"], np.array([0.25, 0.25, 0.25, 0.25]))
        assert np.allclose(by_taxon["t2"], np.array([1.0, 0.0, 0.0, 0.0]))

    def test_calculate_composition_per_taxon_ascii_large_alphabet_uses_bincount(self, args, mocker):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWY"), id="t1"),
                SeqRecord(Seq("YWVTSRQPNMLKIHGFEDCA"), id="t2"),
            ]
        )
        bincount_spy = mocker.spy(composition_per_taxon_module.np, "bincount")

        symbols, rows = svc.calculate_composition_per_taxon(alignment, is_protein=True)
        by_taxon = {taxon: vals for taxon, vals in rows}

        assert bincount_spy.call_count == len(alignment)
        assert symbols == list("ACDEFGHIKLMNPQRSTVWY")
        assert np.allclose(by_taxon["t1"], np.full(20, 0.05))

    def test_calculate_composition_per_taxon_ascii_no_gap_uses_full_lengths(
        self, args, mocker
    ):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWY"), id="t1"),
                SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWA"), id="t2"),
            ]
        )
        full_spy = mocker.spy(composition_per_taxon_module.np, "full")

        symbols, rows = svc.calculate_composition_per_taxon(
            alignment,
            is_protein=True,
        )
        by_taxon = {taxon: vals for taxon, vals in rows}

        full_spy.assert_called_once_with(2, 20, dtype=composition_per_taxon_module.np.intp)
        assert symbols == list("ACDEFGHIKLMNPQRSTVWY")
        assert np.allclose(by_taxon["t1"], np.full(20, 0.05))

    def test_calculate_composition_per_taxon_large_short_protein_uses_symbol_counts(
        self, args, mocker
    ):
        svc = CompositionPerTaxon(args)
        protein = "ACDEFGHIKLMNPQRSTVWY"
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq(protein if idx % 2 == 0 else protein[::-1]), id=f"t{idx}")
                for idx in range(1024)
            ]
        )
        bincount_spy = mocker.spy(composition_per_taxon_module.np, "bincount")

        symbols, rows = svc.calculate_composition_per_taxon(
            alignment,
            is_protein=True,
        )
        by_taxon = {taxon: vals for taxon, vals in rows}

        bincount_spy.assert_not_called()
        assert symbols == list(protein)
        assert np.allclose(by_taxon["t0"], np.full(20, 0.05))

    def test_calculate_composition_per_taxon_identical_sequences_skip_matrix_counts(
        self, args, mocker
    ):
        svc = CompositionPerTaxon(args)
        protein = "ACDEFGHIKLMNPQRSTVWY"
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq(protein.lower()), id="t1"),
                SeqRecord(Seq(protein), id="t2"),
                SeqRecord(Seq(protein.lower()), id="t3"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.composition_per_taxon.np.bincount",
            side_effect=AssertionError(
                "identical sequences should not build per-row histograms"
            ),
        )

        symbols, rows = svc.calculate_composition_per_taxon(
            alignment,
            is_protein=True,
        )
        by_taxon = {taxon: vals for taxon, vals in rows}

        assert symbols == list(protein)
        for taxon in ("t1", "t2", "t3"):
            assert np.allclose(by_taxon[taxon], np.full(20, 0.05))

    def test_identical_composition_rows_return_independent_arrays(self, args):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="t1"),
                SeqRecord(Seq("ACGT"), id="t2"),
            ]
        )

        symbols, rows = svc.calculate_composition_per_taxon(
            alignment,
            is_protein=False,
        )
        rows[0][1][0] = 0.0

        assert symbols == ["A", "C", "G", "T"]
        assert np.allclose(rows[1][1], np.full(4, 0.25))

    def test_ascii_path_uses_invalid_lookup_instead_of_isin(self, args, mocker):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AA--NN??"), id="t1"),
                SeqRecord(Seq("CCGGTTAA"), id="t2"),
            ]
        )
        mocked_isin = mocker.patch.object(composition_per_taxon_module.np, "isin")

        symbols, rows = svc.calculate_composition_per_taxon(alignment, is_protein=False)
        by_taxon = {taxon: vals for taxon, vals in rows}

        mocked_isin.assert_not_called()
        assert symbols == ["A", "C", "G", "T"]
        assert np.allclose(by_taxon["t1"], np.array([1.0, 0.0, 0.0, 0.0]))
        assert np.allclose(by_taxon["t2"], np.array([0.25, 0.25, 0.25, 0.25]))

    def test_invalid_lookup_for_chars_is_case_insensitive(self, args):
        svc = CompositionPerTaxon(args)

        lookup = svc._invalid_lookup_for_chars({"n", "?"})

        assert lookup[ord("N")]
        assert lookup[ord("?")]
        assert not lookup[ord("A")]

    def test_calculate_composition_per_taxon_unicode_fallback(self, args):
        svc = CompositionPerTaxon(args)
        alignment = [
            SimpleNamespace(seq="A\u00d1CC", id="t1"),
            SimpleNamespace(seq="A\u00d1GG", id="t2"),
        ]

        symbols, rows = svc.calculate_composition_per_taxon(alignment, is_protein=False)
        by_taxon = {taxon: vals for taxon, vals in rows}

        assert symbols == ["A", "C", "G", "\u00d1"]
        assert np.allclose(by_taxon["t1"], np.array([0.25, 0.5, 0.0, 0.25]))
        assert np.allclose(by_taxon["t2"], np.array([0.25, 0.0, 0.5, 0.25]))

    def test_calculate_composition_per_taxon_protein(self, args):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACDX"), id="t1"),
                SeqRecord(Seq("AC?A"), id="t2"),
            ]
        )
        symbols, rows = svc.calculate_composition_per_taxon(alignment, is_protein=True)
        assert symbols == ["A", "C", "D"]
        assert len(rows) == 2

    def test_run_json_output(self, mocker):
        svc = CompositionPerTaxon(Namespace(alignment="x.fa", json=True))
        mocker.patch.object(
            CompositionPerTaxon,
            "get_alignment_and_format",
            return_value=(object(), "fasta", False),
        )
        mocker.patch.object(
            CompositionPerTaxon,
            "calculate_composition_per_taxon",
            return_value=(["A", "C"], [("t1", [0.5, 0.5])]),
        )
        mocked_json = mocker.patch.object(composition_per_taxon_module, "print_json")
        svc.run()
        payload = mocked_json.call_args.args[0]
        assert payload["symbols"] == ["A", "C"]
        assert payload["rows"] == payload["taxa"]
        assert payload["rows"][0]["composition"] == {"A": 0.5, "C": 0.5}

    def test_run_text_output(self, mocker, capsys):
        svc = CompositionPerTaxon(Namespace(alignment="x.fa", json=False))
        mocker.patch.object(
            CompositionPerTaxon,
            "get_alignment_and_format",
            return_value=(object(), "fasta", False),
        )
        mocker.patch.object(
            CompositionPerTaxon,
            "calculate_composition_per_taxon",
            return_value=(["A", "C"], [("t1", [0.5, 0.5]), ("t2", [1.0, 0.0])]),
        )
        svc.run()
        out, _ = capsys.readouterr()
        assert "t1\tA:0.5;C:0.5" in out
        assert "t2\tA:1.0;C:0.0" in out

    def test_run_text_output_formats_numpy_rows_exactly(self, mocker, capsys):
        svc = CompositionPerTaxon(Namespace(alignment="x.fa", json=False))
        mocker.patch.object(
            CompositionPerTaxon,
            "get_alignment_and_format",
            return_value=(object(), "fasta", False),
        )
        mocker.patch.object(
            CompositionPerTaxon,
            "calculate_composition_per_taxon",
            return_value=(
                ["A", "C", "G", "T"],
                [
                    ("t1", np.array([0.25, 0.25, 0.25, 0.25])),
                    ("t2", np.array([1.0, 0.0, 0.0, 0.0])),
                ],
            ),
        )

        svc.run()

        out, _ = capsys.readouterr()
        assert out == (
            "t1\tA:0.25;C:0.25;G:0.25;T:0.25\n"
            "t2\tA:1.0;C:0.0;G:0.0;T:0.0\n"
        )

    def test_run_returns_early_when_no_symbols(self, mocker, capsys):
        svc = CompositionPerTaxon(Namespace(alignment="x.fa", json=False))
        mocker.patch.object(
            CompositionPerTaxon,
            "get_alignment_and_format",
            return_value=(object(), "fasta", False),
        )
        mocker.patch.object(
            CompositionPerTaxon,
            "calculate_composition_per_taxon",
            return_value=([], []),
        )
        svc.run()
        out, _ = capsys.readouterr()
        assert out == ""
