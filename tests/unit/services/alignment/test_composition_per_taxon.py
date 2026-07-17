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
    def test_composition_output_rows_index_small_outputs(self, monkeypatch):
        class IndexedRows:
            def __init__(self):
                self.rows = [np.array([1.0]), np.array([0.5])]
                self.indexes = []

            def __getitem__(self, idx):
                self.indexes.append(idx)
                return self.rows[idx]

        monkeypatch.setattr(
            composition_per_taxon_module,
            "_COMPOSITION_ROW_ZIP_MIN_COUNT",
            3,
        )
        freqs = IndexedRows()

        rows = composition_per_taxon_module._composition_output_rows(
            ["t1", "t2"],
            freqs,
        )

        assert [taxon for taxon, _ in rows] == ["t1", "t2"]
        assert freqs.indexes == [0, 1]

    def test_composition_output_rows_zip_large_outputs(self, monkeypatch):
        class IterableRows:
            def __iter__(self):
                return iter([np.array([1.0]), np.array([0.5]), np.array([0.0])])

            def __getitem__(self, _idx):
                raise AssertionError("large composition row output should iterate")

        monkeypatch.setattr(
            composition_per_taxon_module,
            "_COMPOSITION_ROW_ZIP_MIN_COUNT",
            3,
        )

        rows = composition_per_taxon_module._composition_output_rows(
            ["t1", "t2", "t3"],
            IterableRows(),
        )

        assert [taxon for taxon, _ in rows] == ["t1", "t2", "t3"]
        assert [values.tolist() for _, values in rows] == [[1.0], [0.5], [0.0]]

    def test_identical_composition_output_rows_marks_shared_composition(self):
        freqs = np.array([0.5, 0.5])

        rows = composition_per_taxon_module._identical_composition_output_rows(
            ["t1", "t2"],
            freqs,
        )

        assert rows.shared_composition is freqs
        assert [taxon for taxon, _ in rows] == ["t1", "t2"]
        assert [values.tolist() for _, values in rows] == [[0.5, 0.5], [0.5, 0.5]]
        assert rows[0][1] is not rows[1][1]

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

    def test_lazy_numpy_caches_resolved_attributes(self):
        lazy_np = composition_per_taxon_module._LazyNumpy()

        array_attr = lazy_np.array

        assert lazy_np.__dict__["array"] is array_attr
        assert lazy_np.array is array_attr
        assert lazy_np._module is not None

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

    def test_clean_dna_skips_full_symbol_discovery(self, monkeypatch):
        service = CompositionPerTaxon(Namespace(alignment="x.fa"))
        alignment = [
            SeqRecord(Seq("ACGT"), id="t1"),
            SeqRecord(Seq("AGGT"), id="t2"),
            SeqRecord(Seq("TCGA"), id="t3"),
        ]

        def fail_unique(*args, **kwargs):
            raise AssertionError("clean DNA should use known symbol codes")

        monkeypatch.setattr(
            composition_per_taxon_module,
            "_COMPOSITION_SCALAR_MAX_CELLS",
            0,
        )
        monkeypatch.setattr(
            composition_per_taxon_module.np,
            "unique",
            fail_unique,
        )

        symbols, rows = service.calculate_composition_per_taxon(
            alignment,
            is_protein=False,
        )

        assert symbols == ["A", "C", "G", "T"]
        np.testing.assert_allclose(
            [frequencies for _, frequencies in rows],
            [
                [0.25, 0.25, 0.25, 0.25],
                [0.25, 0.0, 0.5, 0.25],
                [0.25, 0.25, 0.25, 0.25],
            ],
        )

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

    def test_calculate_composition_per_taxon_preserves_variable_row_order(self, args):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("aacc"), id="first"),
                SeqRecord(Seq("GGTT"), id="second"),
                SeqRecord(Seq("ACGT"), id="third"),
            ]
        )

        symbols, rows = svc.calculate_composition_per_taxon(alignment, is_protein=False)

        assert symbols == ["A", "C", "G", "T"]
        assert [taxon for taxon, _ in rows] == ["first", "second", "third"]
        assert np.allclose(rows[0][1], np.array([0.5, 0.5, 0.0, 0.0]))
        assert np.allclose(rows[1][1], np.array([0.0, 0.0, 0.5, 0.5]))
        assert np.allclose(rows[2][1], np.full(4, 0.25))

    def test_calculate_composition_per_taxon_ascii_small_alphabet_skips_bincount(
        self, args, mocker, monkeypatch
    ):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AACCGGTT"), id="t1"),
                SeqRecord(Seq("AA--NN??"), id="t2"),
            ]
        )
        monkeypatch.setattr(
            composition_per_taxon_module,
            "_COMPOSITION_SCALAR_MAX_CELLS",
            0,
        )
        bincount_spy = mocker.spy(composition_per_taxon_module.np, "bincount")
        bounded_count_spy = mocker.spy(
            composition_per_taxon_module,
            "_bounded_ascii_row_symbol_counts",
        )

        symbols, rows = svc.calculate_composition_per_taxon(alignment, is_protein=False)
        by_taxon = {taxon: vals for taxon, vals in rows}

        bincount_spy.assert_not_called()
        bounded_count_spy.assert_called_once()
        assert symbols == ["A", "C", "G", "T"]
        assert np.allclose(by_taxon["t1"], np.array([0.25, 0.25, 0.25, 0.25]))
        assert np.allclose(by_taxon["t2"], np.array([1.0, 0.0, 0.0, 0.0]))

    @pytest.mark.parametrize("num_sites", [0xFFFF, 0x10000])
    def test_bounded_ascii_row_symbol_counts_do_not_overflow(self, num_sites):
        alignment_array = np.full((2, num_sites), ord("A"), dtype=np.uint8)
        symbol_values = np.array([ord("A"), ord("C")], dtype=np.uint8)

        observed = composition_per_taxon_module._bounded_ascii_row_symbol_counts(
            alignment_array,
            symbol_values,
        )

        np.testing.assert_array_equal(
            observed,
            np.array(
                [
                    [num_sites, 0],
                    [num_sites, 0],
                ],
                dtype=np.float64,
            ),
        )

    def test_calculate_composition_per_taxon_ascii_large_alphabet_uses_offset_bincount(
        self, args, mocker, monkeypatch
    ):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWY"), id="t1"),
                SeqRecord(Seq("YWVTSRQPNMLKIHGFEDCA"), id="t2"),
            ]
        )
        monkeypatch.setattr(
            composition_per_taxon_module,
            "_COMPOSITION_SCALAR_MAX_CELLS",
            0,
        )
        bincount_spy = mocker.spy(composition_per_taxon_module.np, "bincount")
        offset_spy = mocker.spy(
            composition_per_taxon_module,
            "_row_symbol_counts_from_ascii_codes",
        )

        symbols, rows = svc.calculate_composition_per_taxon(alignment, is_protein=True)
        by_taxon = {taxon: vals for taxon, vals in rows}

        assert bincount_spy.call_count == 1
        offset_spy.assert_called_once()
        assert symbols == list("ACDEFGHIKLMNPQRSTVWY")
        assert np.allclose(by_taxon["t1"], np.full(20, 0.05))

    def test_calculate_composition_per_taxon_ascii_no_gap_uses_full_lengths(
        self, args, mocker, monkeypatch
    ):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWY"), id="t1"),
                SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWA"), id="t2"),
            ]
        )
        monkeypatch.setattr(
            composition_per_taxon_module,
            "_COMPOSITION_SCALAR_MAX_CELLS",
            0,
        )
        full_spy = mocker.spy(composition_per_taxon_module.np, "full")

        symbols, rows = svc.calculate_composition_per_taxon(
            alignment,
            is_protein=True,
        )
        by_taxon = {taxon: vals for taxon, vals in rows}

        full_spy.assert_called_once_with(
            2,
            20,
            dtype=composition_per_taxon_module.np.uint16,
        )
        assert symbols == list("ACDEFGHIKLMNPQRSTVWY")
        assert np.allclose(by_taxon["t1"], np.full(20, 0.05))

    def test_calculate_composition_per_taxon_small_alignment_uses_scalar_path(
        self, args, mocker
    ):
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("aCGN"), id="t1"),
                SeqRecord(Seq("a--t"), id="t2"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.composition_per_taxon.np.frombuffer",
            side_effect=AssertionError(
                "small composition should avoid byte matrix setup"
            ),
        )
        mocker.patch(
            "phykit.services.alignment.composition_per_taxon.np.array",
            side_effect=AssertionError(
                "small composition should avoid NumPy array setup"
            ),
        )

        symbols, rows = svc.calculate_composition_per_taxon(
            alignment,
            is_protein=False,
        )
        by_taxon = {taxon: vals for taxon, vals in rows}

        assert symbols == ["A", "C", "G", "T"]
        assert np.allclose(by_taxon["t1"], [1 / 3, 1 / 3, 1 / 3, 0.0])
        assert np.allclose(by_taxon["t2"], [0.5, 0.0, 0.0, 0.5])

    def test_composition_per_taxon_scalar_rejects_uneven_lengths(self):
        assert composition_per_taxon_module._composition_per_taxon_scalar(
            [("t1", "ACGT"), ("t2", "ACG")],
            {"-", "?", "*", "N", "X"},
        ) is None

    def test_calculate_composition_per_taxon_large_short_protein_uses_offset_bincount(
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
        offset_spy = mocker.spy(
            composition_per_taxon_module,
            "_row_symbol_counts_from_ascii_codes",
        )

        symbols, rows = svc.calculate_composition_per_taxon(
            alignment,
            is_protein=True,
        )
        by_taxon = {taxon: vals for taxon, vals in rows}

        assert bincount_spy.call_count == 1
        offset_spy.assert_called_once()
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
        mocker.patch(
            "phykit.services.alignment.composition_per_taxon.np.sum",
            side_effect=AssertionError(
                "identical sequence frequencies should use counts.sum"
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

    def test_large_indexed_identical_sequences_skip_raw_record_iteration(
        self, args, monkeypatch
    ):
        class IndexedAlignment:
            def __init__(self):
                self.records = [
                    SimpleNamespace(seq="ACGT", id="t1"),
                    SimpleNamespace(seq="acgt", id="t2"),
                    SimpleNamespace(seq="ACGT", id="t3"),
                ]

            def __len__(self):
                return len(self.records)

            def __getitem__(self, idx):
                return self.records[idx]

            def __iter__(self):
                raise AssertionError("indexed identical path should not iterate")

        monkeypatch.setattr(
            composition_per_taxon_module,
            "_IDENTICAL_INDEXED_SCAN_MIN_RECORDS",
            3,
        )
        svc = CompositionPerTaxon(args)

        symbols, rows = svc.calculate_composition_per_taxon(
            IndexedAlignment(),
            is_protein=False,
        )
        rows[0][1][0] = 0.0

        assert symbols == ["A", "C", "G", "T"]
        assert [taxon for taxon, _ in rows] == ["t1", "t2", "t3"]
        assert np.allclose(rows[1][1], np.full(4, 0.25))

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

    def test_large_identical_composition_rows_use_tiled_matrix(
        self, args, monkeypatch, mocker
    ):
        monkeypatch.setattr(
            composition_per_taxon_module,
            "_IDENTICAL_ROW_TILE_MIN_COUNT",
            3,
        )
        tile_spy = mocker.spy(composition_per_taxon_module.np, "tile")
        svc = CompositionPerTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="t1"),
                SeqRecord(Seq("ACGT"), id="t2"),
                SeqRecord(Seq("ACGT"), id="t3"),
            ]
        )

        symbols, rows = svc.calculate_composition_per_taxon(
            alignment,
            is_protein=False,
        )
        rows[0][1][0] = 0.0

        tile_spy.assert_called_once()
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
        assert payload["rows"][0] == {
            "taxon": "t1",
            "composition": {"A": 0.5, "C": 0.5},
        }

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

    def test_run_small_fasta_uses_direct_scalar_path(
        self,
        mocker,
        tmp_path,
        capsys,
    ):
        alignment = tmp_path / "small.fa"
        alignment.write_text(">t1\nACGT\n>t2\nA-GT\n")
        svc = CompositionPerTaxon(Namespace(alignment=str(alignment), json=False))
        get_alignment = mocker.patch.object(
            CompositionPerTaxon,
            "get_alignment_and_format",
            side_effect=AssertionError("small FASTA should use direct scalar path"),
        )

        svc.run()

        out, _ = capsys.readouterr()
        assert out == (
            "t1\tA:0.25;C:0.25;G:0.25;T:0.25\n"
            "t2\tA:0.3333;C:0.0;G:0.3333;T:0.3333\n"
        )
        get_alignment.assert_not_called()

    def test_run_large_fasta_falls_back_to_alignment_reader(self, mocker, tmp_path):
        alignment = tmp_path / "large.fa"
        alignment.write_text(">t1\n" + "A" * 5000 + "\n>t2\n" + "C" * 5000 + "\n")
        svc = CompositionPerTaxon(Namespace(alignment=str(alignment), json=False))
        get_alignment = mocker.patch.object(
            CompositionPerTaxon,
            "get_alignment_and_format",
            return_value=(object(), "fasta", False),
        )
        calculate = mocker.patch.object(
            CompositionPerTaxon,
            "calculate_composition_per_taxon",
            return_value=([], []),
        )

        svc.run()

        get_alignment.assert_called_once_with()
        calculate.assert_called_once()

    def test_simple_fasta_reader_classifies_protein_symbols(self, tmp_path):
        alignment = tmp_path / "protein.fa"
        alignment.write_text(">p1 description\nACD\n>p2\nAXD\n")

        records, is_protein = (
            composition_per_taxon_module._read_small_simple_fasta_records(
                str(alignment)
            )
        )

        assert records == [("p1", "ACD"), ("p2", "AXD")]
        assert is_protein is True

    def test_run_text_output_reuses_shared_identical_composition(
        self,
        mocker,
        capsys,
    ):
        svc = CompositionPerTaxon(Namespace(alignment="x.fa", json=False))
        rows = composition_per_taxon_module._identical_composition_output_rows(
            ["t1", "t2"],
            np.array([0.5, 0.5]),
        )
        mocker.patch.object(
            CompositionPerTaxon,
            "get_alignment_and_format",
            return_value=(object(), "fasta", False),
        )
        mocker.patch.object(
            CompositionPerTaxon,
            "calculate_composition_per_taxon",
            return_value=(["A", "C"], rows),
        )

        svc.run()

        out, _ = capsys.readouterr()
        assert out == "t1\tA:0.5;C:0.5\nt2\tA:0.5;C:0.5\n"

    def test_run_json_output_reuses_shared_identical_composition(self, mocker):
        svc = CompositionPerTaxon(Namespace(alignment="x.fa", json=True))
        rows = composition_per_taxon_module._identical_composition_output_rows(
            ["t1", "t2"],
            np.array([0.5, 0.5]),
        )
        mocker.patch.object(
            CompositionPerTaxon,
            "get_alignment_and_format",
            return_value=(object(), "fasta", False),
        )
        mocker.patch.object(
            CompositionPerTaxon,
            "calculate_composition_per_taxon",
            return_value=(["A", "C"], rows),
        )
        mocked_json = mocker.patch.object(composition_per_taxon_module, "print_json")

        svc.run()

        payload = mocked_json.call_args.args[0]
        assert payload["rows"] == [
            {"taxon": "t1", "composition": {"A": 0.5, "C": 0.5}},
            {"taxon": "t2", "composition": {"A": 0.5, "C": 0.5}},
        ]
        assert (
            payload["rows"][0]["composition"]
            is not payload["rows"][1]["composition"]
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

    def test_scalar_helpers_cover_empty_invalid_and_zero_rows(self):
        assert composition_per_taxon_module._composition_per_taxon_scalar(
            [],
            {"-", "?", "N"},
        ) == ([], [])
        assert composition_per_taxon_module._composition_per_taxon_scalar(
            [("t1", "---"), ("t2", "???")],
            {"-", "?", "N"},
        ) == ([], [])

        symbols, rows = composition_per_taxon_module._composition_per_taxon_scalar(
            [("invalid", "---"), ("valid", "AAA")],
            {"-", "?", "N"},
        )
        assert symbols == ["A"]
        assert rows[0][0] == "invalid"
        assert rows[0][1] == [0.0]
        assert rows[1][1] == [1.0]

    def test_parser_helpers_clean_wrapped_and_unicode_sequences(self):
        assert composition_per_taxon_module._clean_simple_fasta_sequence(
            ["ac ", "g\r"]
        ) == "acg"
        assert composition_per_taxon_module._sequence_has_protein_symbols(
            "ACGTN-?*"
        ) is False
        assert composition_per_taxon_module._sequence_has_protein_symbols(
            "ACΩT"
        ) is True

    @pytest.mark.parametrize(
        "contents",
        [
            "",
            "not fasta\n>t1\nACGT\n",
            ">t1\nACGT\n>t2\nACG\n",
            ">t1\nACGT\n>t2\nACG\n>t3\nACGT\n",
        ],
    )
    def test_simple_fasta_reader_rejects_invalid_inputs(
        self, tmp_path, contents
    ):
        path = tmp_path / "input.fa"
        path.write_text(contents)

        assert composition_per_taxon_module._read_small_simple_fasta_records(
            str(path)
        ) is None

    def test_simple_fasta_reader_handles_missing_and_oversized_inputs(
        self, tmp_path, monkeypatch
    ):
        assert composition_per_taxon_module._read_small_simple_fasta_records(
            str(tmp_path / "missing.fa")
        ) is None

        path = tmp_path / "large.fa"
        path.write_text(">t1\nACGT\n>t2\nACGT\n")
        monkeypatch.setattr(
            composition_per_taxon_module,
            "_COMPOSITION_SCALAR_MAX_CELLS",
            4,
        )
        assert composition_per_taxon_module._read_small_simple_fasta_records(
            str(path)
        ) is None

    def test_bounded_count_dtype_supports_all_integer_ranges(self):
        assert composition_per_taxon_module._bounded_ascii_count_dtype(
            0xFFFF
        ) is np.uint16
        assert composition_per_taxon_module._bounded_ascii_count_dtype(
            0x10000
        ) is np.uint32
        assert composition_per_taxon_module._bounded_ascii_count_dtype(
            0x100000000
        ) is np.uint64

    def test_invalid_lookup_cache_reuses_existing_array(self, args):
        service = CompositionPerTaxon(args)
        service._INVALID_LOOKUP_CACHE.clear()

        first = service._invalid_lookup_for_chars({"n", "?"})
        second = service._invalid_lookup_for_chars({"N", "?"})

        assert second is first

    def test_calculate_composition_accepts_unknown_length_iterables(self, args):
        service = CompositionPerTaxon(args)
        records = iter(
            [
                SimpleNamespace(id="a", seq="ACGT"),
                SimpleNamespace(id="b", seq="AGGT"),
            ]
        )

        symbols, rows = service.calculate_composition_per_taxon(
            records,
            is_protein=False,
        )

        assert symbols == ["A", "C", "G", "T"]
        assert [record_id for record_id, _ in rows] == ["a", "b"]
        assert service.calculate_composition_per_taxon([], False) == ([], [])

    def test_indexed_single_record_skips_zero_sample(self, args, monkeypatch):
        class IndexedAlignment:
            records = [SimpleNamespace(id="only", seq="ACGT")]

            def __len__(self):
                return 1

            def __getitem__(self, index):
                return self.records[index]

            def __iter__(self):
                raise AssertionError("indexed path should not iterate")

        monkeypatch.setattr(
            composition_per_taxon_module,
            "_IDENTICAL_INDEXED_SCAN_MIN_RECORDS",
            1,
        )
        service = CompositionPerTaxon(args)

        symbols, rows = service.calculate_composition_per_taxon(
            IndexedAlignment(),
            is_protein=False,
        )

        assert symbols == ["A", "C", "G", "T"]
        assert [record_id for record_id, _ in rows] == ["only"]

    def test_indexed_identical_unicode_and_invalid_sequences(
        self, args, monkeypatch
    ):
        class IndexedAlignment:
            def __init__(self, sequence):
                self.records = [
                    SimpleNamespace(id=f"t{idx}", seq=sequence)
                    for idx in range(3)
                ]

            def __len__(self):
                return len(self.records)

            def __getitem__(self, index):
                return self.records[index]

            def __iter__(self):
                raise AssertionError("indexed path should not iterate")

        monkeypatch.setattr(
            composition_per_taxon_module,
            "_IDENTICAL_INDEXED_SCAN_MIN_RECORDS",
            3,
        )
        service = CompositionPerTaxon(args)

        symbols, rows = service.calculate_composition_per_taxon(
            IndexedAlignment("AΩA"),
            is_protein=True,
        )
        assert symbols == ["A", "Ω"]
        np.testing.assert_allclose(rows[0][1], [2 / 3, 1 / 3])
        assert service.calculate_composition_per_taxon(
            IndexedAlignment("---"),
            is_protein=False,
        ) == ([], [])

    @pytest.mark.parametrize("mismatch_index", [1, 2])
    def test_indexed_scan_mismatches_fall_back_to_complete_alignment(
        self, args, monkeypatch, mismatch_index
    ):
        class IndexedAlignment:
            def __init__(self):
                self.records = [
                    SimpleNamespace(
                        id=f"t{idx}",
                        seq="TCGT" if idx == mismatch_index else "ACGT",
                    )
                    for idx in range(5)
                ]

            def __len__(self):
                return len(self.records)

            def __getitem__(self, index):
                return self.records[index]

            def __iter__(self):
                return iter(self.records)

        monkeypatch.setattr(
            composition_per_taxon_module,
            "_IDENTICAL_INDEXED_SCAN_MIN_RECORDS",
            5,
        )
        service = CompositionPerTaxon(args)

        symbols, rows = service.calculate_composition_per_taxon(
            IndexedAlignment(),
            is_protein=False,
        )

        assert symbols == ["A", "C", "G", "T"]
        assert len(rows) == 5

    def test_indexed_scan_errors_fall_back_to_iteration(self, args, monkeypatch):
        class FragileIndexedAlignment:
            records = [
                SimpleNamespace(id="t0", seq="ACGT"),
                SimpleNamespace(id="t1", seq="AGGT"),
                SimpleNamespace(id="t2", seq="ATGT"),
            ]

            def __len__(self):
                return len(self.records)

            def __getitem__(self, index):
                if index:
                    raise TypeError("indexing unavailable")
                return self.records[index]

            def __iter__(self):
                return iter(self.records)

        monkeypatch.setattr(
            composition_per_taxon_module,
            "_IDENTICAL_INDEXED_SCAN_MIN_RECORDS",
            3,
        )
        service = CompositionPerTaxon(args)

        symbols, rows = service.calculate_composition_per_taxon(
            FragileIndexedAlignment(),
            is_protein=False,
        )

        assert symbols == ["A", "C", "G", "T"]
        assert len(rows) == 3

    def test_forced_matrix_backends_match_scalar_compositions(
        self, args, monkeypatch
    ):
        def composition_oracle(records, invalid_chars):
            invalid = {char.upper() for char in invalid_chars}
            valid_rows = [
                [char for char in str(record.seq).upper() if char not in invalid]
                for record in records
            ]
            symbols = sorted({char for row in valid_rows for char in row})
            values = [
                [row.count(symbol) / len(row) if row else 0.0 for symbol in symbols]
                for row in valid_rows
            ]
            return symbols, values

        service = CompositionPerTaxon(args)
        monkeypatch.setattr(
            composition_per_taxon_module,
            "_COMPOSITION_SCALAR_MAX_CELLS",
            0,
        )
        monkeypatch.setattr(
            composition_per_taxon_module,
            "_ASCII_OFFSET_BINCOUNT_MAX_ROWS",
            1,
        )

        cases = [
            (
                [
                    SimpleNamespace(id="a", seq="A-A"),
                    SimpleNamespace(id="b", seq="AAA"),
                ],
                False,
            ),
            (
                [
                    SimpleNamespace(id="a", seq="ACDEFGHIKLMNPQRSTVWY"),
                    SimpleNamespace(id="b", seq="YWVTSRQPNMLKIHGFEDCA"),
                ],
                True,
            ),
            (
                [
                    SimpleNamespace(id="a", seq="AΩCC"),
                    SimpleNamespace(id="b", seq="TΩGG"),
                ],
                True,
            ),
        ]

        for records, is_protein in cases:
            symbols, rows = service.calculate_composition_per_taxon(
                records,
                is_protein=is_protein,
            )
            expected_symbols, expected_rows = composition_oracle(
                records,
                service.get_gap_chars(is_protein),
            )
            assert symbols == expected_symbols
            np.testing.assert_allclose(
                [values for _, values in rows],
                expected_rows,
            )

    def test_identical_unicode_and_invalid_matrix_shortcuts(self, args):
        service = CompositionPerTaxon(args)
        unicode_records = [
            SimpleNamespace(id="a", seq="AΩA"),
            SimpleNamespace(id="b", seq="aωa"),
        ]
        invalid_records = [
            SimpleNamespace(id="a", seq="---"),
            SimpleNamespace(id="b", seq="---"),
        ]

        symbols, rows = service.calculate_composition_per_taxon(
            unicode_records,
            is_protein=True,
        )

        assert symbols == ["A", "Ω"]
        np.testing.assert_allclose(rows[0][1], [2 / 3, 1 / 3])
        assert service.calculate_composition_per_taxon(
            invalid_records,
            is_protein=False,
        ) == ([], [])
