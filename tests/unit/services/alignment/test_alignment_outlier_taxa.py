from argparse import Namespace
from pathlib import Path
import subprocess
import sys
import numpy as np
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.alignment_outlier_taxa import AlignmentOutlierTaxa
import phykit.services.alignment.alignment_outlier_taxa as alignment_outlier_taxa_module


here = Path(__file__)


def test_module_import_does_not_import_numpy_or_json_helpers():
    code = """
import sys
import phykit.services.alignment.alignment_outlier_taxa as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_average_site_entropy_common_path_uses_array_sum(monkeypatch):
    site_entropies = np.array([0.5, 1.0, 1.5, 2.0])

    def fail_sum(*_args, **_kwargs):
        raise AssertionError("common site-entropy totals should use ndarray.sum")

    monkeypatch.setattr(
        alignment_outlier_taxa_module.np,
        "sum",
        fail_sum,
        raising=False,
    )

    observed = alignment_outlier_taxa_module._average_site_entropy(
        site_entropies,
        site_entropies.size,
    )

    assert observed == np.mean(site_entropies)


class TestAlignmentOutlierTaxa:
    def _service(self):
        args = Namespace(
            alignment="unused.fa",
            gap_z=3.0,
            composition_z=3.0,
            distance_z=3.0,
            rcvt_z=3.0,
            occupancy_z=3.0,
            entropy_z=3.0,
            json=False,
        )
        return AlignmentOutlierTaxa(args)

    def test_init(self):
        args = Namespace(
            alignment="test.fa",
            gap_z=3.0,
            composition_z=3.0,
            distance_z=3.0,
            rcvt_z=3.0,
            occupancy_z=3.0,
            entropy_z=3.0,
            json=False,
        )
        service = AlignmentOutlierTaxa(args)
        assert service.alignment_file_path == "test.fa"
        assert service.gap_z == 3.0
        assert service.composition_z == 3.0
        assert service.distance_z == 3.0
        assert service.rcvt_z == 3.0
        assert service.occupancy_z == 3.0
        assert service.entropy_z == 3.0

    def test_calculate_outliers_flags_expected_taxa(self):
        args = Namespace(
            alignment="unused.fa",
            gap_z=3.0,
            composition_z=3.0,
            distance_z=3.0,
            rcvt_z=3.0,
            occupancy_z=3.0,
            entropy_z=3.0,
            json=False,
        )
        service = AlignmentOutlierTaxa(args)
        alignment = AlignIO.read(
            f"{here.parent.parent.parent.parent}/sample_files/alignment_outlier_taxa.fa",
            "fasta",
        )

        result = service.calculate_outliers(alignment, is_protein=False)

        flagged = {row["taxon"] for row in result["outliers"]}
        assert "taxon_d" in flagged
        assert "taxon_e" in flagged

        outlier_d = [row for row in result["outliers"] if row["taxon"] == "taxon_d"][0]
        features_d = {reason["feature"] for reason in outlier_d["reasons"]}
        assert "composition_distance" in features_d
        assert "long_branch_proxy" in features_d

        outlier_e = [row for row in result["outliers"] if row["taxon"] == "taxon_e"][0]
        features_e = {reason["feature"] for reason in outlier_e["reasons"]}
        assert "gap_rate" in features_e

    def test_calculate_outliers_handles_lowercase_ambiguity_and_no_overlap(self):
        args = Namespace(
            alignment="unused.fa",
            gap_z=3.0,
            composition_z=3.0,
            distance_z=3.0,
            rcvt_z=3.0,
            occupancy_z=3.0,
            entropy_z=3.0,
            json=False,
        )
        service = AlignmentOutlierTaxa(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACNN"), id="a"),
                SeqRecord(Seq("aG--"), id="b"),
                SeqRecord(Seq("NXTT"), id="c"),
            ]
        )

        dna_result = service.calculate_outliers(alignment, is_protein=False)
        protein_result = service.calculate_outliers(alignment, is_protein=True)

        assert dna_result["rows"][0]["gap_rate"] == 0.5
        assert dna_result["rows"][0]["long_branch_proxy"] == 0.5
        assert dna_result["rows"][2]["long_branch_proxy"] is None
        assert dna_result["rows"][2]["rcvt"] == 0.2222

        assert protein_result["rows"][0]["gap_rate"] == 0.0
        assert protein_result["rows"][2]["long_branch_proxy"] == 1.0
        assert protein_result["rows"][2]["entropy_burden"] == 0.9728

    def test_calculate_outliers_preserves_row_reason_payload(self):
        service = self._service()
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AAAAAA"), id="steady_a"),
                SeqRecord(Seq("AAAAAA"), id="steady_b"),
                SeqRecord(Seq("CCCCCC"), id="composition_outlier"),
                SeqRecord(Seq("------"), id="gap_outlier"),
            ]
        )

        result = service.calculate_outliers(alignment, is_protein=False)
        rows_by_taxon = {row["taxon"]: row for row in result["rows"]}

        gap_row = rows_by_taxon["gap_outlier"]
        assert gap_row["gap_rate"] == 1.0
        assert gap_row["occupancy"] == 0.0
        assert gap_row["long_branch_proxy"] is None
        assert gap_row["flagged"] is True
        assert [
            (reason["feature"], reason["direction"])
            for reason in gap_row["reasons"]
        ] == [
            ("gap_rate", "high"),
            ("occupancy", "low"),
        ]

    def test_calculate_outliers_single_symbol_equal_lengths_shortcut(self, mocker):
        service = self._service()
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("A-A?"), id="a"),
                SeqRecord(Seq("AA??"), id="b"),
                SeqRecord(Seq("-AA-"), id="c"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.alignment_outlier_taxa.np.array",
            side_effect=AssertionError("constant-composition path should skip counts"),
        )

        result = service.calculate_outliers(alignment, is_protein=False)

        assert result["thresholds"] == {
            "gap_rate": None,
            "occupancy": None,
            "composition_distance": None,
            "long_branch_proxy": None,
            "rcvt": None,
            "entropy_burden": None,
        }
        assert result["outliers"] == []
        assert result["rows"] == [
            {
                "taxon": "a",
                "gap_rate": 0.5,
                "occupancy": 0.5,
                "composition_distance": 0.0,
                "long_branch_proxy": 0.0,
                "rcvt": 0.0,
                "entropy_burden": 0.0,
                "flagged": False,
                "reasons": [],
            },
            {
                "taxon": "b",
                "gap_rate": 0.5,
                "occupancy": 0.5,
                "composition_distance": 0.0,
                "long_branch_proxy": 0.0,
                "rcvt": 0.0,
                "entropy_burden": 0.0,
                "flagged": False,
                "reasons": [],
            },
            {
                "taxon": "c",
                "gap_rate": 0.5,
                "occupancy": 0.5,
                "composition_distance": 0.0,
                "long_branch_proxy": 0.0,
                "rcvt": 0.0,
                "entropy_burden": 0.0,
                "flagged": False,
                "reasons": [],
            },
        ]

    def test_calculate_outliers_single_record_skips_matrix_path(self, mocker):
        service = self._service()
        alignment = MultipleSeqAlignment([SeqRecord(Seq("NNNN----"), id="solo")])
        mocker.patch.object(
            alignment_outlier_taxa_module.np,
            "frombuffer",
            side_effect=AssertionError("single-record path should skip matrix work"),
        )

        result = service.calculate_outliers(alignment, is_protein=False)

        assert result["thresholds"] == {
            "gap_rate": None,
            "occupancy": None,
            "composition_distance": None,
            "long_branch_proxy": None,
            "rcvt": None,
            "entropy_burden": None,
        }
        assert result["outliers"] == []
        assert result["rows"] == [
            {
                "taxon": "solo",
                "gap_rate": 1.0,
                "occupancy": 0.0,
                "composition_distance": 0.0,
                "long_branch_proxy": None,
                "rcvt": 0.0,
                "entropy_burden": 0.0,
                "flagged": False,
                "reasons": [],
            }
        ]

    def test_calculate_outliers_identical_multi_symbol_sequences_skip_matrix(
        self, mocker
    ):
        service = self._service()
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTN-?*X"), id="a"),
                SeqRecord(Seq("acgtn-?*x"), id="b"),
                SeqRecord(Seq("ACGTN-?*X"), id="c"),
            ]
        )
        mocker.patch.object(
            alignment_outlier_taxa_module.np,
            "frombuffer",
            side_effect=AssertionError(
                "identical sequences should skip matrix construction"
            ),
        )

        result = service.calculate_outliers(alignment, is_protein=False)

        assert result["thresholds"] == {
            "gap_rate": None,
            "occupancy": None,
            "composition_distance": None,
            "long_branch_proxy": None,
            "rcvt": None,
            "entropy_burden": None,
        }
        assert result["outliers"] == []
        assert result["rows"] == [
            {
                "taxon": "a",
                "gap_rate": 0.5556,
                "occupancy": 0.4444,
                "composition_distance": 0.0,
                "long_branch_proxy": 0.0,
                "rcvt": 0.0,
                "entropy_burden": 0.0,
                "flagged": False,
                "reasons": [],
            },
            {
                "taxon": "b",
                "gap_rate": 0.5556,
                "occupancy": 0.4444,
                "composition_distance": 0.0,
                "long_branch_proxy": 0.0,
                "rcvt": 0.0,
                "entropy_burden": 0.0,
                "flagged": False,
                "reasons": [],
            },
            {
                "taxon": "c",
                "gap_rate": 0.5556,
                "occupancy": 0.4444,
                "composition_distance": 0.0,
                "long_branch_proxy": 0.0,
                "rcvt": 0.0,
                "entropy_burden": 0.0,
                "flagged": False,
                "reasons": [],
            },
        ]

    def test_identical_sequence_helper_does_not_slice_rows(self):
        class NoSliceList(list):
            def __getitem__(self, key):
                if isinstance(key, slice):
                    raise AssertionError("identical-sequence scan should not slice")
                return super().__getitem__(key)

        sequences = NoSliceList(["ACGT", "ACGT", "ACGT"])

        assert alignment_outlier_taxa_module._all_sequences_identical(sequences) is True

    def test_valid_length_for_identical_unicode_sequence_counts_invalid_chars(self):
        sequence = "ACGT\u03a9N-X?"

        assert AlignmentOutlierTaxa._valid_length_for_identical_sequence(
            sequence,
            ["-", "?", "*", "X", "N"],
        ) == 5
        assert AlignmentOutlierTaxa._valid_length_for_identical_sequence(
            sequence,
            ["-", "?", "*", "X"],
        ) == 6

    def test_calculate_outliers_single_symbol_uneven_lengths_uses_normal_path(self):
        service = self._service()
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AAAA"), id="a"),
                SeqRecord(Seq("AA--"), id="b"),
                SeqRecord(Seq("AAAA"), id="c"),
            ]
        )

        result = service.calculate_outliers(alignment, is_protein=False)
        rows = {row["taxon"]: row for row in result["rows"]}

        assert result["thresholds"]["gap_rate"] == 0.0
        assert result["thresholds"]["occupancy"] == 1.0
        assert rows["b"]["gap_rate"] == 0.5
        assert rows["b"]["occupancy"] == 0.5
        assert rows["b"]["flagged"] is True
        assert [
            reason["feature"] for reason in rows["b"]["reasons"]
        ] == ["gap_rate", "occupancy"]

    def test_ascii_path_uses_invalid_lookup_instead_of_isin(self, mocker):
        service = self._service()
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="a"),
                SeqRecord(Seq("A-GN"), id="b"),
                SeqRecord(Seq("TCGT"), id="c"),
            ]
        )
        mocked_isin = mocker.patch.object(alignment_outlier_taxa_module.np, "isin")

        result = service.calculate_outliers(alignment, is_protein=False)

        mocked_isin.assert_not_called()
        assert result["rows"][1]["gap_rate"] == 0.5

    def test_ascii_path_counts_valid_lengths_with_count_nonzero(self, mocker):
        service = self._service()
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="a"),
                SeqRecord(Seq("A-GN"), id="b"),
                SeqRecord(Seq("TCGT"), id="c"),
            ]
        )
        count_nonzero_spy = mocker.spy(
            alignment_outlier_taxa_module.np,
            "count_nonzero",
        )

        result = service.calculate_outliers(alignment, is_protein=False)

        assert result["rows"][1]["gap_rate"] == 0.5
        assert any(
            call.kwargs.get("axis") == 1
            for call in count_nonzero_spy.call_args_list
        )
        assert any(
            call.kwargs.get("axis") == 1 and call.args[0].shape == (3, 3)
            for call in count_nonzero_spy.call_args_list
        )

    def test_all_valid_ascii_path_skips_invalid_lookup(self, mocker):
        service = self._service()
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTACGT"), id="a"),
                SeqRecord(Seq("ACGTTCGT"), id="b"),
                SeqRecord(Seq("TCGTACGA"), id="c"),
            ]
        )
        mocker.patch.object(
            service,
            "_invalid_lookup_for_chars",
            side_effect=AssertionError(
                "all-valid ASCII alignments should not build an invalid lookup"
            ),
        )

        result = service.calculate_outliers(alignment, is_protein=False)

        assert [row["gap_rate"] for row in result["rows"]] == [0.0, 0.0, 0.0]
        assert [row["occupancy"] for row in result["rows"]] == [1.0, 1.0, 1.0]
        assert result["rows"][0]["long_branch_proxy"] == 0.1875

    def test_ascii_symbol_count_helpers_match_equality_reference(self):
        alphabet = b"ACDEFGHIKLMNPQRSTVWY"
        matrix = alignment_outlier_taxa_module.np.frombuffer(
            alphabet + alphabet[::-1] + alphabet[5:] + alphabet[:5],
            dtype=alignment_outlier_taxa_module.np.uint8,
        ).reshape(3, len(alphabet))
        symbols = alignment_outlier_taxa_module.np.unique(matrix)

        row_counts = AlignmentOutlierTaxa._symbol_counts_by_row(matrix, symbols)
        site_counts = AlignmentOutlierTaxa._symbol_counts_by_site(matrix, symbols)
        expected_rows = alignment_outlier_taxa_module.np.array(
            [
                alignment_outlier_taxa_module.np.sum(matrix == symbol, axis=1)
                for symbol in symbols
            ],
            dtype=alignment_outlier_taxa_module.np.float64,
        ).T
        expected_sites = alignment_outlier_taxa_module.np.array(
            [
                alignment_outlier_taxa_module.np.sum(matrix == symbol, axis=0)
                for symbol in symbols
            ],
            dtype=alignment_outlier_taxa_module.np.float64,
        )

        alignment_outlier_taxa_module.np.testing.assert_array_equal(
            row_counts,
            expected_rows,
        )
        alignment_outlier_taxa_module.np.testing.assert_array_equal(
            site_counts,
            expected_sites,
        )

    def test_ascii_symbol_counts_by_row_large_short_uses_single_bincount(
        self,
        monkeypatch,
    ):
        alphabet = b"ACDEFGHIKLMNPQRSTVWY-X?*"
        matrix = alignment_outlier_taxa_module.np.tile(
            alignment_outlier_taxa_module.np.frombuffer(
                alphabet,
                dtype=alignment_outlier_taxa_module.np.uint8,
            ),
            (20_000, 6),
        )
        invalid = alignment_outlier_taxa_module.np.zeros(256, dtype=bool)
        invalid[
            alignment_outlier_taxa_module.np.frombuffer(
                b"-?*X",
                dtype=alignment_outlier_taxa_module.np.uint8,
            )
        ] = True
        symbols = alignment_outlier_taxa_module.np.unique(matrix)
        symbols = symbols[~invalid[symbols]]
        original_bincount = alignment_outlier_taxa_module.np.bincount
        calls = 0

        def count_bincount(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original_bincount(*args, **kwargs)

        monkeypatch.setattr(
            alignment_outlier_taxa_module.np,
            "bincount",
            count_bincount,
        )

        observed = AlignmentOutlierTaxa._symbol_counts_by_row(matrix, symbols)
        expected = alignment_outlier_taxa_module.np.array(
            [
                alignment_outlier_taxa_module.np.sum(matrix == symbol, axis=1)
                for symbol in symbols
            ],
            dtype=alignment_outlier_taxa_module.np.float64,
        ).T

        assert calls == 1
        alignment_outlier_taxa_module.np.testing.assert_array_equal(
            observed,
            expected,
        )

    def test_ascii_symbol_counts_by_site_uses_bincount_above_old_cutoff(
        self,
        monkeypatch,
    ):
        alphabet = b"ACDEFGHIKLMNPQRSTVWY-X?*"
        row = alignment_outlier_taxa_module.np.frombuffer(
            alphabet * 11,
            dtype=alignment_outlier_taxa_module.np.uint8,
        )
        matrix = alignment_outlier_taxa_module.np.tile(row, (32_000, 1))
        invalid = alignment_outlier_taxa_module.np.zeros(256, dtype=bool)
        invalid[
            alignment_outlier_taxa_module.np.frombuffer(
                b"-?*X",
                dtype=alignment_outlier_taxa_module.np.uint8,
            )
        ] = True
        symbols = alignment_outlier_taxa_module.np.unique(matrix)
        symbols = symbols[~invalid[symbols]]
        original_bincount = alignment_outlier_taxa_module.np.bincount
        calls = 0

        def count_bincount(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original_bincount(*args, **kwargs)

        monkeypatch.setattr(
            alignment_outlier_taxa_module.np,
            "bincount",
            count_bincount,
        )

        observed = AlignmentOutlierTaxa._symbol_counts_by_site(matrix, symbols)
        expected = alignment_outlier_taxa_module.np.zeros(
            (symbols.size, matrix.shape[1]),
            dtype=alignment_outlier_taxa_module.np.float64,
        )
        first_row = matrix[0]
        for index, symbol in enumerate(symbols):
            expected[index, first_row == symbol] = matrix.shape[0]

        assert calls == 1
        alignment_outlier_taxa_module.np.testing.assert_array_equal(
            observed,
            expected,
        )

    def test_variable_ascii_outliers_use_symbol_count_helpers(self, monkeypatch):
        service = self._service()
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWY"), id="a"),
                SeqRecord(Seq("YWVTSRQPNMLKIHGFEDCA"), id="b"),
                SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWA"), id="c"),
            ]
        )
        original_row_counts = AlignmentOutlierTaxa._symbol_counts_by_row
        original_site_counts = AlignmentOutlierTaxa._symbol_counts_by_site
        calls = {"rows": 0, "sites": 0}

        def count_rows(alignment_array, symbols):
            calls["rows"] += 1
            return original_row_counts(alignment_array, symbols)

        def count_sites(alignment_array, symbols):
            calls["sites"] += 1
            return original_site_counts(alignment_array, symbols)

        monkeypatch.setattr(
            AlignmentOutlierTaxa,
            "_symbol_counts_by_row",
            staticmethod(count_rows),
        )
        monkeypatch.setattr(
            AlignmentOutlierTaxa,
            "_symbol_counts_by_site",
            staticmethod(count_sites),
        )

        result = service.calculate_outliers(alignment, is_protein=True)

        assert calls == {"rows": 1, "sites": 1}
        assert len(result["rows"]) == 3

    def test_row_l2_norms_matches_linalg_norm(self):
        matrix = np.array(
            [
                [0.1, -0.2, 0.3],
                [0.0, 0.0, 0.0],
                [-0.4, 0.5, -0.6],
            ]
        )

        observed = alignment_outlier_taxa_module._row_l2_norms(matrix)
        expected = np.linalg.norm(matrix, axis=1)

        np.testing.assert_allclose(observed, expected)

    def test_column_dot_matches_explicit_column_sum(self):
        left = np.array(
            [
                [0.2, 0.0, 0.5],
                [0.8, 1.0, 0.5],
            ]
        )
        right = np.array(
            [
                [-2.321928, 0.0, -1.0],
                [-0.321928, 0.0, -1.0],
            ]
        )

        observed = alignment_outlier_taxa_module._column_dot(left, right)
        expected = np.sum(left * right, axis=0)

        np.testing.assert_allclose(observed, expected)

    def test_variable_composition_distance_avoids_linalg_norm(self, monkeypatch):
        service = self._service()
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTACGT"), id="a"),
                SeqRecord(Seq("ACGTTCGT"), id="b"),
                SeqRecord(Seq("TCGTACGA"), id="c"),
            ]
        )

        def fail_norm(*_args, **_kwargs):
            raise AssertionError(
                "composition distances should use the row L2 helper"
            )

        monkeypatch.setattr(
            alignment_outlier_taxa_module.np.linalg,
            "norm",
            fail_norm,
        )

        result = service.calculate_outliers(alignment, is_protein=False)

        assert [row["composition_distance"] for row in result["rows"]] == [
            0.0,
            0.1768,
            0.0,
        ]

    def test_all_valid_long_branch_proxy_matches_pairwise_reference(self):
        matrix = alignment_outlier_taxa_module.np.frombuffer(
            b"ACGT"
            b"ACGA"
            b"TCGA"
            b"TGGA",
            dtype=alignment_outlier_taxa_module.np.uint8,
        ).reshape(4, 4)
        symbols = alignment_outlier_taxa_module.np.unique(matrix)
        site_counts = AlignmentOutlierTaxa._symbol_counts_by_site(matrix, symbols)

        observed = AlignmentOutlierTaxa._all_valid_long_branch_proxy(
            matrix,
            symbols,
            site_counts,
        )
        expected = []
        for row_idx, row in enumerate(matrix):
            distances = []
            for other_idx, other in enumerate(matrix):
                if other_idx != row_idx:
                    distances.append(
                        alignment_outlier_taxa_module.np.mean(row != other)
                    )
            expected.append(sum(distances) / len(distances))

        alignment_outlier_taxa_module.np.testing.assert_allclose(
            observed,
            alignment_outlier_taxa_module.np.array(expected),
        )

    def test_all_valid_ascii_long_branch_uses_site_counts(self, monkeypatch):
        service = self._service()
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTACGT"), id="a"),
                SeqRecord(Seq("ACGTTCGT"), id="b"),
                SeqRecord(Seq("TCGTACGA"), id="c"),
            ]
        )
        original = AlignmentOutlierTaxa._all_valid_long_branch_proxy
        calls = 0

        def count_all_valid(alignment_array, symbols, site_counts):
            nonlocal calls
            calls += 1
            return original(alignment_array, symbols, site_counts)

        monkeypatch.setattr(
            AlignmentOutlierTaxa,
            "_all_valid_long_branch_proxy",
            staticmethod(count_all_valid),
        )

        result = service.calculate_outliers(alignment, is_protein=False)

        assert calls == 1
        assert result["rows"][0]["long_branch_proxy"] == 0.1875

    def test_blocked_long_branch_proxy_matches_pairwise_reference(self):
        service = self._service()
        matrix = alignment_outlier_taxa_module.np.frombuffer(
            b"ACGT"
            b"AG-T"
            b"A--T"
            b"NNNN",
            dtype=alignment_outlier_taxa_module.np.uint8,
        ).reshape(4, 4)
        invalid_lookup = service._invalid_lookup_for_chars(["-", "?", "*", "X", "N"])
        valid_mask = ~invalid_lookup[matrix]
        symbols = alignment_outlier_taxa_module.np.unique(matrix[valid_mask])

        observed = AlignmentOutlierTaxa._blocked_long_branch_proxy(
            matrix,
            valid_mask,
            symbols,
        )

        expected = []
        for row_idx, row in enumerate(matrix):
            distances = []
            for other_idx, other in enumerate(matrix):
                if other_idx == row_idx:
                    continue
                overlap = valid_mask[row_idx] & valid_mask[other_idx]
                if overlap.any():
                    distances.append(
                        alignment_outlier_taxa_module.np.mean(
                            row[overlap] != other[overlap]
                        )
                    )
            expected.append(
                alignment_outlier_taxa_module.np.nan
                if not distances
                else sum(distances) / len(distances)
            )

        alignment_outlier_taxa_module.np.testing.assert_allclose(
            observed,
            alignment_outlier_taxa_module.np.array(expected),
        )

    def test_large_ambiguous_long_branch_uses_blocked_helper(
        self,
        monkeypatch,
    ):
        service = self._service()
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="a"),
                SeqRecord(Seq("AGGT"), id="b"),
                SeqRecord(Seq("A--T"), id="c"),
                SeqRecord(Seq("NNNN"), id="d"),
            ]
        )
        original = AlignmentOutlierTaxa._blocked_long_branch_proxy
        calls = 0

        def count_blocked(alignment_array, valid_mask, symbols):
            nonlocal calls
            calls += 1
            return original(alignment_array, valid_mask, symbols)

        monkeypatch.setattr(AlignmentOutlierTaxa, "_PAIRWISE_MATRIX_MAX_CELLS", 4)
        monkeypatch.setattr(
            AlignmentOutlierTaxa,
            "_blocked_long_branch_proxy",
            staticmethod(count_blocked),
        )

        result = service.calculate_outliers(alignment, is_protein=False)

        assert calls == 1
        assert {
            row["taxon"]: row["long_branch_proxy"] for row in result["rows"]
        } == {
            "a": 0.125,
            "b": 0.125,
            "c": 0.0,
            "d": None,
        }

    def test_long_branch_proxy_matches_pairwise_reference(self):
        service = self._service()
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="a"),
                SeqRecord(Seq("AGGT"), id="b"),
                SeqRecord(Seq("A--T"), id="c"),
                SeqRecord(Seq("NNNN"), id="d"),
            ]
        )
        invalid = {"-", "?", "*", "X", "N"}
        sequences = {record.id: str(record.seq).upper() for record in alignment}

        def pairwise_distance(first, second):
            overlap = [
                (a, b)
                for a, b in zip(first, second)
                if a not in invalid and b not in invalid
            ]
            if not overlap:
                return None
            return sum(a != b for a, b in overlap) / len(overlap)

        expected = {}
        for taxon, sequence in sequences.items():
            distances = []
            for other_taxon, other_sequence in sequences.items():
                if other_taxon == taxon:
                    continue
                distance = pairwise_distance(sequence, other_sequence)
                if distance is not None:
                    distances.append(distance)
            expected[taxon] = (
                None
                if not distances
                else round(sum(distances) / len(distances), 4)
            )

        result = service.calculate_outliers(alignment, is_protein=False)

        assert {
            row["taxon"]: row["long_branch_proxy"] for row in result["rows"]
        } == expected

    def test_invalid_lookup_for_chars_is_case_insensitive(self):
        service = self._service()

        lookup = service._invalid_lookup_for_chars({"n", "?"})

        assert lookup[ord("N")]
        assert lookup[ord("?")]
        assert not lookup[ord("A")]

    def test_unicode_fallback_keeps_non_ascii_symbols(self):
        service = self._service()

        class StringAlignment(list):
            def get_alignment_length(self):
                return len(self[0].seq)

        alignment = StringAlignment(
            [
                type("Record", (), {"id": "a", "seq": "A\u00d1-C"})(),
                type("Record", (), {"id": "b", "seq": "A\u00d1GC"})(),
                type("Record", (), {"id": "c", "seq": "TTGC"})(),
            ]
        )

        result = service.calculate_outliers(alignment, is_protein=False)

        assert result["rows"][0]["gap_rate"] == 0.25
        assert result["rows"][1]["gap_rate"] == 0.0

    def test_run_text_output_batches_report_lines(self, mocker, capsys):
        service = self._service()
        thresholds = {
            "gap_rate": 0.1,
            "occupancy": 0.9,
            "composition_distance": 0.2,
            "long_branch_proxy": 0.3,
            "rcvt": 0.4,
            "entropy_burden": 0.5,
        }
        mocker.patch.object(
            AlignmentOutlierTaxa,
            "get_alignment_and_format",
            return_value=(object(), "fasta", False),
        )
        mocker.patch.object(
            AlignmentOutlierTaxa,
            "calculate_outliers",
            return_value={
                "features": [],
                "thresholds": thresholds,
                "rows": [],
                "outliers": [
                    {
                        "taxon": "taxon_a",
                        "reasons": [
                            {
                                "feature": "gap_rate",
                                "value": 0.25,
                                "threshold": 0.1,
                                "direction": "high",
                                "explanation": "gap rate exceeds threshold",
                            }
                        ],
                    }
                ],
            },
        )

        service.run()

        out, _ = capsys.readouterr()
        assert out == (
            "features_evaluated\t"
            "gap_rate,occupancy,composition_distance,long_branch_proxy,rcvt,entropy_burden\n"
            "thresholds\tgap_rate>0.1;occupancy<0.9;composition_distance>0.2;"
            "long_branch_proxy>0.3;rcvt>0.4;entropy_burden>0.5\n"
            "taxon_a\tgap_rate=0.25>0.1\tgap rate exceeds threshold\n"
        )

    def test_run_text_output_reports_no_outliers(self, mocker, capsys):
        service = self._service()
        thresholds = {
            "gap_rate": 0.1,
            "occupancy": 0.9,
            "composition_distance": 0.2,
            "long_branch_proxy": 0.3,
            "rcvt": 0.4,
            "entropy_burden": 0.5,
        }
        mocker.patch.object(
            AlignmentOutlierTaxa,
            "get_alignment_and_format",
            return_value=(object(), "fasta", False),
        )
        mocker.patch.object(
            AlignmentOutlierTaxa,
            "calculate_outliers",
            return_value={
                "features": [],
                "thresholds": thresholds,
                "rows": [],
                "outliers": [],
            },
        )

        service.run()

        out, _ = capsys.readouterr()
        assert out.endswith("No outlier taxa detected.\n")
