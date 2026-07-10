import pytest
import subprocess
import sys
from argparse import Namespace
from collections import Counter
from math import isclose
from types import SimpleNamespace

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.parsimony_informative_sites import ParsimonyInformative
import phykit.services.alignment.parsimony_informative_sites as pi_module


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestParsimonyInformative(object):
    def test_module_import_does_not_import_numpy_or_biopython_align(self):
        code = """
import sys
import phykit.services.alignment.parsimony_informative_sites as module
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_lazy_numpy_caches_resolved_attributes(self):
        lazy_np = pi_module._LazyNumpy()

        frombuffer_attr = lazy_np.frombuffer

        assert lazy_np.__dict__["frombuffer"] is frombuffer_attr
        assert lazy_np.frombuffer is frombuffer_attr
        assert lazy_np._module is not None

    def test_init_sets_alignment_file_path(self, args):
        pi = ParsimonyInformative(args)
        assert pi.alignment_file_path == args.alignment
        assert pi.output_file_path is None

    def test_parsimony_informative_sites(self, alignment_simple, args):
        pi = ParsimonyInformative(args)
        pi_sites, aln_len, pi_sites_per = pi.calculate_parsimony_informative_sites(
            alignment_simple
        )
        assert isinstance(pi_sites, int)
        assert isinstance(aln_len, int)
        assert isinstance(pi_sites_per, float)
        assert pi_sites == 3
        assert aln_len == 6
        assert isclose(pi_sites_per, 50.0, rel_tol=0.001)

    def test_parsimony_informative_sites_ignores_ambiguous_sites_and_uppercases(self, args):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("AC-N"), id="t1"),
            SeqRecord(Seq("aCGN"), id="t2"),
            SeqRecord(Seq("TGGN"), id="t3"),
            SeqRecord(Seq("tG?N"), id="t4"),
        ])
        pi = ParsimonyInformative(args)

        pi_sites, aln_len, pi_sites_per = pi.calculate_parsimony_informative_sites(
            alignment,
            is_protein=False,
        )

        assert pi_sites == 2
        assert aln_len == 4
        assert isclose(pi_sites_per, 50.0, rel_tol=0.001)

    def test_parsimony_informative_sites_ascii_matrix_path(self, args):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("aa"), id="t1"),
            SeqRecord(Seq("at"), id="t2"),
            SeqRecord(Seq("ta"), id="t3"),
            SeqRecord(Seq("tt"), id="t4"),
        ])
        pi = ParsimonyInformative(args)

        pi_sites, aln_len, pi_sites_per = pi.calculate_parsimony_informative_sites(
            alignment,
            is_protein=False,
        )

        assert pi_sites == 2
        assert aln_len == 2
        assert isclose(pi_sites_per, 100.0, rel_tol=0.001)

    def test_get_number_of_occurrences_counts_records_without_column_slice(self, args):
        class NoColumnSliceAlignment:
            def __init__(self, records):
                self.records = records

            def __iter__(self):
                return iter(self.records)

            def __getitem__(self, key):
                if isinstance(key, tuple):
                    raise AssertionError("column slicing should not be used")
                return self.records[key]

        alignment = NoColumnSliceAlignment(
            [
                SeqRecord(Seq("A"), id="t1"),
                SeqRecord(Seq("C"), id="t2"),
                SeqRecord(Seq("A"), id="t3"),
                SeqRecord(Seq("G"), id="t4"),
                SeqRecord(Seq("C"), id="t5"),
                SeqRecord(Seq("C"), id="t6"),
                SeqRecord(Seq("-"), id="gap"),
            ]
        )
        pi = ParsimonyInformative(args)

        counts = pi.get_number_of_occurrences_per_character(
            alignment, 0, is_protein=False
        )

        assert counts == Counter({"A": 2, "C": 3, "G": 1})

    def test_is_parsimony_informative_stops_after_two_recurrent_states(self, args):
        class RecurrentCounts:
            def values(self):
                yield 2
                yield 2
                raise AssertionError("extra counts should not be scanned")

        pi = ParsimonyInformative(args)

        assert pi.is_parsimony_informative(RecurrentCounts())

    def test_parsimony_informative_sites_ascii_path_avoids_isin(self, mocker, args):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("aa"), id="t1"),
            SeqRecord(Seq("a-"), id="t2"),
            SeqRecord(Seq("ta"), id="t3"),
            SeqRecord(Seq("tt"), id="t4"),
        ])
        pi = ParsimonyInformative(args)
        mocker.patch(
            "phykit.services.alignment.parsimony_informative_sites.np.isin",
            side_effect=AssertionError("ASCII path should use the lookup table"),
        )

        pi_sites, aln_len, pi_sites_per = pi.calculate_parsimony_informative_sites(
            alignment,
            is_protein=False,
        )

        assert pi_sites == 1
        assert aln_len == 2
        assert isclose(pi_sites_per, 50.0, rel_tol=0.001)

    def test_parsimony_informative_sites_ascii_path_uses_gap_code_reduction(
        self, mocker, monkeypatch, args
    ):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("aa"), id="t1"),
            SeqRecord(Seq("a-"), id="t2"),
            SeqRecord(Seq("ta"), id="t3"),
            SeqRecord(Seq("tt"), id="t4"),
        ])
        pi = ParsimonyInformative(args)
        monkeypatch.setattr(
            pi_module,
            "_PARSIMONY_INFORMATIVE_SCALAR_MAX_CELLS",
            0,
        )
        gap_codes_spy = mocker.spy(pi_module, "_get_gap_codes")

        pi_sites, aln_len, pi_sites_per = pi.calculate_parsimony_informative_sites(
            alignment,
            is_protein=False,
        )

        gap_codes_spy.assert_called_once_with(False)
        assert pi_sites == 1
        assert aln_len == 2
        assert isclose(pi_sites_per, 50.0, rel_tol=0.001)

    def test_parsimony_informative_sites_ascii_path_uses_column_histogram(
        self, mocker, monkeypatch, args
    ):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("ACD-"), id="t1"),
            SeqRecord(Seq("ACE-"), id="t2"),
            SeqRecord(Seq("ACF-"), id="t3"),
            SeqRecord(Seq("ACD-"), id="t4"),
            SeqRecord(Seq("ACE-"), id="t5"),
        ])
        pi = ParsimonyInformative(args)
        monkeypatch.setattr(
            pi_module,
            "_PARSIMONY_INFORMATIVE_SCALAR_MAX_CELLS",
            0,
        )
        mocker.patch(
            "phykit.services.alignment.parsimony_informative_sites.np.unique",
            side_effect=AssertionError("ASCII path should not scan each unique symbol"),
        )

        pi_sites, aln_len, pi_sites_per = pi.calculate_parsimony_informative_sites(
            alignment,
            is_protein=True,
        )

        assert pi_sites == 1
        assert aln_len == 4
        assert isclose(pi_sites_per, 25.0, rel_tol=0.001)

    def test_ascii_parsimony_informative_sites_defaults_to_tuned_blocks(
        self,
        monkeypatch,
    ):
        np = pi_module.np
        width = pi_module._ASCII_PI_BLOCK_SIZE + 1
        alignment_array = np.array(
            [
                [ord("A") if (row + col) % 2 else ord("C") for col in range(width)]
                for row in range(4)
            ],
            dtype=np.uint8,
        )
        original_bincount = pi_module.np.bincount
        min_lengths = []

        def bincount_spy(values, *args, **kwargs):
            min_lengths.append(kwargs.get("minlength"))
            return original_bincount(values, *args, **kwargs)

        monkeypatch.setattr(pi_module.np, "bincount", bincount_spy)

        pi_sites = pi_module._count_ascii_parsimony_informative_sites(
            alignment_array,
        )

        assert pi_sites == width
        assert min_lengths == [
            pi_module._ASCII_PI_BLOCK_SIZE * 256,
            256,
        ]

    def test_parsimony_informative_clean_dna_skips_valid_mask(
        self, mocker, monkeypatch, args
    ):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("ACGT"), id="t1"),
            SeqRecord(Seq("ACGA"), id="t2"),
            SeqRecord(Seq("TCGC"), id="t3"),
            SeqRecord(Seq("TCGT"), id="t4"),
        ])
        pi = ParsimonyInformative(args)
        monkeypatch.setattr(
            pi_module,
            "_PARSIMONY_INFORMATIVE_SCALAR_MAX_CELLS",
            0,
        )
        fast_count = mocker.spy(
            pi_module,
            "_count_clean_dna_parsimony_informative_sites",
        )
        generic_count = mocker.spy(
            pi_module,
            "_count_ascii_parsimony_informative_sites",
        )
        mocker.patch(
            "phykit.services.alignment.parsimony_informative_sites.np.sum",
            side_effect=AssertionError("clean DNA helper should use ndarray.sum"),
        )

        pi_sites, aln_len, pi_sites_per = pi.calculate_parsimony_informative_sites(
            alignment,
            is_protein=False,
        )

        fast_count.assert_called_once()
        generic_count.assert_not_called()
        assert pi_sites == 1
        assert aln_len == 4
        assert isclose(pi_sites_per, 25.0, rel_tol=0.001)

    @pytest.mark.parametrize("num_records", [0xFFFF, 0x10000])
    def test_clean_dna_counts_use_non_overflowing_accumulator(self, num_records):
        alignment_array = pi_module.np.full(
            (num_records, 2),
            ord("A"),
            dtype=pi_module.np.uint8,
        )

        assert (
            pi_module._count_clean_dna_parsimony_informative_sites(
                alignment_array,
            )
            == 0
        )

    def test_parsimony_informative_nonstandard_dna_uses_generic_clean_path(
        self, mocker, monkeypatch, args
    ):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("ABGT"), id="t1"),
            SeqRecord(Seq("ABGA"), id="t2"),
            SeqRecord(Seq("TBGG"), id="t3"),
            SeqRecord(Seq("TBGG"), id="t4"),
        ])
        pi = ParsimonyInformative(args)
        monkeypatch.setattr(
            pi_module,
            "_PARSIMONY_INFORMATIVE_SCALAR_MAX_CELLS",
            0,
        )
        generic_count = mocker.spy(
            pi_module,
            "_count_ascii_parsimony_informative_sites",
        )

        pi_sites, aln_len, pi_sites_per = pi.calculate_parsimony_informative_sites(
            alignment,
            is_protein=False,
        )

        assert generic_count.call_args.args[1] is None
        assert pi_sites == 1
        assert aln_len == 4
        assert isclose(pi_sites_per, 25.0, rel_tol=0.001)

    def test_parsimony_informative_clean_ascii_skips_gap_code_setup(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("ACGT"), id="t1"),
            SeqRecord(Seq("ACGA"), id="t2"),
            SeqRecord(Seq("TCGC"), id="t3"),
            SeqRecord(Seq("TCGT"), id="t4"),
        ])
        pi = ParsimonyInformative(args)
        mocker.patch(
            "phykit.services.alignment.parsimony_informative_sites._get_gap_codes",
            side_effect=AssertionError("clean ASCII path should not build gap codes"),
        )

        pi_sites, aln_len, pi_sites_per = pi.calculate_parsimony_informative_sites(
            alignment,
            is_protein=False,
        )

        assert pi_sites == 1
        assert aln_len == 4
        assert isclose(pi_sites_per, 25.0, rel_tol=0.001)

    def test_parsimony_informative_sites_small_ascii_uses_scalar_path(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("AC-N"), id="t1"),
            SeqRecord(Seq("ACGN"), id="t2"),
            SeqRecord(Seq("TGGN"), id="t3"),
            SeqRecord(Seq("TG?N"), id="t4"),
        ])
        pi = ParsimonyInformative(args)
        mocker.patch(
            "phykit.services.alignment.parsimony_informative_sites.np.frombuffer",
            side_effect=AssertionError(
                "small PI-site counts should avoid NumPy setup"
            ),
        )

        pi_sites, aln_len, pi_sites_per = pi.calculate_parsimony_informative_sites(
            alignment,
            is_protein=False,
        )

        assert pi_sites == 2
        assert aln_len == 4
        assert isclose(pi_sites_per, 50.0, rel_tol=0.001)

    def test_parsimony_informative_scalar_path_rejects_uneven_lengths(self):
        assert pi_module._count_parsimony_informative_sites_scalar(
            ["ACGT", "ACG"],
            4,
            False,
        ) is None

    def test_parsimony_informative_identical_sequences_skip_matrix(
        self, mocker, args
    ):
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("AcGt"), id="t1"),
            SeqRecord(Seq("aCgT"), id="t2"),
            SeqRecord(Seq("ACGT"), id="t3"),
            SeqRecord(Seq("acgt"), id="t4"),
        ])
        pi = ParsimonyInformative(args)
        mocker.patch(
            "phykit.services.alignment.parsimony_informative_sites.np.frombuffer",
            side_effect=AssertionError("identical alignments should skip matrix setup"),
        )

        pi_sites, aln_len, pi_sites_per = pi.calculate_parsimony_informative_sites(
            alignment,
            is_protein=False,
        )

        assert pi_sites == 0
        assert aln_len == 4
        assert pi_sites_per == 0.0

    def test_parsimony_informative_single_record_returns_before_sequence_materialization(
        self, args
    ):
        class UnstringableSequence:
            def __str__(self):
                raise AssertionError(
                    "single-record parsimony informative sites should not inspect sequence"
                )

        class DummyAlignment(list):
            def get_alignment_length(self):
                return 7

        alignment = DummyAlignment([SimpleNamespace(seq=UnstringableSequence())])
        pi = ParsimonyInformative(args)

        assert pi.calculate_parsimony_informative_sites(
            alignment,
            is_protein=False,
        ) == (0, 7, 0.0)

    def test_parsimony_informative_sites_unicode_fallback(
        self, args, mocker, monkeypatch
    ):
        class DummyAlignment(list):
            def get_alignment_length(self):
                return 2

        alignment = DummyAlignment([
            SimpleNamespace(seq="A\u00d1", id="t1"),
            SimpleNamespace(seq="A\u00dd", id="t2"),
            SimpleNamespace(seq="T\u00d1", id="t3"),
            SimpleNamespace(seq="T\u00dd", id="t4"),
        ])
        pi = ParsimonyInformative(args)
        monkeypatch.setattr(
            pi_module,
            "_PARSIMONY_INFORMATIVE_SCALAR_MAX_CELLS",
            0,
        )
        count_nonzero_spy = mocker.spy(pi_module.np, "count_nonzero")

        pi_sites, aln_len, pi_sites_per = pi.calculate_parsimony_informative_sites(
            alignment,
            is_protein=False,
        )

        assert pi_sites == 2
        assert aln_len == 2
        assert isclose(pi_sites_per, 100.0, rel_tol=0.001)
        assert any(
            call.args[0].dtype == bool
            for call in count_nonzero_spy.call_args_list
        )
