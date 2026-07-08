from argparse import Namespace
import importlib
import multiprocessing
import subprocess
import sys
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest
import builtins
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.pairwise_identity import PairwiseIdentity
from phykit.services.alignment.pairwise_identity import _identity_for_identical_sequence
import phykit.services.alignment.pairwise_identity as pairwise_identity_module


@pytest.fixture
def args():
    return Namespace(alignment="/some/path/to/file.fa", verbose=False, exclude_gaps=False)


def _simple_alignment():
    return MultipleSeqAlignment(
        [
            SeqRecord(Seq("ACGT"), id="a"),
            SeqRecord(Seq("AC-T"), id="b"),
            SeqRecord(Seq("TCGT"), id="c"),
        ]
    )


def test_identity_for_identical_sequence_counts_dna_gap_characters():
    assert _identity_for_identical_sequence("ACGTN-?*X", False, True) == pytest.approx(
        4 / 9
    )


def test_identity_for_identical_sequence_dna_ascii_skips_count_loop():
    class NoCountStr(str):
        def count(self, *args, **kwargs):
            raise AssertionError("ASCII DNA gap counting should use byte deletion")

    sequence = NoCountStr("ACGTN-?*X")

    assert _identity_for_identical_sequence(sequence, False, True) == pytest.approx(
        4 / 9
    )


def test_identity_for_identical_sequence_keeps_protein_n_as_nongap():
    assert _identity_for_identical_sequence("ACGTN-?*X", True, True) == pytest.approx(
        5 / 9
    )


def test_all_sequences_identical_does_not_slice():
    class NoSliceList(list):
        def __getitem__(self, key):
            if isinstance(key, slice):
                raise AssertionError("identical sequence scan should not slice")
            return super().__getitem__(key)

    assert pairwise_identity_module._all_sequences_identical(
        NoSliceList(["ACGT", "ACGT", "ACGT"])
    )
    assert not pairwise_identity_module._all_sequences_identical(
        NoSliceList(["ACGT", "TGCA", "ACGT"])
    )
    assert not pairwise_identity_module._all_sequences_identical(
        NoSliceList(["ACGT", "ACGT", "TGCA"])
    )


def test_module_import_does_not_import_scipy_clustering(monkeypatch):
    module_name = "phykit.services.alignment.pairwise_identity"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if (
            name == "scipy.cluster.hierarchy"
            or name.startswith("scipy.cluster.hierarchy.")
            or name == "scipy.spatial.distance"
            or name.startswith("scipy.spatial.distance.")
        ):
            raise AssertionError(
                "pairwise_identity module import should not import SciPy clustering"
            )
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", guarded_import)
    try:
        importlib.import_module(module_name)
    finally:
        imported = sys.modules.pop(module_name, None)
        if previous is not None:
            sys.modules[module_name] = previous
        parent_name, _, child_name = module_name.rpartition(".")
        parent = sys.modules.get(parent_name)
        if parent is not None:
            if previous is not None:
                setattr(parent, child_name, previous)
            elif getattr(parent, child_name, None) is imported:
                delattr(parent, child_name)


def test_module_import_does_not_import_biopython_align_or_tqdm():
    code = """
import sys
import phykit.services.alignment.pairwise_identity as module
assert hasattr(module.np, "__getattr__")
assert callable(module.mp.cpu_count)
assert callable(module.mp.Pool)
assert callable(module.print_json)
assert callable(module.print_summary_statistics)
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "multiprocessing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.stats_summary" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "tqdm" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_multiprocessing_caches_module_and_keeps_cpu_count_patchable():
    lazy_mp = pairwise_identity_module._LazyMultiprocessing()

    with patch.object(multiprocessing, "cpu_count", return_value=11):
        assert lazy_mp.cpu_count() == 11
        assert lazy_mp._module is multiprocessing

    with patch.object(multiprocessing, "cpu_count", return_value=13) as cpu_count:
        assert lazy_mp.cpu_count() == 13

    cpu_count.assert_called_once_with()


def test_lazy_numpy_caches_module_and_attributes():
    lazy_np = pairwise_identity_module._LazyNumpy()

    first_fromiter = lazy_np.fromiter
    second_fromiter = lazy_np.fromiter

    assert lazy_np._module is not None
    assert first_fromiter is second_fromiter
    assert lazy_np.__dict__["fromiter"] is first_fromiter


class TestPairwiseIdentity:
    def test_init_sets_expected_attrs(self, args):
        service = PairwiseIdentity(args)
        assert service.alignment_file_path == args.alignment
        assert service.verbose is False
        assert service.exclude_gaps is False
        assert service.json_output is False
        assert service.plot is False

    def test_calculate_identity_vectorized(self, args):
        service = PairwiseIdentity(args)
        seq1 = np.array(list("ACGT"), dtype="U1")
        seq2 = np.array(list("AC-T"), dtype="U1")
        identity = service._calculate_identity_vectorized(seq1, seq2)
        assert identity == 0.75

    def test_calculate_identity_vectorized_excluding_gaps(self, args):
        service = PairwiseIdentity(args)
        seq1 = np.array(list("ACGT"), dtype="U1")
        seq2 = np.array(list("AC-T"), dtype="U1")
        gap_mask1 = np.isin(seq1, ["-"])
        gap_mask2 = np.isin(seq2, ["-"])
        identity = service._calculate_identity_vectorized(
            seq1, seq2, (gap_mask1, gap_mask2), exclude_gaps=True
        )
        # denominator is still full length by implementation
        assert identity == 0.75

    def test_sequence_to_array_uppercases_ascii_bytes(self, args):
        service = PairwiseIdentity(args)

        seq_array = service._sequence_to_array("ac-g")

        assert seq_array.dtype.kind == "S"
        assert seq_array.tolist() == [b"A", b"C", b"-", b"G"]

    def test_sequence_to_array_falls_back_for_non_ascii(self, args):
        service = PairwiseIdentity(args)

        seq_array = service._sequence_to_array("A\u00d1")

        assert seq_array.dtype.kind == "U"
        assert seq_array.tolist() == ["A", "\u00d1"]

    def test_gap_chars_array_matches_byte_sequences(self, args):
        service = PairwiseIdentity(args)
        seq_array = service._sequence_to_array("AC-T")

        gap_chars_array = service._gap_chars_array_for_sequence(seq_array, {"-"})

        assert gap_chars_array.dtype.kind == "S"
        assert np.isin(seq_array, gap_chars_array).tolist() == [
            False,
            False,
            True,
            False,
        ]

    def test_calculate_identity_vectorized_counts_masked_matches(self, args):
        service = PairwiseIdentity(args)
        seq_one = np.array(list("AC-T"), dtype="U1")
        seq_two = np.array(list("ACGT"), dtype="U1")
        gap_mask_one = np.array([False, False, True, False])
        gap_mask_two = np.array([False, False, False, False])

        identity = service._calculate_identity_vectorized(
            seq_one,
            seq_two,
            gap_mask=(gap_mask_one, gap_mask_two),
            exclude_gaps=True,
        )

        assert identity == 0.75

    def test_process_pair_batch(self, args):
        service = PairwiseIdentity(args)
        alignment_data = [
            {"id": "a", "seq": np.array(list("AAAA"), dtype="U1")},
            {"id": "b", "seq": np.array(list("AAAT"), dtype="U1")},
            {"id": "c", "seq": np.array(list("AATT"), dtype="U1")},
        ]
        results = service._process_pair_batch(
            alignment_data,
            pair_indices=[(0, 1), (0, 2), (1, 2)],
            exclude_gaps=False,
            gap_chars=set(["-"]),
        )
        assert [result["pair_id"] for result in results] == [
            ["a", "b"],
            ["a", "c"],
            ["b", "c"],
        ]
        assert [round(result["identity"], 4) for result in results] == [
            0.75,
            0.5,
            0.75,
        ]

    def test_process_pair_batch_uses_precomputed_gap_masks(self, args, mocker):
        service = PairwiseIdentity(args)
        alignment_data = [
            {
                "id": "a",
                "seq": np.array(list("AC-T"), dtype="U1"),
                "gap_mask": np.array([False, False, True, False]),
            },
            {
                "id": "b",
                "seq": np.array(list("ACGT"), dtype="U1"),
                "gap_mask": np.array([False, False, False, False]),
            },
        ]
        mocked_isin = mocker.patch(
            "phykit.services.alignment.pairwise_identity.np.isin"
        )

        results = service._process_pair_batch(
            alignment_data,
            pair_indices=[(0, 1)],
            exclude_gaps=True,
            gap_chars=set(["-"]),
        )

        mocked_isin.assert_not_called()
        assert results[0]["pair_id"] == ["a", "b"]
        assert round(results[0]["identity"], 4) == 0.75

    def test_calculate_pairwise_identities_sequential(self, mocker, args):
        service = PairwiseIdentity(args)
        mocker.patch.object(PairwiseIdentity, "_should_use_multiprocessing", return_value=False)
        pair_ids, identities, stats = service.calculate_pairwise_identities(
            _simple_alignment(), exclude_gaps=False, is_protein=False
        )
        assert len(pair_ids) == 3
        assert "a-b" in identities
        assert round(identities["a-b"], 4) == 0.75
        assert "mean" in stats

    def test_single_record_pairwise_identities_return_before_sequence_materialization(
        self, args
    ):
        class UnstringableSequence:
            def __str__(self):
                raise AssertionError(
                    "single-record pairwise identity should not inspect sequence"
                )

        alignment = [SimpleNamespace(id="solo", seq=UnstringableSequence())]
        service = PairwiseIdentity(args)

        pair_ids, identities, stats = service.calculate_pairwise_identities(
            alignment,
            exclude_gaps=False,
            is_protein=False,
        )

        assert pair_ids == []
        assert identities == {}
        assert stats is None

    def test_single_record_pairwise_identity_stats_return_before_sequence_materialization(
        self, args
    ):
        class UnstringableSequence:
            def __str__(self):
                raise AssertionError(
                    "single-record pairwise identity stats should not inspect sequence"
                )

        alignment = [SimpleNamespace(id="solo", seq=UnstringableSequence())]
        service = PairwiseIdentity(args)

        stats = service.calculate_pairwise_identity_stats(
            alignment,
            exclude_gaps=False,
            is_protein=False,
        )

        assert stats is None

    def test_calculate_pairwise_identities_sequential_streams_pairs(
        self, mocker, args
    ):
        service = PairwiseIdentity(args)
        mocker.patch.object(
            PairwiseIdentity,
            "_calculate_pairwise_identities_matrix",
            return_value=None,
        )
        should_use_mp = mocker.patch.object(
            PairwiseIdentity, "_should_use_multiprocessing", return_value=False
        )

        def process_batch(_alignment_data, pair_indices, _exclude_gaps, _gap_chars):
            assert not isinstance(pair_indices, list)
            return [
                {"pair_id": ["a", "b"], "identity": 0.75},
                {"pair_id": ["a", "c"], "identity": 0.5},
                {"pair_id": ["b", "c"], "identity": 0.25},
            ]

        process_spy = mocker.patch.object(
            PairwiseIdentity,
            "_process_pair_batch",
            side_effect=process_batch,
        )

        pair_ids, identities, stats = service.calculate_pairwise_identities(
            _simple_alignment(), exclude_gaps=False, is_protein=False
        )

        should_use_mp.assert_called_once_with(3)
        process_spy.assert_called_once()
        assert pair_ids == [["a", "b"], ["a", "c"], ["b", "c"]]
        assert identities == {"a-b": 0.75, "a-c": 0.5, "b-c": 0.25}
        assert stats["mean"] == pytest.approx(0.5)

    def test_batched_pair_indices_streams_all_pairs(self, args):
        service = PairwiseIdentity(args)

        batches = list(service._batched_pair_indices(5, 3))
        observed_pairs = [pair for batch in batches for pair in batch]

        assert [len(batch) for batch in batches] == [3, 3, 3, 1]
        assert observed_pairs == list(
            pairwise_identity_module.itertools.combinations(range(5), 2)
        )

    def test_calculate_pairwise_identities_exclude_gaps_matrix_path(self, mocker, args):
        service = PairwiseIdentity(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AC-T"), id="a"),
                SeqRecord(Seq("ACGT"), id="b"),
                SeqRecord(Seq("TCGT"), id="c"),
            ]
        )
        mocker.patch.object(
            PairwiseIdentity,
            "_process_pair_batch",
            side_effect=AssertionError("matrix path should avoid pair batches"),
        )

        pair_ids, identities, stats = service.calculate_pairwise_identities(
            alignment,
            exclude_gaps=True,
            is_protein=False,
        )

        assert pair_ids == [["a", "b"], ["a", "c"], ["b", "c"]]
        assert identities == {"a-b": 0.75, "a-c": 0.5, "b-c": 0.75}
        assert stats["mean"] == pytest.approx(2 / 3)

    def test_calculate_pairwise_identities_default_matrix_path(self, mocker, args):
        service = PairwiseIdentity(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("A--T"), id="a"),
                SeqRecord(Seq("A--G"), id="b"),
                SeqRecord(Seq("N--G"), id="c"),
            ]
        )
        mocker.patch.object(
            PairwiseIdentity,
            "_process_pair_batch",
            side_effect=AssertionError("matrix path should avoid pair batches"),
        )

        pair_ids, identities, stats = service.calculate_pairwise_identities(
            alignment,
            exclude_gaps=False,
            is_protein=False,
        )

        assert pair_ids == [["a", "b"], ["a", "c"], ["b", "c"]]
        assert identities == {"a-b": 0.75, "a-c": 0.5, "b-c": 0.75}
        assert stats["mean"] == pytest.approx(2 / 3)

    def test_calculate_pairwise_identities_identical_sequences_skip_matrix(
        self, mocker, args
    ):
        service = PairwiseIdentity(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("A--T"), id="a"),
                SeqRecord(Seq("a--t"), id="b"),
                SeqRecord(Seq("A--T"), id="c"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.pairwise_identity.np.frombuffer",
            side_effect=AssertionError("identical alignments should skip matrix setup"),
        )

        pair_ids, identities, stats = service.calculate_pairwise_identities(
            alignment,
            exclude_gaps=True,
            is_protein=False,
        )

        assert pair_ids == [["a", "b"], ["a", "c"], ["b", "c"]]
        assert identities == {"a-b": 0.5, "a-c": 0.5, "b-c": 0.5}
        assert stats == {
            "mean": 0.5,
            "median": 0.5,
            "twenty_fifth": 0.5,
            "seventy_fifth": 0.5,
            "minimum": 0.5,
            "maximum": 0.5,
            "standard_deviation": 0.0,
            "variance": 0.0,
        }

    def test_calculate_pairwise_identity_stats_identical_sequences_skip_matrix(
        self, mocker, args
    ):
        service = PairwiseIdentity(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("A--T"), id="a"),
                SeqRecord(Seq("a--t"), id="b"),
                SeqRecord(Seq("A--T"), id="c"),
            ]
        )
        mocker.patch(
            "phykit.services.alignment.pairwise_identity.np.frombuffer",
            side_effect=AssertionError("identical alignments should skip matrix setup"),
        )

        stats = service.calculate_pairwise_identity_stats(
            alignment,
            exclude_gaps=False,
            is_protein=False,
        )

        assert stats == {
            "mean": 1.0,
            "median": 1.0,
            "twenty_fifth": 1.0,
            "seventy_fifth": 1.0,
            "minimum": 1.0,
            "maximum": 1.0,
            "standard_deviation": 0.0,
            "variance": 0.0,
        }

    def test_exclude_gaps_clean_ascii_matrix_path_skips_gap_lookup(self, mocker, args):
        service = PairwiseIdentity(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="a"),
                SeqRecord(Seq("ACGA"), id="b"),
                SeqRecord(Seq("TCGT"), id="c"),
            ]
        )
        mocker.patch.object(
            PairwiseIdentity,
            "_get_gap_lookup",
            side_effect=AssertionError(
                "clean ASCII exclude-gaps matrix path should skip gap lookup"
            ),
        )
        mocker.patch.object(
            PairwiseIdentity,
            "_process_pair_batch",
            side_effect=AssertionError("matrix path should avoid pair batches"),
        )

        pair_ids, identities, stats = service.calculate_pairwise_identities(
            alignment,
            exclude_gaps=True,
            is_protein=False,
        )
        stats_only = service.calculate_pairwise_identity_stats(
            alignment,
            exclude_gaps=True,
            is_protein=False,
        )

        assert pair_ids == [["a", "b"], ["a", "c"], ["b", "c"]]
        assert identities == {"a-b": 0.75, "a-c": 0.75, "b-c": 0.5}
        assert stats_only == pytest.approx(stats)

    def test_calculate_pairwise_identity_stats_matrix_matches_full_result(self, args):
        service = PairwiseIdentity(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("A--T"), id="a"),
                SeqRecord(Seq("A--G"), id="b"),
                SeqRecord(Seq("N--G"), id="c"),
                SeqRecord(Seq("A--T"), id="d"),
            ]
        )

        full_stats = service.calculate_pairwise_identities(
            alignment,
            exclude_gaps=False,
            is_protein=False,
        )[2]
        stats_only = service.calculate_pairwise_identity_stats(
            alignment,
            exclude_gaps=False,
            is_protein=False,
        )

        assert stats_only == pytest.approx(full_stats)

    def test_calculate_pairwise_identities_matrix_summarizes_condensed_values(
        self, args, mocker
    ):
        service = PairwiseIdentity(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="a"),
                SeqRecord(Seq("ACGA"), id="b"),
                SeqRecord(Seq("TCGA"), id="c"),
            ]
        )
        mocked_arr_stats = mocker.spy(
            pairwise_identity_module,
            "calculate_summary_statistics_from_arr",
        )
        mocker.patch(
            "phykit.services.alignment.pairwise_identity.calculate_summary_statistics_from_dict",
            side_effect=AssertionError(
                "matrix output should summarize condensed identity values"
            ),
        )

        pair_ids, identities, stats = service.calculate_pairwise_identities(
            alignment,
            exclude_gaps=False,
            is_protein=False,
        )

        assert pair_ids == [["a", "b"], ["a", "c"], ["b", "c"]]
        assert identities == {"a-b": 0.75, "a-c": 0.5, "b-c": 0.75}
        mocked_arr_stats.assert_called_once()
        assert stats["mean"] == pytest.approx(2 / 3)

    def test_calculate_pairwise_identity_stats_gappy_matrix_uses_squareform(
        self, args, mocker
    ):
        service = PairwiseIdentity(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("A--TN"), id="a"),
                SeqRecord(Seq("A--GN"), id="b"),
                SeqRecord(Seq("N--GA"), id="c"),
                SeqRecord(Seq("A--TA"), id="d"),
            ]
        )
        full_stats = service.calculate_pairwise_identities(
            alignment,
            exclude_gaps=True,
            is_protein=False,
        )[2]
        mocked_squareform = mocker.spy(pairwise_identity_module, "squareform")
        mocker.patch(
            "phykit.services.alignment.pairwise_identity.np.triu_indices",
            side_effect=AssertionError(
                "gappy stats-only matrix path should use condensed squareform"
            ),
        )

        stats_only = service.calculate_pairwise_identity_stats(
            alignment,
            exclude_gaps=True,
            is_protein=False,
        )

        mocked_squareform.assert_called()
        assert stats_only == pytest.approx(full_stats)

    def test_calculate_pairwise_identity_stats_matrix_matches_full_result_multiblock(
        self, args
    ):
        service = PairwiseIdentity(args)
        alphabet = "ACGT"
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(
                    Seq(
                        "".join(
                            alphabet[(idx + pos) % len(alphabet)]
                            for pos in range(16)
                        )
                    ),
                    id=f"taxon{idx}",
                )
                for idx in range(70)
            ]
        )

        full_stats = service.calculate_pairwise_identities(
            alignment,
            exclude_gaps=False,
            is_protein=False,
        )[2]
        stats_only = service.calculate_pairwise_identity_stats(
            alignment,
            exclude_gaps=False,
            is_protein=False,
        )

        assert stats_only == pytest.approx(full_stats)

    def test_run_summary_uses_stats_only_path(self, mocker):
        args = Namespace(
            alignment="/some/path/to/file.fa",
            verbose=False,
            exclude_gaps=False,
            json=False,
            plot=False,
        )
        service = PairwiseIdentity(args)
        alignment = _simple_alignment()
        mocker.patch.object(
            PairwiseIdentity,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        stats_only = mocker.patch.object(
            service,
            "calculate_pairwise_identity_stats",
            return_value={"mean": 0.75},
        )
        full_result = mocker.patch.object(
            service,
            "calculate_pairwise_identities",
            side_effect=AssertionError("summary run should not build pair rows"),
        )
        mocked_summary = mocker.patch(
            "phykit.services.alignment.pairwise_identity.print_summary_statistics"
        )

        service.run()

        stats_only.assert_called_once_with(alignment, False, False)
        full_result.assert_not_called()
        mocked_summary.assert_called_once_with({"mean": 0.75})

    def test_run_summary_only_does_not_materialize_taxa_ids(self, mocker):
        class RecordWithExpensiveId:
            @property
            def id(self):
                raise AssertionError("summary-only output should not read taxa IDs")

        args = Namespace(
            alignment="/some/path/to/file.fa",
            verbose=False,
            exclude_gaps=False,
            json=False,
            plot=False,
        )
        service = PairwiseIdentity(args)
        alignment = [RecordWithExpensiveId()]
        mocker.patch.object(
            PairwiseIdentity,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocker.patch.object(
            service,
            "calculate_pairwise_identity_stats",
            return_value={"mean": 0.75},
        )
        mocker.patch(
            "phykit.services.alignment.pairwise_identity.print_summary_statistics"
        )

        service.run()

    def test_matrix_pairwise_identity_matches_batch_reference(self, mocker, args):
        service = PairwiseIdentity(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("AnXT"), id="a"),
                SeqRecord(Seq("ACxT"), id="b"),
                SeqRecord(Seq("AG-T"), id="c"),
            ]
        )
        mocker.patch.object(PairwiseIdentity, "_should_use_multiprocessing", return_value=False)

        matrix_result = service._calculate_pairwise_identities_matrix(
            alignment,
            is_protein=False,
            exclude_gaps=True,
        )
        mocker.patch.object(
            PairwiseIdentity,
            "_calculate_pairwise_identities_matrix",
            return_value=None,
        )
        batch_result = service.calculate_pairwise_identities(
            alignment,
            exclude_gaps=True,
            is_protein=False,
        )

        assert matrix_result == batch_result

    def test_run_json_verbose_and_plot(self, mocker):
        args = Namespace(
            alignment="/some/path/to/file.fa",
            verbose=True,
            exclude_gaps=False,
            json=True,
            plot=True,
            plot_output="out.png",
        )
        service = PairwiseIdentity(args)
        alignment = _simple_alignment()
        mocker.patch.object(
            PairwiseIdentity,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocker.patch.object(
            PairwiseIdentity,
            "calculate_pairwise_identities",
            return_value=(
                [["a", "b"]],
                {"a-b": 0.75},
                {"mean": 0.75},
            ),
        )
        mocked_plot = mocker.patch.object(PairwiseIdentity, "_plot_pairwise_identity_heatmap")
        mocked_json = mocker.patch(
            "phykit.services.alignment.pairwise_identity.print_json"
        )

        service.run()

        mocked_plot.assert_called_once()
        payload = mocked_json.call_args.args[0]
        assert payload["verbose"] is True
        assert payload["pairs"] == payload["rows"]
        assert payload["plot_output"] == "out.png"
        assert payload["rows"][0] == {"taxon_a": "a", "taxon_b": "b", "identity": 0.75}

    def test_run_summary_prints_and_plot_message(self, mocker):
        args = Namespace(
            alignment="/some/path/to/file.fa",
            verbose=False,
            exclude_gaps=False,
            json=False,
            plot=True,
            plot_output="out.png",
        )
        service = PairwiseIdentity(args)
        alignment = _simple_alignment()
        mocker.patch.object(
            PairwiseIdentity,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocker.patch.object(
            PairwiseIdentity,
            "calculate_pairwise_identities",
            return_value=(
                [["a", "b"]],
                {"a-b": 0.75},
                {"mean": 0.75},
            ),
        )
        mocker.patch.object(PairwiseIdentity, "_plot_pairwise_identity_heatmap")
        mocked_summary = mocker.patch(
            "phykit.services.alignment.pairwise_identity.print_summary_statistics"
        )
        mocked_print = mocker.patch("builtins.print")

        service.run()
        mocked_summary.assert_called_once_with({"mean": 0.75})
        mocked_print.assert_called_once_with("Saved pairwise identity heatmap: out.png")

    def test_run_verbose_batches_text_output(self, mocker, capsys):
        args = Namespace(
            alignment="/some/path/to/file.fa",
            verbose=True,
            exclude_gaps=False,
            json=False,
            plot=True,
            plot_output="out.png",
        )
        service = PairwiseIdentity(args)
        alignment = _simple_alignment()
        mocker.patch.object(
            PairwiseIdentity,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocker.patch.object(
            PairwiseIdentity,
            "calculate_pairwise_identities",
            return_value=(
                [["a", "b"], ["a", "c"]],
                {"a-b": 0.75, "a-c": 0.5},
                {"mean": 0.625},
            ),
        )
        mocker.patch.object(PairwiseIdentity, "_plot_pairwise_identity_heatmap")

        service.run()

        out, _ = capsys.readouterr()
        assert out == "a\tb\t0.75\na\tc\t0.5\nSaved pairwise identity heatmap: out.png\n"

    def test_should_use_multiprocessing_respects_env(self, args, monkeypatch):
        service = PairwiseIdentity(args)
        monkeypatch.setenv("PHYKIT_DISABLE_MP", "1")
        assert service._should_use_multiprocessing(99999) is False
        monkeypatch.delenv("PHYKIT_DISABLE_MP")
        monkeypatch.setenv("PHYKIT_FORCE_MP", "1")
        assert service._should_use_multiprocessing(1) is True
        monkeypatch.delenv("PHYKIT_FORCE_MP")
        assert service._should_use_multiprocessing(service.MP_MIN_PAIRS - 1) is False
        assert service._should_use_multiprocessing(service.MP_MIN_PAIRS) is True

    def test_process_args_defaults(self):
        parsed = PairwiseIdentity(
            Namespace(alignment="x.fa", verbose=False, exclude_gaps=False)
        ).process_args(
            Namespace(alignment="x.fa", verbose=False, exclude_gaps=False)
        )
        assert parsed["json_output"] is False
        assert parsed["plot"] is False
        assert parsed["plot_output"] == "pairwise_identity_heatmap.png"

    def test_print_json_output_non_verbose(self, mocker):
        service = PairwiseIdentity(
            Namespace(alignment="x.fa", verbose=False, exclude_gaps=True, json=True, plot=False)
        )
        mocked_json = mocker.patch("phykit.services.alignment.pairwise_identity.print_json")
        service._print_json_output(pair_ids=[["a", "b"]], pairwise_identities={"a-b": 0.5}, stats={"mean": 0.5})
        payload = mocked_json.call_args.args[0]
        assert payload["verbose"] is False
        assert payload["exclude_gaps"] is True
        assert payload["summary"] == {"mean": 0.5}

    def test_pairwise_identity_matrix_from_canonical_pairs_uses_squareform(self, mocker):
        taxa = ["a", "b", "c"]
        pair_ids = [["a", "b"], ["a", "c"], ["b", "c"]]
        identities = {"a-b": 0.8, "a-c": 0.5, "b-c": 0.6}
        mocked_squareform = mocker.spy(pairwise_identity_module, "squareform")
        mocker.patch(
            "phykit.services.alignment.pairwise_identity.np.triu_indices",
            side_effect=AssertionError(
                "canonical pair order should use condensed squareform fill"
            ),
        )

        matrix = pairwise_identity_module._pairwise_identity_matrix_from_pairs(
            taxa,
            pair_ids,
            identities,
        )

        mocked_squareform.assert_called_once()
        np.testing.assert_allclose(
            matrix,
            np.array(
                [
                    [1.0, 0.8, 0.5],
                    [0.8, 1.0, 0.6],
                    [0.5, 0.6, 1.0],
                ],
                dtype=np.float32,
            ),
        )

    def test_pairwise_identity_matrix_from_arbitrary_pairs_uses_fallback_order(self, mocker):
        taxa = ["a", "b", "c"]
        pair_ids = [["b", "c"], ["a", "b"], ["a", "c"]]
        identities = {"b-c": 0.6, "a-b": 0.8, "a-c": 0.5}
        mocker.patch(
            "phykit.services.alignment.pairwise_identity.np.triu_indices",
            side_effect=AssertionError("arbitrary pair order should use fallback fill"),
        )

        matrix = pairwise_identity_module._pairwise_identity_matrix_from_pairs(
            taxa,
            pair_ids,
            identities,
        )

        np.testing.assert_allclose(
            matrix,
            np.array(
                [
                    [1.0, 0.8, 0.5],
                    [0.8, 1.0, 0.6],
                    [0.5, 0.6, 1.0],
                ],
                dtype=np.float32,
            ),
        )

    def test_plot_pairwise_identity_heatmap_creates_file(self, tmp_path):
        pytest.importorskip("matplotlib")
        output = tmp_path / "pairwise.png"
        service = PairwiseIdentity(
            Namespace(alignment="x.fa", verbose=False, exclude_gaps=False, plot=True, plot_output=str(output))
        )
        service._plot_pairwise_identity_heatmap(
            taxa=["a", "b", "c"],
            pair_ids=[["a", "b"], ["a", "c"], ["b", "c"]],
            pairwise_identities={"a-b": 0.8, "a-c": 0.5, "b-c": 0.6},
        )
        assert output.exists()

    def test_plot_pairwise_identity_skips_redundant_tight_layout(
        self, tmp_path, monkeypatch
    ):
        pytest.importorskip("matplotlib")
        from matplotlib.figure import Figure

        output = tmp_path / "pairwise.png"
        service = PairwiseIdentity(
            Namespace(alignment="x.fa", verbose=False, exclude_gaps=False, plot=True, plot_output=str(output))
        )

        def fail_tight_layout(self, *args, **kwargs):
            raise AssertionError("bbox_inches='tight' handles saved bounds")

        monkeypatch.setattr(Figure, "tight_layout", fail_tight_layout)
        service._plot_pairwise_identity_heatmap(
            taxa=["a", "b", "c"],
            pair_ids=[["a", "b"], ["a", "c"], ["b", "c"]],
            pairwise_identities={"a-b": 0.8, "a-c": 0.5, "b-c": 0.6},
        )

        assert output.exists()

    def test_plot_pairwise_identity_heatmap_importerror(self, monkeypatch, capsys):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name.startswith("matplotlib"):
                raise ImportError("no matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        service = PairwiseIdentity(
            Namespace(alignment="x.fa", verbose=False, exclude_gaps=False)
        )
        with pytest.raises(SystemExit) as exc:
            service._plot_pairwise_identity_heatmap(["a"], [], {})
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "matplotlib is required for --plot in pairwise_identity" in out

    def test_plot_pairwise_identity_heatmap_empty_taxa_returns(self):
        pytest.importorskip("matplotlib")
        service = PairwiseIdentity(
            Namespace(alignment="x.fa", verbose=False, exclude_gaps=False)
        )
        service._plot_pairwise_identity_heatmap([], [], {})

    def test_run_verbose_handles_broken_pipe(self, mocker):
        service = PairwiseIdentity(
            Namespace(alignment="x.fa", verbose=True, exclude_gaps=False, json=False, plot=False)
        )
        alignment = _simple_alignment()
        mocker.patch.object(PairwiseIdentity, "get_alignment_and_format", return_value=(alignment, "fasta", False))
        mocker.patch.object(
            PairwiseIdentity,
            "calculate_pairwise_identities",
            return_value=([["a", "b"]], {"a-b": 0.75}, {"mean": 0.75}),
        )
        mocker.patch("builtins.print", side_effect=BrokenPipeError)
        service.run()

    def test_calculate_pairwise_identities_multiprocessing_branch(self, mocker, args, monkeypatch):
        class FakePool:
            def __init__(self, processes):
                self.processes = processes

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

            def map(self, process_func, chunks):
                return [process_func(chunk) for chunk in chunks]

            def imap(self, process_func, chunks):
                for chunk in chunks:
                    yield process_func(chunk)

        service = PairwiseIdentity(args)
        mocker.patch.object(PairwiseIdentity, "_should_use_multiprocessing", return_value=True)
        mocker.patch("phykit.services.alignment.pairwise_identity.mp.Pool", FakePool)
        mocker.patch("phykit.services.alignment.pairwise_identity.mp.cpu_count", return_value=2)
        monkeypatch.setattr(pairwise_identity_module.sys.stderr, "isatty", lambda: False)

        pair_ids, identities, stats = service.calculate_pairwise_identities(
            _simple_alignment(),
            exclude_gaps=False,
            is_protein=False,
        )
        assert len(pair_ids) == 3
        assert "a-b" in identities
        assert "mean" in stats
