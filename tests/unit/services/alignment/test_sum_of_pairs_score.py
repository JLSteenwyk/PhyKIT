import itertools
import subprocess
import sys
import pytest
from argparse import Namespace
from types import SimpleNamespace
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.sum_of_pairs_score import SumOfPairsScore
import phykit.services.alignment.sum_of_pairs_score as module


def test_module_import_does_not_import_biopython_fasta_parser():
    code = """
import sys
import phykit.services.alignment.sum_of_pairs_score as module
assert hasattr(module.np, "__getattr__")
assert hasattr(module.mp, "Pool")
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "multiprocessing" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def args():
    return Namespace(fasta="query.fasta", reference="reference.fasta")


def _make_records(seq_map):
    return {key: SeqRecord(Seq(seq), id=key) for key, seq in seq_map.items()}


class TestSumOfPairsScore:
    def test_process_args(self, args):
        sop = SumOfPairsScore(args)
        assert sop.fasta == args.fasta
        assert sop.reference == args.reference

    def test_read_fasta_uses_first_header_token_and_preserves_case(self, tmp_path):
        path = tmp_path / "alignment.fa"
        path.write_text(">id1 description\nAa -\nC\n>id2\ncg TA\n")

        records = SumOfPairsScore._read_fasta(str(path))

        assert records == {"id1": "Aa-C", "id2": "cgTA"}

    def test_read_fasta_raises_for_duplicate_ids(self, tmp_path):
        path = tmp_path / "alignment.fa"
        path.write_text(">id1 description\nAAAA\n>id1 other\nCCCC\n")

        with pytest.raises(ValueError, match="Duplicate key 'id1'"):
            SumOfPairsScore._read_fasta(str(path))

    def test_run_complete_equal_lengths_skips_pair_list_setup(
        self, mocker, args, capsys
    ):
        sop = SumOfPairsScore(args)
        reference_records = _make_records({
            "id1": "AAAA",
            "id2": "AAAA",
            "id3": "AAAA",
        })
        query_records = _make_records({
            "id1": "AAAA",
            "id2": "AAAA",
            "id3": "AAAA",
        })
        mocker.patch.object(
            SumOfPairsScore,
            "_read_fasta",
            side_effect=[query_records, reference_records],
        )
        mocked_determine = mocker.patch.object(
            sop,
            "determine_number_of_matches_and_total_pairs",
            side_effect=AssertionError("complete run path should not build pair list"),
        )

        sop.run()

        assert capsys.readouterr().out == "1.0\n"
        mocked_determine.assert_not_called()

    def test_run_mixed_lengths_uses_pair_list_fallback(self, mocker, args, capsys):
        sop = SumOfPairsScore(args)
        reference_records = _make_records({
            "id1": "AAAA",
            "id2": "AAAA",
            "id3": "AAAA",
        })
        query_records = _make_records({
            "id1": "AAA",
            "id2": "AAAA",
            "id3": "AAAA",
        })
        mocker.patch.object(
            SumOfPairsScore,
            "_read_fasta",
            side_effect=[query_records, reference_records],
        )
        mocked_determine = mocker.patch.object(
            sop,
            "determine_number_of_matches_and_total_pairs",
            return_value=(1, 2),
        )

        sop.run()

        assert capsys.readouterr().out == "0.5\n"
        mocked_determine.assert_called_once_with(
            list(itertools.combinations(reference_records.keys(), 2)),
            reference_records,
            query_records,
        )

    def test_determine_matches_sequential_equal_lengths(self, args):
        sop = SumOfPairsScore(args)
        ids = ["id1", "id2", "id3"]
        reference_records = _make_records({key: "AAAA" for key in ids})
        query_records = _make_records({key: "AAAA" for key in ids})
        record_id_pairs = list(itertools.combinations(ids, 2))

        matches, pairs = sop.determine_number_of_matches_and_total_pairs(
            record_id_pairs, reference_records, query_records
        )

        assert matches == 4 * len(record_id_pairs)
        assert pairs == 4 * len(record_id_pairs)

    def test_determine_matches_sequential_equal_lengths_caches_arrays(
        self, mocker, args
    ):
        sop = SumOfPairsScore(args)
        ids = [f"id{i}" for i in range(8)]
        reference_records = _make_records({key: "AAAAA" for key in ids})
        query_records = _make_records({
            key: "AAAAT" if key == "id0" else "AAAAA"
            for key in ids
        })
        record_id_pairs = list(itertools.combinations(ids, 2))[:-1]
        array_spy = mocker.spy(SumOfPairsScore, "_sequence_arrays")

        matches, pairs = sop.determine_number_of_matches_and_total_pairs(
            record_id_pairs, reference_records, query_records
        )

        pairs_with_changed_taxon = sum("id0" in pair for pair in record_id_pairs)
        assert matches == (
            5 * (len(record_id_pairs) - pairs_with_changed_taxon)
            + 4 * pairs_with_changed_taxon
        )
        assert pairs == 5 * len(record_id_pairs)
        assert array_spy.call_count == len(reference_records)

    def test_complete_pair_set_ordered_path_does_not_slice_record_ids(self):
        class NoSliceIds:
            def __init__(self, values):
                self.values = values

            def __len__(self):
                return len(self.values)

            def __iter__(self):
                return iter(self.values)

            def __getitem__(self, item):
                if isinstance(item, slice):
                    raise AssertionError("ordered pair validation should not slice ids")
                return self.values[item]

        ids = [f"id{i}" for i in range(8)]
        record_id_pairs = list(itertools.combinations(ids, 2))

        assert (
            SumOfPairsScore._has_complete_pair_set(record_id_pairs, NoSliceIds(ids))
            is True
        )

    def test_complete_pair_set_unordered_fallback(self):
        ids = [f"id{i}" for i in range(6)]
        record_id_pairs = list(reversed(list(itertools.combinations(ids, 2))))

        assert SumOfPairsScore._has_complete_pair_set(record_id_pairs, ids) is True
        assert (
            SumOfPairsScore._has_complete_pair_set(record_id_pairs[:-1], ids)
            is False
        )

    def test_determine_matches_complete_equal_lengths_uses_fast_path(self, mocker, args):
        sop = SumOfPairsScore(args)
        ids = [f"id{i}" for i in range(12)]
        reference_records = _make_records({
            key: "AAAA" if idx % 2 == 0 else "CCCC"
            for idx, key in enumerate(ids)
        })
        query_records = _make_records({
            key: "AAAA" if idx % 3 == 0 else str(reference_records[key].seq)
            for idx, key in enumerate(ids)
        })
        record_id_pairs = list(itertools.combinations(ids, 2))
        mocked_pool = mocker.patch(
            "phykit.services.alignment.sum_of_pairs_score.mp.Pool"
        )

        matches, pairs = sop.determine_number_of_matches_and_total_pairs(
            record_id_pairs, reference_records, query_records
        )

        matching_taxa = sum(
            str(reference_records[key].seq) == str(query_records[key].seq)
            for key in ids
        )
        assert matches == 4 * (matching_taxa * (matching_taxa - 1) // 2)
        assert pairs == 4 * len(record_id_pairs)
        mocked_pool.assert_not_called()

    def test_complete_equal_lengths_counts_site_matches_with_count_nonzero(
        self, monkeypatch
    ):
        reference_records = _make_records({
            "id1": "AAAA",
            "id2": "ACAC",
            "id3": "GCAT",
        })
        query_records = _make_records({
            "id1": "AATA",
            "id2": "ATAC",
            "id3": "GCAA",
        })
        original_count_nonzero = module.np.count_nonzero
        count_nonzero_axes = []

        def counting_count_nonzero(values, axis=None):
            count_nonzero_axes.append(axis)
            return original_count_nonzero(values, axis=axis)

        monkeypatch.setattr(module.np, "count_nonzero", counting_count_nonzero)
        monkeypatch.setattr(
            module.np,
            "sum",
            lambda *args, **kwargs: pytest.fail(
                "complete equal-length fast path should use ndarray.sum"
            ),
        )

        matches, pairs = SumOfPairsScore._calculate_equal_length_complete_records(
            reference_records,
            query_records,
        )

        assert matches == 6
        assert pairs == 12
        assert count_nonzero_axes == [0]

    def test_complete_equal_lengths_identical_records_skip_matrix_stack(
        self, mocker, args
    ):
        sop = SumOfPairsScore(args)
        ids = [f"id{i}" for i in range(12)]
        reference_records = _make_records({key: "ACGT" for key in ids})
        query_records = _make_records({key: "ACGT" for key in ids})
        record_id_pairs = list(itertools.combinations(ids, 2))
        mocked_stack = mocker.patch.object(
            SumOfPairsScore,
            "_stack_equal_length_sequence_pairs",
            side_effect=AssertionError("identical records should not build matrices"),
        )

        matches, pairs = sop.determine_number_of_matches_and_total_pairs(
            record_id_pairs, reference_records, query_records
        )

        assert matches == 4 * len(record_id_pairs)
        assert pairs == 4 * len(record_id_pairs)
        mocked_stack.assert_not_called()

    def test_incomplete_unchanged_pair_set_skips_sequence_arrays(
        self, mocker, args
    ):
        sop = SumOfPairsScore(args)
        reference_records = _make_records({
            "id1": "AAAA",
            "id2": "AAA",
            "id3": "AA",
        })
        query_records = _make_records({
            "id1": "AAAA",
            "id2": "AAA",
            "id3": "AA",
        })
        record_id_pairs = [("id1", "id2"), ("id1", "id3")]
        mocked_arrays = mocker.patch.object(
            SumOfPairsScore,
            "_sequence_arrays",
            side_effect=AssertionError("unchanged pair set should skip arrays"),
        )
        mocked_pool = mocker.patch(
            "phykit.services.alignment.sum_of_pairs_score.mp.Pool",
            side_effect=AssertionError("unchanged pair set should skip pool"),
        )

        matches, pairs = sop.determine_number_of_matches_and_total_pairs(
            record_id_pairs, reference_records, query_records
        )

        assert matches == 5
        assert pairs == 5
        mocked_arrays.assert_not_called()
        mocked_pool.assert_not_called()

    def test_determine_matches_sequential_mismatched_lengths(self, args):
        sop = SumOfPairsScore(args)
        reference_records = _make_records({
            "id1": "ABCD",
            "id2": "WXYZ",
        })
        query_records = _make_records({
            "id1": "ABX",
            "id2": "WXY",
        })
        record_id_pairs = [("id1", "id2")]

        matches, pairs = sop.determine_number_of_matches_and_total_pairs(
            record_id_pairs, reference_records, query_records
        )

        assert matches == 2
        assert pairs == 3

    def test_determine_matches_sequential_mismatched_lengths_caches_arrays(
        self, mocker, args
    ):
        sop = SumOfPairsScore(args)
        reference_records = _make_records({
            "id1": "AAAAA",
            "id2": "AAAAC",
            "id3": "GGGGG",
        })
        query_records = _make_records({
            "id1": "AAAA",
            "id2": "AAAT",
            "id3": "GGGA",
        })
        record_id_pairs = list(itertools.combinations(reference_records, 2))
        array_spy = mocker.spy(SumOfPairsScore, "_sequence_arrays")

        matches, pairs = sop.determine_number_of_matches_and_total_pairs(
            record_id_pairs, reference_records, query_records
        )

        assert matches == 9
        assert pairs == 12
        assert array_spy.call_count == len(reference_records)

    def test_determine_matches_sequential_mismatched_lengths_non_ascii(self, args):
        sop = SumOfPairsScore(args)
        reference_records = {
            "id1": SimpleNamespace(seq="AΩCT"),
            "id2": SimpleNamespace(seq="AΩGT"),
        }
        query_records = {
            "id1": SimpleNamespace(seq="AΩA"),
            "id2": SimpleNamespace(seq="ATG"),
        }
        record_id_pairs = [("id1", "id2")]

        matches, pairs = sop.determine_number_of_matches_and_total_pairs(
            record_id_pairs, reference_records, query_records
        )

        assert matches == 1
        assert pairs == 3

    def test_process_pair_batch_mixed_lengths(self, args):
        reference_records = _make_records({
            "id1": "AAAA",
            "id2": "AAAA",
            "id3": "ABCD",
            "id4": "WXYZ",
        })
        query_records = _make_records({
            "id1": "AAAA",
            "id2": "AAAA",
            "id3": "ABX",
            "id4": "WXY",
        })
        pair_batch = [("id1", "id2"), ("id3", "id4")]

        matches, pairs = SumOfPairsScore._process_pair_batch(
            pair_batch, reference_records, query_records
        )

        assert matches == 6
        assert pairs == 7

    def test_process_pair_batch_matches_scalar_reference(self, args):
        reference_records = _make_records({
            "id1": "ACGT",
            "id2": "ACG",
            "id3": "A-GT",
        })
        query_records = _make_records({
            "id1": "ACGA",
            "id2": "ATG",
            "id3": "A-G",
        })
        pair_batch = [("id1", "id2"), ("id1", "id3"), ("id2", "id3")]

        expected_matches = 0
        expected_pairs = 0
        for first, second in pair_batch:
            ref1 = str(reference_records[first].seq)
            ref2 = str(reference_records[second].seq)
            query1 = str(query_records[first].seq)
            query2 = str(query_records[second].seq)
            for i in range(min(len(ref1), len(ref2), len(query1), len(query2))):
                if ref1[i] == query1[i] and ref2[i] == query2[i]:
                    expected_matches += 1
                expected_pairs += 1

        matches, pairs = SumOfPairsScore._process_pair_batch(
            pair_batch, reference_records, query_records
        )

        assert matches == expected_matches
        assert pairs == expected_pairs

    def test_process_pair_batch_non_ascii_fallback_matches_scalar_reference(self, args):
        reference_records = {
            "id1": SimpleNamespace(seq="AΩC"),
            "id2": SimpleNamespace(seq="ATC"),
        }
        query_records = {
            "id1": SimpleNamespace(seq="AΩA"),
            "id2": SimpleNamespace(seq="AAC"),
        }
        pair_batch = [("id1", "id2")]

        matches, pairs = SumOfPairsScore._process_pair_batch(
            pair_batch, reference_records, query_records
        )

        assert matches == 1
        assert pairs == 3

    def test_determine_matches_parallel_path(self, mocker, args):
        sop = SumOfPairsScore(args)
        ids = [f"id{i}" for i in range(13)]
        reference_records = _make_records({key: "AAAAA" for key in ids})
        query_records = _make_records({
            key: "AAAAT" if key == "id0" else "AAAAA"
            for key in ids
        })
        record_id_pairs = list(itertools.combinations(ids, 2))[:-1]

        created_pools = []
        observed_batches = []

        class DummyPool:
            def __init__(self, processes):
                self.processes = processes
                created_pools.append(self)

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

            def map(self, func, batches):
                observed_batches.extend(batches)
                return [func(batch) for batch in batches]

        mocker.patch("phykit.services.alignment.sum_of_pairs_score.mp.Pool", DummyPool)
        mocker.patch("phykit.services.alignment.sum_of_pairs_score.mp.cpu_count", return_value=4)

        matches, pairs = sop.determine_number_of_matches_and_total_pairs(
            record_id_pairs, reference_records, query_records
        )

        assert created_pools
        assert all(batch for batch in observed_batches)
        pairs_with_changed_taxon = sum("id0" in pair for pair in record_id_pairs)
        assert matches == (
            5 * (len(record_id_pairs) - pairs_with_changed_taxon)
            + 4 * pairs_with_changed_taxon
        )
        assert pairs == 5 * len(record_id_pairs)
