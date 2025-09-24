import itertools
import pytest
from argparse import Namespace
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.sum_of_pairs_score import SumOfPairsScore


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

    def test_determine_matches_parallel_path(self, mocker, args):
        sop = SumOfPairsScore(args)
        ids = [f"id{i}" for i in range(12)]
        reference_records = _make_records({key: "AAAAA" for key in ids})
        query_records = _make_records({key: "AAAAA" for key in ids})
        record_id_pairs = list(itertools.combinations(ids, 2))

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
        assert matches == 5 * len(record_id_pairs)
        assert pairs == 5 * len(record_id_pairs)
