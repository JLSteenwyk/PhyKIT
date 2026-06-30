from __future__ import annotations

import itertools

from ._fasta import read_unique_fasta_first_token
from .base import Alignment


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyMultiprocessing:
    def cpu_count(self):
        import multiprocessing as _mp

        return _mp.cpu_count()

    def Pool(self, *args, **kwargs):
        import multiprocessing as _mp

        return _mp.Pool(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


mp = _LazyMultiprocessing()
np = _LazyNumpy()


class SumOfPairsScore(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(fasta=parsed["fasta"], reference=parsed["reference"])
        self.json_output = parsed["json_output"]

    def run(self):
        query_records = self._read_fasta(self.fasta)
        reference_records = self._read_fasta(self.reference)

        fast_result = self._calculate_equal_length_complete_records(
            reference_records, query_records
        )
        if fast_result is None:
            record_id_pairs = list(
                itertools.combinations(reference_records.keys(), 2)
            )
            number_of_matches, number_of_total_pairs = \
                self.determine_number_of_matches_and_total_pairs(
                    record_id_pairs, reference_records, query_records
                )
        else:
            number_of_matches, number_of_total_pairs = fast_result

        score = round(number_of_matches / number_of_total_pairs, 4)

        if self.json_output:
            print_json(dict(sum_of_pairs_score=score))
            return

        print(score)

    def process_args(self, args) -> dict[str, str]:
        return dict(
            fasta=args.fasta,
            reference=args.reference,
            json_output=getattr(args, "json", False),
        )

    @staticmethod
    def _read_fasta(path: str) -> dict[str, str]:
        return read_unique_fasta_first_token(path)

    @staticmethod
    def _record_seq(record: object) -> str:
        return str(getattr(record, "seq", record))

    @staticmethod
    def _has_complete_pair_set(
        record_id_pairs: list[tuple[str, str]],
        record_ids: list[str],
    ) -> bool:
        expected_len = (len(record_ids) * (len(record_ids) - 1)) // 2
        if len(record_id_pairs) != expected_len:
            return False

        for observed_pair, expected_pair in zip(
            record_id_pairs, itertools.combinations(record_ids, 2)
        ):
            if observed_pair != expected_pair:
                return set(record_id_pairs) == set(
                    itertools.combinations(record_ids, 2)
                )

        return True

    @staticmethod
    def _sequence_array(seq: str):
        return np.frombuffer(seq.encode("ascii"), dtype=np.uint8)

    @staticmethod
    def _sequence_arrays(ref_seq: str, query_seq: str):
        try:
            return (
                SumOfPairsScore._sequence_array(ref_seq),
                SumOfPairsScore._sequence_array(query_seq),
            )
        except UnicodeEncodeError:
            return (
                np.array(list(ref_seq), dtype='U1'),
                np.array(list(query_seq), dtype='U1'),
            )

    @staticmethod
    def _stack_equal_length_sequence_pairs(
        ref_seqs: list[str],
        query_seqs: list[str],
        seq_len: int,
    ):
        try:
            ref_array = np.frombuffer(
                "".join(ref_seqs).encode("ascii"),
                dtype=np.uint8,
            ).reshape(len(ref_seqs), seq_len)
            query_array = np.frombuffer(
                "".join(query_seqs).encode("ascii"),
                dtype=np.uint8,
            ).reshape(len(query_seqs), seq_len)
            return ref_array, query_array
        except UnicodeEncodeError:
            return (
                np.array([list(seq) for seq in ref_seqs], dtype='U1'),
                np.array([list(seq) for seq in query_seqs], dtype='U1'),
            )

    @staticmethod
    def _calculate_equal_length_complete_pairs(
        record_id_pairs: list[tuple[str, str]],
        reference_records: dict[str, object],
        query_records: dict[str, object],
    ) -> tuple[int, int] | None:
        record_ids = list(reference_records.keys())
        if not SumOfPairsScore._has_complete_pair_set(record_id_pairs, record_ids):
            return None
        return SumOfPairsScore._calculate_equal_length_complete_records(
            reference_records,
            query_records,
        )

    @staticmethod
    def _calculate_equal_length_complete_records(
        reference_records: dict[str, object],
        query_records: dict[str, object],
    ) -> tuple[int, int] | None:
        record_ids = list(reference_records.keys())
        ref_seqs = [
            SumOfPairsScore._record_seq(reference_records[seq_id])
            for seq_id in record_ids
        ]
        query_seqs = [
            SumOfPairsScore._record_seq(query_records[seq_id])
            for seq_id in record_ids
        ]
        lengths = {len(seq) for seq in ref_seqs}
        lengths.update(len(seq) for seq in query_seqs)
        if len(lengths) != 1:
            return None

        seq_len = lengths.pop()
        if seq_len == 0:
            return 0, 0

        total_pairs = (len(record_ids) * (len(record_ids) - 1) // 2) * seq_len
        if ref_seqs == query_seqs:
            return total_pairs, total_pairs

        ref_array, query_array = SumOfPairsScore._stack_equal_length_sequence_pairs(
            ref_seqs,
            query_seqs,
            seq_len,
        )
        matching_taxa_per_site = np.count_nonzero(ref_array == query_array, axis=0)
        matches_per_site = (
            matching_taxa_per_site * (matching_taxa_per_site - 1)
        ) // 2
        return int(np.sum(matches_per_site)), total_pairs

    @staticmethod
    def _calculate_unchanged_record_pairs(
        record_id_pairs: list[tuple[str, str]],
        reference_records: dict[str, object],
        query_records: dict[str, object],
    ) -> tuple[int, int] | None:
        seq_cache = {}
        record_seq = SumOfPairsScore._record_seq
        total_pairs = 0

        for first_in_pair, second_in_pair in record_id_pairs:
            for seq_id in (first_in_pair, second_in_pair):
                if seq_id in seq_cache:
                    continue
                ref_seq = record_seq(reference_records[seq_id])
                query_seq = record_seq(query_records[seq_id])
                if ref_seq != query_seq:
                    return None
                seq_cache[seq_id] = ref_seq

            total_pairs += min(
                len(seq_cache[first_in_pair]),
                len(seq_cache[second_in_pair]),
            )

        return total_pairs, total_pairs

    @staticmethod
    def _process_pair_batch(
        pair_batch: list[tuple[str, str]],
        reference_records: dict[str, object],
        query_records: dict[str, object],
    ) -> tuple[int, int]:
        """Process a batch of sequence pairs in parallel."""
        batch_matches = 0
        batch_total = 0
        count_nonzero = np.count_nonzero

        # Pre-convert sequences to numpy arrays for the batch
        seq_arrays = {}
        for pair in pair_batch:
            for seq_id in pair:
                if seq_id not in seq_arrays:
                    ref_seq = SumOfPairsScore._record_seq(reference_records[seq_id])
                    query_seq = SumOfPairsScore._record_seq(query_records[seq_id])
                    ref_array, query_array = SumOfPairsScore._sequence_arrays(
                        ref_seq,
                        query_seq,
                    )
                    seq_arrays[seq_id] = (ref_seq, query_seq, ref_array, query_array)

        for first_in_pair, second_in_pair in pair_batch:
            ref_seq1_str, query_seq1_str, ref_seq1, query_seq1 = seq_arrays[
                first_in_pair
            ]
            ref_seq2_str, query_seq2_str, ref_seq2, query_seq2 = seq_arrays[
                second_in_pair
            ]

            # Check if all sequences have the same length
            if (len(ref_seq1) == len(query_seq1) and
                len(ref_seq2) == len(query_seq2) and
                len(ref_seq1) == len(ref_seq2)):
                # Use vectorized comparison
                matches = (ref_seq1 == query_seq1) & (ref_seq2 == query_seq2)
                batch_matches += count_nonzero(matches)
                batch_total += len(ref_seq1)
            else:
                min_len = min(len(ref_seq1_str), len(ref_seq2_str),
                             len(query_seq1_str), len(query_seq2_str))

                # Vectorize the mismatched comparison when possible
                if min_len > 0:
                    ref1_trimmed = ref_seq1[:min_len]
                    ref2_trimmed = ref_seq2[:min_len]
                    query1_trimmed = query_seq1[:min_len]
                    query2_trimmed = query_seq2[:min_len]

                    matches = (ref1_trimmed == query1_trimmed) & (
                        ref2_trimmed == query2_trimmed
                    )
                    batch_matches += count_nonzero(matches)
                    batch_total += min_len

        return int(batch_matches), batch_total

    def determine_number_of_matches_and_total_pairs(
        self,
        record_id_pairs: list[tuple[str, str]],
        reference_records: dict[str, object],
        query_records: dict[str, object],
    ) -> tuple[int, int]:
        fast_result = self._calculate_equal_length_complete_pairs(
            record_id_pairs,
            reference_records,
            query_records,
        )
        if fast_result is not None:
            return fast_result

        unchanged_result = self._calculate_unchanged_record_pairs(
            record_id_pairs,
            reference_records,
            query_records,
        )
        if unchanged_result is not None:
            return unchanged_result

        # For small datasets, use sequential processing
        if len(record_id_pairs) < 50:
            number_of_matches = 0
            number_of_total_pairs = 0
            seq_arrays = {}
            count_nonzero = np.count_nonzero
            record_seq = self._record_seq
            sequence_arrays = self._sequence_arrays

            for first_in_pair, second_in_pair in record_id_pairs:
                if first_in_pair not in seq_arrays:
                    ref_seq = record_seq(reference_records[first_in_pair])
                    query_seq = record_seq(query_records[first_in_pair])
                    ref_array, query_array = sequence_arrays(ref_seq, query_seq)
                    seq_arrays[first_in_pair] = (
                        ref_seq,
                        query_seq,
                        ref_array,
                        query_array,
                    )
                if second_in_pair not in seq_arrays:
                    ref_seq = record_seq(reference_records[second_in_pair])
                    query_seq = record_seq(query_records[second_in_pair])
                    ref_array, query_array = sequence_arrays(ref_seq, query_seq)
                    seq_arrays[second_in_pair] = (
                        ref_seq,
                        query_seq,
                        ref_array,
                        query_array,
                    )

                ref_seq1_str, query_seq1_str, ref_seq1, query_seq1 = seq_arrays[
                    first_in_pair
                ]
                ref_seq2_str, query_seq2_str, ref_seq2, query_seq2 = seq_arrays[
                    second_in_pair
                ]
                min_len = min(
                    len(ref_seq1_str),
                    len(ref_seq2_str),
                    len(query_seq1_str),
                    len(query_seq2_str),
                )

                if min_len > 0:
                    matches = (
                        (ref_seq1[:min_len] == query_seq1[:min_len])
                        & (ref_seq2[:min_len] == query_seq2[:min_len])
                    )
                    number_of_matches += int(count_nonzero(matches))
                    number_of_total_pairs += min_len

            return int(number_of_matches), number_of_total_pairs

        # Use multiprocessing for larger datasets
        num_workers = min(mp.cpu_count(), 8)
        batch_size = max(10, len(record_id_pairs) // (num_workers * 4))

        # Create batches
        pair_batches = [record_id_pairs[i:i + batch_size]
                       for i in range(0, len(record_id_pairs), batch_size)]

        # Process batches in parallel
        from functools import partial

        process_func = partial(self._process_pair_batch,
                              reference_records=reference_records,
                              query_records=query_records)

        with mp.Pool(processes=num_workers) as pool:
            batch_results = pool.map(process_func, pair_batches)

        # Aggregate results
        total_matches = sum(matches for matches, _ in batch_results)
        total_pairs = sum(pairs for _, pairs in batch_results)

        return int(total_matches), total_pairs
