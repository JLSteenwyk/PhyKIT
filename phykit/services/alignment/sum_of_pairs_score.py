from __future__ import annotations

import itertools

from ._fasta import read_unique_fasta_first_token
from .base import Alignment


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyMultiprocessing:
    _module = None

    def _load(self):
        module = self._module
        if module is None:
            import multiprocessing as module

            self._module = module
        return module

    def cpu_count(self):
        return self._load().cpu_count()

    def Pool(self, *args, **kwargs):
        return self._load().Pool(*args, **kwargs)


class _LazyNumpy:
    _module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = self._module = _np

        value = getattr(module, name)
        setattr(self, name, value)
        return value


mp = _LazyMultiprocessing()
np = _LazyNumpy()
_SOPS_SCALAR_PAIR_MAX_CELLS = 8192


def _bounded_match_count_dtype(max_count: int):
    if max_count <= 0xFFFF:
        return np.uint16
    if max_count <= 0xFFFFFFFF:
        return np.uint32
    return np.uint64


def _bounded_matches_per_site(match_mask):
    count_dtype = _bounded_match_count_dtype(match_mask.shape[0])
    return match_mask.sum(axis=0, dtype=count_dtype).astype(np.intp)


class SumOfPairsScore(Alignment):
    MP_MIN_PAIRS = 500_000

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(fasta=parsed["fasta"], reference=parsed["reference"])
        self.json_output = parsed["json_output"]

    def run(self):
        query_records = self._read_fasta(self.fasta)
        reference_records = (
            query_records
            if self.reference == self.fasta
            else self._read_fasta(self.reference)
        )

        fast_result = self._calculate_equal_length_complete_records(
            reference_records, query_records
        )
        if fast_result is None:
            fast_result = self._calculate_unchanged_complete_records(
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
        equal_length_result = SumOfPairsScore._calculate_equal_length_complete_records(
            reference_records,
            query_records,
        )
        if equal_length_result is not None:
            return equal_length_result
        return SumOfPairsScore._calculate_unchanged_complete_records(
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
        if not ref_seqs:
            return None

        seq_len = len(ref_seqs[0])
        number_of_records = len(record_ids)
        if reference_records is query_records:
            for ref_seq in ref_seqs:
                if len(ref_seq) != seq_len:
                    return None
            if seq_len == 0:
                return 0, 0
            total_pairs = (number_of_records * (number_of_records - 1) // 2) * seq_len
            return total_pairs, total_pairs

        query_seqs = [
            SumOfPairsScore._record_seq(query_records[seq_id])
            for seq_id in record_ids
        ]
        if not query_seqs:
            return None

        changed_count = 0
        changed_ref_seqs = []
        changed_query_seqs = []
        for ref_seq, query_seq in zip(ref_seqs, query_seqs):
            if len(ref_seq) != seq_len or len(query_seq) != seq_len:
                return None
            if ref_seq != query_seq:
                changed_count += 1
                changed_ref_seqs.append(ref_seq)
                changed_query_seqs.append(query_seq)
        if seq_len == 0:
            return 0, 0

        total_pairs = (len(record_ids) * (len(record_ids) - 1) // 2) * seq_len
        if changed_count == 0:
            return total_pairs, total_pairs

        ref_array, query_array = SumOfPairsScore._stack_equal_length_sequence_pairs(
            changed_ref_seqs,
            changed_query_seqs,
            seq_len,
        )
        matching_taxa_per_site = (
            _bounded_matches_per_site(ref_array == query_array)
            + number_of_records
            - changed_count
        )
        matches_per_site = (
            matching_taxa_per_site * (matching_taxa_per_site - 1)
        ) // 2
        return int(matches_per_site.sum()), total_pairs

    @staticmethod
    def _calculate_unchanged_complete_records(
        reference_records: dict[str, object],
        query_records: dict[str, object],
    ) -> tuple[int, int] | None:
        if not reference_records:
            return None

        record_seq = SumOfPairsScore._record_seq
        lengths = []
        if reference_records is query_records:
            for record in reference_records.values():
                lengths.append(len(record_seq(record)))
        else:
            for seq_id, reference_record in reference_records.items():
                ref_seq = record_seq(reference_record)
                if ref_seq != record_seq(query_records[seq_id]):
                    return None
                lengths.append(len(ref_seq))

        lengths.sort()
        total_pairs = 0
        remaining = len(lengths) - 1
        for length in lengths:
            total_pairs += length * remaining
            remaining -= 1
        return total_pairs, total_pairs

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
    def _calculate_pair_matches_scalar(
        record_id_pairs: list[tuple[str, str]],
        reference_records: dict[str, object],
        query_records: dict[str, object],
    ) -> tuple[int, int] | None:
        seq_cache = {}
        record_seq = SumOfPairsScore._record_seq
        total_cells = 0

        for first_in_pair, second_in_pair in record_id_pairs:
            for seq_id in (first_in_pair, second_in_pair):
                if seq_id in seq_cache:
                    continue
                ref_seq = record_seq(reference_records[seq_id])
                query_seq = record_seq(query_records[seq_id])
                comparable_len = min(len(ref_seq), len(query_seq))
                total_cells += comparable_len
                if total_cells > _SOPS_SCALAR_PAIR_MAX_CELLS:
                    return None
                seq_cache[seq_id] = (comparable_len, ref_seq, query_seq)

        number_of_matches = 0
        number_of_total_pairs = 0
        for first_in_pair, second_in_pair in record_id_pairs:
            len1, ref_seq1, query_seq1 = seq_cache[first_in_pair]
            len2, ref_seq2, query_seq2 = seq_cache[second_in_pair]
            min_len = min(len1, len2)
            for idx in range(min_len):
                if (
                    ref_seq1[idx] == query_seq1[idx]
                    and ref_seq2[idx] == query_seq2[idx]
                ):
                    number_of_matches += 1
            number_of_total_pairs += min_len

        return number_of_matches, number_of_total_pairs

    @staticmethod
    def _sequence_match_masks_for_pairs(
        record_id_pairs: list[tuple[str, str]],
        reference_records: dict[str, object],
        query_records: dict[str, object],
    ) -> dict[str, tuple[int, object]]:
        seq_masks = {}
        record_seq = SumOfPairsScore._record_seq
        sequence_arrays = SumOfPairsScore._sequence_arrays

        for pair in record_id_pairs:
            for seq_id in pair:
                if seq_id in seq_masks:
                    continue
                ref_seq = record_seq(reference_records[seq_id])
                query_seq = record_seq(query_records[seq_id])
                ref_array, query_array = sequence_arrays(ref_seq, query_seq)
                comparable_len = min(len(ref_array), len(query_array))
                if comparable_len > 0:
                    match_mask = (
                        ref_array[:comparable_len]
                        == query_array[:comparable_len]
                    )
                else:
                    match_mask = None
                seq_masks[seq_id] = (comparable_len, match_mask)

        return seq_masks

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
        logical_and = np.logical_and
        seq_masks = SumOfPairsScore._sequence_match_masks_for_pairs(
            pair_batch,
            reference_records,
            query_records,
        )
        scratch_len = max((length for length, _ in seq_masks.values()), default=0)
        scratch = np.empty(scratch_len, dtype=np.bool_)

        for first_in_pair, second_in_pair in pair_batch:
            len1, matches1 = seq_masks[first_in_pair]
            len2, matches2 = seq_masks[second_in_pair]
            min_len = min(len1, len2)
            if min_len > 0:
                out = scratch[:min_len]
                logical_and(matches1[:min_len], matches2[:min_len], out=out)
                batch_matches += int(count_nonzero(out))
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

        scalar_result = self._calculate_pair_matches_scalar(
            record_id_pairs,
            reference_records,
            query_records,
        )
        if scalar_result is not None:
            return scalar_result

        # For small and medium fallback datasets, multiprocessing overhead
        # dominates the cached sequential pair kernel.
        if len(record_id_pairs) < self.MP_MIN_PAIRS:
            number_of_matches = 0
            number_of_total_pairs = 0
            count_nonzero = np.count_nonzero
            logical_and = np.logical_and
            seq_masks = self._sequence_match_masks_for_pairs(
                record_id_pairs,
                reference_records,
                query_records,
            )
            scratch_len = max((length for length, _ in seq_masks.values()), default=0)
            scratch = np.empty(scratch_len, dtype=np.bool_)

            for first_in_pair, second_in_pair in record_id_pairs:
                len1, matches1 = seq_masks[first_in_pair]
                len2, matches2 = seq_masks[second_in_pair]
                min_len = min(len1, len2)
                if min_len > 0:
                    out = scratch[:min_len]
                    logical_and(matches1[:min_len], matches2[:min_len], out=out)
                    number_of_matches += int(count_nonzero(out))
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
