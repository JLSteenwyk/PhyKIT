import itertools
from typing import Dict, List, Tuple
import numpy as np
import multiprocessing as mp
from functools import partial

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .base import Alignment


class SumOfPairsScore(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        query_records = SeqIO.to_dict(SeqIO.parse(self.fasta, "fasta"))
        reference_records = SeqIO.to_dict(SeqIO.parse(self.reference, "fasta"))

        record_id_pairs = list(
            itertools.combinations(reference_records.keys(), 2)
        )

        number_of_matches, number_of_total_pairs = \
            self.determine_number_of_matches_and_total_pairs(
                record_id_pairs, reference_records, query_records
            )

        print(round(number_of_matches / number_of_total_pairs, 4))

    def process_args(self, args) -> Dict[str, str]:
        return dict(fasta=args.fasta, reference=args.reference)

    @staticmethod
    def _process_pair_batch(
        pair_batch: List[Tuple[str, str]],
        reference_records: Dict[str, SeqRecord],
        query_records: Dict[str, SeqRecord],
    ) -> Tuple[int, int]:
        """Process a batch of sequence pairs in parallel."""
        batch_matches = 0
        batch_total = 0

        # Pre-convert sequences to numpy arrays for the batch
        seq_arrays = {}
        for pair in pair_batch:
            for seq_id in [pair[0], pair[1]]:
                if seq_id not in seq_arrays:
                    ref_seq = str(reference_records[seq_id].seq)
                    query_seq = str(query_records[seq_id].seq)
                    seq_arrays[seq_id] = {
                        'ref': ref_seq,
                        'query': query_seq,
                        'ref_array': np.array(list(ref_seq), dtype='U1'),
                        'query_array': np.array(list(query_seq), dtype='U1')
                    }

        for first_in_pair, second_in_pair in pair_batch:
            ref_seq1 = seq_arrays[first_in_pair]['ref_array']
            ref_seq2 = seq_arrays[second_in_pair]['ref_array']
            query_seq1 = seq_arrays[first_in_pair]['query_array']
            query_seq2 = seq_arrays[second_in_pair]['query_array']

            # Check if all sequences have the same length
            if (len(ref_seq1) == len(query_seq1) and
                len(ref_seq2) == len(query_seq2) and
                len(ref_seq1) == len(ref_seq2)):
                # Use vectorized comparison
                matches = (ref_seq1 == query_seq1) & (ref_seq2 == query_seq2)
                batch_matches += np.sum(matches)
                batch_total += len(ref_seq1)
            else:
                # Fall back to optimized comparison for mismatched lengths
                ref_seq1_str = seq_arrays[first_in_pair]['ref']
                ref_seq2_str = seq_arrays[second_in_pair]['ref']
                query_seq1_str = seq_arrays[first_in_pair]['query']
                query_seq2_str = seq_arrays[second_in_pair]['query']

                min_len = min(len(ref_seq1_str), len(ref_seq2_str),
                             len(query_seq1_str), len(query_seq2_str))

                # Vectorize the mismatched comparison when possible
                if min_len > 0:
                    ref1_trimmed = np.array(list(ref_seq1_str[:min_len]), dtype='U1')
                    ref2_trimmed = np.array(list(ref_seq2_str[:min_len]), dtype='U1')
                    query1_trimmed = np.array(list(query_seq1_str[:min_len]), dtype='U1')
                    query2_trimmed = np.array(list(query_seq2_str[:min_len]), dtype='U1')

                    matches = (ref1_trimmed == query1_trimmed) & (ref2_trimmed == query2_trimmed)
                    batch_matches += np.sum(matches)
                    batch_total += min_len

        return int(batch_matches), batch_total

    def determine_number_of_matches_and_total_pairs(
        self,
        record_id_pairs: List[Tuple[str, str]],
        reference_records: Dict[str, SeqRecord],
        query_records: Dict[str, SeqRecord],
    ) -> Tuple[int, int]:
        # For small datasets, use sequential processing
        if len(record_id_pairs) < 50:
            number_of_matches = 0
            number_of_total_pairs = 0

            for first_in_pair, second_in_pair in record_id_pairs:
                ref_seq1_str = str(reference_records[first_in_pair].seq)
                ref_seq2_str = str(reference_records[second_in_pair].seq)
                query_seq1_str = str(query_records[first_in_pair].seq)
                query_seq2_str = str(query_records[second_in_pair].seq)

                if (len(ref_seq1_str) == len(query_seq1_str) and
                    len(ref_seq2_str) == len(query_seq2_str) and
                    len(ref_seq1_str) == len(ref_seq2_str)):
                    # Use vectorized comparison
                    ref_seq1 = np.array(list(ref_seq1_str), dtype='U1')
                    ref_seq2 = np.array(list(ref_seq2_str), dtype='U1')
                    query_seq1 = np.array(list(query_seq1_str), dtype='U1')
                    query_seq2 = np.array(list(query_seq2_str), dtype='U1')

                    matches = (ref_seq1 == query_seq1) & (ref_seq2 == query_seq2)
                    number_of_matches += np.sum(matches)
                    number_of_total_pairs += len(ref_seq1)
                else:
                    min_len = min(len(ref_seq1_str), len(ref_seq2_str),
                                 len(query_seq1_str), len(query_seq2_str))

                    for i in range(min_len):
                        if (ref_seq1_str[i] == query_seq1_str[i] and
                            ref_seq2_str[i] == query_seq2_str[i]):
                            number_of_matches += 1
                        number_of_total_pairs += 1

            return int(number_of_matches), number_of_total_pairs

        # Use multiprocessing for larger datasets
        num_workers = min(mp.cpu_count(), 8)
        batch_size = max(10, len(record_id_pairs) // (num_workers * 4))

        # Create batches
        pair_batches = [record_id_pairs[i:i + batch_size]
                       for i in range(0, len(record_id_pairs), batch_size)]

        # Process batches in parallel
        process_func = partial(self._process_pair_batch,
                              reference_records=reference_records,
                              query_records=query_records)

        with mp.Pool(processes=num_workers) as pool:
            batch_results = pool.map(process_func, pair_batches)

        # Aggregate results
        total_matches = sum(matches for matches, _ in batch_results)
        total_pairs = sum(pairs for _, pairs in batch_results)

        return int(total_matches), total_pairs
