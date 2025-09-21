import itertools
from typing import Dict, List, Tuple
import numpy as np
import multiprocessing as mp
from functools import partial
import sys

from Bio.Align import MultipleSeqAlignment
try:
    from tqdm import tqdm
except ImportError:
    # Fallback if tqdm is not installed
    def tqdm(iterable, *args, **kwargs):
        return iterable

from .base import Alignment
from ...helpers.stats_summary import (
    calculate_summary_statistics_from_dict,
    print_summary_statistics,
)


class PairwiseIdentity(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _, is_protein = self.get_alignment_and_format()

        pair_ids, pairwise_identities, stats = \
            self.calculate_pairwise_identities(
                alignment, self.exclude_gaps
            )

        if self.verbose:
            try:
                for pair, identity in zip(
                    pair_ids, pairwise_identities.values()
                ):
                    print(f"{pair[0]}\t{pair[1]}\t{round(identity, 4)}")
            except BrokenPipeError:
                pass
        else:
            print_summary_statistics(stats)

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            verbose=args.verbose,
            exclude_gaps=args.exclude_gaps,
        )

    def _calculate_identity_vectorized(self, seq_arr1, seq_arr2, gap_mask=None, exclude_gaps=False):
        """Vectorized calculation of sequence identity."""
        matches = (seq_arr1 == seq_arr2)

        if exclude_gaps and gap_mask is not None:
            # Match original behavior: count identities when at least one doesn't have a gap
            # This matches the original "res_one not in gap_chars or res_two not in gap_chars"
            valid_for_identity = ~gap_mask[0] | ~gap_mask[1]
            identities = np.sum(matches & valid_for_identity)
        else:
            identities = np.sum(matches)

        # Total compared is always the full length (matching original behavior)
        total_compared = len(seq_arr1)

        return identities / total_compared if total_compared > 0 else 0

    def _process_pair_batch(self, alignment_data, pair_indices, exclude_gaps, gap_chars):
        """Process a batch of sequence pairs."""
        results = []
        for idx1, idx2 in pair_indices:
            seq_one = alignment_data[idx1]['seq']
            seq_two = alignment_data[idx2]['seq']

            if exclude_gaps:
                # Create boolean masks for gap positions
                gap_mask1 = np.isin(seq_one, list(gap_chars))
                gap_mask2 = np.isin(seq_two, list(gap_chars))
                identity = self._calculate_identity_vectorized(
                    seq_one, seq_two, (gap_mask1, gap_mask2), exclude_gaps
                )
            else:
                identity = self._calculate_identity_vectorized(seq_one, seq_two)

            results.append({
                'pair_id': [alignment_data[idx1]['id'], alignment_data[idx2]['id']],
                'identity': identity
            })
        return results

    def calculate_pairwise_identities(
        self,
        alignment: MultipleSeqAlignment,
        exclude_gaps: bool,
    ) -> Tuple[List[List[str]], Dict[str, float], Dict[str, float]]:
        gap_chars = self.get_gap_chars()

        # Convert sequences to numpy arrays for faster comparison
        alignment_data = []
        for record in alignment:
            seq_array = np.array([c.upper() for c in str(record.seq)], dtype='U1')
            alignment_data.append({
                'id': record.id,
                'seq': seq_array
            })

        # Generate all pairwise combinations
        all_pairs = list(itertools.combinations(range(len(alignment)), 2))

        pairwise_identities = {}
        pair_ids = []

        # For small datasets or when not using multiprocessing
        if len(all_pairs) < 50:
            # Process all pairs without multiprocessing
            results = self._process_pair_batch(alignment_data, all_pairs, exclude_gaps, gap_chars)
            for result in results:
                pair_id = result['pair_id']
                pair_ids.append(pair_id)
                pairwise_identities["-".join(pair_id)] = result['identity']
        else:
            # Use multiprocessing for larger datasets
            num_workers = min(mp.cpu_count(), 8)
            chunk_size = max(1, len(all_pairs) // (num_workers * 4))
            pair_chunks = [all_pairs[i:i + chunk_size] for i in range(0, len(all_pairs), chunk_size)]

            # Create partial function
            process_func = partial(
                self._process_pair_batch,
                alignment_data,
                exclude_gaps=exclude_gaps,
                gap_chars=gap_chars
            )

            # Process in parallel with progress bar
            with mp.Pool(processes=num_workers) as pool:
                # Only show progress bar if stderr is a tty (not redirected)
                if sys.stderr.isatty():
                    chunk_results = list(tqdm(
                        pool.imap(process_func, pair_chunks),
                        total=len(pair_chunks),
                        desc="Calculating pairwise identities",
                        unit="batch"
                    ))
                else:
                    chunk_results = pool.map(process_func, pair_chunks)

            # Combine results
            for chunk_result in chunk_results:
                for result in chunk_result:
                    pair_id = result['pair_id']
                    pair_ids.append(pair_id)
                    pairwise_identities["-".join(pair_id)] = result['identity']

        stats = calculate_summary_statistics_from_dict(pairwise_identities)

        return pair_ids, pairwise_identities, stats
