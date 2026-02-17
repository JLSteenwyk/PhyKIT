import itertools
from typing import Dict, List, Tuple
import numpy as np
import multiprocessing as mp
from functools import partial
import sys
import os
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

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
from ...helpers.json_output import print_json


class PairwiseIdentity(Alignment):
    MP_MIN_PAIRS = 2000
    MAX_MP_WORKERS = 8

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            alignment_file_path=parsed["alignment_file_path"],
            verbose=parsed["verbose"],
            exclude_gaps=parsed["exclude_gaps"],
        )
        self.json_output = parsed["json_output"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]

    def _should_use_multiprocessing(self, n_pairs: int) -> bool:
        if os.environ.get("PHYKIT_DISABLE_MP", "0") == "1":
            return False
        if os.environ.get("PHYKIT_FORCE_MP", "0") == "1":
            return True
        return n_pairs >= self.MP_MIN_PAIRS

    def run(self):
        alignment, _, is_protein = self.get_alignment_and_format()
        taxa = [record.id for record in alignment]

        pair_ids, pairwise_identities, stats = \
            self.calculate_pairwise_identities(
                alignment, self.exclude_gaps, is_protein
            )

        if self.plot:
            self._plot_pairwise_identity_heatmap(taxa, pair_ids, pairwise_identities)

        if self.json_output:
            self._print_json_output(pair_ids, pairwise_identities, stats)
            return

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

        if self.plot:
            print(f"Saved pairwise identity heatmap: {self.plot_output}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            verbose=args.verbose,
            exclude_gaps=args.exclude_gaps,
            json_output=getattr(args, "json", False),
            plot=getattr(args, "plot", False),
            plot_output=getattr(args, "plot_output", "pairwise_identity_heatmap.png"),
        )

    def _print_json_output(
        self,
        pair_ids: List[List[str]],
        pairwise_identities: Dict[str, float],
        stats: Dict[str, float],
    ) -> None:
        payload = dict(verbose=self.verbose, exclude_gaps=self.exclude_gaps)
        if self.verbose:
            rows = [
                dict(
                    taxon_a=pair[0],
                    taxon_b=pair[1],
                    identity=round(identity, 4),
                )
                for pair, identity in zip(pair_ids, pairwise_identities.values())
            ]
            payload["rows"] = rows
            payload["pairs"] = rows
        else:
            payload["summary"] = stats
        if self.plot:
            payload["plot_output"] = self.plot_output
        print_json(payload)

    def _plot_pairwise_identity_heatmap(
        self,
        taxa: List[str],
        pair_ids: List[List[str]],
        pairwise_identities: Dict[str, float],
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in pairwise_identity. Install matplotlib and retry.")
            raise SystemExit(2)

        n_taxa = len(taxa)
        if n_taxa == 0:
            return

        taxon_to_index = {taxon: idx for idx, taxon in enumerate(taxa)}
        matrix = np.ones((n_taxa, n_taxa), dtype=np.float32)

        for pair, identity in zip(pair_ids, pairwise_identities.values()):
            i = taxon_to_index[pair[0]]
            j = taxon_to_index[pair[1]]
            matrix[i, j] = identity
            matrix[j, i] = identity

        if n_taxa >= 3:
            distance_matrix = np.clip(1.0 - matrix, 0.0, 1.0)
            np.fill_diagonal(distance_matrix, 0.0)
            condensed = squareform(distance_matrix, checks=False)
            order = leaves_list(linkage(condensed, method="average"))
        else:
            order = np.arange(n_taxa)

        ordered_matrix = matrix[np.ix_(order, order)]
        ordered_taxa = [taxa[idx] for idx in order]

        fig_size = max(6, min(20, n_taxa * 0.35))
        fig, ax = plt.subplots(figsize=(fig_size, fig_size))
        image = ax.imshow(ordered_matrix, cmap="viridis", vmin=0, vmax=1, interpolation="nearest")
        ax.set_title("Pairwise Identity Heatmap")
        ax.set_xticks(np.arange(n_taxa))
        ax.set_yticks(np.arange(n_taxa))
        ax.set_xticklabels(ordered_taxa, rotation=90, fontsize=7)
        ax.set_yticklabels(ordered_taxa, fontsize=7)
        ax.set_xlabel("Taxa (clustered)")
        ax.set_ylabel("Taxa (clustered)")
        colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        colorbar.set_label("Identity")
        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=300, bbox_inches="tight")
        plt.close(fig)

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
        is_protein: bool = False,
    ) -> Tuple[List[List[str]], Dict[str, float], Dict[str, float]]:
        gap_chars = self.get_gap_chars(is_protein)

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

        # For small/medium workloads, multiprocessing overhead dominates.
        if not self._should_use_multiprocessing(len(all_pairs)):
            # Process all pairs without multiprocessing
            results = self._process_pair_batch(alignment_data, all_pairs, exclude_gaps, gap_chars)
            for result in results:
                pair_id = result['pair_id']
                pair_ids.append(pair_id)
                pairwise_identities["-".join(pair_id)] = result['identity']
        else:
            # Use multiprocessing for larger datasets
            num_workers = min(mp.cpu_count(), self.MAX_MP_WORKERS)
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
