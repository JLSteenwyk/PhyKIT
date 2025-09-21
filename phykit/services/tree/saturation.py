from enum import Enum
import itertools
import sys
from typing import Dict, List, Tuple
import multiprocessing as mp
from functools import partial

from Bio import Align
from Bio.Phylo import Newick
import numpy as np
from sklearn.linear_model import LinearRegression

from .base import Tree
from ...helpers.files import (
    get_alignment_and_format as get_alignment_and_format_helper
)


class FileFormat(Enum):
    fasta = "fasta"
    clustal = "clustal"
    maf = "maf"
    mauve = "mauve"
    phylip = "phylip"
    phylip_seq = "phylip-sequential"
    phylip_rel = "phylip-relaxed"
    stockholm = "stockholm"


class Saturation(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        alignment, _, is_protein = get_alignment_and_format_helper(
            self.alignment_file_path
        )

        tree = self.read_tree_file()

        tips = self.get_tip_names_from_tree(tree)
        combos = list(itertools.combinations(tips, 2))

        (
            patristic_distances,
            uncorrected_distances,
        ) = self.loop_through_combos_and_calculate_pds_and_pis(
            combos, alignment, tree, self.exclude_gaps
        )

        # calculate slope and fit the y-intercept to zero
        # Fitting the y-intercept to zero follows Jeffroy et al.
        # See fig 2 https://www.cell.com/trends/genetics/fulltext/S0168-9525(06)00051-5
        model = LinearRegression(fit_intercept=False)
        model.fit(
            np.array(patristic_distances).reshape(-1, 1),
            np.array(uncorrected_distances)
        )

        self.print_res(
            self.verbose, combos, uncorrected_distances, patristic_distances, model.coef_[0]
        )

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            alignment_file_path=args.alignment,
            exclude_gaps=args.exclude_gaps,
            verbose=args.verbose,
        )

    def _process_combo_batch(self, tree, seq_arrays, gap_mask, exclude_gaps, combo_batch):
        """Process a batch of combinations in parallel."""
        results = []
        for combo in combo_batch:
            # Calculate patristic distance
            pd = tree.distance(combo[0], combo[1])

            # Calculate uncorrected distance using numpy operations
            seq1_arr = seq_arrays[combo[0]]
            seq2_arr = seq_arrays[combo[1]]

            if exclude_gaps:
                # Use pre-computed gap masks
                gap_mask1 = gap_mask[combo[0]]
                gap_mask2 = gap_mask[combo[1]]
                valid_positions = ~(gap_mask1 | gap_mask2)

                if np.any(valid_positions):
                    matches = seq1_arr[valid_positions] == seq2_arr[valid_positions]
                    identities = np.sum(matches)
                    adjusted_len = np.sum(valid_positions)
                    ud = 1 - (identities / adjusted_len)
                else:
                    ud = float('nan')
            else:
                matches = seq1_arr == seq2_arr
                identities = np.sum(matches)
                ud = 1 - (identities / len(seq1_arr))

            results.append((pd, ud))
        return results

    def loop_through_combos_and_calculate_pds_and_pis(
        self,
        combos: List[Tuple[str, str]],
        alignment: Align.MultipleSeqAlignment,
        tree: Newick.Tree,
        exclude_gaps: bool,
    ) -> Tuple[
        List[float],
        List[float]
    ]:
        """
        loop through all taxon combinations and determine
        their patristic distance and pairwise identity
        """
        gap_chars = self.get_gap_chars()

        # Convert sequences to numpy arrays for vectorized operations
        seq_arrays = {}
        gap_mask = {}
        for record in alignment:
            seq_arr = np.array([c.upper() for c in str(record.seq)], dtype='U1')
            seq_arrays[record.name] = seq_arr
            if exclude_gaps:
                gap_mask[record.name] = np.isin(seq_arr, list(gap_chars))

        # For small datasets, process sequentially
        if len(combos) < 50:
            patristic_distances = []
            uncorrected_distances = []
            for combo in combos:
                pd = tree.distance(combo[0], combo[1])
                patristic_distances.append(pd)

                seq1_arr = seq_arrays[combo[0]]
                seq2_arr = seq_arrays[combo[1]]

                if exclude_gaps:
                    gap_mask1 = gap_mask[combo[0]]
                    gap_mask2 = gap_mask[combo[1]]
                    valid_positions = ~(gap_mask1 | gap_mask2)

                    if np.any(valid_positions):
                        matches = seq1_arr[valid_positions] == seq2_arr[valid_positions]
                        identities = np.sum(matches)
                        adjusted_len = np.sum(valid_positions)
                        ud = 1 - (identities / adjusted_len)
                    else:
                        ud = float('nan')
                else:
                    matches = seq1_arr == seq2_arr
                    identities = np.sum(matches)
                    ud = 1 - (identities / len(seq1_arr))

                uncorrected_distances.append(ud)
        else:
            # Use multiprocessing for larger datasets
            num_workers = min(mp.cpu_count(), 8)
            chunk_size = max(1, len(combos) // (num_workers * 4))
            combo_chunks = [combos[i:i + chunk_size] for i in range(0, len(combos), chunk_size)]

            # Create partial function
            process_func = partial(
                self._process_combo_batch,
                tree,
                seq_arrays,
                gap_mask,
                exclude_gaps
            )

            # Process in parallel
            with mp.Pool(processes=num_workers) as pool:
                chunk_results = pool.map(process_func, combo_chunks)

            # Flatten results
            patristic_distances = []
            uncorrected_distances = []
            for chunk_result in chunk_results:
                for pd, ud in chunk_result:
                    patristic_distances.append(pd)
                    uncorrected_distances.append(ud)

        return patristic_distances, uncorrected_distances

    def print_res(
        self,
        verbose: bool,
        combos: List[str],
        uncorrected_distances: List[float],
        patristic_distances: List[float],
        slope: float,
    ) -> None:
        """
        print results to stdout
        """
        try:
            if verbose:
                for cbo, dist, pd in zip(
                    combos, uncorrected_distances, patristic_distances
                ):
                    print(
                        f"{cbo[0]}\t{cbo[1]}\t{round(dist,4)}\t{round(pd, 4)}"
                    )
            else:
                print(f"{round(slope, 4)}\t{abs(round(1-slope, 4))}")
        except BrokenPipeError:
            pass
