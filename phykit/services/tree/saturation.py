import itertools
from typing import Dict, List, Tuple
import multiprocessing as mp
from functools import partial
import os

from Bio import Align
from Bio.Phylo import Newick
import numpy as np

from .base import Tree
from ...helpers.files import (
    get_alignment_and_format as get_alignment_and_format_helper
)
from ...helpers.json_output import print_json


class Saturation(Tree):
    MP_MIN_COMBOS = 2000
    MAX_MP_WORKERS = 8

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            alignment_file_path=parsed["alignment_file_path"],
            exclude_gaps=parsed["exclude_gaps"],
            verbose=parsed["verbose"],
        )
        self.json_output = parsed["json_output"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]

    def _should_use_multiprocessing(self, n_combos: int) -> bool:
        if os.environ.get("PHYKIT_DISABLE_MP", "0") == "1":
            return False
        if os.environ.get("PHYKIT_FORCE_MP", "0") == "1":
            return True
        return n_combos >= self.MP_MIN_COMBOS

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
            combos, alignment, tree, self.exclude_gaps, is_protein
        )

        # calculate slope while fitting the y-intercept to zero.
        # This follows Jeffroy et al. (fig 2):
        # https://www.cell.com/trends/genetics/fulltext/S0168-9525(06)00051-5
        x = np.asarray(patristic_distances, dtype=float)
        y = np.asarray(uncorrected_distances, dtype=float)
        finite_mask = np.isfinite(x) & np.isfinite(y)
        x = x[finite_mask]
        y = y[finite_mask]

        denom = float(np.dot(x, x))
        slope = float(np.dot(x, y) / denom) if denom != 0.0 else 0.0

        if self.plot:
            self._plot_saturation_scatter(x, y, slope)

        self.print_res(
            self.verbose, combos, uncorrected_distances, patristic_distances, slope
        )

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            alignment_file_path=args.alignment,
            exclude_gaps=args.exclude_gaps,
            verbose=args.verbose,
            json_output=getattr(args, "json", False),
            plot=getattr(args, "plot", False),
            plot_output=getattr(args, "plot_output", "saturation_plot.png"),
        )

    def _plot_saturation_scatter(
        self,
        patristic_distances: np.ndarray,
        uncorrected_distances: np.ndarray,
        slope: float,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in saturation. Install matplotlib and retry.")
            raise SystemExit(2)

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.scatter(
            patristic_distances,
            uncorrected_distances,
            s=14,
            alpha=0.6,
            color="#2b8cbe",
            edgecolors="none",
        )

        if patristic_distances.size > 0:
            x_line = np.linspace(0.0, float(np.max(patristic_distances)), 200)
            y_line = slope * x_line
            ax.plot(
                x_line,
                y_line,
                color="#000000",
                linestyle="--",
                linewidth=2.0,
                label=f"Fit through origin (slope={slope:.4f})",
            )
            ax.legend(loc="best", frameon=False)

        ax.set_title("Saturation: Patristic vs Uncorrected Distance")
        ax.set_xlabel("Patristic distance")
        ax.set_ylabel("Uncorrected distance")
        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=300, bbox_inches="tight")
        plt.close(fig)

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
        is_protein: bool = False,
    ) -> Tuple[
        List[float],
        List[float]
    ]:
        """
        loop through all taxon combinations and determine
        their patristic distance and pairwise identity
        """
        gap_chars = self.get_gap_chars(is_protein)

        # Convert sequences to numpy arrays for vectorized operations
        seq_arrays = {}
        gap_mask = {}
        for record in alignment:
            seq_arr = np.array([c.upper() for c in str(record.seq)], dtype='U1')
            seq_arrays[record.name] = seq_arr
            if exclude_gaps:
                gap_mask[record.name] = np.isin(seq_arr, list(gap_chars))

        # For small/medium workloads, multiprocessing overhead dominates.
        if not self._should_use_multiprocessing(len(combos)):
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
            num_workers = min(mp.cpu_count(), self.MAX_MP_WORKERS)
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
            if self.json_output:
                payload = dict(verbose=verbose, exclude_gaps=self.exclude_gaps)
                if verbose:
                    rows = [
                        dict(
                            taxon_a=cbo[0],
                            taxon_b=cbo[1],
                            uncorrected_distance=round(dist, 4),
                            patristic_distance=round(pd, 4),
                        )
                        for cbo, dist, pd in zip(
                            combos, uncorrected_distances, patristic_distances
                        )
                    ]
                    payload["rows"] = rows
                    payload["pairs"] = rows
                else:
                    payload["summary"] = dict(
                        slope=round(slope, 4),
                        one_minus_slope_abs=abs(round(1 - slope, 4)),
                    )
                if self.plot:
                    payload["plot_output"] = self.plot_output
                print_json(payload)
                return

            if verbose:
                for cbo, dist, pd in zip(
                    combos, uncorrected_distances, patristic_distances
                ):
                    print(
                        f"{cbo[0]}\t{cbo[1]}\t{round(dist,4)}\t{round(pd, 4)}"
                    )
            else:
                print(f"{round(slope, 4)}\t{abs(round(1-slope, 4))}")
            if self.plot:
                print(f"Saved saturation plot: {self.plot_output}")
        except BrokenPipeError:
            pass
