import sys
import itertools
from typing import Dict, List, Tuple

import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

from Bio.Phylo.BaseTree import TreeMixin
from Bio.Phylo import Newick

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig


class TipToTipDistance(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            tip_1=parsed["tip_1"],
            tip_2=parsed["tip_2"],
        )
        self.json_output = parsed["json_output"]
        self.all_pairs = parsed["all_pairs"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def run(self):
        tree_zero = self.read_tree_file()

        if self.plot and not self.all_pairs:
            print("--plot requires --all-pairs for tip_to_tip_distance.")
            sys.exit(2)

        if self.all_pairs:
            rows = self.calculate_all_pairwise_distances(tree_zero)
            if self.plot:
                self._plot_tip_distance_heatmap(rows)

            if self.json_output:
                payload = dict(all_pairs=True, rows=rows, pairs=rows)
                if self.plot:
                    payload["plot_output"] = self.plot_output
                print_json(payload)
                return

            for row in rows:
                print(f"{row['taxon_a']}\t{row['taxon_b']}\t{row['tip_to_tip_distance']}")
            if self.plot:
                print(f"Saved tip-to-tip heatmap: {self.plot_output}")
            return

        if not self.tip_1 or not self.tip_2:
            print("tip_1 and tip_2 are required unless --all-pairs is provided.\nExiting...")
            sys.exit(2)

        self.check_leaves(tree_zero, self.tip_1, self.tip_2)
        distance = round(TreeMixin.distance(tree_zero, self.tip_1, self.tip_2), 4)
        if self.json_output:
            print_json(
                dict(
                    taxon_a=self.tip_1,
                    taxon_b=self.tip_2,
                    tip_to_tip_distance=distance,
                )
            )
            return
        print(distance)

    def check_leaves(
        self,
        tree_zero: Newick.Tree,
        tip_1: str,
        tip_2: str,
    ) -> None:
        leaf1 = TreeMixin.find_any(tree_zero, tip_1)
        if not bool(leaf1):
            print(tip_1, "not on tree\nExiting...")
            sys.exit(2)
        leaf2 = TreeMixin.find_any(tree_zero, tip_2)
        if not bool(leaf2):
            print(tip_2, "not on tree\nExiting...")
            sys.exit(2)

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree_zero,
            tip_1=getattr(args, "tip_1", None),
            tip_2=getattr(args, "tip_2", None),
            json_output=getattr(args, "json", False),
            all_pairs=getattr(args, "all_pairs", False),
            plot=getattr(args, "plot", False),
            plot_output=getattr(args, "plot_output", "tip_to_tip_distance_heatmap.png"),
            plot_config=PlotConfig.from_args(args),
        )

    def calculate_all_pairwise_distances(self, tree_zero: Newick.Tree) -> List[Dict[str, object]]:
        tips = [tip.name for tip in tree_zero.get_terminals()]
        rows = []
        for taxon_a, taxon_b in itertools.combinations(tips, 2):
            rows.append(
                dict(
                    taxon_a=taxon_a,
                    taxon_b=taxon_b,
                    tip_to_tip_distance=round(float(TreeMixin.distance(tree_zero, taxon_a, taxon_b)), 4),
                )
            )
        return rows

    def _build_distance_matrix(
        self,
        rows: List[Dict[str, object]],
    ) -> Tuple[List[str], np.ndarray]:
        taxa = sorted({str(row["taxon_a"]) for row in rows} | {str(row["taxon_b"]) for row in rows})
        n_taxa = len(taxa)
        taxon_to_index = {taxon: idx for idx, taxon in enumerate(taxa)}
        matrix = np.zeros((n_taxa, n_taxa), dtype=np.float64)

        for row in rows:
            i = taxon_to_index[str(row["taxon_a"])]
            j = taxon_to_index[str(row["taxon_b"])]
            distance = float(row["tip_to_tip_distance"])
            matrix[i, j] = distance
            matrix[j, i] = distance
        return taxa, matrix

    def _plot_tip_distance_heatmap(self, rows: List[Dict[str, object]]) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in tip_to_tip_distance. Install matplotlib and retry.")
            raise SystemExit(2)

        if not rows:
            return

        taxa, matrix = self._build_distance_matrix(rows)
        n_taxa = len(taxa)
        if n_taxa >= 3:
            condensed = squareform(matrix, checks=False)
            order = leaves_list(linkage(condensed, method="average"))
        else:
            order = np.arange(n_taxa)

        ordered_matrix = matrix[np.ix_(order, order)]
        ordered_taxa = [taxa[idx] for idx in order]

        config = self.plot_config
        config.resolve(n_rows=n_taxa, n_cols=n_taxa)

        fig_w = config.fig_width or max(6, min(20, n_taxa * 0.35))
        fig_h = config.fig_height or fig_w
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        image = ax.imshow(ordered_matrix, cmap="viridis", interpolation="nearest")

        if config.ylabel_fontsize and config.ylabel_fontsize > 0:
            ax.set_xticks(np.arange(n_taxa))
            ax.set_yticks(np.arange(n_taxa))
            ax.set_xticklabels(ordered_taxa, rotation=90, fontsize=config.xlabel_fontsize or config.ylabel_fontsize)
            ax.set_yticklabels(ordered_taxa, fontsize=config.ylabel_fontsize)
        else:
            ax.set_xticks([])
            ax.set_yticks([])

        ax.set_xlabel("Taxa (clustered)")
        ax.set_ylabel("Taxa (clustered)")
        colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        colorbar.set_label("Distance")

        if config.show_title:
            ax.set_title(config.title or "Tip-to-Tip Distance Heatmap", fontsize=config.title_fontsize)
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
