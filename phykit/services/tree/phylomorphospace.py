import os
import pickle
from typing import Dict, List, Tuple

import numpy as np

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig
from ...helpers.trait_parsing import parse_multi_trait_file
from ...errors import PhykitUserError


class Phylomorphospace(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.trait_x = parsed["trait_x"]
        self.trait_y = parsed["trait_y"]
        self.color_by = parsed["color_by"]
        self.plot_output = parsed["plot_output"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            trait_x=getattr(args, "trait_x", None),
            trait_y=getattr(args, "trait_y", None),
            color_by=getattr(args, "color_by", None),
            plot_output=getattr(args, "plot_output", "phylomorphospace_plot.png"),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        tree = self.read_tree_file()
        self.validate_tree(tree, min_tips=3, require_branch_lengths=True, context="phylomorphospace")

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(
            self.trait_data_path, tree_tips
        )

        # Resolve trait_x and trait_y
        trait_x_name, trait_y_name, x_idx, y_idx = self._resolve_trait_axes(
            trait_names
        )

        ordered_names = sorted(traits.keys())

        # Build 2-column data matrix for the selected traits
        data = np.array(
            [[traits[name][x_idx], traits[name][y_idx]] for name in ordered_names]
        )

        # Full data matrix for color-by column lookup
        Y = np.array(
            [[traits[name][j] for j in range(len(trait_names))] for name in ordered_names]
        )

        # Ancestral reconstruction on the 2-column matrix
        node_estimates, node_distances, tree_pruned = (
            self._reconstruct_ancestral_scores(tree, data, ordered_names)
        )

        self._plot_phylomorphospace(
            tree_pruned, data, ordered_names,
            node_estimates, node_distances,
            trait_x_name, trait_y_name,
            trait_names, Y,
        )

        if self.json_output:
            result = {
                "trait_x": trait_x_name,
                "trait_y": trait_y_name,
                "tip_data": {
                    name: {
                        trait_x_name: float(data[i, 0]),
                        trait_y_name: float(data[i, 1]),
                    }
                    for i, name in enumerate(ordered_names)
                },
                "plot_output": self.plot_output,
            }
            print_json(result)
        else:
            print(f"Saved phylomorphospace plot: {self.plot_output}")

    def _resolve_trait_axes(
        self, trait_names: List[str]
    ) -> Tuple[str, str, int, int]:
        if self.trait_x is None and self.trait_y is None:
            if len(trait_names) == 2:
                return trait_names[0], trait_names[1], 0, 1
            else:
                raise PhykitUserError(
                    [
                        f"Trait file has {len(trait_names)} traits: "
                        f"{', '.join(trait_names)}.",
                        "Use --trait-x and --trait-y to specify which two to plot.",
                    ],
                    code=2,
                )

        if self.trait_x is None or self.trait_y is None:
            raise PhykitUserError(
                ["Both --trait-x and --trait-y must be specified together."],
                code=2,
            )

        if self.trait_x not in trait_names:
            raise PhykitUserError(
                [
                    f"--trait-x '{self.trait_x}' not found in trait file columns: "
                    f"{', '.join(trait_names)}.",
                ],
                code=2,
            )
        if self.trait_y not in trait_names:
            raise PhykitUserError(
                [
                    f"--trait-y '{self.trait_y}' not found in trait file columns: "
                    f"{', '.join(trait_names)}.",
                ],
                code=2,
            )

        x_idx = trait_names.index(self.trait_x)
        y_idx = trait_names.index(self.trait_y)
        return self.trait_x, self.trait_y, x_idx, y_idx

    def _reconstruct_ancestral_scores(
        self, tree, scores: np.ndarray, ordered_names: List[str]
    ) -> Tuple[Dict, Dict, object]:
        tree_copy = pickle.loads(pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL))

        tip_names_in_tree = [t.name for t in tree_copy.get_terminals()]
        tips_to_prune = [t for t in tip_names_in_tree if t not in ordered_names]
        if tips_to_prune:
            tree_copy = self.prune_tree_using_taxa_list(tree_copy, tips_to_prune)

        name_to_idx = {name: i for i, name in enumerate(ordered_names)}

        node_estimates = {}
        node_variances = {}
        node_distances = {}

        root = tree_copy.root
        for clade in tree_copy.find_clades(order="postorder"):
            if clade == root:
                node_distances[id(clade)] = 0.0
            else:
                node_distances[id(clade)] = tree_copy.distance(root, clade)

        for clade in tree_copy.find_clades(order="postorder"):
            if clade.is_terminal():
                name = clade.name
                if name in name_to_idx:
                    node_estimates[id(clade)] = scores[name_to_idx[name]]
                    node_variances[id(clade)] = 0.0
            else:
                children = clade.clades
                precisions = []
                weighted_scores = []
                for child in children:
                    child_id = id(child)
                    if child_id not in node_estimates:
                        continue
                    v_i = child.branch_length if child.branch_length else 0.0
                    child_var = node_variances[child_id]
                    denom = v_i + child_var
                    if denom == 0:
                        denom = 1e-10
                    prec = 1.0 / denom
                    precisions.append(prec)
                    weighted_scores.append(prec * node_estimates[child_id])

                if precisions:
                    total_prec = sum(precisions)
                    est = sum(weighted_scores) / total_prec
                    node_estimates[id(clade)] = est
                    node_variances[id(clade)] = 1.0 / total_prec

        return node_estimates, node_distances, tree_copy

    def _parse_color_by(
        self,
        color_by: str,
        trait_names: List[str],
        Y: np.ndarray,
        ordered_names: List[str],
    ) -> Tuple[np.ndarray, List[str], str]:
        if color_by in trait_names:
            col_idx = trait_names.index(color_by)
            return Y[:, col_idx], [], "continuous"

        if not os.path.isfile(color_by):
            raise PhykitUserError(
                [
                    f"--color-by '{color_by}' is not a trait column name "
                    f"({', '.join(trait_names)}) and is not a valid file path.",
                ],
                code=2,
            )

        name_to_idx = {name: i for i, name in enumerate(ordered_names)}
        values = [None] * len(ordered_names)

        with open(color_by) as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                parts = stripped.split("\t")
                if len(parts) < 2:
                    continue
                taxon = parts[0]
                if taxon in name_to_idx:
                    values[name_to_idx[taxon]] = parts[1]

        missing = [ordered_names[i] for i, v in enumerate(values) if v is None]
        if missing:
            raise PhykitUserError(
                [
                    f"--color-by file missing values for taxa: {', '.join(missing)}",
                ],
                code=2,
            )

        is_numeric = True
        for v in values:
            try:
                float(v)
            except (ValueError, TypeError):
                is_numeric = False
                break

        if is_numeric:
            return np.array([float(v) for v in values]), [], "continuous"
        else:
            categories = sorted(set(values))
            return np.array(values), categories, "discrete"

    def _plot_phylomorphospace(
        self,
        tree,
        data: np.ndarray,
        ordered_names: List[str],
        node_estimates: Dict,
        node_distances: Dict,
        trait_x_name: str,
        trait_y_name: str,
        trait_names: List[str],
        Y: np.ndarray,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.collections import LineCollection
            from matplotlib.colors import Normalize
        except ImportError:
            print("matplotlib is required for phylomorphospace. Install matplotlib and retry.")
            raise SystemExit(2)

        config = self.plot_config
        config.resolve(n_rows=None, n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        # Collect all distances for normalization
        all_dists = [d for d in node_distances.values()]
        max_dist = max(all_dists) if all_dists else 1.0

        # Draw tree edges colored by distance from root
        segments = []
        colors = []
        norm = Normalize(vmin=0, vmax=max_dist)
        cmap = plt.get_cmap("coolwarm")

        for clade in tree.find_clades(order="preorder"):
            parent_id = id(clade)
            if parent_id not in node_estimates:
                continue
            parent_vals = node_estimates[parent_id]

            for child in clade.clades:
                child_id = id(child)
                if child_id not in node_estimates:
                    continue
                child_vals = node_estimates[child_id]

                x0, y0 = parent_vals[0], parent_vals[1]
                x1, y1 = child_vals[0], child_vals[1]
                segments.append([(x0, y0), (x1, y1)])

                parent_dist = node_distances.get(parent_id, 0)
                child_dist = node_distances.get(child_id, 0)
                mid_dist = (parent_dist + child_dist) / 2.0
                colors.append(mid_dist)

        if segments:
            lc = LineCollection(
                segments,
                array=np.array(colors),
                cmap=cmap,
                norm=norm,
                linewidths=1.0,
                alpha=0.7,
                zorder=2,
            )
            ax.add_collection(lc)
            cbar = fig.colorbar(lc, ax=ax, pad=0.02, fraction=0.046)
            cbar.set_label("Distance from root")

        # Determine tip point colors
        color_values = None
        color_categories = None
        color_kind = None
        if self.color_by:
            color_values, color_categories, color_kind = self._parse_color_by(
                self.color_by, trait_names, Y, ordered_names
            )

        if color_values is not None and color_kind == "continuous":
            scatter = ax.scatter(
                data[:, 0],
                data[:, 1],
                s=40,
                alpha=0.8,
                c=color_values,
                cmap="viridis",
                edgecolors="white",
                linewidth=0.5,
                zorder=3,
            )
            cbar = fig.colorbar(scatter, ax=ax, pad=0.02, fraction=0.046)
            cbar.set_label(
                self.color_by if self.color_by in trait_names else "color value"
            )
        elif color_values is not None and color_kind == "discrete":
            unique_cats = color_categories
            cat_colors = plt.get_cmap("tab10")
            cat_map = {
                cat: cat_colors(i / max(len(unique_cats) - 1, 1))
                for i, cat in enumerate(unique_cats)
            }
            point_colors = [cat_map[v] for v in color_values]
            ax.scatter(
                data[:, 0],
                data[:, 1],
                s=40,
                alpha=0.8,
                c=point_colors,
                edgecolors="white",
                linewidth=0.5,
                zorder=3,
            )
            from matplotlib.lines import Line2D
            handles = [
                Line2D(
                    [0], [0], marker="o", color="w",
                    markerfacecolor=cat_map[cat], markersize=7, label=cat,
                )
                for cat in unique_cats
            ]
            ax.legend(
                handles=handles, title=self.color_by, fontsize=7, title_fontsize=8
            )
        else:
            ax.scatter(
                data[:, 0],
                data[:, 1],
                s=40,
                alpha=0.8,
                color="#2b8cbe",
                edgecolors="white",
                linewidth=0.5,
                zorder=3,
            )

        for k, name in enumerate(ordered_names):
            ax.annotate(
                name,
                (data[k, 0], data[k, 1]),
                textcoords="offset points",
                xytext=(5, 5),
                fontsize=8,
            )

        ax.axhline(0, color="#cccccc", linewidth=0.8, linestyle="--", zorder=1)
        ax.axvline(0, color="#cccccc", linewidth=0.8, linestyle="--", zorder=1)

        ax.set_xlabel(trait_x_name)
        ax.set_ylabel(trait_y_name)
        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
