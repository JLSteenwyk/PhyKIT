from __future__ import annotations

import os

from .base import Tree
from ...errors import PhykitUserError


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def parse_multi_trait_file(*args, **kwargs):
    from ...helpers.trait_parsing import (
        parse_multi_trait_file as _parse_multi_trait_file,
    )

    return _parse_multi_trait_file(*args, **kwargs)


def trait_matrix_from_rows(*args, **kwargs):
    from ...helpers.trait_parsing import (
        trait_matrix_from_rows as _trait_matrix_from_rows,
    )

    return _trait_matrix_from_rows(*args, **kwargs)


def _root_distance_max(values):
    return max(values, default=1.0)


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

    def process_args(self, args) -> dict[str, str]:
        from ...helpers.plot_config import PlotConfig

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
        tree = self.read_tree_file_unmodified()
        self.validate_tree(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="phylomorphospace",
        )

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(
            self.trait_data_path, tree_tips
        )

        # Resolve trait_x and trait_y
        trait_x_name, trait_y_name, x_idx, y_idx = self._resolve_trait_axes(
            trait_names
        )

        ordered_names = sorted(traits.keys())

        # Full data matrix for color-by column lookup, plus selected plot axes.
        Y = trait_matrix_from_rows(traits, ordered_names)
        data = Y[:, [x_idx, y_idx]]

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
        self, trait_names: list[str]
    ) -> tuple[str, str, int, int]:
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
        self, tree, scores: np.ndarray, ordered_names: list[str]
    ) -> tuple[dict, dict, object]:
        ordered_name_set = set(ordered_names)
        tip_names_in_tree = self.get_tip_names_from_tree(tree)
        tips_to_prune = [t for t in tip_names_in_tree if t not in ordered_name_set]
        tree_for_analysis = self._fast_copy(tree) if tips_to_prune else tree
        if tips_to_prune:
            tree_for_analysis = self.prune_tree_using_taxa_list(
                tree_for_analysis, tips_to_prune
            )

        name_to_idx = {name: i for i, name in enumerate(ordered_names)}

        root = tree_for_analysis.root
        try:
            preorder_clades = []
            node_distances = {}
            stack = [(root, 0.0)]
            while stack:
                clade, distance = stack.pop()
                preorder_clades.append(clade)
                node_distances[id(clade)] = distance
                children = clade.clades
                if children:
                    for child in reversed(children):
                        branch_length = (
                            child.branch_length if child.branch_length else 0.0
                        )
                        stack.append((child, distance + branch_length))
        except AttributeError:
            node_distances = {}
            for clade in tree_for_analysis.find_clades(order="postorder"):
                if clade == root:
                    node_distances[id(clade)] = 0.0
                else:
                    node_distances[id(clade)] = tree_for_analysis.distance(root, clade)
            preorder_clades = list(tree_for_analysis.find_clades(order="preorder"))

        node_estimates = {}
        node_variances = {}

        for clade in reversed(preorder_clades):
            children = clade.clades
            if not children:
                name = clade.name
                if name in name_to_idx:
                    node_estimates[id(clade)] = scores[name_to_idx[name]]
                    node_variances[id(clade)] = 0.0
            else:
                total_prec = 0.0
                weighted_score = None
                for child in children:
                    child_id = id(child)
                    child_estimate = node_estimates.get(child_id)
                    if child_estimate is None:
                        continue
                    v_i = child.branch_length if child.branch_length else 0.0
                    child_var = node_variances[child_id]
                    denom = v_i + child_var
                    if denom == 0:
                        denom = 1e-10
                    prec = 1.0 / denom
                    if weighted_score is None:
                        weighted_score = prec * child_estimate
                    else:
                        weighted_score += prec * child_estimate
                    total_prec += prec

                if total_prec:
                    node_estimates[id(clade)] = weighted_score / total_prec
                    node_variances[id(clade)] = 1.0 / total_prec

        return node_estimates, node_distances, tree_for_analysis

    def _parse_color_by(
        self,
        color_by: str,
        trait_names: list[str],
        Y: np.ndarray,
        ordered_names: list[str],
    ) -> tuple[np.ndarray, list[str], str]:
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
                if not stripped or stripped[0] == "#":
                    continue
                taxon, sep, rest = stripped.partition("\t")
                if not sep:
                    continue
                idx = name_to_idx.get(taxon)
                if idx is not None:
                    values[idx] = rest.partition("\t")[0]

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
            return (
                np.fromiter(
                    (float(v) for v in values),
                    dtype=np.float64,
                    count=len(values),
                ),
                [],
                "continuous",
            )
        else:
            categories = sorted(set(values))
            return np.array(values), categories, "discrete"

    def _plot_phylomorphospace(
        self,
        tree,
        data: np.ndarray,
        ordered_names: list[str],
        node_estimates: dict,
        node_distances: dict,
        trait_x_name: str,
        trait_y_name: str,
        trait_names: list[str],
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
        max_dist = _root_distance_max(node_distances.values())

        # Draw tree edges colored by distance from root
        segments = []
        colors = []
        norm = Normalize(vmin=0, vmax=max_dist)
        cmap = plt.get_cmap("coolwarm")

        preorder_clades = self._preorder_clades_direct(tree)
        if preorder_clades is None:
            preorder_clades = tree.find_clades(order="preorder")

        for clade in preorder_clades:
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

    @staticmethod
    def _preorder_clades_direct(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        clades = []
        stack = [root]
        try:
            pop = stack.pop
            append = stack.append
            append_clade = clades.append
            while stack:
                clade = pop()
                append_clade(clade)
                children = clade.clades
                child_count = len(children)
                if child_count == 2:
                    append(children[1])
                    append(children[0])
                elif child_count:
                    for idx in range(child_count - 1, -1, -1):
                        append(children[idx])
        except AttributeError:
            return None
        return clades
