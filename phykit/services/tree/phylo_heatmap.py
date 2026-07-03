"""
Phylogenetic heatmap: display a phylogeny alongside a matrix of trait values.

Draws a phylogram on the left and a color-coded heatmap of numeric values
on the right, with rows aligned to tree tips. Analogous to R's
phytools::phylo.heatmap().
"""
from __future__ import annotations

import sys

from .base import Tree
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def build_parent_map(*args, **kwargs):
    from ...helpers.plot_config import build_parent_map as _build_parent_map

    return _build_parent_map(*args, **kwargs)


def compute_node_positions(*args, **kwargs):
    from ...helpers.plot_config import compute_node_positions as _compute_node_positions

    return _compute_node_positions(*args, **kwargs)


def draw_tree_branches(*args, **kwargs):
    from ...helpers.plot_config import draw_tree_branches as _draw_tree_branches

    return _draw_tree_branches(*args, **kwargs)


def parse_color_file(*args, **kwargs):
    from ...helpers.color_annotations import parse_color_file as _parse_color_file

    return _parse_color_file(*args, **kwargs)


def resolve_mrca(*args, **kwargs):
    from ...helpers.color_annotations import resolve_mrca as _resolve_mrca

    return _resolve_mrca(*args, **kwargs)


def draw_range_rect(*args, **kwargs):
    from ...helpers.color_annotations import draw_range_rect as _draw_range_rect

    return _draw_range_rect(*args, **kwargs)


def draw_range_wedge(*args, **kwargs):
    from ...helpers.color_annotations import draw_range_wedge as _draw_range_wedge

    return _draw_range_wedge(*args, **kwargs)


def get_clade_branch_ids(*args, **kwargs):
    from ...helpers.color_annotations import get_clade_branch_ids as _get_clade_branch_ids

    return _get_clade_branch_ids(*args, **kwargs)


def build_color_legend_handles(*args, **kwargs):
    from ...helpers.color_annotations import (
        build_color_legend_handles as _build_color_legend_handles,
    )

    return _build_color_legend_handles(*args, **kwargs)


def apply_label_colors(*args, **kwargs):
    from ...helpers.color_annotations import apply_label_colors as _apply_label_colors

    return _apply_label_colors(*args, **kwargs)


def compute_circular_coords(*args, **kwargs):
    from ...helpers.circular_layout import compute_circular_coords as _compute_circular_coords

    return _compute_circular_coords(*args, **kwargs)


def draw_circular_branches(*args, **kwargs):
    from ...helpers.circular_layout import draw_circular_branches as _draw_circular_branches

    return _draw_circular_branches(*args, **kwargs)


def draw_circular_tip_labels(*args, **kwargs):
    from ...helpers.circular_layout import draw_circular_tip_labels as _draw_circular_tip_labels

    return _draw_circular_tip_labels(*args, **kwargs)


class _LazyNumpy:
    _module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = _np
            self._module = module

        value = getattr(module, name)
        setattr(self, name, value)
        return value


np = _LazyNumpy()


class PhyloHeatmap(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.data_path = parsed["data_path"]
        self.output_path = parsed["output_path"]
        self.split = parsed["split"]
        self.standardize = parsed["standardize"]
        self.cmap_name = parsed["cmap"]
        self.cluster_columns = parsed["cluster_columns"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        copied_tree = False
        if self._needs_default_branch_lengths(tree):
            tree = self._fast_copy(tree)
            copied_tree = True
        self.validate_tree(tree, min_tips=3, assign_default_branch_length=1e-8, context="phylo heatmap")

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_names, trait_data = self._parse_trait_matrix(
            self.data_path, tree_tips
        )

        # Prune tree to shared taxa
        tips_to_prune = self._tips_to_prune_for_ordered_mapping(
            tree_tips, trait_data
        )
        if tips_to_prune:
            if not copied_tree:
                tree = self._fast_copy(tree)
                copied_tree = True
            tree = self.prune_tree_using_taxa_list(tree, tips_to_prune)

        if self.plot_config.ladderize:
            if not copied_tree:
                tree = self._fast_copy(tree)
            tree.ladderize()

        # Get tip order from tree traversal
        if tips_to_prune or self.plot_config.ladderize:
            tip_order = self.get_tip_names_from_tree(tree)
        else:
            tip_order = tree_tips

        self._plot_phylo_heatmap(
            tree, tip_order, trait_names, trait_data, self.output_path
        )

        if self.json_output:
            self._print_json(tip_order, trait_names, trait_data)
        else:
            print(f"Phylogenetic heatmap saved: {self.output_path}")

    @staticmethod
    def _needs_default_branch_lengths(tree) -> bool:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return True

        stack = [root]
        pop = stack.pop
        extend = stack.extend
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if clade is not root and clade.branch_length is None:
                    return True
                if children:
                    extend(children)
        except AttributeError:
            return True

        return False

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        split = getattr(args, "split", 0.3)
        if split <= 0 or split >= 1:
            raise PhykitUserError(
                ["--split must be between 0 and 1 (exclusive)."], code=2
            )
        return dict(
            tree_file_path=args.tree,
            data_path=args.data,
            output_path=args.output,
            split=split,
            standardize=getattr(args, "standardize", False),
            cmap=getattr(args, "cmap", "viridis"),
            cluster_columns=getattr(args, "cluster_columns", False),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def _parse_trait_matrix(
        self, path: str, tree_tips: list[str]
    ) -> tuple[list[str], dict[str, list[float]]]:
        """Parse a multi-column numeric TSV file.

        Format: taxon<tab>col1<tab>col2<tab>...
        First non-comment line is the header with trait names.
        Returns (trait_names, {taxon: [val1, val2, ...]}).
        """
        try:
            with open(path) as f:
                header = None
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped[0] == "#":
                        continue
                    header = stripped.split("\t")
                    break

                if header is None:
                    raise PhykitUserError(
                        ["Data file must have a header row and at least one data row."],
                        code=2,
                    )

                if len(header) < 2:
                    raise PhykitUserError(
                        ["Header must have at least 2 columns (taxon + 1 trait)."],
                        code=2,
                    )
                trait_names = header[1:]
                n_cols = len(header)

                trait_data = {}
                data_line_idx = 2
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped[0] == "#":
                        continue
                    parts = stripped.split("\t")
                    if len(parts) != n_cols:
                        raise PhykitUserError(
                            [
                                f"Line {data_line_idx}: expected {n_cols} columns, got {len(parts)}."
                            ],
                            code=2,
                        )
                    taxon = parts[0]
                    try:
                        values = list(map(float, parts[1:]))
                    except ValueError:
                        raise PhykitUserError(
                            [
                                f"Line {data_line_idx}: non-numeric value for taxon '{taxon}'."
                            ],
                            code=2,
                        )
                    trait_data[taxon] = values
                    data_line_idx += 1
        except FileNotFoundError:
            raise PhykitUserError(
                [f"{path} not found. Check filename and path."], code=2
            )

        if not trait_data:
            raise PhykitUserError(
                ["Data file must have a header row and at least one data row."],
                code=2,
            )

        # Validate shared taxa
        tree_tip_set = set(tree_tips)
        if (
            len(tree_tip_set) >= 3
            and len(tree_tip_set) == len(trait_data)
            and tree_tip_set == trait_data.keys()
        ):
            return trait_names, trait_data

        data_taxa_set = set(trait_data)
        shared = tree_tip_set & data_taxa_set

        if len(shared) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa between tree and data.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        tree_only = tree_tip_set - data_taxa_set
        data_only = data_taxa_set - tree_tip_set
        if tree_only:
            print(
                f"Warning: {len(tree_only)} taxa in tree but not in data: "
                f"{', '.join(sorted(tree_only))}",
                file=sys.stderr,
            )
        if data_only:
            print(
                f"Warning: {len(data_only)} taxa in data but not in tree: "
                f"{', '.join(sorted(data_only))}",
                file=sys.stderr,
            )

        return trait_names, {t: trait_data[t] for t in shared}

    @staticmethod
    def _build_heatmap_matrix(
        tip_order: list[str],
        trait_data: dict[str, list[float]],
    ) -> np.ndarray:
        return np.asarray([trait_data[taxon] for taxon in tip_order], dtype=float)

    @staticmethod
    def _standardize_heatmap_matrix(matrix: np.ndarray) -> np.ndarray:
        if np.isfinite(matrix).all():
            col_means = matrix.mean(axis=0)
            col_stds = matrix.std(axis=0)
        else:
            col_means = np.nanmean(matrix, axis=0)
            col_stds = np.nanstd(matrix, axis=0)
        col_stds[col_stds == 0] = 1.0
        return (matrix - col_means) / col_stds

    def _plot_phylo_heatmap(
        self,
        tree,
        tip_order: list[str],
        trait_names: list[str],
        trait_data: dict[str, list[float]],
        output_path: str,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib import gridspec
        except ImportError:
            print("matplotlib is required. Install matplotlib and retry.")
            raise SystemExit(2)

        config = self.plot_config
        n_tips = len(tip_order)
        n_traits = len(trait_names)
        config.resolve(n_rows=n_tips, n_cols=n_traits)

        # Build the data matrix in tip_order
        matrix = self._build_heatmap_matrix(tip_order, trait_data)
        if self.standardize:
            matrix = self._standardize_heatmap_matrix(matrix)

        if config.circular:
            self._plot_phylo_heatmap_circular(
                tree, tip_order, trait_names, trait_data, matrix,
                config, n_tips, n_traits, output_path, plt,
            )
        else:
            self._plot_phylo_heatmap_rect(
                tree, tip_order, trait_names, trait_data, matrix,
                config, n_tips, n_traits, output_path, plt, gridspec,
            )

    def _plot_phylo_heatmap_rect(
        self, tree, tip_order, trait_names, trait_data, matrix,
        config, n_tips, n_traits, output_path, plt, gridspec,
    ) -> None:
        fig = plt.figure(figsize=(config.fig_width, config.fig_height))

        # Cluster columns if requested
        col_order = list(range(n_traits))
        col_dendro_Z = None
        if self.cluster_columns and n_traits >= 2:
            from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
            from scipy.spatial.distance import pdist

            # Cluster traits by their column values (transpose, compute distance)
            col_dist = pdist(matrix.T, metric="euclidean")
            col_dendro_Z = linkage(col_dist, method="average")
            col_order = list(leaves_list(col_dendro_Z))
            matrix = matrix[:, col_order]
            trait_names = [trait_names[i] for i in col_order]

        if col_dendro_Z is not None:
            # 2x3 grid: tree | heatmap | colorbar, with dendrogram row on top
            gs = gridspec.GridSpec(
                2, 3,
                width_ratios=[self.split, (1 - self.split) * 0.92, (1 - self.split) * 0.08],
                height_ratios=[0.10, 0.90],
                wspace=0.02, hspace=0.0,
            )
            ax_tree = fig.add_subplot(gs[1, 0])
            ax_heat = fig.add_subplot(gs[1, 1])
            ax_cbar = fig.add_subplot(gs[1, 2])
            ax_col_dendro = fig.add_subplot(gs[0, 1])

            # Hide unused cells
            ax_empty_tl = fig.add_subplot(gs[0, 0])
            ax_empty_tl.axis("off")
            ax_empty_tr = fig.add_subplot(gs[0, 2])
            ax_empty_tr.axis("off")
        else:
            gs = gridspec.GridSpec(
                1, 2,
                width_ratios=[self.split, 1 - self.split],
                wspace=0.02,
            )
            ax_tree = fig.add_subplot(gs[0])
            ax_heat = fig.add_subplot(gs[1])

        # --- Draw phylogram on left panel ---
        parent_map = build_parent_map(tree)
        preorder_clades = list(self._iter_preorder(tree.root))
        node_x, node_y = compute_node_positions(
            tree,
            parent_map,
            cladogram=self.plot_config.cladogram,
            preorder_clades=preorder_clades,
        )
        draw_tree_branches(ax_tree, tree, node_x, node_y, parent_map)

        # Apply color annotations to tree panel
        if self.plot_config.color_file:
            color_data = parse_color_file(self.plot_config.color_file)
            for taxa_list, clr, lbl in color_data["ranges"]:
                mrca = resolve_mrca(tree, taxa_list)
                if mrca is not None:
                    draw_range_rect(ax_tree, tree, mrca, clr, node_x, node_y)
            for taxa_list, clade_color, lbl in color_data["clades"]:
                mrca = resolve_mrca(tree, taxa_list)
                if mrca is not None:
                    from matplotlib.collections import LineCollection

                    clade_ids = get_clade_branch_ids(tree, mrca, parent_map)
                    horizontal_segments = []
                    vertical_segments = []
                    for cl in preorder_clades:
                        if cl == tree.root:
                            continue
                        if id(cl) in clade_ids and id(cl) in parent_map:
                            pid_val = id(parent_map[id(cl)])
                            cid_val = id(cl)
                            x0, x1 = node_x[pid_val], node_x[cid_val]
                            y0 = node_y.get(pid_val, 0)
                            y1 = node_y.get(cid_val, 0)
                            horizontal_segments.append([(x0, y1), (x1, y1)])
                            vertical_segments.append([(x0, y0), (x0, y1)])
                    if vertical_segments:
                        ax_tree.add_collection(
                            LineCollection(
                                vertical_segments,
                                colors=clade_color,
                                linewidths=1.5,
                                zorder=2,
                            )
                        )
                    if horizontal_segments:
                        ax_tree.add_collection(
                            LineCollection(
                                horizontal_segments,
                                colors=clade_color,
                                linewidths=1.5,
                                zorder=2,
                            )
                        )

        ax_tree.set_ylim(-0.5, n_tips - 0.5)
        ax_tree.invert_yaxis()
        ax_tree.axis("off")

        # --- Draw heatmap on right panel ---
        im = ax_heat.imshow(
            matrix, aspect="auto", cmap=self.cmap_name,
            interpolation="nearest",
        )

        # Draw column dendrogram now (after imshow sets the x-axis range)
        if col_dendro_Z is not None:
            # imshow places columns at 0, 1, ..., n_traits-1
            # dendrogram places leaves at 5, 15, 25, ... (10*i + 5)
            # We need to draw dendrogram bars manually in heatmap x-coords
            # by converting: heatmap_x = (dendro_x - 5) / 10
            dendro_data = dendrogram(
                col_dendro_Z, ax=ax_col_dendro,
                no_labels=True, no_plot=True,
            )
            # Draw the dendrogram lines manually in the correct coordinates
            for xs, ys in zip(dendro_data["icoord"], dendro_data["dcoord"]):
                # Convert dendrogram x-coords to heatmap column indices
                hx = [(x - 5.0) / 10.0 for x in xs]
                ax_col_dendro.plot(hx, ys, color="#333", lw=1.0)
            ax_col_dendro.set_xlim(-0.5, n_traits - 0.5)
            ax_col_dendro.axis("off")

        # Column labels: bottom if dendrogram on top, otherwise top
        ax_heat.set_xticks(np.arange(n_traits))
        xlabel_fs = config.xlabel_fontsize if config.xlabel_fontsize and config.xlabel_fontsize > 0 else 8
        ax_heat.set_xticklabels(trait_names, rotation=90, fontsize=xlabel_fs)
        if col_dendro_Z is not None:
            ax_heat.xaxis.set_ticks_position("bottom")
            ax_heat.xaxis.set_label_position("bottom")
        else:
            ax_heat.xaxis.set_ticks_position("top")
            ax_heat.xaxis.set_label_position("top")

        # Row labels (taxon names) on right side of heatmap
        ax_heat.set_yticks(np.arange(n_tips))
        ylabel_fs = config.ylabel_fontsize if config.ylabel_fontsize and config.ylabel_fontsize > 0 else 8
        if ylabel_fs > 0:
            ax_heat.set_yticklabels(tip_order, fontsize=ylabel_fs)
            ax_heat.yaxis.set_ticks_position("right")
            ax_heat.yaxis.set_label_position("right")
        else:
            ax_heat.set_yticklabels([])

        # Colorbar
        if col_dendro_Z is not None:
            # Use dedicated colorbar axes (no space stealing)
            cbar = fig.colorbar(im, cax=ax_cbar)
        else:
            cbar = fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.12)
        if self.standardize:
            cbar.set_label("Z-score")
        else:
            cbar.set_label("Value")

        # Apply label color annotations to heatmap y-tick labels
        if self.plot_config.color_file:
            color_data = parse_color_file(self.plot_config.color_file)
            for label_obj in ax_heat.get_yticklabels():
                taxon = label_obj.get_text()
                if taxon in color_data["labels"]:
                    label_obj.set_color(color_data["labels"][taxon])
            color_legend = build_color_legend_handles(color_data)
            if color_legend:
                ax_tree.legend(handles=color_legend, loc="upper right", fontsize=8, frameon=True)

        # Title
        if config.show_title:
            fig.suptitle(
                config.title or "Phylogenetic Heatmap",
                fontsize=config.title_fontsize,
                y=1.02,
            )

        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def _plot_phylo_heatmap_circular(
        self, tree, tip_order, trait_names, trait_data, matrix,
        config, n_tips, n_traits, output_path, plt,
    ) -> None:
        import math
        from matplotlib.collections import LineCollection, PatchCollection
        from matplotlib.patches import Wedge
        import matplotlib.colors as mcolors

        # Build parent_map and node_x
        parent_map = build_parent_map(tree)
        preorder_clades = list(self._iter_preorder(tree.root))
        tips = [clade for clade in preorder_clades if not clade.clades]
        node_x, _node_y = compute_node_positions(
            tree,
            parent_map,
            cladogram=config.cladogram,
            preorder_clades=preorder_clades,
        )

        coords = compute_circular_coords(
            tree,
            node_x,
            parent_map,
            preorder_clades=preorder_clades,
            terminal_clades=tips,
        )

        # Build tip name -> id lookup
        tip_name_to_id = {
            tip.name: id(tip)
            for tip in preorder_clades
            if not tip.clades
        }

        # Compute max_radius from tips
        max_radius = max(
            coords[tip_name_to_id[name]]["radius"]
            for name in tip_order
            if name in tip_name_to_id
        )

        # Ring parameters
        gap = max_radius * 0.15
        ring_width = max_radius * 0.08

        # Single axes figure
        fig, ax = plt.subplots(
            1, 1, figsize=(config.fig_width, config.fig_height)
        )
        ax.set_aspect("equal")
        ax.axis("off")

        # Draw circular tree
        draw_circular_branches(ax, tree, coords, parent_map)
        draw_circular_tip_labels(ax, tree, coords, fontsize=7,
                                 offset=max_radius * 0.03)

        # Apply color annotations
        if self.plot_config.color_file:
            color_data = parse_color_file(self.plot_config.color_file)
            for taxa_list, clr, lbl in color_data["ranges"]:
                mrca = resolve_mrca(tree, taxa_list)
                if mrca is not None:
                    draw_range_wedge(ax, tree, mrca, clr, coords)
            for taxa_list, clade_color, lbl in color_data["clades"]:
                mrca = resolve_mrca(tree, taxa_list)
                if mrca is not None:
                    clade_ids = get_clade_branch_ids(tree, mrca, parent_map)
                    radial_segments = []
                    for cl in preorder_clades:
                        if cl == tree.root:
                            continue
                        if id(cl) in clade_ids and id(cl) in parent_map:
                            parent_coords = coords[id(parent_map[id(cl)])]
                            child_coords = coords[id(cl)]
                            angle = child_coords["angle"]
                            cos_angle = math.cos(angle)
                            sin_angle = math.sin(angle)
                            r_parent = parent_coords["radius"]
                            r_child = child_coords["radius"]
                            radial_segments.append(
                                [
                                    (r_parent * cos_angle, r_parent * sin_angle),
                                    (r_child * cos_angle, r_child * sin_angle),
                                ]
                            )
                    if radial_segments:
                        ax.add_collection(
                            LineCollection(
                                radial_segments,
                                colors=clade_color,
                                linewidths=1.5,
                                capstyle="round",
                                zorder=2,
                            )
                        )
            apply_label_colors(ax, color_data["labels"])
            color_legend = build_color_legend_handles(color_data)
            if color_legend:
                ax.legend(handles=color_legend, loc="upper right", fontsize=8, frameon=True)

        # Build colormap for heatmap
        cmap = plt.get_cmap(self.cmap_name)
        vmin = float(np.nanmin(matrix))
        vmax = float(np.nanmax(matrix))
        if vmin == vmax:
            vmax = vmin + 1.0
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

        # Draw heatmap rings
        wedge_half = (2.0 * math.pi / n_tips) / 2.0 * 0.9  # 90% of slot
        heatmap_wedges = []
        heatmap_values = []

        for col_idx in range(n_traits):
            r_inner = max_radius + gap + col_idx * ring_width
            r_outer = r_inner + ring_width
            for row_idx, taxon in enumerate(tip_order):
                if taxon not in tip_name_to_id:
                    continue
                tid = tip_name_to_id[taxon]
                tip_angle = coords[tid]["angle"]
                theta1 = math.degrees(tip_angle - wedge_half)
                theta2 = math.degrees(tip_angle + wedge_half)
                value = matrix[row_idx, col_idx]
                wedge = Wedge(
                    (0, 0), r_outer, theta1, theta2,
                    width=ring_width,
                )
                heatmap_wedges.append(wedge)
                heatmap_values.append(value)
        if heatmap_wedges:
            heatmap_collection = PatchCollection(
                heatmap_wedges,
                cmap=cmap,
                norm=norm,
                edgecolors="none",
                match_original=False,
            )
            heatmap_collection.set_array(np.asarray(heatmap_values, dtype=float))
            ax.add_collection(heatmap_collection)

        # Trait name labels at outside of each ring
        label_fontsize = 6
        for col_idx, trait_name in enumerate(trait_names):
            r_label = max_radius + gap + col_idx * ring_width + ring_width / 2.0
            ax.text(
                0, r_label, trait_name,
                fontsize=label_fontsize, ha="center", va="bottom",
                rotation=0,
            )

        # Colorbar
        import matplotlib.cm as mcm
        sm = mcm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04, shrink=0.6)
        if self.standardize:
            cbar.set_label("Z-score")
        else:
            cbar.set_label("Value")

        # Title
        if config.show_title:
            fig.suptitle(
                config.title or "Phylogenetic Heatmap",
                fontsize=config.title_fontsize,
                y=1.02,
            )

        ax.autoscale_view()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    @staticmethod
    def _iter_preorder(root):
        stack = [root]
        pop = stack.pop
        append = stack.append
        while stack:
            clade = pop()
            yield clade
            children = clade.clades
            if children:
                append(children[-1])
                if len(children) == 2:
                    append(children[0])
                else:
                    for idx in range(len(children) - 2, -1, -1):
                        append(children[idx])

    def _print_json(self, tip_order, trait_names, trait_data):
        payload = {
            "n_taxa": len(tip_order),
            "n_traits": len(trait_names),
            "trait_names": trait_names,
            "tip_order": tip_order,
            "standardized": self.standardize,
            "output_file": self.output_path,
        }
        print_json(payload)
