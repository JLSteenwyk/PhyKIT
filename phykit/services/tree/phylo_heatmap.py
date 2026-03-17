"""
Phylogenetic heatmap: display a phylogeny alongside a matrix of trait values.

Draws a phylogram on the left and a color-coded heatmap of numeric values
on the right, with rows aligned to tree tips. Analogous to R's
phytools::phylo.heatmap().
"""
import sys
from typing import Dict, List, Tuple

import numpy as np

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig, compute_node_x_cladogram
from ...helpers.color_annotations import (
    parse_color_file,
    resolve_mrca,
    draw_range_rect,
    draw_range_wedge,
    get_clade_branch_ids,
    build_color_legend_handles,
)
from ...errors import PhykitUserError


class PhyloHeatmap(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.data_path = parsed["data_path"]
        self.output_path = parsed["output_path"]
        self.split = parsed["split"]
        self.standardize = parsed["standardize"]
        self.cmap_name = parsed["cmap"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        tree = self.read_tree_file()
        self._validate_tree(tree)

        tree_tips = [t.name for t in tree.get_terminals()]
        trait_names, trait_data = self._parse_trait_matrix(
            self.data_path, tree_tips
        )

        # Prune tree to shared taxa
        shared_taxa = set(trait_data.keys())
        tips_to_prune = [t for t in tree_tips if t not in shared_taxa]
        if tips_to_prune:
            tree = self.prune_tree_using_taxa_list(tree, tips_to_prune)

        if self.plot_config.ladderize:
            tree.ladderize()

        # Get tip order from tree traversal
        tip_order = [t.name for t in tree.get_terminals()]

        self._plot_phylo_heatmap(
            tree, tip_order, trait_names, trait_data, self.output_path
        )

        if self.json_output:
            self._print_json(tip_order, trait_names, trait_data)
        else:
            print(f"Phylogenetic heatmap saved: {self.output_path}")

    def process_args(self, args) -> Dict:
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
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips."], code=2
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                clade.branch_length = 1e-8

    def _parse_trait_matrix(
        self, path: str, tree_tips: List[str]
    ) -> Tuple[List[str], Dict[str, List[float]]]:
        """Parse a multi-column numeric TSV file.

        Format: taxon<tab>col1<tab>col2<tab>...
        First non-comment line is the header with trait names.
        Returns (trait_names, {taxon: [val1, val2, ...]}).
        """
        try:
            with open(path) as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise PhykitUserError(
                [f"{path} not found. Check filename and path."], code=2
            )

        data_lines = []
        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            data_lines.append(stripped)

        if len(data_lines) < 2:
            raise PhykitUserError(
                ["Data file must have a header row and at least one data row."],
                code=2,
            )

        header = data_lines[0].split("\t")
        if len(header) < 2:
            raise PhykitUserError(
                ["Header must have at least 2 columns (taxon + 1 trait)."],
                code=2,
            )
        trait_names = header[1:]
        n_cols = len(header)

        trait_data = {}
        for line_idx, line in enumerate(data_lines[1:], 2):
            parts = line.split("\t")
            if len(parts) != n_cols:
                raise PhykitUserError(
                    [f"Line {line_idx}: expected {n_cols} columns, got {len(parts)}."],
                    code=2,
                )
            taxon = parts[0]
            try:
                values = [float(v) for v in parts[1:]]
            except ValueError:
                raise PhykitUserError(
                    [f"Line {line_idx}: non-numeric value for taxon '{taxon}'."],
                    code=2,
                )
            trait_data[taxon] = values

        # Validate shared taxa
        tree_tip_set = set(tree_tips)
        data_taxa_set = set(trait_data.keys())
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

    def _plot_phylo_heatmap(
        self,
        tree,
        tip_order: List[str],
        trait_names: List[str],
        trait_data: Dict[str, List[float]],
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
        matrix = np.array([trait_data[taxon] for taxon in tip_order])
        if self.standardize:
            col_means = np.nanmean(matrix, axis=0)
            col_stds = np.nanstd(matrix, axis=0)
            col_stds[col_stds == 0] = 1.0
            matrix = (matrix - col_means) / col_stds

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
        gs = gridspec.GridSpec(
            1, 2,
            width_ratios=[self.split, 1 - self.split],
            wspace=0.02,
        )

        ax_tree = fig.add_subplot(gs[0])
        ax_heat = fig.add_subplot(gs[1])

        # --- Draw phylogram on left panel ---
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade

        node_x = {}
        node_y = {}
        for i, tip_name in enumerate(tip_order):
            for tip in tree.get_terminals():
                if tip.name == tip_name:
                    node_y[id(tip)] = i
                    break

        root = tree.root
        if self.plot_config.cladogram:
            node_x = compute_node_x_cladogram(tree, parent_map)
        else:
            for clade in tree.find_clades(order="preorder"):
                if clade == root:
                    node_x[id(clade)] = 0.0
                elif id(clade) in parent_map:
                    parent = parent_map[id(clade)]
                    t = clade.branch_length if clade.branch_length else 0.0
                    node_x[id(clade)] = node_x.get(id(parent), 0.0) + t

        for clade in tree.find_clades(order="postorder"):
            if not clade.is_terminal() and id(clade) not in node_y:
                child_ys = [
                    node_y[id(c)] for c in clade.clades if id(c) in node_y
                ]
                if child_ys:
                    node_y[id(clade)] = np.mean(child_ys)
                else:
                    node_y[id(clade)] = 0.0

        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            if id(clade) not in parent_map:
                continue
            parent = parent_map[id(clade)]
            if id(parent) not in node_x or id(clade) not in node_x:
                continue
            x0 = node_x[id(parent)]
            x1 = node_x[id(clade)]
            y0 = node_y.get(id(parent), 0)
            y1 = node_y.get(id(clade), 0)
            ax_tree.plot([x0, x1], [y1, y1], color="black", lw=1.5)
            ax_tree.plot([x0, x0], [y0, y1], color="black", lw=1.5)

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
                    clade_ids = get_clade_branch_ids(tree, mrca, parent_map)
                    for cl in tree.find_clades(order="preorder"):
                        if cl == tree.root:
                            continue
                        if id(cl) in clade_ids and id(cl) in parent_map:
                            pid_val = id(parent_map[id(cl)])
                            cid_val = id(cl)
                            x0, x1 = node_x[pid_val], node_x[cid_val]
                            y0 = node_y.get(pid_val, 0)
                            y1 = node_y.get(cid_val, 0)
                            ax_tree.plot([x0, x1], [y1, y1], color=clade_color, lw=1.5, zorder=2)
                            ax_tree.plot([x0, x0], [y0, y1], color=clade_color, lw=1.5, zorder=2)

        ax_tree.set_ylim(-0.5, n_tips - 0.5)
        ax_tree.invert_yaxis()
        ax_tree.axis("off")

        # --- Draw heatmap on right panel ---
        im = ax_heat.imshow(
            matrix, aspect="auto", cmap=self.cmap_name,
            interpolation="nearest",
        )

        # Column labels on top
        ax_heat.set_xticks(np.arange(n_traits))
        xlabel_fs = config.xlabel_fontsize if config.xlabel_fontsize and config.xlabel_fontsize > 0 else 8
        ax_heat.set_xticklabels(trait_names, rotation=90, fontsize=xlabel_fs)
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
        from matplotlib.patches import Wedge
        import matplotlib.colors as mcolors
        from ...helpers.circular_layout import (
            compute_circular_coords,
            draw_circular_branches,
            draw_circular_tip_labels,
        )

        # Build parent_map and node_x
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade

        root = tree.root
        if config.cladogram:
            node_x = compute_node_x_cladogram(tree, parent_map)
        else:
            node_x = {}
            for clade in tree.find_clades(order="preorder"):
                if clade == root:
                    node_x[id(clade)] = 0.0
                elif id(clade) in parent_map:
                    parent = parent_map[id(clade)]
                    bl = clade.branch_length if clade.branch_length else 0.0
                    node_x[id(clade)] = node_x.get(id(parent), 0.0) + bl

        coords = compute_circular_coords(tree, node_x, parent_map)

        # Build tip name -> id lookup
        tip_name_to_id = {}
        for tip in tree.get_terminals():
            tip_name_to_id[tip.name] = id(tip)

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
            from ...helpers.circular_layout import draw_circular_colored_branch
            color_data = parse_color_file(self.plot_config.color_file)
            for taxa_list, clr, lbl in color_data["ranges"]:
                mrca = resolve_mrca(tree, taxa_list)
                if mrca is not None:
                    draw_range_wedge(ax, tree, mrca, clr, coords)
            for taxa_list, clade_color, lbl in color_data["clades"]:
                mrca = resolve_mrca(tree, taxa_list)
                if mrca is not None:
                    clade_ids = get_clade_branch_ids(tree, mrca, parent_map)
                    for cl in tree.find_clades(order="preorder"):
                        if cl == tree.root:
                            continue
                        if id(cl) in clade_ids and id(cl) in parent_map:
                            draw_circular_colored_branch(ax, coords[id(parent_map[id(cl)])], coords[id(cl)], clade_color, lw=1.5)
            for taxon, lbl_color in color_data["labels"].items():
                for text_obj in ax.texts:
                    if text_obj.get_text() == taxon:
                        text_obj.set_color(lbl_color)
                        break
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
                color = cmap(norm(value))
                wedge = Wedge(
                    (0, 0), r_outer, theta1, theta2,
                    width=ring_width,
                    facecolor=color, edgecolor="none",
                )
                ax.add_patch(wedge)

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
