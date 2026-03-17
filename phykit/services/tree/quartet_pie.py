"""
Quartet pie chart visualization of gene tree concordance.

Draws a phylogram with pie charts at internal nodes showing the
proportion of gene trees supporting the species tree topology (gCF)
versus the two NNI alternative topologies (gDF1, gDF2). Supports
both native computation from gene trees and parsing of ASTRAL -t 2
or wASTRAL --support 3 annotations.
"""
import sys
from typing import Dict, List, Tuple

import numpy as np

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig, compute_node_x_cladogram
from ...helpers.quartet_utils import (
    compute_gcf_per_node,
    parse_astral_annotations,
)
from ...helpers.circular_layout import (
    compute_circular_coords,
    draw_circular_branches,
    draw_circular_tip_labels,
    draw_circular_colored_branch,
    draw_circular_colored_arc,
    circular_branch_points,
    radial_offset,
)
from ...helpers.color_annotations import (
    parse_color_file,
    resolve_mrca,
    draw_range_rect,
    draw_range_wedge,
    get_clade_branch_ids,
    build_color_legend_handles,
)
from ...errors import PhykitUserError


class QuartetPie(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.gene_trees_path = parsed["gene_trees_path"]
        self.output_path = parsed["output_path"]
        self.annotate = parsed["annotate"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        tree = self.read_tree_file()
        self._validate_tree(tree)

        if self.gene_trees_path:
            # Native mode: compute gCF from gene trees
            gene_trees = self._parse_gene_trees(self.gene_trees_path)
            proportions = compute_gcf_per_node(tree, gene_trees)
            input_mode = "native"
            n_gene_trees = len(gene_trees)
        else:
            # ASTRAL mode: parse q1/q2/q3 from node labels
            astral_props = parse_astral_annotations(tree)
            if not astral_props:
                raise PhykitUserError(
                    [
                        "No ASTRAL q1/q2/q3 annotations found in the tree.",
                        "Either provide gene trees with -g, or use an ASTRAL",
                        "-t 2 or wASTRAL --support 3 output tree with quartet",
                        "annotations.",
                    ],
                    code=2,
                )
            # Convert to same format: (q1, q2, q3, 0, 0, 0) — no raw counts
            proportions = {
                cid: (q1, q2, q3, 0, 0, 0)
                for cid, (q1, q2, q3) in astral_props.items()
            }
            input_mode = "astral"
            n_gene_trees = 0

        if self.plot_config.ladderize:
            tree.ladderize()

        self._plot_quartet_pie(tree, proportions, self.output_path)

        if self.json_output:
            self._print_json(tree, proportions, input_mode, n_gene_trees)
        else:
            print(f"Quartet pie chart saved: {self.output_path}")

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            gene_trees_path=getattr(args, "gene_trees", None),
            output_path=args.output,
            annotate=getattr(args, "annotate", False),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 4:
            raise PhykitUserError(
                ["Tree must have at least 4 tips for quartet analysis."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                clade.branch_length = 1e-8

    def _parse_gene_trees(self, path: str) -> list:
        from Bio import Phylo
        try:
            return list(Phylo.parse(path, "newick"))
        except Exception:
            raise PhykitUserError(
                [
                    f"Could not parse gene trees from {path}.",
                    "File should contain one Newick tree per line.",
                ],
                code=2,
            )

    def _plot_quartet_pie(
        self,
        tree,
        proportions: Dict[int, Tuple],
        output_path: str,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.patches import Patch
        except ImportError:
            print("matplotlib is required for quartet_pie. Install matplotlib and retry.")
            raise SystemExit(2)

        config = self.plot_config
        tips = list(tree.get_terminals())
        config.resolve(n_rows=len(tips), n_cols=None)

        default_colors = ["#2b8cbe", "#d62728", "#969696"]
        colors = config.merge_colors(default_colors)

        # Build parent map
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade

        # Compute node positions
        node_x = {}
        node_y = {}

        for i, tip in enumerate(tips):
            node_y[id(tip)] = i

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

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        if self.plot_config.circular:
            # --- Circular mode ---
            coords = compute_circular_coords(tree, node_x, parent_map)
            ax.set_aspect("equal")
            ax.axis("off")

            draw_circular_branches(ax, tree, coords, parent_map)
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize and config.ylabel_fontsize > 0 else 9
            max_x = max(node_x.values()) if node_x else 1.0
            draw_circular_tip_labels(ax, tree, coords, fontsize=label_fontsize, offset=max_x * 0.03)

            # Apply color annotations
            if self.plot_config.color_file:
                color_data = parse_color_file(self.plot_config.color_file)
                for taxa_list, color, label in color_data["ranges"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        draw_range_wedge(ax, tree, mrca, color, coords)
                for taxa_list, clade_color, label in color_data["clades"]:
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

            # Legend
            legend_handles = [
                Patch(facecolor=colors[0], edgecolor="black", linewidth=0.5,
                      label="Concordant (gCF / q1)"),
                Patch(facecolor=colors[1], edgecolor="black", linewidth=0.5,
                      label="Discordant alt 1 (gDF1 / q2)"),
                Patch(facecolor=colors[2], edgecolor="black", linewidth=0.5,
                      label="Discordant alt 2 (gDF2 / q3)"),
            ]
            if self.plot_config.color_file:
                color_legend = build_color_legend_handles(color_data)
                legend_handles.extend(color_legend)
            legend_loc = config.legend_position or "upper right"
            if legend_loc != "none":
                ax.legend(handles=legend_handles, loc=legend_loc, fontsize=8, frameon=True)

            if config.show_title:
                ax.set_title(
                    config.title or "Quartet Concordance Pie Chart",
                    fontsize=config.title_fontsize,
                )

            # Finalize layout BEFORE placing pie insets, so transData is stable
            fig.subplots_adjust(left=0.05, right=0.85, top=0.92, bottom=0.12)
            fig.canvas.draw()

            # Pie charts at internal nodes
            # Circular mode: scale with n_tips but stay larger than rectangular
            # since radial spacing gives more room between nodes
            n_tips = len(tips)
            pie_size = min(0.05, 0.6 / max(n_tips, 1))

            for clade in tree.find_clades(order="preorder"):
                if clade.is_terminal() or clade == root:
                    continue
                cid = id(clade)
                if cid not in proportions:
                    continue

                props = proportions[cid]
                gcf, gdf1, gdf2 = props[0], props[1], props[2]
                cx = coords[cid]["x"]
                cy = coords[cid]["y"]

                # Convert data coords to figure-fraction coords for the inset
                disp = ax.transData.transform((cx, cy))
                fig_coord = fig.transFigure.inverted().transform(disp)
                fx, fy = fig_coord

                # Create a small inset axes centered on the node
                inset = fig.add_axes(
                    [fx - pie_size / 2, fy - pie_size / 2, pie_size, pie_size],
                    zorder=10,
                )
                wedge_vals = [gcf, gdf1, gdf2]
                wedge_colors = [c for v, c in zip(wedge_vals, colors) if v > 1e-6]
                wedge_vals = [v for v in wedge_vals if v > 1e-6]
                if wedge_vals:
                    inset.pie(
                        wedge_vals, colors=wedge_colors, startangle=90,
                        wedgeprops={"edgecolor": "black", "linewidth": 0.5},
                    )
                inset.set_aspect("equal")
                inset.axis("off")

                # Annotate with values if requested
                if self.annotate:
                    ax.annotate(
                        f"{gcf:.2f}/{gdf1:.2f}/{gdf2:.2f}",
                        (cx, cy),
                        textcoords="offset points",
                        xytext=(8, 8),
                        fontsize=6,
                        color="black",
                        zorder=11,
                    )

            # Save without bbox_inches="tight" to preserve inset positions
            fig.savefig(output_path, dpi=config.dpi)
            plt.close(fig)

        else:
            # --- Rectangular mode ---
            # Draw branches
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

                ax.plot([x0, x1], [y1, y1], color="black", lw=1.5)
                ax.plot([x0, x0], [y0, y1], color="black", lw=1.5)

            # Tip labels
            max_x = max(node_x.values()) if node_x else 1.0
            offset = max_x * 0.03
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize and config.ylabel_fontsize > 0 else 9
            for tip in tips:
                ax.text(
                    node_x[id(tip)] + offset, node_y[id(tip)],
                    tip.name, va="center", fontsize=label_fontsize,
                )

            # Apply color annotations
            if self.plot_config.color_file:
                color_data = parse_color_file(self.plot_config.color_file)
                for taxa_list, color, label in color_data["ranges"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        draw_range_rect(ax, tree, mrca, color, node_x, node_y)
                for taxa_list, clade_color, label in color_data["clades"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        clade_ids = get_clade_branch_ids(tree, mrca, parent_map)
                        for cl in tree.find_clades(order="preorder"):
                            if cl == tree.root:
                                continue
                            if id(cl) in clade_ids and id(cl) in parent_map:
                                pid = id(parent_map[id(cl)])
                                cid = id(cl)
                                x0, x1 = node_x[pid], node_x[cid]
                                y0 = node_y.get(pid, 0)
                                y1 = node_y.get(cid, 0)
                                ax.plot([x0, x1], [y1, y1], color=clade_color, lw=1.5, zorder=2)
                                ax.plot([x0, x0], [y0, y1], color=clade_color, lw=1.5, zorder=2)
                for taxon, lbl_color in color_data["labels"].items():
                    for text_obj in ax.texts:
                        if text_obj.get_text() == taxon:
                            text_obj.set_color(lbl_color)
                            break

            # Legend
            legend_handles = [
                Patch(facecolor=colors[0], edgecolor="black", linewidth=0.5,
                      label="Concordant (gCF / q1)"),
                Patch(facecolor=colors[1], edgecolor="black", linewidth=0.5,
                      label="Discordant alt 1 (gDF1 / q2)"),
                Patch(facecolor=colors[2], edgecolor="black", linewidth=0.5,
                      label="Discordant alt 2 (gDF2 / q3)"),
            ]
            if self.plot_config.color_file:
                color_legend = build_color_legend_handles(color_data)
                legend_handles.extend(color_legend)
            legend_loc = config.legend_position or "upper right"
            if legend_loc != "none":
                ax.legend(handles=legend_handles, loc=legend_loc, fontsize=8, frameon=True)

            ax.set_xlabel("Branch length (subs/site)")
            ax.set_yticks([])
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)

            if config.show_title:
                ax.set_title(
                    config.title or "Quartet Concordance Pie Chart",
                    fontsize=config.title_fontsize,
                )
            if config.axis_fontsize:
                ax.xaxis.label.set_fontsize(config.axis_fontsize)

            # Finalize layout BEFORE placing pie insets, so transData is stable
            fig.subplots_adjust(left=0.05, right=0.85, top=0.92, bottom=0.12)
            fig.canvas.draw()

            # Pie charts at internal nodes — rendered as inset axes so they
            # appear as perfect circles regardless of axis scaling, and are
            # drawn above the phylogeny branches.
            n_tips = len(tips)
            pie_size = min(0.06, 0.8 / max(n_tips, 1))

            for clade in tree.find_clades(order="preorder"):
                if clade.is_terminal() or clade == root:
                    continue
                cid = id(clade)
                if cid not in proportions:
                    continue

                props = proportions[cid]
                gcf, gdf1, gdf2 = props[0], props[1], props[2]
                cx = node_x.get(cid, 0)
                cy = node_y.get(cid, 0)

                # Convert data coords to figure-fraction coords for the inset
                disp = ax.transData.transform((cx, cy))
                fig_coord = fig.transFigure.inverted().transform(disp)
                fx, fy = fig_coord

                # Create a small inset axes centered on the node
                inset = fig.add_axes(
                    [fx - pie_size / 2, fy - pie_size / 2, pie_size, pie_size],
                    zorder=10,
                )
                wedge_vals = [gcf, gdf1, gdf2]
                wedge_colors = [c for v, c in zip(wedge_vals, colors) if v > 1e-6]
                wedge_vals = [v for v in wedge_vals if v > 1e-6]
                if wedge_vals:
                    inset.pie(
                        wedge_vals, colors=wedge_colors, startangle=90,
                        wedgeprops={"edgecolor": "black", "linewidth": 0.5},
                    )
                inset.set_aspect("equal")
                inset.axis("off")

                # Annotate with values if requested
                if self.annotate:
                    ax.annotate(
                        f"{gcf:.2f}/{gdf1:.2f}/{gdf2:.2f}",
                        (cx, cy),
                        textcoords="offset points",
                        xytext=(8, 8),
                        fontsize=6,
                        color="black",
                        zorder=11,
                    )

            # Save without bbox_inches="tight" to preserve inset positions
            fig.savefig(output_path, dpi=config.dpi)
            plt.close(fig)

    def _print_json(self, tree, proportions, input_mode, n_gene_trees):
        nodes = []
        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            cid = id(clade)
            if cid not in proportions:
                continue
            props = proportions[cid]
            tip_names = sorted(t.name for t in clade.get_terminals())
            nodes.append({
                "node_tips": tip_names,
                "gCF": round(props[0], 4),
                "gDF1": round(props[1], 4),
                "gDF2": round(props[2], 4),
                "concordant_count": props[3],
                "disc1_count": props[4],
                "disc2_count": props[5],
            })

        payload = {
            "n_taxa": tree.count_terminals(),
            "n_gene_trees": n_gene_trees,
            "input_mode": input_mode,
            "nodes": nodes,
            "output_file": self.output_path,
        }
        print_json(payload)
