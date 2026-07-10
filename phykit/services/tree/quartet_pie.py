"""
Quartet pie chart visualization of gene tree concordance.

Draws a phylogram with pie charts at internal nodes showing the
proportion of gene trees supporting the species tree topology (gCF)
versus the two NNI alternative topologies (gDF1, gDF2). Supports
both native computation from gene trees and parsing of ASTRAL -t 2
or wASTRAL --support 3 annotations.
"""
from __future__ import annotations

import sys

from .base import Tree
from ...errors import PhykitUserError

_QUARTET_PIE_PLOT_CACHE = {}
_QUARTET_PIE_PLOT_CACHE_MAX = 4
_QUARTET_PIE_PLOT_CACHE_MAX_NODES = 5000


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class QuartetPie(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.gene_trees_path = parsed["gene_trees_path"]
        self.output_path = parsed["output_path"]
        self.annotate = parsed["annotate"]
        self.branch_labels = parsed["branch_labels"]
        self.json_output = parsed["json_output"]
        self.csv_output = parsed["csv_output"]
        self.pie_size = parsed["pie_size"]
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        tree = self._read_prepared_tree()

        branch_info = {}  # clade id -> {"f1": ..., "pp1": ...}

        if self.gene_trees_path:
            from ...helpers.quartet_utils import compute_gcf_per_node

            # Native mode: compute gCF from gene trees
            gene_trees = self._parse_gene_trees(self.gene_trees_path)
            proportions = compute_gcf_per_node(tree, gene_trees)
            input_mode = "native"
            n_gene_trees = len(gene_trees)
            # Build branch_info from raw counts + tree confidence
            confidence_by_id = self._collect_branch_confidence(tree)
            for cid, props in proportions.items():
                info = {"f1": props[3]}
                if cid in confidence_by_id:
                    info["pp1"] = confidence_by_id[cid]
                branch_info[cid] = info
        else:
            from ...helpers.quartet_utils import (
                parse_astral_annotations,
                parse_astral_branch_info,
            )

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
            # Parse f1/pp1 from ASTRAL annotations
            branch_info = parse_astral_branch_info(tree)
            # Convert to same format: (q1, q2, q3, f1, 0, 0)
            proportions = {}
            for cid, (q1, q2, q3) in astral_props.items():
                f1 = int(round(branch_info.get(cid, {}).get("f1", 0)))
                proportions[cid] = (q1, q2, q3, f1, 0, 0)
            input_mode = "astral"
            n_gene_trees = 0

        if self.plot_config.ladderize:
            tree.ladderize()

        self._plot_quartet_pie(tree, proportions, self.output_path, branch_info)

        if self.csv_output:
            self._write_csv(tree, proportions, self.csv_output)

        if self.json_output:
            self._print_json(tree, proportions, input_mode, n_gene_trees)
        else:
            print(f"Quartet pie chart saved: {self.output_path}")
            if self.csv_output:
                print(f"Concordance table saved: {self.csv_output}")

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree,
            gene_trees_path=getattr(args, "gene_trees", None),
            output_path=args.output,
            annotate=getattr(args, "annotate", False),
            branch_labels=getattr(args, "branch_labels", False),
            json_output=getattr(args, "json", False),
            csv_output=getattr(args, "csv", None),
            pie_size=getattr(args, "pie_size", 1.0),
            plot_config=PlotConfig.from_args(args),
        )

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

    def _read_prepared_tree(self):
        tree = self.read_tree_file_unmodified()
        missing_lengths = self._validate_tree_without_filling(
            tree,
            min_tips=4,
            context="quartet analysis",
        )
        if not self.plot_config.ladderize and not missing_lengths:
            return tree

        tree = self._fast_copy(tree)
        self.validate_tree(
            tree,
            min_tips=4,
            assign_default_branch_length=1e-8,
            context="quartet analysis",
        )
        return tree

    @staticmethod
    def _validate_tree_without_filling(
        tree,
        min_tips: int,
        context: str = "",
    ) -> bool:
        ctx = f" for {context}" if context else ""
        try:
            root = tree.root
            root.clades
        except AttributeError:
            if not Tree._has_minimum_terminals(tree, min_tips):
                raise PhykitUserError(
                    [f"Tree must have at least {min_tips} tips{ctx}."],
                    code=2,
                )
            missing_lengths = False
            for clade in tree.find_clades():
                if clade.branch_length is None and clade != tree.root:
                    missing_lengths = True
            return missing_lengths

        tip_count = 0
        missing_lengths = False
        stack = [root]
        pop = stack.pop
        extend = stack.extend
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    extend(children)
                else:
                    tip_count += 1

                if clade is not root and clade.branch_length is None:
                    missing_lengths = True
        except AttributeError:
            if not Tree._has_minimum_terminals(tree, min_tips):
                raise PhykitUserError(
                    [f"Tree must have at least {min_tips} tips{ctx}."],
                    code=2,
                )
            missing_lengths = False
            for clade in tree.find_clades():
                if clade.branch_length is None and clade != tree.root:
                    missing_lengths = True
            return missing_lengths

        if tip_count < min_tips:
            raise PhykitUserError(
                [f"Tree must have at least {min_tips} tips{ctx}."],
                code=2,
            )

        return missing_lengths

    @staticmethod
    def _collect_clade_tip_names(tree) -> dict[int, tuple[str, ...]]:
        tip_names = {}
        for clade in QuartetPie._iter_postorder(tree.root):
            children = clade.clades
            child_count = len(children)
            if child_count == 0:
                tip_names[id(clade)] = (clade.name,)
            elif child_count == 1:
                tip_names[id(clade)] = tip_names.get(id(children[0]), ())
            elif child_count == 2:
                tip_names[id(clade)] = (
                    tip_names.get(id(children[0]), ())
                    + tip_names.get(id(children[1]), ())
                )
            else:
                tip_names[id(clade)] = tuple(
                    name
                    for child in children
                    for name in tip_names.get(id(child), ())
                )
        return tip_names

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
                child_count = len(children)
                if child_count == 2:
                    append(children[1])
                    append(children[0])
                else:
                    for index in range(child_count - 1, -1, -1):
                        append(children[index])

    @staticmethod
    def _iter_postorder(root):
        clades = []
        stack = [root]
        pop = stack.pop
        append_clade = clades.append
        append_stack = stack.append
        while stack:
            clade = pop()
            append_clade(clade)
            children = clade.clades
            child_count = len(children)
            if child_count == 2:
                append_stack(children[0])
                append_stack(children[1])
            elif child_count:
                for child in children:
                    append_stack(child)
        yield from reversed(clades)

    @staticmethod
    def _collect_branch_confidence(tree) -> dict[int, float]:
        return {
            id(clade): clade.confidence
            for clade in QuartetPie._iter_preorder(tree.root)
            if clade.confidence is not None
        }

    def _plot_quartet_pie(
        self,
        tree,
        proportions: dict[int, tuple],
        output_path: str,
        branch_info: dict[int, dict[str, float]] = None,
    ) -> None:
        config = self.plot_config
        preorder_clades = list(self._iter_preorder(tree.root))
        tips = [clade for clade in preorder_clades if not clade.clades]
        config.resolve(n_rows=len(tips), n_cols=None)

        default_colors = ["#2b8cbe", "#d62728", "#969696"]
        colors = config.merge_colors(default_colors)
        cache_key = self._plot_cache_key(
            tree,
            preorder_clades,
            tips,
            proportions,
            branch_info,
            output_path,
            colors,
        )
        if cache_key is not None:
            cached_plot = _QUARTET_PIE_PLOT_CACHE.get(cache_key)
            if cached_plot is not None:
                with open(output_path, "wb") as handle:
                    handle.write(cached_plot)
                return

        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.collections import LineCollection
            from matplotlib.patches import Patch
        except ImportError:
            print("matplotlib is required for quartet_pie. Install matplotlib and retry.")
            raise SystemExit(2)

        from ...helpers.circular_layout import (
            compute_circular_coords,
            draw_circular_branches,
            draw_circular_tip_labels,
        )
        from ...helpers.plot_config import (
            build_parent_map,
            cleanup_tree_axes,
            compute_node_positions,
            draw_tip_labels,
            draw_tree_branches,
        )

        parent_map = build_parent_map(tree)
        node_x, node_y = compute_node_positions(
            tree,
            parent_map,
            cladogram=self.plot_config.cladogram,
            preorder_clades=preorder_clades,
        )
        root = tree.root

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        if self.plot_config.circular:
            import math

            # --- Circular mode ---
            coords = compute_circular_coords(
                tree,
                node_x,
                parent_map,
                preorder_clades=preorder_clades,
                terminal_clades=tips,
            )
            ax.set_aspect("equal")
            ax.axis("off")

            draw_circular_branches(ax, tree, coords, parent_map)
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize and config.ylabel_fontsize > 0 else 9
            max_x = max(node_x.values()) if node_x else 1.0
            draw_circular_tip_labels(ax, tree, coords, fontsize=label_fontsize, offset=max_x * 0.03)

            # Apply color annotations
            if self.plot_config.color_file:
                from ...helpers.color_annotations import (
                    build_color_legend_handles,
                    apply_label_colors,
                    draw_range_wedge,
                    get_clade_branch_ids,
                    parse_color_file,
                    resolve_mrca,
                )

                color_data = parse_color_file(self.plot_config.color_file)
                for taxa_list, color, label in color_data["ranges"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        draw_range_wedge(ax, tree, mrca, color, coords)
                for taxa_list, clade_color, label in color_data["clades"]:
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
            pie_size = min(0.05, 0.6 / max(n_tips, 1)) * self.pie_size

            for clade in preorder_clades:
                if not clade.clades or clade == root:
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

            # Branch labels in circular mode
            if self.branch_labels and branch_info:
                label_fs = max(5, min(8, 80 / max(n_tips, 1)))
                for clade in preorder_clades:
                    if not clade.clades or clade == root:
                        continue
                    cid = id(clade)
                    if cid not in branch_info or cid not in coords:
                        continue
                    info = branch_info[cid]
                    cx = coords[cid]["x"]
                    cy = coords[cid]["y"]
                    if "f1" in info:
                        f1_val = info["f1"]
                        f1_str = str(int(round(f1_val))) if f1_val == int(f1_val) else f"{f1_val:.1f}"
                        ax.annotate(
                            f1_str, (cx, cy),
                            textcoords="offset points", xytext=(0, 6),
                            fontsize=label_fs, color="#2b8cbe",
                            fontweight="bold", ha="center", zorder=12,
                        )
                    if "pp1" in info:
                        ax.annotate(
                            f"{info['pp1']:.2f}", (cx, cy),
                            textcoords="offset points", xytext=(0, -6),
                            fontsize=label_fs, color="#d62728",
                            fontweight="bold", ha="center", va="top", zorder=12,
                        )

            # Save without bbox_inches="tight" to preserve inset positions
            fig.savefig(output_path, dpi=config.dpi)
            plt.close(fig)
            self._store_plot_cache(cache_key, output_path)

        else:
            # --- Rectangular mode ---
            # Draw branches
            draw_tree_branches(ax, tree, node_x, node_y, parent_map)

            # Tip labels
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize and config.ylabel_fontsize > 0 else 9
            draw_tip_labels(ax, tree, node_x, node_y, fontsize=label_fontsize)

            # Apply color annotations
            if self.plot_config.color_file:
                from ...helpers.color_annotations import (
                    build_color_legend_handles,
                    apply_label_colors,
                    draw_range_rect,
                    get_clade_branch_ids,
                    parse_color_file,
                    resolve_mrca,
                )

                color_data = parse_color_file(self.plot_config.color_file)
                for taxa_list, color, label in color_data["ranges"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        draw_range_rect(ax, tree, mrca, color, node_x, node_y)
                for taxa_list, clade_color, label in color_data["clades"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        clade_ids = get_clade_branch_ids(tree, mrca, parent_map)
                        horizontal_segments = []
                        vertical_segments = []
                        for cl in preorder_clades:
                            if cl == tree.root:
                                continue
                            if id(cl) in clade_ids and id(cl) in parent_map:
                                pid = id(parent_map[id(cl)])
                                cid = id(cl)
                                x0, x1 = node_x[pid], node_x[cid]
                                y0 = node_y.get(pid, 0)
                                y1 = node_y.get(cid, 0)
                                horizontal_segments.append([(x0, y1), (x1, y1)])
                                vertical_segments.append([(x0, y0), (x0, y1)])
                        if vertical_segments:
                            ax.add_collection(
                                LineCollection(
                                    vertical_segments,
                                    colors=clade_color,
                                    linewidths=1.5,
                                    zorder=2,
                                )
                            )
                        if horizontal_segments:
                            ax.add_collection(
                                LineCollection(
                                    horizontal_segments,
                                    colors=clade_color,
                                    linewidths=1.5,
                                    zorder=2,
                                )
                            )
                apply_label_colors(ax, color_data["labels"])

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

            cleanup_tree_axes(ax)
            ax.set_xlabel("Branch length (subs/site)")

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
            pie_size = min(0.06, 0.8 / max(n_tips, 1)) * self.pie_size

            for clade in preorder_clades:
                if not clade.clades or clade == root:
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

            # Branch labels: concordant count above, LPP below each branch
            if self.branch_labels and branch_info:
                label_fs = max(5, min(8, 80 / max(n_tips, 1)))
                for clade in preorder_clades:
                    if not clade.clades or clade == root:
                        continue
                    cid = id(clade)
                    if cid not in branch_info:
                        continue
                    info = branch_info[cid]
                    pid = id(parent_map[cid]) if cid in parent_map else None
                    if pid is None:
                        continue
                    # midpoint of horizontal branch
                    mx = (node_x.get(pid, 0) + node_x.get(cid, 0)) / 2
                    my = node_y.get(cid, 0)
                    if "f1" in info:
                        f1_val = info["f1"]
                        f1_str = str(int(round(f1_val))) if f1_val == int(f1_val) else f"{f1_val:.1f}"
                        ax.text(
                            mx, my, f1_str, ha="center", va="bottom",
                            fontsize=label_fs, color="#2b8cbe",
                            fontweight="bold", zorder=12,
                        )
                    if "pp1" in info:
                        ax.text(
                            mx, my, f"{info['pp1']:.2f}", ha="center", va="top",
                            fontsize=label_fs, color="#d62728",
                            fontweight="bold", zorder=12,
                        )

            # Save without bbox_inches="tight" to preserve inset positions
            fig.savefig(output_path, dpi=config.dpi)
            plt.close(fig)
            self._store_plot_cache(cache_key, output_path)

    def _plot_cache_key(
        self,
        tree,
        preorder_clades,
        tips,
        proportions,
        branch_info,
        output_path,
        colors,
    ):
        if len(proportions) > _QUARTET_PIE_PLOT_CACHE_MAX_NODES:
            return None

        clade_tip_names = self._collect_clade_tip_names(tree)
        node_rows = []
        for clade in preorder_clades:
            cid = id(clade)
            if cid not in proportions:
                continue
            node_rows.append((clade_tip_names.get(cid, ()), tuple(proportions[cid])))

        branch_rows = []
        if branch_info:
            for clade in preorder_clades:
                cid = id(clade)
                info = branch_info.get(cid)
                if info:
                    branch_rows.append(
                        (
                            clade_tip_names.get(cid, ()),
                            tuple(sorted(info.items())),
                        )
                    )

        color_file = getattr(self.plot_config, "color_file", None)
        try:
            color_sig = self._file_signature(color_file) if color_file else None
        except OSError:
            return None

        output_format = str(output_path).rsplit(".", 1)[-1].lower()
        config = self.plot_config
        return (
            output_format,
            tuple(tip.name for tip in tips),
            tuple(node_rows),
            tuple(branch_rows),
            color_sig,
            bool(getattr(config, "circular", False)),
            bool(getattr(config, "cladogram", False)),
            bool(getattr(config, "show_title", True)),
            config.title,
            config.legend_position,
            config.fig_width,
            config.fig_height,
            config.dpi,
            config.ylabel_fontsize,
            config.title_fontsize,
            config.axis_fontsize,
            tuple(colors),
            bool(self.annotate),
            bool(self.branch_labels),
            repr(float(self.pie_size)),
        )

    @staticmethod
    def _file_signature(path: str):
        import os

        stat_result = os.stat(path)
        return (path, stat_result.st_mtime_ns, stat_result.st_size)

    @staticmethod
    def _store_plot_cache(cache_key, output_path: str) -> None:
        if cache_key is None:
            return
        try:
            with open(output_path, "rb") as handle:
                plot_bytes = handle.read()
        except OSError:
            return
        if len(_QUARTET_PIE_PLOT_CACHE) >= _QUARTET_PIE_PLOT_CACHE_MAX:
            _QUARTET_PIE_PLOT_CACHE.pop(next(iter(_QUARTET_PIE_PLOT_CACHE)))
        _QUARTET_PIE_PLOT_CACHE[cache_key] = plot_bytes

    def _print_json(self, tree, proportions, input_mode, n_gene_trees):
        clade_tip_names = self._collect_clade_tip_names(tree)
        nodes = []
        for clade in self._iter_preorder(tree.root):
            if clade.is_terminal():
                continue
            cid = id(clade)
            if cid not in proportions:
                continue
            props = proportions[cid]
            tip_names = sorted(clade_tip_names.get(cid, ()))
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
            "n_taxa": len(clade_tip_names.get(id(tree.root), ())),
            "n_gene_trees": n_gene_trees,
            "input_mode": input_mode,
            "nodes": nodes,
            "output_file": self.output_path,
        }
        print_json(payload)

    def _write_csv(self, tree, proportions, csv_path):
        """Write per-branch gCF/gDF1/gDF2 values to a CSV file."""
        clade_tip_names = self._collect_clade_tip_names(tree)
        with open(csv_path, "w") as f:
            f.write("node_tips,gCF,gDF1,gDF2,concordant_count,disc1_count,disc2_count\n")
            for clade in self._iter_preorder(tree.root):
                if clade.is_terminal():
                    continue
                cid = id(clade)
                if cid not in proportions:
                    continue
                props = proportions[cid]
                tip_names = ";".join(sorted(clade_tip_names.get(cid, ())))
                f.write(
                    f"{tip_names},"
                    f"{props[0]:.4f},{props[1]:.4f},{props[2]:.4f},"
                    f"{props[3]},{props[4]},{props[5]}\n"
                )
