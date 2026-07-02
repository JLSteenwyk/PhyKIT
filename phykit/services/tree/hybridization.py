from __future__ import annotations

import os
from io import StringIO
from pathlib import Path

from .base import Tree
from ...errors import PhykitUserError

_path_isabs = os.path.isabs


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()
_FDR_VECTOR_MIN_LENGTH = 2048


def _binomial_two_sided_p_value(successes: int, total: int) -> float:
    from scipy.special import bdtr

    tail_count = min(successes, total - successes)
    return min(1.0, 2.0 * float(bdtr(tail_count, total, 0.5)))


class Hybridization(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.gene_trees_path = parsed["gene_trees_path"]
        self.support_threshold = parsed["support_threshold"]
        self.alpha = parsed["alpha"]
        self.json_output = parsed["json_output"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree,
            gene_trees_path=args.gene_trees,
            support_threshold=getattr(args, "support", None),
            alpha=getattr(args, "alpha", 0.05),
            json_output=getattr(args, "json", False),
            plot_output=getattr(args, "plot_output", None),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        species_tree = self.read_tree_file_unmodified()
        gene_trees = self._parse_gene_trees(self.gene_trees_path)

        topology_counts, shared_taxa = self._count_topologies(
            species_tree, gene_trees, self.support_threshold
        )

        # Test each branch and collect results
        branch_results = []
        for branch_key in sorted(topology_counts.keys()):
            data = topology_counts[branch_key]
            total = data["n_concordant"] + data["n_alt1"] + data["n_alt2"]
            gcf = data["n_concordant"] / total if total > 0 else 1.0

            test_result = self._test_asymmetry(data["n_alt1"], data["n_alt2"])

            entry = dict(
                taxa=data["split"],
                n_concordant=data["n_concordant"],
                n_alt1=data["n_alt1"],
                n_alt2=data["n_alt2"],
                gcf=round(gcf, 4),
            )
            entry.update(test_result)
            branch_results.append(entry)

        # FDR correction across testable p-values
        testable_indices = []
        testable_pvals = []
        for i, entry in enumerate(branch_results):
            if entry["p_value"] is not None:
                testable_indices.append(i)
                testable_pvals.append(entry["p_value"])

        fdr_corrected = self._fdr(testable_pvals)
        for idx, fdr_p in zip(testable_indices, fdr_corrected):
            branch_results[idx]["fdr_p"] = fdr_p

        # Set fdr_p to None for untestable branches
        for entry in branch_results:
            if "fdr_p" not in entry:
                entry["fdr_p"] = None

        # Compute hybrid score and significance per branch
        alpha = self.alpha
        for entry in branch_results:
            if entry["fdr_p"] is not None and entry["fdr_p"] < alpha:
                entry["significant"] = True
                entry["hybrid_score"] = entry["asymmetry_ratio"]
            else:
                entry["significant"] = False
                entry["hybrid_score"] = 0.0

        n_reticulations = sum(1 for e in branch_results if e["significant"])

        # Summary
        summary = dict(
            n_branches=len(branch_results),
            n_gene_trees=len(gene_trees),
            support_threshold=self.support_threshold,
            alpha=alpha,
            n_reticulations=n_reticulations,
        )

        # Output
        if self.json_output:
            self._output_json(branch_results, summary)
        else:
            self._output_text(branch_results, summary)

        if self.plot_output:
            plot_tree = species_tree
            if self.plot_config.ladderize:
                plot_tree = self.read_tree_file()
                plot_tree.ladderize()
            self._plot(plot_tree, branch_results, self.plot_output,
                       shared_taxa=shared_taxa)

    # ------------------------------------------------------------------
    # Output methods
    # ------------------------------------------------------------------

    def _output_text(self, branch_results, summary) -> None:
        lines = [
            "Hybridization Analysis",
            "======================",
            f"Species tree branches: {summary['n_branches']}",
            f"Gene trees: {summary['n_gene_trees']}",
        ]
        if summary["support_threshold"] is not None:
            lines.append(f"Support threshold: {summary['support_threshold']}")
        lines.extend([
            "",
            f"Estimated reticulation events: {summary['n_reticulations']}",
            "",
        ])

        row_format = (
            "  {:<35}{:>8.2f}{:>8.2f}{:>8.2f}{:>8}{:>10}{:>12}"
        ).format
        append_line = lines.append
        header_written = False
        nonsig_count = 0

        for entry in branch_results:
            if not entry["significant"]:
                nonsig_count += 1
                continue
            if not header_written:
                append_line(f"Significant branches (FDR < {summary['alpha']}):")
                append_line(
                    f"  {'Branch (taxa)':<35}"
                    f"{'gCF':>8}"
                    f"{'gDF1':>8}"
                    f"{'gDF2':>8}"
                    f"{'Ratio':>8}"
                    f"{'FDR_p':>10}"
                    f"{'Favored':>12}"
                )
                header_written = True
            branch_label = ",".join(entry["taxa"])
            total = entry["n_concordant"] + entry["n_alt1"] + entry["n_alt2"]
            gdf1 = entry["n_alt1"] / total if total > 0 else 0.0
            gdf2 = entry["n_alt2"] / total if total > 0 else 0.0
            ratio_str = (
                f"{entry['asymmetry_ratio']:.2f}"
                if entry["asymmetry_ratio"] is not None
                else "NA"
            )
            fdr_str = (
                f"{entry['fdr_p']:.4f}"
                if entry["fdr_p"] is not None
                else "NA"
            )
            favored = entry.get("favored_alt", "-") or "-"
            append_line(
                row_format(
                    branch_label,
                    entry["gcf"],
                    gdf1,
                    gdf2,
                    ratio_str,
                    fdr_str,
                    favored,
                )
            )
        if header_written:
            append_line("")

        append_line(f"Non-significant branches: {nonsig_count}")

        try:
            print("\n".join(lines))
        except BrokenPipeError:
            pass

    def _output_json(self, branch_results, summary) -> None:
        result = dict(
            n_branches=summary["n_branches"],
            n_gene_trees=summary["n_gene_trees"],
            support_threshold=summary["support_threshold"],
            alpha=summary["alpha"],
            n_reticulations=summary["n_reticulations"],
            branches=branch_results,
        )
        print_json(result)

    def _plot(self, species_tree, branch_results, output_path,
              shared_taxa=None) -> None:
        """Phylogram colored by hybridization score at each branch."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        import numpy as np

        # Build lookup from taxa label -> branch result
        branch_lookup = {}
        for entry in branch_results:
            key = ",".join(entry["taxa"])
            branch_lookup[key] = entry

        parent_map = self._build_parent_map(species_tree)
        tips = self._get_terminal_clades(species_tree)
        preorder_clades = self._preorder_clades(species_tree)
        postorder_clades = self._postorder_clades(species_tree)
        all_taxa_fs = (shared_taxa if shared_taxa is not None
                       else frozenset(t.name for t in tips))

        # Compute node positions
        node_x = {}
        node_y = {}

        for i, tip in enumerate(tips):
            node_y[id(tip)] = i

        root = species_tree.root
        if self.plot_config.cladogram:
            from ...helpers.plot_config import compute_node_x_cladogram

            node_x = compute_node_x_cladogram(species_tree, parent_map)
        else:
            for clade in preorder_clades:
                if clade == root:
                    node_x[id(clade)] = 0.0
                else:
                    if id(clade) in parent_map:
                        parent = parent_map[id(clade)]
                        t = clade.branch_length if clade.branch_length else 0.0
                        node_x[id(clade)] = node_x.get(id(parent), 0.0) + t

        for clade in postorder_clades:
            if not clade.is_terminal() and id(clade) not in node_y:
                child_ys = [
                    node_y[id(c)] for c in clade.clades if id(c) in node_y
                ]
                if child_ys:
                    node_y[id(clade)] = sum(child_ys) / len(child_ys)
                else:
                    node_y[id(clade)] = 0.0

        # Map internal nodes to their branch result
        clade_taxa = self._collect_clade_taxa(species_tree)
        node_to_result = {}
        for clade in preorder_clades:
            if clade.is_terminal():
                continue
            node_tips = clade_taxa.get(id(clade), frozenset()) & all_taxa_fs
            split_label = (
                sorted(node_tips)
                if len(node_tips) <= len(all_taxa_fs) - len(node_tips)
                else sorted(all_taxa_fs - node_tips)
            )
            key = ",".join(split_label)
            if key in branch_lookup:
                node_to_result[id(clade)] = branch_lookup[key]

        # Color setup: gray for 0, warm gradient for significant branches
        # Use OrRd colormap for significant branches (orange to red)
        cmap = plt.cm.OrRd
        norm = Normalize(vmin=0.5, vmax=1.0)

        config = self.plot_config
        config.resolve(n_rows=len(tips), n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        if self.plot_config.circular:
            self._plot_circular(ax, fig, species_tree, root, tips, parent_map,
                                node_x, node_y, node_to_result, cmap, norm, config,
                                preorder_clades)
        else:
            self._plot_rectangular(ax, fig, species_tree, root, tips, parent_map,
                                   node_x, node_y, node_to_result, cmap, norm, config,
                                   preorder_clades)

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        import matplotlib.pyplot as plt2
        plt.close(fig)

    def _plot_circular(self, ax, fig, species_tree, root, tips, parent_map,
                       node_x, node_y, node_to_result, cmap, norm, config,
                       preorder_clades=None):
        import math
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.collections import LineCollection
        from ...helpers.circular_layout import (
            compute_circular_coords,
            draw_circular_tip_labels,
        )

        coords = compute_circular_coords(species_tree, node_x, parent_map)
        ax.set_aspect("equal")
        ax.axis("off")

        radial_segments = []
        radial_colors = []
        radial_linewidths = []
        score_color_cache = {}
        if preorder_clades is None:
            preorder_clades = self._preorder_clades(species_tree)

        for clade in preorder_clades:
            if clade == root:
                continue
            cid = id(clade)
            if cid not in parent_map:
                continue
            parent = parent_map[cid]
            if id(parent) not in coords or cid not in coords:
                continue

            color = "gray"
            lw = 2
            if cid in node_to_result:
                entry = node_to_result[cid]
                score = entry["hybrid_score"]
                if score > 0:
                    color = score_color_cache.get(score)
                    if color is None:
                        color = cmap(norm(score))
                        score_color_cache[score] = color
                    lw = 3

            parent_coords = coords[id(parent)]
            child_coords = coords[cid]
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
            radial_colors.append(color)
            radial_linewidths.append(lw)

        if radial_segments:
            ax.add_collection(
                LineCollection(
                    radial_segments,
                    colors=radial_colors,
                    linewidths=radial_linewidths,
                    capstyle="round",
                    zorder=2,
                )
            )

        arc_segments = []
        arc_colors = []
        arc_fractions = [idx / 60 for idx in range(61)]
        for clade in preorder_clades:
            if clade.is_terminal() or not clade.clades:
                continue
            cid = id(clade)
            pc = coords[cid]
            child_angles = [coords[id(ch)]["angle"] for ch in clade.clades]
            if len(child_angles) < 2:
                continue
            start_a = min(child_angles)
            end_a = max(child_angles)
            span = (end_a - start_a) % (2.0 * math.pi)
            if span > math.pi:
                start_a, end_a = end_a, start_a

            arc_color = "gray"
            if cid in node_to_result:
                entry = node_to_result[cid]
                score = entry["hybrid_score"]
                if score > 0:
                    arc_color = score_color_cache.get(score)
                    if arc_color is None:
                        arc_color = cmap(norm(score))
                        score_color_cache[score] = arc_color
            start = start_a % (2.0 * math.pi)
            end = end_a % (2.0 * math.pi)
            diff = (end - start) % (2.0 * math.pi)
            if diff > math.pi:
                diff -= 2.0 * math.pi
            arc_segments.append(
                [
                    (
                        pc["radius"] * math.cos(start + diff * fraction),
                        pc["radius"] * math.sin(start + diff * fraction),
                    )
                    for fraction in arc_fractions
                ]
            )
            arc_colors.append(arc_color)

        if arc_segments:
            ax.add_collection(
                LineCollection(
                    arc_segments,
                    colors=arc_colors,
                    linewidths=1.5,
                    capstyle="round",
                    zorder=1,
                )
            )

        if radial_segments or arc_segments:
            ax.autoscale_view()

        # Tip labels
        max_x = max(node_x.values()) if node_x else 1.0
        label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize is not None else 9
        if label_fontsize > 0:
            draw_circular_tip_labels(ax, species_tree, coords, fontsize=label_fontsize, offset=max_x * 0.02)

        # Color annotations
        if self.plot_config.color_file:
            from ...helpers.color_annotations import (
                build_color_legend_handles,
                apply_label_colors,
                draw_range_wedge,
                parse_color_file,
                resolve_mrca,
            )

            color_data = parse_color_file(self.plot_config.color_file)
            for taxa_list, clr, lbl in color_data["ranges"]:
                mrca = resolve_mrca(species_tree, taxa_list)
                if mrca is not None:
                    draw_range_wedge(ax, species_tree, mrca, clr, coords)
            apply_label_colors(ax, color_data["labels"])
            color_legend = build_color_legend_handles(color_data)
            if color_legend:
                ax.legend(handles=color_legend, loc="upper right", fontsize=8, frameon=True)

        # Mark significant branches with stars
        star_x = []
        star_y = []
        for clade in preorder_clades:
            if clade.is_terminal():
                continue
            cid = id(clade)
            if cid not in node_to_result:
                continue
            entry = node_to_result[cid]
            if cid not in coords:
                continue
            cx = coords[cid]["x"]
            cy = coords[cid]["y"]

            if entry["significant"]:
                star_x.append(cx)
                star_y.append(cy)
        if star_x:
            ax.scatter(star_x, star_y, s=100, c="red", marker="*", zorder=5)

        # Colorbar
        legend_loc = config.legend_position or "upper right"
        if legend_loc != "none":
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.15)
            cbar_fontsize = config.axis_fontsize if config.axis_fontsize else 10
            cbar.set_label("Hybridization score", fontsize=cbar_fontsize)

        if config.show_title:
            ax.set_title(config.title or "Hybridization Analysis", fontsize=config.title_fontsize)

    def _plot_rectangular(self, ax, fig, species_tree, root, tips, parent_map,
                          node_x, node_y, node_to_result, cmap, norm, config,
                          preorder_clades=None):
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.collections import LineCollection

        horizontal_segments = []
        horizontal_colors = []
        horizontal_linewidths = []
        vertical_segments = []
        score_color_cache = {}
        if preorder_clades is None:
            preorder_clades = self._preorder_clades(species_tree)

        for clade in preorder_clades:
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

            color = "gray"
            lw = 2
            if id(clade) in node_to_result:
                entry = node_to_result[id(clade)]
                score = entry["hybrid_score"]
                if score > 0:
                    color = score_color_cache.get(score)
                    if color is None:
                        color = cmap(norm(score))
                        score_color_cache[score] = color
                    lw = 3

            horizontal_segments.append([(x0, y1), (x1, y1)])
            horizontal_colors.append(color)
            horizontal_linewidths.append(lw)
            vertical_segments.append([(x0, y0), (x0, y1)])

        if vertical_segments:
            ax.add_collection(
                LineCollection(
                    vertical_segments,
                    colors="gray",
                    linewidths=1.5,
                    zorder=1,
                )
            )
        if horizontal_segments:
            ax.add_collection(
                LineCollection(
                    horizontal_segments,
                    colors=horizontal_colors,
                    linewidths=horizontal_linewidths,
                    zorder=2,
                )
            )
        if vertical_segments or horizontal_segments:
            ax.autoscale_view()

        # Mark significant branches with stars
        star_x = []
        star_y = []
        for clade in preorder_clades:
            if clade.is_terminal():
                continue
            if id(clade) not in node_to_result:
                continue
            entry = node_to_result[id(clade)]
            x = node_x.get(id(clade), 0)
            y = node_y.get(id(clade), 0)

            if entry["significant"]:
                star_x.append(x)
                star_y.append(y)
        if star_x:
            ax.scatter(star_x, star_y, s=100, c="red", marker="*", zorder=5)

        # Tip labels
        label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize is not None else 9
        if label_fontsize > 0:
            max_x = max(node_x.values()) if node_x else 0
            offset = max_x * 0.02
            for tip in tips:
                ax.text(
                    node_x[id(tip)] + offset, node_y[id(tip)],
                    tip.name, va="center", fontsize=label_fontsize,
                )

        # Color annotations
        if self.plot_config.color_file:
            from ...helpers.color_annotations import (
                build_color_legend_handles,
                apply_label_colors,
                draw_range_rect,
                parse_color_file,
                resolve_mrca,
            )

            color_data = parse_color_file(self.plot_config.color_file)
            for taxa_list, clr, lbl in color_data["ranges"]:
                mrca = resolve_mrca(species_tree, taxa_list)
                if mrca is not None:
                    draw_range_rect(ax, species_tree, mrca, clr, node_x, node_y)
            apply_label_colors(ax, color_data["labels"])
            color_legend = build_color_legend_handles(color_data)
            if color_legend:
                ax.legend(handles=color_legend, loc="upper right", fontsize=8, frameon=True)

        # Colorbar
        legend_loc = config.legend_position or "upper right"
        if legend_loc != "none":
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.15)
            cbar_fontsize = config.axis_fontsize if config.axis_fontsize else 10
            cbar.set_label("Hybridization score", fontsize=cbar_fontsize)

        ax.set_xlabel("Branch length (subs/site)")
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        if config.show_title:
            ax.set_title(config.title or "Hybridization Analysis", fontsize=config.title_fontsize)
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)

    # ------------------------------------------------------------------
    # Gene tree parsing
    # ------------------------------------------------------------------

    def _parse_gene_trees(self, path: str) -> list:
        from Bio import Phylo

        source = Path(path)
        try:
            with source.open() as handle:
                cleaned = [
                    stripped
                    for line in handle
                    if (stripped := line.strip())
                    and stripped[0] != "#"
                ]
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        trees = []
        parent_str = str(source.parent)
        parent_prefix = "" if parent_str == "." else parent_str + os.sep
        for line in cleaned:
            if line.startswith("("):
                trees.append(Phylo.read(StringIO(line), "newick"))
            else:
                tree_path = line if _path_isabs(line) else parent_prefix + line
                trees.append(Phylo.read(tree_path, "newick"))
        return trees

    # ------------------------------------------------------------------
    # Bipartition extraction and topology counting
    # ------------------------------------------------------------------

    @staticmethod
    def _canonical_split(taxa_side, all_taxa):
        taxa_side_len = len(taxa_side)
        all_taxa_len = len(all_taxa)
        if taxa_side_len * 2 < all_taxa_len:
            return frozenset(taxa_side)
        complement = all_taxa - taxa_side
        if taxa_side_len * 2 > all_taxa_len:
            return frozenset(complement)
        if not taxa_side:
            return frozenset(taxa_side)
        if min(taxa_side) <= min(complement):
            return frozenset(taxa_side)
        return frozenset(complement)

    @staticmethod
    def _build_parent_map(tree) -> dict:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            root = None

        parent_map = {}
        if root is not None:
            stack = [root]
            pop = stack.pop
            extend = stack.extend
            try:
                while stack:
                    clade = pop()
                    children = clade.clades
                    for child in children:
                        parent_map[id(child)] = clade
                    if children:
                        extend(children)
            except AttributeError:
                parent_map = {}
            else:
                return parent_map

        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    @staticmethod
    def _preorder_clades(tree) -> list:
        direct_clades = Hybridization._preorder_clades_direct(tree)
        if direct_clades is not None:
            return direct_clades
        return list(tree.find_clades(order="preorder"))

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
            while stack:
                clade = stack.pop()
                clades.append(clade)
                children = clade.clades
                if children:
                    stack.extend(reversed(children))
        except AttributeError:
            return None
        return clades

    @staticmethod
    def _postorder_clades(tree) -> list:
        direct_clades = Hybridization._postorder_clades_direct(tree)
        if direct_clades is not None:
            return direct_clades
        return list(tree.find_clades(order="postorder"))

    @staticmethod
    def _postorder_clades_direct(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        clades = []
        stack = [root]
        try:
            while stack:
                clade = stack.pop()
                clades.append(clade)
                children = clade.clades
                if children:
                    stack.extend(children)
        except AttributeError:
            return None
        clades.reverse()
        return clades

    @staticmethod
    def _collect_clade_taxa(tree) -> dict[int, frozenset]:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            root = None

        clade_taxa: dict[int, frozenset] = {}
        if root is not None:
            preorder = []
            stack = [root]
            try:
                while stack:
                    clade = stack.pop()
                    preorder.append(clade)
                    children = clade.clades
                    if children:
                        stack.extend(children)

                for clade in reversed(preorder):
                    children = clade.clades
                    if children:
                        child_count = len(children)
                        if child_count == 2:
                            taxa = (
                                clade_taxa[id(children[0])]
                                | clade_taxa[id(children[1])]
                            )
                        elif child_count == 1:
                            taxa = clade_taxa[id(children[0])]
                        else:
                            taxa = frozenset().union(
                                *(clade_taxa[id(child)] for child in children)
                            )
                    else:
                        taxa = frozenset({clade.name})
                    clade_taxa[id(clade)] = taxa
            except (AttributeError, KeyError, TypeError):
                clade_taxa = {}
            else:
                return clade_taxa

        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                clade_taxa[id(clade)] = frozenset({clade.name})
            else:
                taxa = frozenset()
                for child in clade.clades:
                    taxa = taxa | clade_taxa.get(id(child), frozenset())
                clade_taxa[id(clade)] = taxa
        return clade_taxa

    @staticmethod
    def _get_terminal_clades(tree) -> list:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return list(tree.get_terminals())

        terminals = []
        stack = [root]
        pop = stack.pop
        append = stack.append
        append_terminal = terminals.append
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    if len(children) == 2:
                        left, right = children
                        append(right)
                        append(left)
                    else:
                        for index in range(len(children) - 1, -1, -1):
                            append(children[index])
                else:
                    append_terminal(clade)
        except AttributeError:
            return list(tree.get_terminals())
        return terminals

    @staticmethod
    def _collect_taxa(trees) -> set:
        taxa = set()
        for tree in trees:
            stack = [tree.root]
            while stack:
                node = stack.pop()
                if not node.clades:
                    taxa.add(node.name)
                else:
                    stack.extend(node.clades)
        return taxa

    def _extract_splits(self, tree, all_taxa_fs, support_threshold=None) -> set:
        """Extract canonical bipartitions from a single tree.

        If support_threshold is provided, skip any internal node whose
        clade.confidence is not None and is below the threshold.
        """
        splits = set()
        tip_sets = {}
        stack = [(tree.root, False)]
        while stack:
            node, children_done = stack[-1]
            if not node.clades:
                stack.pop()
                name = node.name
                tip_sets[id(node)] = frozenset((name,)) if name in all_taxa_fs else frozenset()
            elif children_done:
                stack.pop()
                merged = frozenset().union(*(tip_sets[id(c)] for c in node.clades))
                tip_sets[id(node)] = merged
                if len(merged) > 1 and merged != all_taxa_fs:
                    # Check support threshold
                    if support_threshold is not None:
                        conf = getattr(node, 'confidence', None)
                        if conf is not None and conf < support_threshold:
                            # Skip this bipartition (collapse it)
                            continue
                    splits.add(self._canonical_split(merged, all_taxa_fs))
            else:
                stack[-1] = (node, True)
                for child in reversed(node.clades):
                    stack.append((child, False))
        return splits

    def _get_four_groups(self, tree, node, parent_map, all_taxa_fs, clade_taxa=None):
        if node.is_terminal() or len(node.clades) < 2:
            return None

        if clade_taxa is None:
            clade_taxa = self._collect_clade_taxa(tree)

        empty = frozenset()
        children = node.clades
        C1 = clade_taxa.get(id(children[0]), empty) & all_taxa_fs
        C2 = clade_taxa.get(id(children[1]), empty) & all_taxa_fs
        if len(children) > 2:
            C2 = C2.union(
                *(
                    clade_taxa.get(id(extra_child), empty) & all_taxa_fs
                    for extra_child in children[2:]
                )
            )

        parent = parent_map.get(id(node))
        if parent is None:
            return None

        siblings = [c for c in parent.clades if id(c) != id(node)]
        if not siblings:
            return None

        chosen_sib = None
        S = empty
        for sib in siblings:
            candidate = clade_taxa.get(id(sib), empty) & all_taxa_fs
            if candidate:
                S = candidate
                chosen_sib = sib
                break

        if not S:
            return None

        D = all_taxa_fs - C1 - C2 - S

        if not D:
            if chosen_sib.is_terminal() or len(chosen_sib.clades) < 2:
                return None
            sibling_children = chosen_sib.clades
            S = clade_taxa.get(id(sibling_children[0]), empty) & all_taxa_fs
            D = clade_taxa.get(id(sibling_children[1]), empty) & all_taxa_fs
            if len(sibling_children) > 2:
                D = D.union(
                    *(
                        clade_taxa.get(id(extra), empty) & all_taxa_fs
                        for extra in sibling_children[2:]
                    )
                )

        if not C1 or not C2 or not S:
            return None

        return C1, C2, S, D

    def _count_topologies(self, species_tree, gene_trees,
                          support_threshold=None) -> dict:
        """Count concordant and two NNI-alternative topologies for each
        internal branch of the species tree across gene trees.

        If support_threshold is provided, gene tree bipartitions with
        support below the threshold are collapsed before counting.
        """
        species_taxa = set(self.get_tip_names_from_tree(species_tree))
        gene_taxa = self._collect_taxa(gene_trees)
        all_taxa_fs = frozenset(sorted(species_taxa & gene_taxa))
        parent_map = self._build_parent_map(species_tree)
        species_clade_taxa = self._collect_clade_taxa(species_tree)

        gene_tree_splits = []
        for gt in gene_trees:
            splits = self._extract_splits(gt, all_taxa_fs, support_threshold)
            gene_tree_splits.append(splits)

        result = {}
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            groups = self._get_four_groups(
                species_tree, clade, parent_map, all_taxa_fs, species_clade_taxa
            )
            if groups is None:
                continue
            C1, C2, S, D = groups

            concordant_bp = self._canonical_split(C1 | C2, all_taxa_fs)
            nni_alt1_bp = self._canonical_split(S | C2, all_taxa_fs)
            nni_alt2_bp = self._canonical_split(C1 | S, all_taxa_fs)

            n_concordant, n_alt1, n_alt2 = self._count_split_matches(
                gene_tree_splits,
                concordant_bp,
                nni_alt1_bp,
                nni_alt2_bp,
            )

            node_tips = species_clade_taxa.get(id(clade), frozenset()) & all_taxa_fs
            split_label = (
                sorted(node_tips)
                if len(node_tips) <= len(all_taxa_fs) - len(node_tips)
                else sorted(all_taxa_fs - node_tips)
            )
            branch_key = ",".join(split_label)
            result[branch_key] = dict(
                split=split_label,
                n_concordant=n_concordant,
                n_alt1=n_alt1,
                n_alt2=n_alt2,
            )
        return result, all_taxa_fs

    @staticmethod
    def _count_split_matches(
        gene_tree_splits,
        concordant_bp,
        nni_alt1_bp,
        nni_alt2_bp,
    ) -> tuple[int, int, int]:
        n_concordant = n_alt1 = n_alt2 = 0
        for splits in gene_tree_splits:
            if concordant_bp in splits:
                n_concordant += 1
            if nni_alt1_bp in splits:
                n_alt1 += 1
            if nni_alt2_bp in splits:
                n_alt2 += 1
        return n_concordant, n_alt1, n_alt2

    # ------------------------------------------------------------------
    # Statistical testing
    # ------------------------------------------------------------------

    def _test_asymmetry(self, n_alt1: int, n_alt2: int) -> dict:
        total = n_alt1 + n_alt2
        if total == 0:
            return dict(
                asymmetry_ratio=None,
                p_value=None,
                favored_alt=None,
            )

        asymmetry_ratio = max(n_alt1, n_alt2) / total
        p_value = _binomial_two_sided_p_value(n_alt1, total)

        if n_alt1 > n_alt2:
            favored_alt = "alt1"
        elif n_alt2 > n_alt1:
            favored_alt = "alt2"
        else:
            favored_alt = None

        return dict(
            asymmetry_ratio=asymmetry_ratio,
            p_value=p_value,
            favored_alt=favored_alt,
        )

    @staticmethod
    def _fdr(p_values: list[float]) -> list[float]:
        """Benjamini-Hochberg FDR correction."""
        n = len(p_values)
        if n == 0:
            return []
        if n < _FDR_VECTOR_MIN_LENGTH:
            indexed = sorted(enumerate(p_values), key=lambda item: item[1])
            corrected = [0.0] * n
            previous = 1.0
            for rank_index in range(n - 1, -1, -1):
                original_index, p_value = indexed[rank_index]
                rank = rank_index + 1
                adjusted = min(p_value * n / rank, previous)
                adjusted = min(adjusted, 1.0)
                corrected[original_index] = adjusted
                previous = adjusted
            return corrected
        p_arr = np.asarray(p_values, dtype=float)
        order = np.argsort(p_arr)
        sorted_p = p_arr[order]
        ranks = np.arange(1, n + 1, dtype=float)
        adjusted = np.minimum.accumulate((sorted_p * n / ranks)[::-1])[::-1]
        adjusted = np.minimum(adjusted, 1.0)
        corrected = np.empty(n, dtype=float)
        corrected[order] = adjusted
        return corrected.tolist()
