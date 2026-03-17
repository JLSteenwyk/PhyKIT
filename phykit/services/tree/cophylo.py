import sys
from typing import Dict, List, Optional, Tuple

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


class Cophylo(Tree):
    """Cophylogenetic (tanglegram) visualization of two phylogenies.

    Draws two trees facing each other with connecting lines between
    matching taxa, analogous to R's phytools::cophylo().
    """

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree1_file_path"])
        self.tree2_file_path = parsed["tree2_file_path"]
        self.mapping_file_path = parsed["mapping_file_path"]
        self.output_path = parsed["output_path"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        tree1 = self.read_tree_file()
        self._validate_tree(tree1, "tree1")

        # Read tree2 manually via Phylo
        from Bio import Phylo
        try:
            tree2 = Phylo.read(self.tree2_file_path, "newick")
        except Exception:
            raise PhykitUserError(
                [
                    f"{self.tree2_file_path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )
        self._validate_tree(tree2, "tree2")

        tips1 = [t.name for t in tree1.get_terminals()]
        tips2 = [t.name for t in tree2.get_terminals()]

        # Build mapping
        if self.mapping_file_path:
            mapping = self._parse_mapping_file(self.mapping_file_path)
        else:
            # Default: match by identical names
            shared = set(tips1) & set(tips2)
            if len(shared) < 2:
                raise PhykitUserError(
                    [
                        f"Only {len(shared)} shared taxa between trees.",
                        "At least 2 shared taxa are required.",
                        "Use -m/--mapping to specify a mapping file if tip names differ.",
                    ],
                    code=2,
                )
            mapping = {name: name for name in shared}

        # Validate mapping
        mapping_valid = {}
        tips1_set = set(tips1)
        tips2_set = set(tips2)
        for t1, t2 in mapping.items():
            if t1 in tips1_set and t2 in tips2_set:
                mapping_valid[t1] = t2

        if len(mapping_valid) < 2:
            raise PhykitUserError(
                [
                    f"Only {len(mapping_valid)} valid mapped taxa found.",
                    "At least 2 matched taxa are required.",
                ],
                code=2,
            )

        n_matched = len(mapping_valid)
        n_tips1 = len(tips1)
        n_tips2 = len(tips2)

        if self.plot_config.ladderize:
            tree1.ladderize()
            tree2.ladderize()

        # Optimize tip ordering to reduce line crossings
        tree1_order, tree2_order = self._optimize_tip_order(
            tree1, tree2, mapping_valid
        )

        # Plot
        self._plot_cophylo(
            tree1, tree2, mapping_valid,
            tree1_order, tree2_order, self.output_path,
        )

        if self.json_output:
            result = {
                "n_tips_tree1": n_tips1,
                "n_tips_tree2": n_tips2,
                "n_matched": n_matched,
                "matched_taxa": {k: v for k, v in sorted(mapping_valid.items())},
                "plot_output": self.output_path,
            }
            print_json(result)
        else:
            print("Cophylogenetic Plot (Tanglegram)")
            print(f"\nTree 1 tips: {n_tips1}")
            print(f"Tree 2 tips: {n_tips2}")
            print(f"Matched taxa: {n_matched}")
            print(f"Saved cophylo plot: {self.output_path}")

    def process_args(self, args) -> Dict:
        return dict(
            tree1_file_path=args.tree1,
            tree2_file_path=args.tree2,
            mapping_file_path=getattr(args, "mapping", None),
            output_path=args.output,
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def _validate_tree(self, tree, label: str) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 2:
            raise PhykitUserError(
                [f"{label} must have at least 2 tips."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                clade.branch_length = 1.0

    def _parse_mapping_file(self, path: str) -> Dict[str, str]:
        try:
            with open(path) as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        mapping = {}
        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise PhykitUserError(
                    [
                        f"Line {line_num} in mapping file has {len(parts)} columns; expected 2.",
                        "Each line should be: taxon_in_tree1<tab>taxon_in_tree2",
                    ],
                    code=2,
                )
            mapping[parts[0]] = parts[1]

        return mapping

    def _build_parent_map(self, tree) -> Dict:
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    def _get_tip_order(self, tree) -> List[str]:
        """Get tip names in tree traversal order (postorder gives
        bottom-to-top in a standard phylogram layout)."""
        return [tip.name for tip in tree.get_terminals()]

    def _optimize_tip_order(
        self, tree1, tree2, mapping: Dict[str, str],
    ) -> Tuple[Dict[str, int], Dict[str, int]]:
        """Compute tip y-positions for both trees.

        Uses the natural tree traversal order for tree1, then
        attempts to reorder tree2 tips to minimize line crossings
        by rotating internal nodes.
        """
        # Tree 1: use natural traversal order
        tips1 = self._get_tip_order(tree1)
        tree1_order = {name: i for i, name in enumerate(tips1)}

        # Tree 2: try to match tree1's ordering for connected taxa
        # Build target y-positions for tree2 tips based on their
        # connected partner's position in tree1
        reverse_mapping = {v: k for k, v in mapping.items()}

        # Rotate tree2 nodes to minimize crossings
        self._rotate_tree(tree2, tree1_order, reverse_mapping)

        tips2 = self._get_tip_order(tree2)
        tree2_order = {name: i for i, name in enumerate(tips2)}

        return tree1_order, tree2_order

    def _rotate_tree(
        self, tree, target_order: Dict[str, int],
        reverse_mapping: Dict[str, str],
    ) -> None:
        """Rotate internal nodes of tree to minimize crossing lines.

        For each internal node, check if swapping children improves
        the alignment with target_order (via the mapping).
        """
        def _get_mean_target_pos(clade):
            tips = [t.name for t in clade.get_terminals()]
            positions = []
            for t in tips:
                partner = reverse_mapping.get(t)
                if partner and partner in target_order:
                    positions.append(target_order[partner])
            if positions:
                return np.mean(positions)
            return 0.0

        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                continue
            if len(clade.clades) == 2:
                pos0 = _get_mean_target_pos(clade.clades[0])
                pos1 = _get_mean_target_pos(clade.clades[1])
                if pos0 > pos1:
                    clade.clades[0], clade.clades[1] = (
                        clade.clades[1], clade.clades[0]
                    )

    def _plot_cophylo(
        self, tree1, tree2, mapping: Dict[str, str],
        tree1_order: Dict[str, int], tree2_order: Dict[str, int],
        output_path: str,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.lines import Line2D
        except ImportError:
            print(
                "matplotlib is required for cophylo plotting. "
                "Install matplotlib and retry."
            )
            raise SystemExit(2)

        tips1 = list(tree1.get_terminals())
        tips2 = list(tree2.get_terminals())
        n_max = max(len(tips1), len(tips2))
        config = self.plot_config
        config.resolve(n_rows=n_max, n_cols=None)

        if config.circular:
            self._plot_cophylo_circular(
                tree1, tree2, mapping, output_path, config, plt,
            )
        else:
            self._plot_cophylo_rect(
                tree1, tree2, mapping, tree1_order, tree2_order,
                output_path, config, n_max, plt,
            )

    def _plot_cophylo_rect(
        self, tree1, tree2, mapping, tree1_order, tree2_order,
        output_path, config, n_max, plt,
    ) -> None:
        fig, (ax1, ax_mid, ax2) = plt.subplots(
            1, 3,
            figsize=(config.fig_width, config.fig_height),
            gridspec_kw={"width_ratios": [4, 2, 4]},
        )

        # Draw tree1 (left, normal orientation: root on left, tips on right)
        self._draw_phylogram(
            ax1, tree1, tree1_order, direction="right",
            color="#2171b5",
        )

        # Draw tree2 (right, mirrored: root on right, tips on left)
        self._draw_phylogram(
            ax2, tree2, tree2_order, direction="left",
            color="#cb181d",
        )

        # Draw connecting lines in the middle panel
        ax_mid.set_xlim(0, 1)
        ax_mid.set_ylim(-0.5, n_max - 0.5)
        ax_mid.axis("off")

        for t1_name, t2_name in mapping.items():
            if t1_name in tree1_order and t2_name in tree2_order:
                y1 = tree1_order[t1_name]
                y2 = tree2_order[t2_name]
                ax_mid.plot(
                    [0, 1], [y1, y2],
                    color="gray", alpha=0.5, linewidth=0.8,
                    zorder=1,
                )

        if config.show_title:
            ax1.set_title("Tree 1", fontsize=config.title_fontsize or 11, fontweight="bold")
            ax2.set_title("Tree 2", fontsize=config.title_fontsize or 11, fontweight="bold")

        fig.suptitle("Cophylogenetic Plot (Tanglegram)", fontsize=13)
        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved cophylo plot: {output_path}")

    def _plot_cophylo_circular(
        self, tree1, tree2, mapping, output_path, config, plt,
    ) -> None:
        from matplotlib.patches import ConnectionPatch
        from ...helpers.circular_layout import (
            compute_circular_coords,
            draw_circular_branches,
            draw_circular_tip_labels,
        )

        fig, (ax1, ax2) = plt.subplots(
            1, 2, figsize=(config.fig_width, config.fig_height),
        )

        def _build_circular(tree, ax, color):
            parent_map = self._build_parent_map(tree)
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
            draw_circular_branches(ax, tree, coords, parent_map, color=color)
            draw_circular_tip_labels(ax, tree, coords, fontsize=7)
            ax.set_aspect("equal")
            ax.axis("off")
            return coords, parent_map

        coords1, pmap1 = _build_circular(tree1, ax1, color="#2171b5")
        coords2, pmap2 = _build_circular(tree2, ax2, color="#cb181d")

        # Apply color annotations to both circular trees
        if self.plot_config.color_file:
            from ...helpers.circular_layout import draw_circular_colored_branch
            color_data = parse_color_file(self.plot_config.color_file)
            for tree_obj, coords_obj, pmap_obj, ax_obj in [
                (tree1, coords1, pmap1, ax1), (tree2, coords2, pmap2, ax2)
            ]:
                for taxa_list, clr, lbl in color_data["ranges"]:
                    mrca = resolve_mrca(tree_obj, taxa_list)
                    if mrca is not None:
                        draw_range_wedge(ax_obj, tree_obj, mrca, clr, coords_obj)
                for taxa_list, clade_color, lbl in color_data["clades"]:
                    mrca = resolve_mrca(tree_obj, taxa_list)
                    if mrca is not None:
                        clade_ids = get_clade_branch_ids(tree_obj, mrca, pmap_obj)
                        for cl in tree_obj.find_clades(order="preorder"):
                            if cl == tree_obj.root:
                                continue
                            if id(cl) in clade_ids and id(cl) in pmap_obj:
                                draw_circular_colored_branch(ax_obj, coords_obj[id(pmap_obj[id(cl)])], coords_obj[id(cl)], clade_color, lw=1.5)
                for taxon, lbl_color in color_data["labels"].items():
                    for text_obj in ax_obj.texts:
                        if text_obj.get_text() == taxon:
                            text_obj.set_color(lbl_color)
                            break
            color_legend = build_color_legend_handles(color_data)
            if color_legend:
                ax1.legend(handles=color_legend, loc="upper right", fontsize=8, frameon=True)

        # Build tip name -> id lookups
        tip_name_to_id1 = {t.name: id(t) for t in tree1.get_terminals()}
        tip_name_to_id2 = {t.name: id(t) for t in tree2.get_terminals()}

        # Draw connecting lines between matched taxa
        for t1_name, t2_name in mapping.items():
            tid1 = tip_name_to_id1.get(t1_name)
            tid2 = tip_name_to_id2.get(t2_name)
            if tid1 is None or tid2 is None:
                continue
            if tid1 not in coords1 or tid2 not in coords2:
                continue
            con = ConnectionPatch(
                xyA=(coords1[tid1]["x"], coords1[tid1]["y"]),
                xyB=(coords2[tid2]["x"], coords2[tid2]["y"]),
                coordsA="data", coordsB="data",
                axesA=ax1, axesB=ax2,
                color="gray", alpha=0.5, lw=0.5,
            )
            fig.add_artist(con)

        if config.show_title:
            ax1.set_title("Tree 1", fontsize=config.title_fontsize or 11, fontweight="bold")
            ax2.set_title("Tree 2", fontsize=config.title_fontsize or 11, fontweight="bold")

        fig.suptitle("Cophylogenetic Plot (Tanglegram)", fontsize=13)
        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved cophylo plot: {output_path}")

    def _draw_phylogram(
        self, ax, tree, tip_order: Dict[str, int],
        direction: str = "right", color: str = "#333333",
    ) -> None:
        """Draw a phylogram on the given axes.

        direction='right': root on left, tips on right (normal)
        direction='left': root on right, tips on left (mirrored)
        """
        parent_map = self._build_parent_map(tree)
        tips = list(tree.get_terminals())

        node_x = {}
        node_y = {}

        # Assign y-positions from tip_order
        for tip in tips:
            if tip.name in tip_order:
                node_y[id(tip)] = tip_order[tip.name]

        # Assign x-positions (cumulative branch length from root)
        root = tree.root
        if self.plot_config.cladogram:
            node_x = compute_node_x_cladogram(tree, parent_map)
        else:
            node_x[id(root)] = 0.0
            for clade in tree.find_clades(order="preorder"):
                if clade == root:
                    continue
                if id(clade) in parent_map:
                    parent = parent_map[id(clade)]
                    bl = clade.branch_length if clade.branch_length else 0.0
                    node_x[id(clade)] = node_x[id(parent)] + bl

        # Internal node y-positions (mean of children)
        for clade in tree.find_clades(order="postorder"):
            if not clade.is_terminal() and id(clade) not in node_y:
                child_ys = [
                    node_y[id(c)] for c in clade.clades if id(c) in node_y
                ]
                if child_ys:
                    node_y[id(clade)] = np.mean(child_ys)
                else:
                    node_y[id(clade)] = 0.0

        # For mirrored trees, flip x so root is on the right
        max_x = max(node_x.values()) if node_x else 1.0
        if direction == "left":
            for k in node_x:
                node_x[k] = max_x - node_x[k]

        # Draw branches
        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            if id(clade) not in parent_map:
                continue
            parent = parent_map[id(clade)]
            if id(parent) not in node_x or id(clade) not in node_x:
                continue

            px = node_x[id(parent)]
            py = node_y.get(id(parent), 0)
            cx = node_x[id(clade)]
            cy = node_y.get(id(clade), 0)

            # Vertical connector at parent x
            ax.plot(
                [px, px], [py, cy],
                color=color, linewidth=1.0, solid_capstyle="butt",
            )
            # Horizontal branch
            ax.plot(
                [px, cx], [cy, cy],
                color=color, linewidth=1.5, solid_capstyle="butt",
            )

        # Tip labels
        label_offset = max_x * 0.03
        for tip in tips:
            if id(tip) not in node_x or tip.name not in tip_order:
                continue
            tx = node_x[id(tip)]
            ty = node_y[id(tip)]

            if direction == "right":
                ax.text(
                    tx + label_offset, ty, tip.name,
                    va="center", ha="left", fontsize=8,
                )
            else:
                ax.text(
                    tx - label_offset, ty, tip.name,
                    va="center", ha="right", fontsize=8,
                )

        # Apply color annotations
        if self.plot_config.color_file:
            color_data = parse_color_file(self.plot_config.color_file)
            for taxa_list, clr, lbl in color_data["ranges"]:
                mrca = resolve_mrca(tree, taxa_list)
                if mrca is not None:
                    draw_range_rect(ax, tree, mrca, clr, node_x, node_y)
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
                            ax.plot([x0, x1], [y1, y1], color=clade_color, lw=1.5, zorder=2)
                            ax.plot([x0, x0], [y0, y1], color=clade_color, lw=1.5, zorder=2)
            for taxon, lbl_color in color_data["labels"].items():
                for text_obj in ax.texts:
                    if text_obj.get_text() == taxon:
                        text_obj.set_color(lbl_color)
                        break
            color_legend = build_color_legend_handles(color_data)
            if color_legend:
                ax.legend(handles=color_legend, loc="upper right", fontsize=8, frameon=True)

        n_max = max(tip_order.values()) + 1 if tip_order else 1
        ax.set_ylim(-0.5, n_max - 0.5)
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.set_xlabel("Branch length")
