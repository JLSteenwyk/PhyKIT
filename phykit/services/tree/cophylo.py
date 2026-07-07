from __future__ import annotations

import sys

from .base import Tree
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def compute_circular_coords(*args, **kwargs):
    from ...helpers.circular_layout import compute_circular_coords as _compute_circular_coords

    return _compute_circular_coords(*args, **kwargs)


def draw_circular_branches(*args, **kwargs):
    from ...helpers.circular_layout import draw_circular_branches as _draw_circular_branches

    return _draw_circular_branches(*args, **kwargs)


def draw_circular_tip_labels(*args, **kwargs):
    from ...helpers.circular_layout import draw_circular_tip_labels as _draw_circular_tip_labels

    return _draw_circular_tip_labels(*args, **kwargs)


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

        try:
            tree2 = self._read_tree_with_error(
                self.tree2_file_path,
                "tree2_file_path",
            )
        except Exception:
            raise PhykitUserError(
                [
                    f"{self.tree2_file_path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )
        self._validate_tree(tree2, "tree2")

        tips1 = self.get_tip_names_from_tree(tree1)
        tips2 = self.get_tip_names_from_tree(tree2)
        tips1_set = set(tips1)
        tips2_set = set(tips2)

        # Build mapping
        if self.mapping_file_path:
            mapping = self._parse_mapping_file(self.mapping_file_path)
            mapping_valid = {}
            for t1, t2 in mapping.items():
                if t1 in tips1_set and t2 in tips2_set:
                    mapping_valid[t1] = t2
        else:
            # Default: match by identical names
            shared = tips1_set & tips2_set
            if len(shared) < 2:
                raise PhykitUserError(
                    [
                        f"Only {len(shared)} shared taxa between trees.",
                        "At least 2 shared taxa are required.",
                        "Use -m/--mapping to specify a mapping file if tip names differ.",
                    ],
                    code=2,
                )
            mapping_valid = {name: name for name in shared}

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
            self._print_text_output(n_tips1, n_tips2, n_matched)

    def _print_text_output(
        self,
        n_tips1: int,
        n_tips2: int,
        n_matched: int,
    ) -> None:
        print(
            "\n".join(
                [
                    "Cophylogenetic Plot (Tanglegram)",
                    f"\nTree 1 tips: {n_tips1}",
                    f"Tree 2 tips: {n_tips2}",
                    f"Matched taxa: {n_matched}",
                    f"Saved cophylo plot: {self.output_path}",
                ]
            )
        )

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree1_file_path=args.tree1,
            tree2_file_path=args.tree2,
            mapping_file_path=getattr(args, "mapping", None),
            output_path=args.output,
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def _validate_tree(self, tree, label: str) -> None:
        tip_count = self._validate_standard_tree(tree)
        if tip_count is None:
            tips = list(tree.get_terminals())
            tip_count = len(tips)
            for clade in tree.find_clades():
                if clade.branch_length is None and clade != tree.root:
                    clade.branch_length = 1.0

        if tip_count < 2:
            raise PhykitUserError(
                [f"{label} must have at least 2 tips."],
                code=2,
            )

    @staticmethod
    def _validate_standard_tree(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        tip_count = 0
        stack = [root]
        pop = stack.pop
        extend = stack.extend
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if clade is not root and clade.branch_length is None:
                    clade.branch_length = 1.0
                if children:
                    extend(children)
                else:
                    tip_count += 1
        except AttributeError:
            return None

        return tip_count

    def _parse_mapping_file(self, path: str) -> dict[str, str]:
        try:
            mapping = {}
            with open(path) as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line[0] == "#":
                        continue
                    taxon1, sep, taxon2 = line.partition("\t")
                    if not sep or "\t" in taxon2:
                        column_count = line.count("\t") + 1
                        raise PhykitUserError(
                            [
                                f"Line {line_num} in mapping file has {column_count} columns; expected 2.",
                                "Each line should be: taxon_in_tree1<tab>taxon_in_tree2",
                            ],
                            code=2,
                        )
                    mapping[taxon1] = taxon2
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        return mapping

    def _build_parent_map(self, tree) -> dict:
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
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return list(tree.find_clades(order="preorder"))

        clades = []
        stack = [root]
        try:
            pop = stack.pop
            append_stack = stack.append
            append_clade = clades.append
            while stack:
                clade = pop()
                append_clade(clade)
                children = clade.clades
                if children:
                    child_count = len(children)
                    if child_count == 2:
                        append_stack(children[1])
                        append_stack(children[0])
                    else:
                        for index in range(child_count - 1, -1, -1):
                            append_stack(children[index])
        except AttributeError:
            return list(tree.find_clades(order="preorder"))
        return clades

    @staticmethod
    def _postorder_clades(tree) -> list:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return list(tree.find_clades(order="postorder"))

        clades = []
        stack = [root]
        pop = stack.pop
        extend = stack.extend
        append_clade = clades.append
        try:
            while stack:
                clade = pop()
                append_clade(clade)
                children = clade.clades
                if children:
                    extend(children)
        except AttributeError:
            return list(tree.find_clades(order="postorder"))
        clades.reverse()
        return clades

    @staticmethod
    def _terminal_clades(tree) -> list:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return list(tree.get_terminals())

        terminals = []
        stack = [root]
        try:
            pop = stack.pop
            append_stack = stack.append
            append_terminal = terminals.append
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    child_count = len(children)
                    if child_count == 2:
                        append_stack(children[1])
                        append_stack(children[0])
                    else:
                        for index in range(child_count - 1, -1, -1):
                            append_stack(children[index])
                else:
                    append_terminal(clade)
        except AttributeError:
            return list(tree.get_terminals())
        return terminals

    def _get_tip_order(self, tree) -> list[str]:
        """Get tip names in tree traversal order (postorder gives
        bottom-to-top in a standard phylogram layout)."""
        return self.get_tip_names_from_tree(tree)

    def _optimize_tip_order(
        self, tree1, tree2, mapping: dict[str, str],
    ) -> tuple[dict[str, int], dict[str, int]]:
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
        self, tree, target_order: dict[str, int],
        reverse_mapping: dict[str, str],
    ) -> None:
        """Rotate internal nodes of tree to minimize crossing lines.

        For each internal node, check if swapping children improves
        the alignment with target_order (via the mapping).
        """
        position_sums = {}
        position_counts = {}

        for clade in self._postorder_clades(tree):
            children = clade.clades
            if not children:
                partner = reverse_mapping.get(clade.name)
                if partner and partner in target_order:
                    position_sum = target_order[partner]
                    position_count = 1
                else:
                    position_sum = 0.0
                    position_count = 0
            else:
                position_sum = 0.0
                position_count = 0
                for child in children:
                    child_id = id(child)
                    position_sum += position_sums[child_id]
                    position_count += position_counts[child_id]

                if len(children) == 2:
                    child0, child1 = children
                    child0_id = id(child0)
                    child1_id = id(child1)
                    child0_count = position_counts[child0_id]
                    child1_count = position_counts[child1_id]
                    pos0 = (
                        position_sums[child0_id] / child0_count
                        if child0_count
                        else 0.0
                    )
                    pos1 = (
                        position_sums[child1_id] / child1_count
                        if child1_count
                        else 0.0
                    )
                    if pos0 > pos1:
                        children[0], children[1] = child1, child0

            position_sums[id(clade)] = position_sum
            position_counts[id(clade)] = position_count

    def _plot_cophylo(
        self, tree1, tree2, mapping: dict[str, str],
        tree1_order: dict[str, int], tree2_order: dict[str, int],
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

        tips1 = self._terminal_clades(tree1)
        tips2 = self._terminal_clades(tree2)
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
        from matplotlib.collections import LineCollection

        fig, (ax1, ax_mid, ax2) = plt.subplots(
            1, 3,
            figsize=(config.fig_width, config.fig_height),
            gridspec_kw={"width_ratios": [4, 2, 4]},
        )
        color_data = None
        if self.plot_config.color_file:
            from ...helpers.color_annotations import parse_color_file

            color_data = parse_color_file(self.plot_config.color_file)

        # Draw tree1 (left, normal orientation: root on left, tips on right)
        self._draw_phylogram(
            ax1, tree1, tree1_order, direction="right",
            color="#2171b5", color_data=color_data,
        )

        # Draw tree2 (right, mirrored: root on right, tips on left)
        self._draw_phylogram(
            ax2, tree2, tree2_order, direction="left",
            color="#cb181d", color_data=color_data,
        )

        # Draw connecting lines in the middle panel
        ax_mid.set_xlim(0, 1)
        ax_mid.set_ylim(-0.5, n_max - 0.5)
        ax_mid.axis("off")

        association_segments = []
        for t1_name, t2_name in mapping.items():
            if t1_name in tree1_order and t2_name in tree2_order:
                y1 = tree1_order[t1_name]
                y2 = tree2_order[t2_name]
                association_segments.append(((0, y1), (1, y2)))
        if association_segments:
            ax_mid.add_collection(
                LineCollection(
                    association_segments,
                    colors="gray",
                    alpha=0.5,
                    linewidths=0.8,
                    zorder=1,
                ),
                autolim=True,
            )

        if config.show_title:
            ax1.set_title("Tree 1", fontsize=config.title_fontsize or 11, fontweight="bold")
            ax2.set_title("Tree 2", fontsize=config.title_fontsize or 11, fontweight="bold")

        fig.suptitle("Cophylogenetic Plot (Tanglegram)", fontsize=13)
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved cophylo plot: {output_path}")

    def _plot_cophylo_circular(
        self, tree1, tree2, mapping, output_path, config, plt,
    ) -> None:
        from matplotlib.collections import LineCollection

        fig, (ax1, ax2) = plt.subplots(
            1, 2, figsize=(config.fig_width, config.fig_height),
        )

        def _build_circular(tree, ax, color):
            parent_map = self._build_parent_map(tree)
            preorder_clades = self._preorder_clades(tree)
            terminal_clades = [
                clade for clade in preorder_clades if not clade.clades
            ]
            root = tree.root
            if config.cladogram:
                from ...helpers.plot_config import compute_node_x_cladogram

                node_x = compute_node_x_cladogram(tree, parent_map)
            else:
                node_x = {}
                for clade in preorder_clades:
                    if clade == root:
                        node_x[id(clade)] = 0.0
                    elif id(clade) in parent_map:
                        parent = parent_map[id(clade)]
                        bl = clade.branch_length if clade.branch_length else 0.0
                        node_x[id(clade)] = node_x.get(id(parent), 0.0) + bl

            coords = compute_circular_coords(
                tree,
                node_x,
                parent_map,
                preorder_clades=preorder_clades,
                terminal_clades=terminal_clades,
            )
            draw_circular_branches(ax, tree, coords, parent_map, color=color)
            draw_circular_tip_labels(ax, tree, coords, fontsize=7)
            ax.set_aspect("equal")
            ax.axis("off")
            return coords, parent_map, preorder_clades, terminal_clades

        coords1, pmap1, preorder1, terminals1 = _build_circular(
            tree1, ax1, color="#2171b5"
        )
        coords2, pmap2, preorder2, terminals2 = _build_circular(
            tree2, ax2, color="#cb181d"
        )

        # Apply color annotations to both circular trees
        if self.plot_config.color_file:
            import math
            from matplotlib.collections import LineCollection
            from ...helpers.color_annotations import (
                parse_color_file,
                resolve_mrca,
                draw_range_wedge,
                get_clade_branch_ids,
                build_color_legend_handles,
                apply_label_colors,
            )
            color_data = parse_color_file(self.plot_config.color_file)
            for tree_obj, coords_obj, pmap_obj, ax_obj, preorder_obj in [
                (tree1, coords1, pmap1, ax1, preorder1),
                (tree2, coords2, pmap2, ax2, preorder2),
            ]:
                for taxa_list, clr, lbl in color_data["ranges"]:
                    mrca = resolve_mrca(tree_obj, taxa_list)
                    if mrca is not None:
                        draw_range_wedge(ax_obj, tree_obj, mrca, clr, coords_obj)
                for taxa_list, clade_color, lbl in color_data["clades"]:
                    mrca = resolve_mrca(tree_obj, taxa_list)
                    if mrca is not None:
                        clade_ids = get_clade_branch_ids(tree_obj, mrca, pmap_obj)
                        radial_segments = []
                        for cl in preorder_obj:
                            if cl == tree_obj.root:
                                continue
                            if id(cl) in clade_ids and id(cl) in pmap_obj:
                                parent_coords = coords_obj[id(pmap_obj[id(cl)])]
                                child_coords = coords_obj[id(cl)]
                                angle = child_coords["angle"]
                                r_parent = parent_coords["radius"]
                                r_child = child_coords["radius"]
                                cos_angle = math.cos(angle)
                                sin_angle = math.sin(angle)
                                radial_segments.append(
                                    [
                                        (r_parent * cos_angle, r_parent * sin_angle),
                                        (r_child * cos_angle, r_child * sin_angle),
                                    ]
                                )
                        if radial_segments:
                            ax_obj.add_collection(
                                LineCollection(
                                    radial_segments,
                                    colors=clade_color,
                                    linewidths=1.5,
                                    capstyle="round",
                                    zorder=2,
                                )
                            )
                apply_label_colors(ax_obj, color_data["labels"])
            color_legend = build_color_legend_handles(color_data)
            if color_legend:
                ax1.legend(handles=color_legend, loc="upper right", fontsize=8, frameon=True)

        # Build tip name -> id lookups
        tip_name_to_id1 = {t.name: id(t) for t in terminals1}
        tip_name_to_id2 = {t.name: id(t) for t in terminals2}

        if config.show_title:
            ax1.set_title("Tree 1", fontsize=config.title_fontsize or 11, fontweight="bold")
            ax2.set_title("Tree 2", fontsize=config.title_fontsize or 11, fontweight="bold")

        fig.suptitle("Cophylogenetic Plot (Tanglegram)", fontsize=13)
        fig.tight_layout()

        # Draw connecting lines between matched taxa in one figure-level collection
        figure_segments = []
        to_figure = fig.transFigure.inverted().transform
        for t1_name, t2_name in mapping.items():
            tid1 = tip_name_to_id1.get(t1_name)
            tid2 = tip_name_to_id2.get(t2_name)
            if tid1 is None or tid2 is None:
                continue
            if tid1 not in coords1 or tid2 not in coords2:
                continue
            point1 = to_figure(
                ax1.transData.transform((coords1[tid1]["x"], coords1[tid1]["y"]))
            )
            point2 = to_figure(
                ax2.transData.transform((coords2[tid2]["x"], coords2[tid2]["y"]))
            )
            figure_segments.append((point1, point2))
        if figure_segments:
            fig.add_artist(
                LineCollection(
                    figure_segments,
                    colors="gray",
                    alpha=0.5,
                    linewidths=0.5,
                    transform=fig.transFigure,
                    clip_on=False,
                    zorder=1,
                )
            )
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved cophylo plot: {output_path}")

    @staticmethod
    def _assign_internal_y_positions(preorder_clades, node_y) -> None:
        for clade in reversed(preorder_clades):
            children = clade.clades
            if not children or id(clade) in node_y:
                continue

            child_count = len(children)
            if child_count == 2:
                left_y = node_y.get(id(children[0]))
                right_y = node_y.get(id(children[1]))
                if left_y is not None and right_y is not None:
                    node_y[id(clade)] = (left_y + right_y) * 0.5
                elif left_y is not None:
                    node_y[id(clade)] = left_y
                elif right_y is not None:
                    node_y[id(clade)] = right_y
                else:
                    node_y[id(clade)] = 0.0
                continue

            if child_count == 1:
                node_y[id(clade)] = node_y.get(id(children[0]), 0.0)
                continue

            total_y = 0.0
            seen = 0
            for child in children:
                child_y = node_y.get(id(child))
                if child_y is not None:
                    total_y += child_y
                    seen += 1
            node_y[id(clade)] = total_y / seen if seen else 0.0

    def _draw_phylogram(
        self, ax, tree, tip_order: dict[str, int],
        direction: str = "right", color: str = "#333333",
        color_data=None,
    ) -> None:
        """Draw a phylogram on the given axes.

        direction='right': root on left, tips on right (normal)
        direction='left': root on right, tips on left (mirrored)
        """
        parent_map = self._build_parent_map(tree)
        tips = self._terminal_clades(tree)
        preorder_clades = self._preorder_clades(tree)

        node_x = {}
        node_y = {}

        # Assign y-positions from tip_order
        for tip in tips:
            if tip.name in tip_order:
                node_y[id(tip)] = tip_order[tip.name]

        # Assign x-positions (cumulative branch length from root)
        root = tree.root
        if self.plot_config.cladogram:
            from ...helpers.plot_config import compute_node_x_cladogram

            node_x = compute_node_x_cladogram(tree, parent_map)
        else:
            node_x[id(root)] = 0.0
            for clade in preorder_clades:
                if clade == root:
                    continue
                if id(clade) in parent_map:
                    parent = parent_map[id(clade)]
                    bl = clade.branch_length if clade.branch_length else 0.0
                    node_x[id(clade)] = node_x[id(parent)] + bl

        # Internal node y-positions (mean of children)
        self._assign_internal_y_positions(preorder_clades, node_y)

        # For mirrored trees, flip x so root is on the right
        max_x = max(node_x.values()) if node_x else 1.0
        if direction == "left":
            for k in node_x:
                node_x[k] = max_x - node_x[k]

        # Draw branches
        if hasattr(ax, "add_collection"):
            from matplotlib.collections import LineCollection

            vertical_segments = []
            horizontal_segments = []
            for clade in preorder_clades:
                if clade == root:
                    continue
                cid = id(clade)
                if cid not in parent_map:
                    continue
                parent = parent_map[cid]
                pid = id(parent)
                if pid not in node_x or cid not in node_x:
                    continue

                px = node_x[pid]
                py = node_y.get(pid, 0)
                cx = node_x[cid]
                cy = node_y.get(cid, 0)

                vertical_segments.append(((px, py), (px, cy)))
                horizontal_segments.append(((px, cy), (cx, cy)))

            if vertical_segments:
                ax.add_collection(
                    LineCollection(
                        vertical_segments,
                        colors=color,
                        linewidths=1.0,
                        capstyle="butt",
                    ),
                    autolim=True,
                )
            if horizontal_segments:
                ax.add_collection(
                    LineCollection(
                        horizontal_segments,
                        colors=color,
                        linewidths=1.5,
                        capstyle="butt",
                    ),
                    autolim=True,
                )
            ax.autoscale_view()
        else:
            for clade in preorder_clades:
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
            from ...helpers.color_annotations import (
                parse_color_file,
                resolve_mrca,
                draw_range_rect,
                get_clade_branch_ids,
                build_color_legend_handles,
                apply_label_colors,
            )
            if color_data is None:
                color_data = parse_color_file(self.plot_config.color_file)
            for taxa_list, clr, lbl in color_data["ranges"]:
                mrca = resolve_mrca(tree, taxa_list)
                if mrca is not None:
                    draw_range_rect(ax, tree, mrca, clr, node_x, node_y)
            for taxa_list, clade_color, lbl in color_data["clades"]:
                mrca = resolve_mrca(tree, taxa_list)
                if mrca is not None:
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
                    if hasattr(ax, "add_collection"):
                        from matplotlib.collections import LineCollection

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
                    else:
                        for (x0, y0), (_, y1) in vertical_segments:
                            ax.plot([x0, x0], [y0, y1], color=clade_color, lw=1.5, zorder=2)
                        for (x0, y1), (x1, _) in horizontal_segments:
                            ax.plot([x0, x1], [y1, y1], color=clade_color, lw=1.5, zorder=2)
            apply_label_colors(ax, color_data["labels"])
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
