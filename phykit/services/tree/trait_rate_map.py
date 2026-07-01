from __future__ import annotations

import math
import sys

from .base import Tree
from ...errors import PhykitUserError


_BL_FLOOR = 1e-8


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class TraitRateMap(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.output_path = parsed["output_path"]
        self.trait_column = parsed["trait_column"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        self.validate_tree(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="trait rate map visualization",
        )

        tree_tips = self.get_tip_names_from_tree(tree)

        if self.trait_column:
            trait_values = self._parse_multi_trait_data(
                self.trait_data_path, tree_tips, self.trait_column
            )
            trait_name = self.trait_column
        else:
            trait_values = self._parse_single_trait_data(
                self.trait_data_path, tree_tips
            )
            trait_name = "trait"

        ordered_names = sorted(trait_values.keys())
        n = len(ordered_names)

        # Prune tree to shared taxa
        tips_to_prune = [t for t in tree_tips if t not in trait_values]
        tree_for_analysis = tree
        if tips_to_prune or self.plot_config.ladderize:
            tree_for_analysis = self._fast_copy(tree)
            if tips_to_prune:
                tree_for_analysis = self.prune_tree_using_taxa_list(
                    tree_for_analysis,
                    tips_to_prune,
                )

        if self.plot_config.ladderize:
            tree_for_analysis.ladderize()

        preorder_clades = list(self._iter_preorder(tree_for_analysis.root))
        postorder_clades = list(reversed(preorder_clades))

        # Label internal nodes
        node_labels = self._label_internal_nodes(tree_for_analysis, preorder_clades)

        # Reconstruct ancestral states (Felsenstein weighted average)
        node_values = self._ancestral_reconstruction(
            tree_for_analysis, trait_values, postorder_clades
        )

        # Compute per-branch rates
        parent_map = self._build_parent_map(tree_for_analysis, preorder_clades)
        branch_rates = self._compute_branch_rates(
            tree_for_analysis,
            node_values,
            parent_map,
            node_labels,
            preorder_clades,
        )

        # Plot
        self._plot_rate_map(
            tree_for_analysis,
            node_values,
            branch_rates,
            parent_map,
            node_labels,
            trait_values,
            trait_name,
            self.output_path,
        )

        # Compute summary statistics
        all_rates = [entry["rate"] for entry in branch_rates]
        mean_rate = sum(all_rates) / len(all_rates) if all_rates else 0.0

        min_entry = min(branch_rates, key=lambda e: e["rate"]) if branch_rates else None
        max_entry = max(branch_rates, key=lambda e: e["rate"]) if branch_rates else None

        # Text output
        self._print_text_output(n, trait_name, mean_rate, min_entry, max_entry)

        if self.json_output:
            result = {
                "n_taxa": n,
                "trait": trait_name,
                "mean_rate": float(mean_rate),
                "min_rate": {
                    "value": float(min_entry["rate"]),
                    "branch": f"{min_entry['parent']} -> {min_entry['child']}",
                } if min_entry else None,
                "max_rate": {
                    "value": float(max_entry["rate"]),
                    "branch": f"{max_entry['parent']} -> {max_entry['child']}",
                } if max_entry else None,
                "branches": [
                    {
                        "parent": e["parent"],
                        "child": e["child"],
                        "rate": float(e["rate"]),
                        "change": float(e["change"]),
                        "branch_length": float(e["branch_length"]),
                    }
                    for e in branch_rates
                ],
                "output_file": self.output_path,
            }
            print_json(result)

    def _print_text_output(
        self,
        n: int,
        trait_name: str,
        mean_rate: float,
        min_entry,
        max_entry,
    ) -> None:
        lines = [
            "Trait Rate Map",
            f"Taxa: {n}",
            f"Trait: {trait_name}",
            f"Mean rate: {mean_rate:.4f}",
        ]
        if min_entry:
            lines.append(
                f"Min rate: {min_entry['rate']:.4f} "
                f"(branch to {min_entry['child']})"
            )
        if max_entry:
            lines.append(
                f"Max rate: {max_entry['rate']:.4f} "
                f"(branch to {max_entry['child']})"
            )
        lines.append(f"Output: {self.output_path}")
        print("\n".join(lines))

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            output_path=args.output,
            trait_column=getattr(args, "trait", None),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    # ------------------------------------------------------------------
    # Trait data parsing
    # ------------------------------------------------------------------

    def _parse_single_trait_data(
        self, path: str, tree_tips: list[str]
    ) -> dict[str, float]:
        """Parse two-column TSV: taxon<tab>value, no header."""
        try:
            traits = {}
            with open(path) as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line[0] == "#":
                        continue
                    parts = line.split("\t", 2)
                    if len(parts) != 2:
                        column_count = line.count("\t") + 1
                        raise PhykitUserError(
                            [
                                f"Line {line_num} in trait file has {column_count} columns; expected 2.",
                                "Each line should be: taxon_name<tab>trait_value",
                            ],
                            code=2,
                        )
                    taxon, value_str = parts
                    try:
                        traits[taxon] = float(value_str)
                    except ValueError:
                        raise PhykitUserError(
                            [
                                f"Non-numeric trait value '{value_str}' for taxon '{taxon}' on line {line_num}.",
                            ],
                            code=2,
                        )
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        return self._intersect_traits(traits, tree_tips)

    def _parse_multi_trait_data(
        self, path: str, tree_tips: list[str], column: str
    ) -> dict[str, float]:
        """Parse multi-column TSV with header row, extracting one column."""
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
                        ["Multi-trait file must have a header row and at least one data row."],
                        code=2,
                    )

                if column not in header:
                    raise PhykitUserError(
                        [
                            f"Column '{column}' not found in trait file header.",
                            f"Available columns: {', '.join(header[1:])}",
                        ],
                        code=2,
                    )
                col_idx = header.index(column)
                split_limit = col_idx + 1

                traits = {}
                saw_data = False
                line_num = 2
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped[0] == "#":
                        continue
                    saw_data = True
                    parts = stripped.split("\t", split_limit)
                    if len(parts) <= col_idx:
                        raise PhykitUserError(
                            [f"Line {line_num} in trait file has too few columns."],
                            code=2,
                        )
                    taxon = parts[0]
                    value_str = parts[col_idx]
                    try:
                        traits[taxon] = float(value_str)
                    except ValueError:
                        raise PhykitUserError(
                            [
                                f"Non-numeric trait value '{value_str}' for taxon '{taxon}' on line {line_num}.",
                            ],
                            code=2,
                        )
                    line_num += 1

                if not saw_data:
                    raise PhykitUserError(
                        ["Multi-trait file must have a header row and at least one data row."],
                        code=2,
                    )
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        return self._intersect_traits(traits, tree_tips)

    def _intersect_traits(
        self, traits: dict[str, float], tree_tips: list[str]
    ) -> dict[str, float]:
        """Return shared taxa between traits and tree; warn about mismatches."""
        tree_tip_set = set(tree_tips)
        if (
            len(tree_tip_set) >= 3
            and len(tree_tip_set) == len(traits)
            and tree_tip_set == traits.keys()
        ):
            return traits

        trait_taxa_set = set(traits)
        shared = tree_tip_set & trait_taxa_set

        tree_only = tree_tip_set - trait_taxa_set
        trait_only = trait_taxa_set - tree_tip_set

        if tree_only:
            print(
                f"Warning: {len(tree_only)} taxa in tree but not in trait file: "
                f"{', '.join(sorted(tree_only))}",
                file=sys.stderr,
            )
        if trait_only:
            print(
                f"Warning: {len(trait_only)} taxa in trait file but not in tree: "
                f"{', '.join(sorted(trait_only))}",
                file=sys.stderr,
            )

        if len(shared) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa between tree and trait file.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        return {taxon: traits[taxon] for taxon in shared}

    # ------------------------------------------------------------------
    # Node labelling and parent map
    # ------------------------------------------------------------------

    def _label_internal_nodes(self, tree, preorder_clades=None) -> dict:
        """Assign N1, N2, ... labels to unnamed internal nodes."""
        labels = {}
        counter = 1
        if preorder_clades is None:
            preorder_clades = self._iter_preorder(tree.root)
        for clade in preorder_clades:
            if not clade.is_terminal():
                if clade.name:
                    labels[id(clade)] = clade.name
                else:
                    labels[id(clade)] = f"N{counter}"
                    counter += 1
        return labels

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

    @staticmethod
    def _iter_postorder(root):
        clades = []
        stack = [root]
        while stack:
            clade = stack.pop()
            clades.append(clade)
            stack.extend(clade.clades)
        yield from reversed(clades)

    def _build_parent_map(self, tree, preorder_clades=None) -> dict:
        parent_map = {}
        if preorder_clades is not None:
            for clade in preorder_clades:
                for child in clade.clades:
                    parent_map[id(child)] = clade
            return parent_map

        stack = [tree.root]
        pop = stack.pop
        extend = stack.extend
        while stack:
            clade = pop()
            children = clade.clades
            for child in children:
                parent_map[id(child)] = clade
            if children:
                extend(children)
        return parent_map

    # ------------------------------------------------------------------
    # Ancestral reconstruction (Felsenstein weighted average)
    # ------------------------------------------------------------------

    def _ancestral_reconstruction(
        self, tree, trait_values: dict[str, float], postorder_clades=None
    ) -> dict[int, float]:
        """Postorder weighted-average ancestral reconstruction.

        For each internal node, the reconstructed value is the
        inverse-branch-length-weighted average of its children's values.
        """
        node_values: dict[int, float] = {}

        if postorder_clades is None:
            postorder_clades = self._iter_postorder(tree.root)

        for clade in postorder_clades:
            if not clade.clades:
                if clade.name in trait_values:
                    node_values[id(clade)] = trait_values[clade.name]
                continue
            total_w = 0.0
            weighted_sum = 0.0
            for child in clade.clades:
                child_value = node_values.get(id(child))
                if child_value is None:
                    continue
                bl = child.branch_length if child.branch_length else _BL_FLOOR
                bl = max(bl, _BL_FLOOR)
                w = 1.0 / bl
                total_w += w
                weighted_sum += w * child_value
            if total_w:
                node_values[id(clade)] = weighted_sum / total_w

        return node_values

    # ------------------------------------------------------------------
    # Branch rate computation
    # ------------------------------------------------------------------

    def _compute_branch_rates(
        self, tree, node_values: dict[int, float],
        parent_map: dict, node_labels: dict, preorder_clades=None,
    ) -> list[dict]:
        """Compute per-branch evolutionary rate.

        Rate = (child_val - parent_val)^2 / branch_length
        (squared standardized contrast, proportional to local BM rate).
        """
        branch_rates = []
        root = tree.root
        root_id = id(root)
        if preorder_clades is None:
            preorder_clades = self._iter_preorder(root)

        for clade in preorder_clades:
            clade_id = id(clade)
            if clade_id == root_id:
                continue
            parent = parent_map.get(clade_id)
            if parent is None:
                continue
            parent_id = id(parent)
            parent_val = node_values.get(parent_id)
            if parent_val is None:
                continue
            child_val = node_values.get(clade_id)
            if child_val is None:
                continue

            bl = clade.branch_length if clade.branch_length else _BL_FLOOR
            bl = max(bl, _BL_FLOOR)

            change = child_val - parent_val
            rate = change * change / bl

            # Determine labels
            parent_label = node_labels.get(parent_id)
            if parent_label is None:
                parent_label = parent.name if parent.name else "root"

            if not clade.clades:
                child_label = clade.name
            else:
                child_label = node_labels.get(clade_id, "unknown")

            branch_rates.append({
                "parent": parent_label,
                "child": child_label,
                "rate": rate,
                "change": change,
                "branch_length": bl,
                "clade_id": clade_id,
            })

        return branch_rates

    def _prepare_rate_map_plot_data(self, tree, branch_rates):
        rate_by_clade = {}
        all_rates = []
        for entry in branch_rates:
            rate = entry["rate"]
            rate_by_clade[entry["clade_id"]] = rate
            all_rates.append(rate)

        preorder_clades = []
        tips = []
        stack = [tree.root]
        while stack:
            clade = stack.pop()
            preorder_clades.append(clade)
            children = clade.clades
            if children:
                stack.extend(reversed(children))
            else:
                tips.append(clade)

        return rate_by_clade, preorder_clades, tips, all_rates

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def _plot_rate_map(
        self, tree, node_values, branch_rates, parent_map,
        node_labels, trait_values, trait_name, output_path,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.colors import Normalize
        except ImportError:
            print(
                "matplotlib is required for trait rate map plotting. "
                "Install matplotlib and retry."
            )
            raise SystemExit(2)

        root = tree.root
        rate_by_clade, preorder_clades, tips, all_rates = (
            self._prepare_rate_map_plot_data(tree, branch_rates)
        )
        from ...helpers.plot_config import compute_node_positions

        node_x, node_y = compute_node_positions(
            tree,
            parent_map,
            cladogram=self.plot_config.cladogram,
            preorder_clades=preorder_clades,
        )

        # Normalize rates for colormap
        if not all_rates:
            return
        vmin = min(all_rates)
        vmax = max(all_rates)
        if vmin == vmax:
            vmax = vmin + 1.0
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap("inferno")

        config = self.plot_config
        config.resolve(n_rows=len(tips), n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        if self.plot_config.circular:
            from matplotlib.collections import LineCollection
            from ...helpers.circular_layout import (
                compute_circular_coords,
                draw_circular_tip_labels,
            )

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

            radial_segments = []
            radial_colors = []
            radial_widths = []
            rate_color_cache = {}
            for clade in preorder_clades:
                if clade == root:
                    continue
                cid = id(clade)
                if cid not in parent_map:
                    continue
                parent = parent_map[cid]
                pid = id(parent)
                if pid not in coords or cid not in coords:
                    continue

                color = "gray"
                lw = 2
                if cid in rate_by_clade:
                    rate = rate_by_clade[cid]
                    color = rate_color_cache.get(rate)
                    if color is None:
                        color = cmap(norm(rate))
                        rate_color_cache[rate] = color
                    lw = 3

                angle = coords[cid]["angle"]
                r_p = coords[pid]["radius"]
                r_c = coords[cid]["radius"]
                x0 = r_p * math.cos(angle)
                y0 = r_p * math.sin(angle)
                x1 = r_c * math.cos(angle)
                y1 = r_c * math.sin(angle)
                radial_segments.append(((x0, y0), (x1, y1)))
                radial_colors.append(color)
                radial_widths.append(lw)

            if radial_segments:
                ax.add_collection(
                    LineCollection(
                        radial_segments,
                        colors=radial_colors,
                        linewidths=radial_widths,
                        capstyle="round",
                        zorder=2,
                    ),
                    autolim=True,
                )

            # Arcs at internal nodes
            arc_segments = []
            arc_fractions = [idx / 60 for idx in range(61)]
            for clade in preorder_clades:
                if not clade.clades:
                    continue
                cid = id(clade)
                if cid not in coords:
                    continue
                child_angles = [
                    coords[id(ch)]["angle"]
                    for ch in clade.clades if id(ch) in coords
                ]
                if len(child_angles) < 2:
                    continue
                start_a = min(child_angles)
                end_a = max(child_angles)
                span = (end_a - start_a) % (2.0 * math.pi)
                if span > math.pi:
                    start_a, end_a = end_a, start_a

                start = start_a % (2.0 * math.pi)
                end = end_a % (2.0 * math.pi)
                diff = (end - start) % (2.0 * math.pi)
                if diff > math.pi:
                    diff = diff - 2.0 * math.pi
                radius = coords[cid]["radius"]
                arc_segments.append([
                    (
                        radius * math.cos(start + diff * fraction),
                        radius * math.sin(start + diff * fraction),
                    )
                    for fraction in arc_fractions
                ])

            if arc_segments:
                ax.add_collection(
                    LineCollection(
                        arc_segments,
                        colors="gray",
                        linewidths=1.5,
                        capstyle="round",
                        zorder=1,
                    ),
                    autolim=True,
                )
            if radial_segments or arc_segments:
                ax.autoscale_view()

            # Tip labels
            max_x = max(node_x.values()) if node_x else 1.0
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize is not None else 9
            if label_fontsize > 0:
                draw_circular_tip_labels(
                    ax, tree, coords, fontsize=label_fontsize, offset=max_x * 0.03
                )

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
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        draw_range_wedge(ax, tree, mrca, clr, coords)
                apply_label_colors(ax, color_data["labels"])
                color_legend = build_color_legend_handles(color_data)
                if color_legend:
                    ax.legend(
                        handles=color_legend, loc="upper right",
                        fontsize=8, frameon=True,
                    )

            # Colorbar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.15)
            cbar.set_label("Evolutionary rate")

            if config.show_title:
                ax.set_title(
                    config.title or "Trait Rate Map",
                    fontsize=config.title_fontsize,
                )

        else:
            # --- Rectangular mode ---
            from matplotlib.collections import LineCollection

            vertical_segments = []
            horizontal_segments = []
            horizontal_colors = []
            horizontal_widths = []
            rate_color_cache = {}
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

                x0 = node_x[pid]
                x1 = node_x[cid]
                y0 = node_y.get(pid, 0)
                y1 = node_y.get(cid, 0)

                # Horizontal branch colored by rate
                color = "gray"
                lw = 2
                if cid in rate_by_clade:
                    rate = rate_by_clade[cid]
                    color = rate_color_cache.get(rate)
                    if color is None:
                        color = cmap(norm(rate))
                        rate_color_cache[rate] = color
                    lw = 3

                horizontal_segments.append(((x0, y1), (x1, y1)))
                horizontal_colors.append(color)
                horizontal_widths.append(lw)
                vertical_segments.append(((x0, y0), (x0, y1)))

            if vertical_segments:
                ax.add_collection(
                    LineCollection(
                        vertical_segments,
                        colors="gray",
                        linewidths=1.5,
                    ),
                    autolim=True,
                )
            if horizontal_segments:
                ax.add_collection(
                    LineCollection(
                        horizontal_segments,
                        colors=horizontal_colors,
                        linewidths=horizontal_widths,
                    ),
                    autolim=True,
                )
            ax.autoscale_view()

            # Tip labels
            max_x = max(node_x.values()) if node_x else 0
            offset = max_x * 0.02
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize is not None else 9
            if label_fontsize > 0:
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
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        draw_range_rect(ax, tree, mrca, clr, node_x, node_y)
                apply_label_colors(ax, color_data["labels"])
                color_legend = build_color_legend_handles(color_data)
                if color_legend:
                    ax.legend(
                        handles=color_legend, loc="upper right",
                        fontsize=8, frameon=True,
                    )

            # Colorbar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.15)
            cbar.set_label("Evolutionary rate")

            ax.set_xlabel("Branch length")
            ax.set_yticks([])
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            if config.show_title:
                ax.set_title(
                    config.title or "Trait Rate Map",
                    fontsize=config.title_fontsize,
                )

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
