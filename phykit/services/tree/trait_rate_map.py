import copy
import math
import sys
from typing import Dict, List, Tuple

import numpy as np

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig, compute_node_x_cladogram
from ...helpers.circular_layout import (
    compute_circular_coords,
    draw_circular_tip_labels,
    draw_circular_colored_branch,
    draw_circular_colored_arc,
)
from ...helpers.color_annotations import (
    parse_color_file,
    resolve_mrca,
    draw_range_rect,
    draw_range_wedge,
    build_color_legend_handles,
)
from ...errors import PhykitUserError


_BL_FLOOR = 1e-8


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
        tree = self.read_tree_file()
        self._validate_tree(tree)

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
        tree_copy = copy.deepcopy(tree)
        tip_names_in_tree = [t.name for t in tree_copy.get_terminals()]
        tips_to_prune = [t for t in tip_names_in_tree if t not in trait_values]
        if tips_to_prune:
            tree_copy = self.prune_tree_using_taxa_list(tree_copy, tips_to_prune)

        if self.plot_config.ladderize:
            tree_copy.ladderize()

        # Label internal nodes
        node_labels = self._label_internal_nodes(tree_copy)

        # Reconstruct ancestral states (Felsenstein weighted average)
        node_values = self._ancestral_reconstruction(tree_copy, trait_values)

        # Compute per-branch rates
        parent_map = self._build_parent_map(tree_copy)
        branch_rates = self._compute_branch_rates(
            tree_copy, node_values, parent_map, node_labels
        )

        # Plot
        self._plot_rate_map(
            tree_copy, node_values, branch_rates, parent_map,
            node_labels, trait_values, trait_name, self.output_path,
        )

        # Compute summary statistics
        all_rates = [entry["rate"] for entry in branch_rates]
        mean_rate = np.mean(all_rates) if all_rates else 0.0

        min_entry = min(branch_rates, key=lambda e: e["rate"]) if branch_rates else None
        max_entry = max(branch_rates, key=lambda e: e["rate"]) if branch_rates else None

        # Text output
        print("Trait Rate Map")
        print(f"Taxa: {n}")
        print(f"Trait: {trait_name}")
        print(f"Mean rate: {mean_rate:.4f}")
        if min_entry:
            print(f"Min rate: {min_entry['rate']:.4f} (branch to {min_entry['child']})")
        if max_entry:
            print(f"Max rate: {max_entry['rate']:.4f} (branch to {max_entry['child']})")
        print(f"Output: {self.output_path}")

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

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            output_path=args.output,
            trait_column=getattr(args, "trait", None),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for trait rate map visualization."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                raise PhykitUserError(
                    ["All branches in the tree must have lengths."],
                    code=2,
                )

    # ------------------------------------------------------------------
    # Trait data parsing
    # ------------------------------------------------------------------

    def _parse_single_trait_data(
        self, path: str, tree_tips: List[str]
    ) -> Dict[str, float]:
        """Parse two-column TSV: taxon<tab>value, no header."""
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

        traits = {}
        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise PhykitUserError(
                    [
                        f"Line {line_num} in trait file has {len(parts)} columns; expected 2.",
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

        return self._intersect_traits(traits, tree_tips)

    def _parse_multi_trait_data(
        self, path: str, tree_tips: List[str], column: str
    ) -> Dict[str, float]:
        """Parse multi-column TSV with header row, extracting one column."""
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

        clean_lines = [l.strip() for l in lines if l.strip() and not l.startswith("#")]
        if len(clean_lines) < 2:
            raise PhykitUserError(
                ["Multi-trait file must have a header row and at least one data row."],
                code=2,
            )

        header = clean_lines[0].split("\t")
        if column not in header:
            raise PhykitUserError(
                [
                    f"Column '{column}' not found in trait file header.",
                    f"Available columns: {', '.join(header[1:])}",
                ],
                code=2,
            )
        col_idx = header.index(column)

        traits = {}
        for line_num, line in enumerate(clean_lines[1:], 2):
            parts = line.split("\t")
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

        return self._intersect_traits(traits, tree_tips)

    def _intersect_traits(
        self, traits: Dict[str, float], tree_tips: List[str]
    ) -> Dict[str, float]:
        """Return shared taxa between traits and tree; warn about mismatches."""
        tree_tip_set = set(tree_tips)
        trait_taxa_set = set(traits.keys())
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

    def _label_internal_nodes(self, tree) -> Dict:
        """Assign N1, N2, ... labels to unnamed internal nodes."""
        labels = {}
        counter = 1
        for clade in tree.find_clades(order="preorder"):
            if not clade.is_terminal():
                if clade.name:
                    labels[id(clade)] = clade.name
                else:
                    labels[id(clade)] = f"N{counter}"
                    counter += 1
        return labels

    def _build_parent_map(self, tree) -> Dict:
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    # ------------------------------------------------------------------
    # Ancestral reconstruction (Felsenstein weighted average)
    # ------------------------------------------------------------------

    def _ancestral_reconstruction(
        self, tree, trait_values: Dict[str, float]
    ) -> Dict[int, float]:
        """Postorder weighted-average ancestral reconstruction.

        For each internal node, the reconstructed value is the
        inverse-branch-length-weighted average of its children's values.
        """
        node_values: Dict[int, float] = {}

        # Tips: observed values
        for tip in tree.get_terminals():
            if tip.name in trait_values:
                node_values[id(tip)] = trait_values[tip.name]

        # Internal nodes (postorder)
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                continue
            weights = []
            vals = []
            for child in clade.clades:
                if id(child) not in node_values:
                    continue
                bl = child.branch_length if child.branch_length else _BL_FLOOR
                bl = max(bl, _BL_FLOOR)
                w = 1.0 / bl
                weights.append(w)
                vals.append(node_values[id(child)])
            if weights:
                total_w = sum(weights)
                node_values[id(clade)] = sum(
                    w * v for w, v in zip(weights, vals)
                ) / total_w

        return node_values

    # ------------------------------------------------------------------
    # Branch rate computation
    # ------------------------------------------------------------------

    def _compute_branch_rates(
        self, tree, node_values: Dict[int, float],
        parent_map: Dict, node_labels: Dict,
    ) -> List[Dict]:
        """Compute per-branch evolutionary rate.

        Rate = (child_val - parent_val)^2 / branch_length
        (squared standardized contrast, proportional to local BM rate).
        """
        branch_rates = []
        root = tree.root

        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            if id(clade) not in parent_map:
                continue
            parent = parent_map[id(clade)]
            if id(parent) not in node_values or id(clade) not in node_values:
                continue

            parent_val = node_values[id(parent)]
            child_val = node_values[id(clade)]
            bl = clade.branch_length if clade.branch_length else _BL_FLOOR
            bl = max(bl, _BL_FLOOR)

            change = child_val - parent_val
            rate = change ** 2 / bl

            # Determine labels
            if id(parent) in node_labels:
                parent_label = node_labels[id(parent)]
            else:
                parent_label = parent.name if parent.name else "root"

            if clade.is_terminal():
                child_label = clade.name
            elif id(clade) in node_labels:
                child_label = node_labels[id(clade)]
            else:
                child_label = "unknown"

            branch_rates.append({
                "parent": parent_label,
                "child": child_label,
                "rate": rate,
                "change": change,
                "branch_length": bl,
                "clade_id": id(clade),
            })

        return branch_rates

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

        # Build rate lookup: id(clade) -> rate
        rate_by_clade = {}
        for entry in branch_rates:
            rate_by_clade[entry["clade_id"]] = entry["rate"]

        tips = list(tree.get_terminals())
        node_x = {}
        node_y = {}

        # Assign tip y-positions
        for i, tip in enumerate(tips):
            node_y[id(tip)] = i

        # Assign x-positions
        root = tree.root
        if self.plot_config.cladogram:
            node_x = compute_node_x_cladogram(tree, parent_map)
        else:
            for clade in tree.find_clades(order="preorder"):
                if clade == root:
                    node_x[id(clade)] = 0.0
                else:
                    if id(clade) in parent_map:
                        parent = parent_map[id(clade)]
                        t = clade.branch_length if clade.branch_length else 0.0
                        node_x[id(clade)] = node_x[id(parent)] + t

        # Assign internal y-positions
        for clade in tree.find_clades(order="postorder"):
            if not clade.is_terminal() and id(clade) not in node_y:
                child_ys = [
                    node_y[id(c)] for c in clade.clades if id(c) in node_y
                ]
                if child_ys:
                    node_y[id(clade)] = np.mean(child_ys)
                else:
                    node_y[id(clade)] = 0.0

        # Normalize rates for colormap
        all_rates = [entry["rate"] for entry in branch_rates]
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
            # --- Circular mode ---
            coords = compute_circular_coords(tree, node_x, parent_map)
            ax.set_aspect("equal")
            ax.axis("off")

            for clade in tree.find_clades(order="preorder"):
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
                    color = cmap(norm(rate_by_clade[cid]))
                    lw = 3

                draw_circular_colored_branch(
                    ax, coords[pid], coords[cid], color, lw=lw
                )

            # Arcs at internal nodes
            for clade in tree.find_clades(order="preorder"):
                if clade.is_terminal() or not clade.clades:
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

                draw_circular_colored_arc(
                    ax, 0, 0, coords[cid]["radius"],
                    start_a, end_a, color="gray", lw=1.5,
                )

            # Tip labels
            max_x = max(node_x.values()) if node_x else 1.0
            draw_circular_tip_labels(
                ax, tree, coords, fontsize=9, offset=max_x * 0.03
            )

            # Color annotations
            if self.plot_config.color_file:
                color_data = parse_color_file(self.plot_config.color_file)
                for taxa_list, clr, lbl in color_data["ranges"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        draw_range_wedge(ax, tree, mrca, clr, coords)
                for taxon, lbl_color in color_data["labels"].items():
                    for text_obj in ax.texts:
                        if text_obj.get_text() == taxon:
                            text_obj.set_color(lbl_color)
                            break
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

                # Horizontal branch colored by rate
                color = "gray"
                lw = 2
                if id(clade) in rate_by_clade:
                    color = cmap(norm(rate_by_clade[id(clade)]))
                    lw = 3

                ax.plot([x0, x1], [y1, y1], color=color, lw=lw)

                # Vertical connector in gray
                ax.plot([x0, x0], [y0, y1], color="gray", lw=1.5)

            # Tip labels
            max_x = max(node_x.values()) if node_x else 0
            offset = max_x * 0.02
            for tip in tips:
                ax.text(
                    node_x[id(tip)] + offset, node_y[id(tip)],
                    tip.name, va="center", fontsize=9,
                )

            # Color annotations
            if self.plot_config.color_file:
                color_data = parse_color_file(self.plot_config.color_file)
                for taxa_list, clr, lbl in color_data["ranges"]:
                    mrca = resolve_mrca(tree, taxa_list)
                    if mrca is not None:
                        draw_range_rect(ax, tree, mrca, clr, node_x, node_y)
                for taxon, lbl_color in color_data["labels"].items():
                    for text_obj in ax.texts:
                        if text_obj.get_text() == taxon:
                            text_obj.set_color(lbl_color)
                            break
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
