"""
Chronogram: time-calibrated phylogeny with geological timescale.

Plots an ultrametric tree with geological epoch/period/era bands
and a time axis in millions of years ago (Ma).
"""
import math
from typing import Dict, List

import numpy as np

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import (
    PlotConfig,
    build_parent_map,
    compute_node_positions,
    draw_tip_labels,
    cleanup_tree_axes,
)
from ...helpers.geological_timescale import get_timescale_for_range
from ...helpers.circular_layout import (
    compute_circular_coords,
    draw_circular_tip_labels,
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


class Chronogram(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.root_age = parsed["root_age"]
        self.timescale = parsed["timescale"]
        self.node_ages = parsed["node_ages"]
        self.plot_output = parsed["plot_output"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            root_age=args.root_age,
            timescale=getattr(args, "timescale", "auto"),
            node_ages=getattr(args, "node_ages", False),
            plot_output=args.plot_output,
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        tree = self.read_tree_file()
        self.validate_tree(
            tree, min_tips=3, require_branch_lengths=True,
            context="chronogram",
        )

        if self.plot_config.ladderize:
            tree.ladderize()

        parent_map = build_parent_map(tree)

        # Compute node ages from root
        root_to_tip = self._compute_root_to_tip(tree)
        max_height = max(root_to_tip.values()) if root_to_tip else 1.0

        # Scale factor: convert tree units to Ma
        scale = self.root_age / max_height if max_height > 0 else 1.0

        # Parse HPD confidence intervals (auto-detected from annotations)
        hpd_intervals = self._parse_hpd_intervals(tree, self.root_age, scale)

        if self.plot_config.circular:
            self._plot_circular(tree, parent_map, scale, root_to_tip, hpd_intervals)
        else:
            self._plot_rectangular(tree, parent_map, scale, root_to_tip, hpd_intervals)

        if self.json_output:
            self._print_json(tree, parent_map, scale, root_to_tip, hpd_intervals)

    def _compute_root_to_tip(self, tree) -> Dict[int, float]:
        """Compute root-to-node distance for every node."""
        distances = {}
        for clade in tree.find_clades(order="preorder"):
            if clade == tree.root:
                distances[id(clade)] = 0.0
            else:
                parent_dist = 0.0
                path = tree.get_path(clade)
                if path:
                    for node in path:
                        if node.branch_length:
                            parent_dist += node.branch_length
                distances[id(clade)] = parent_dist
        return distances

    @staticmethod
    def _parse_hpd_intervals(tree, root_age, scale) -> Dict[int, tuple]:
        """Parse 95% HPD intervals from BEAST or MCMCTree annotations.

        Looks for these patterns in node comments:
          BEAST:    [&height_95%_HPD={lower,upper}]
          MCMCTree: [&95%HPD={lower,upper}]

        If annotations use 'height' (distance from tips), they are
        converted to age. If they already represent age, they are
        used directly.

        Returns dict mapping node id -> (age_lower, age_upper).
        """
        import re
        intervals = {}

        # Detect which format: height-based or age-based
        # BEAST uses height (distance from tips); MCMCTree may use age
        height_pattern = re.compile(
            r'height_95%_HPD\s*=\s*\{?\s*([\d.eE+-]+)\s*,\s*([\d.eE+-]+)\s*\}?'
        )
        hpd_pattern = re.compile(
            r'95%HPD\s*=\s*\{?\s*([\d.eE+-]+)\s*,\s*([\d.eE+-]+)\s*\}?'
        )
        # Also try height_range or CI
        ci_pattern = re.compile(
            r'(?:CI|HPD|hpd_range)\s*=\s*\{?\s*([\d.eE+-]+)\s*,\s*([\d.eE+-]+)\s*\}?'
        )

        for clade in tree.find_clades(order="preorder"):
            comment = clade.comment or ""
            if not comment:
                continue

            match = height_pattern.search(comment)
            if match:
                # BEAST height_95%_HPD: values are heights from tips
                h_lo = float(match.group(1))
                h_hi = float(match.group(2))
                # Convert height to age: these are already ages in BEAST
                # (height = time from present)
                intervals[id(clade)] = (
                    min(h_lo, h_hi), max(h_lo, h_hi)
                )
                continue

            match = hpd_pattern.search(comment)
            if match:
                lo = float(match.group(1))
                hi = float(match.group(2))
                intervals[id(clade)] = (min(lo, hi), max(lo, hi))
                continue

            match = ci_pattern.search(comment)
            if match:
                lo = float(match.group(1))
                hi = float(match.group(2))
                intervals[id(clade)] = (min(lo, hi), max(lo, hi))
                continue

        return intervals

    # ---- Rectangular mode ----

    def _plot_rectangular(self, tree, parent_map, scale, root_to_tip, hpd_intervals=None):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.patches import Rectangle
        except ImportError:
            print("matplotlib is required for chronogram plotting.")
            raise SystemExit(2)

        config = self.plot_config
        tips = list(tree.get_terminals())
        n_tips = len(tips)
        config.resolve(n_rows=n_tips, n_cols=None)

        fig, (ax_tree, ax_scale) = plt.subplots(
            2, 1, figsize=(config.fig_width, config.fig_height),
            gridspec_kw={"height_ratios": [20, 1], "hspace": 0.02},
            sharex=True,
        )

        # Compute positions
        node_y = {}
        for i, tip in enumerate(tips):
            node_y[id(tip)] = i

        for clade in tree.find_clades(order="postorder"):
            if not clade.is_terminal() and id(clade) not in node_y:
                child_ys = [
                    node_y[id(c)] for c in clade.clades if id(c) in node_y
                ]
                node_y[id(clade)] = float(np.mean(child_ys)) if child_ys else 0.0

        node_x = {}
        for cid, dist in root_to_tip.items():
            node_x[cid] = self.root_age - (dist * scale)

        # Geological timescale bands + scale bar
        intervals, colors = get_timescale_for_range(
            self.root_age, self.timescale
        )
        y_min = -0.5
        y_max = n_tips - 0.5

        for name, start_ma, end_ma in intervals:
            if start_ma <= 0 or end_ma >= self.root_age * 1.05:
                continue
            x_left = max(end_ma, 0)
            x_right = min(start_ma, self.root_age * 1.05)
            color = colors.get(name, "#f0f0f0")

            # Bands on tree panel
            rect = Rectangle(
                (x_left, y_min), x_right - x_left, y_max - y_min,
                facecolor=color, alpha=0.20, edgecolor="none", zorder=0,
            )
            ax_tree.add_patch(rect)

            # Faint vertical boundary lines at epoch borders
            ax_tree.axvline(
                x_left, color="#cccccc", lw=0.3, zorder=0, linestyle="-",
            )

            # Scale bar below the tree with labels inside
            scale_rect = Rectangle(
                (x_left, 0), x_right - x_left, 1,
                facecolor=color, alpha=0.5, edgecolor="#999999",
                linewidth=0.3, zorder=1,
            )
            ax_scale.add_patch(scale_rect)

            # Label inside scale bar (horizontal)
            band_width = x_right - x_left
            mid_x = (x_left + x_right) / 2
            if band_width > self.root_age * 0.06:
                ax_scale.text(
                    mid_x, 0.5, name,
                    ha="center", va="center", fontsize=5.5,
                    color="#333333", zorder=2, clip_on=True,
                )
            elif band_width > self.root_age * 0.025:
                ax_scale.text(
                    mid_x, 0.5, name,
                    ha="center", va="center", fontsize=4,
                    color="#333333", rotation=90, zorder=2, clip_on=True,
                )

        # Draw branches
        root = tree.root
        branch_color = "#2c2c2c"
        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            pid = id(parent_map.get(id(clade), root))
            cid = id(clade)
            if pid not in node_x or cid not in node_x:
                continue

            x0 = node_x[pid]
            x1 = node_x[cid]
            y0 = node_y.get(pid, 0)
            y1 = node_y.get(cid, 0)

            ax_tree.plot(
                [x0, x1], [y1, y1], color=branch_color,
                lw=1.2, solid_capstyle="round", zorder=2,
            )
            ax_tree.plot(
                [x0, x0], [y0, y1], color=branch_color,
                lw=0.8, solid_capstyle="round", zorder=2,
            )

        # Color annotations
        if config.color_file:
            color_data = parse_color_file(config.color_file)
            for taxa_list, clr, lbl in color_data["ranges"]:
                mrca = resolve_mrca(tree, taxa_list)
                if mrca is not None:
                    draw_range_rect(ax_tree, tree, mrca, clr, node_x, node_y)
            for taxa_list, clade_color, lbl in color_data["clades"]:
                mrca = resolve_mrca(tree, taxa_list)
                if mrca is not None:
                    clade_ids = get_clade_branch_ids(tree, mrca, parent_map)
                    for cl in tree.find_clades(order="preorder"):
                        if cl == root:
                            continue
                        if id(cl) in clade_ids and id(cl) in parent_map:
                            p = id(parent_map[id(cl)])
                            c = id(cl)
                            ax_tree.plot(
                                [node_x[p], node_x[c]],
                                [node_y.get(c, 0), node_y.get(c, 0)],
                                color=clade_color, lw=1.5, zorder=3,
                            )
                            ax_tree.plot(
                                [node_x[p], node_x[p]],
                                [node_y.get(p, 0), node_y.get(c, 0)],
                                color=clade_color, lw=1.5, zorder=3,
                            )
            for taxon, lbl_color in color_data["labels"].items():
                for text_obj in ax_tree.texts:
                    if text_obj.get_text() == taxon:
                        text_obj.set_color(lbl_color)
                        break

        # Tip labels (italic for species names)
        label_fs = config.ylabel_fontsize if config.ylabel_fontsize else 9
        if label_fs > 0:
            offset = self.root_age * 0.008
            for tip in tips:
                ax_tree.text(
                    node_x[id(tip)] - offset, node_y[id(tip)],
                    tip.name, va="center", ha="right",
                    fontsize=label_fs, fontstyle="italic",
                    zorder=4,
                )

        # 95% HPD confidence interval bars
        if hpd_intervals:
            bar_height = max(0.15, 0.6 / max(n_tips, 1))
            for clade in tree.find_clades(order="preorder"):
                if clade.is_terminal() or clade == root:
                    continue
                cid = id(clade)
                if cid not in hpd_intervals:
                    continue
                lo, hi = hpd_intervals[cid]
                cy = node_y.get(cid, 0)
                ax_tree.barh(
                    cy, hi - lo, left=lo, height=bar_height,
                    color="#2b8cbe", alpha=0.25, edgecolor="none",
                    zorder=1,
                )

        # Node age labels
        if self.node_ages:
            for clade in tree.find_clades(order="preorder"):
                if clade.is_terminal() or clade == root:
                    continue
                cid = id(clade)
                age = node_x.get(cid, 0)
                ax_tree.text(
                    age, node_y.get(cid, 0) + 0.15,
                    f"{age:.1f}", ha="center", va="bottom",
                    fontsize=5, color="#555555",
                    bbox=dict(facecolor="white", edgecolor="none",
                              alpha=0.8, pad=0.5),
                    zorder=5,
                )

        # Tree panel axes
        ax_tree.set_xlim(self.root_age * 1.03, -self.root_age * 0.12)
        ax_tree.set_ylim(y_min - 0.3, y_max + 0.3)
        ax_tree.set_yticks([])
        for spine in ax_tree.spines.values():
            spine.set_visible(False)
        ax_tree.tick_params(axis="x", which="both", length=0, labelbottom=False)

        # Scale bar axes
        ax_scale.set_xlim(self.root_age * 1.03, -self.root_age * 0.12)
        ax_scale.set_ylim(0, 1)
        ax_scale.set_yticks([])
        for spine in ax_scale.spines.values():
            spine.set_visible(False)
        ax_scale.set_xlabel("Millions of years ago (Ma)", fontsize=9)
        ax_scale.tick_params(axis="x", labelsize=7)

        if config.axis_fontsize:
            ax_scale.xaxis.label.set_fontsize(config.axis_fontsize)

        if config.show_title:
            ax_tree.set_title(
                config.title or "Chronogram",
                fontsize=config.title_fontsize, pad=10,
            )

        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        if not self.json_output:
            print(f"Chronogram saved: {self.plot_output}")

    # ---- Circular mode ----

    def _plot_circular(self, tree, parent_map, scale, root_to_tip, hpd_intervals=None):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.patches import Wedge
        except ImportError:
            print("matplotlib is required for chronogram plotting.")
            raise SystemExit(2)

        config = self.plot_config
        tips = list(tree.get_terminals())
        n_tips = len(tips)
        config.resolve(n_rows=n_tips, n_cols=None)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        ax.set_aspect("equal")
        ax.axis("off")

        # Compute node_x (distance from root) for circular layout
        node_x = {}
        for clade in tree.find_clades(order="preorder"):
            cid = id(clade)
            if clade == tree.root:
                node_x[cid] = 0.0
            else:
                parent = parent_map.get(id(clade))
                if parent is not None:
                    t = clade.branch_length if clade.branch_length else 0.0
                    node_x[cid] = node_x.get(id(parent), 0.0) + t

        coords = compute_circular_coords(tree, node_x, parent_map)

        max_radius = max(
            c["radius"] for c in coords.values()
        ) if coords else 1.0
        radius_scale = self.root_age / max_radius if max_radius > 0 else 1.0

        # Draw geological timescale as concentric rings
        intervals, colors = get_timescale_for_range(
            self.root_age, self.timescale
        )
        for name, start_ma, end_ma in intervals:
            if start_ma <= 0:
                continue
            r_inner = (self.root_age - start_ma) / self.root_age * max_radius * radius_scale
            r_outer = (self.root_age - end_ma) / self.root_age * max_radius * radius_scale
            if r_inner < 0:
                r_inner = 0
            if r_outer <= r_inner:
                continue

            color = colors.get(name, "#f0f0f0")
            ring = Wedge(
                (0, 0), r_outer, 0, 360,
                width=r_outer - r_inner,
                facecolor=color, alpha=0.20,
                edgecolor="none", zorder=0,
            )
            ax.add_patch(ring)

            # Faint concentric boundary circle
            circle_pts = np.linspace(0, 2 * np.pi, 200)
            ax.plot(
                r_inner * np.cos(circle_pts),
                r_inner * np.sin(circle_pts),
                color="#cccccc", lw=0.3, zorder=0,
            )

        # Radial age tick marks
        tick_interval = self._nice_tick_interval(self.root_age)
        tick_angle = -np.pi / 2 - 0.15  # just below bottom
        for age in np.arange(0, self.root_age + tick_interval, tick_interval):
            r = (self.root_age - age) / self.root_age * max_radius * radius_scale
            if r < 0:
                continue
            # Faint radial tick line
            tick_len = max_radius * radius_scale * 0.015
            tx = r * math.cos(tick_angle)
            ty = r * math.sin(tick_angle)
            tx2 = (r + tick_len) * math.cos(tick_angle)
            ty2 = (r + tick_len) * math.sin(tick_angle)
            ax.plot([tx, tx2], [ty, ty2], color="#888888", lw=0.5, zorder=1)
            # Label
            lx = (r + tick_len * 2.5) * math.cos(tick_angle)
            ly = (r + tick_len * 2.5) * math.sin(tick_angle)
            label = str(int(age)) if age == int(age) else f"{age:.0f}"
            ax.text(
                lx, ly, label, ha="center", va="top",
                fontsize=5, color="#666666", zorder=1,
            )

        # "Ma" label
        r_label = (max_radius * radius_scale + max_radius * radius_scale * 0.08)
        lx = r_label * math.cos(tick_angle)
        ly = r_label * math.sin(tick_angle)
        ax.text(
            lx, ly, "Ma", ha="center", va="top",
            fontsize=6, color="#666666", fontstyle="italic", zorder=1,
        )

        # Draw branches
        branch_color = "#2c2c2c"
        root = tree.root
        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            parent = parent_map.get(id(clade))
            if parent is None:
                continue
            pid = id(parent)
            cid = id(clade)
            if pid not in coords or cid not in coords:
                continue

            angle = coords[cid]["angle"]
            r_p = coords[pid]["radius"]
            r_c = coords[cid]["radius"]

            x0 = r_p * math.cos(angle)
            y0 = r_p * math.sin(angle)
            x1 = r_c * math.cos(angle)
            y1 = r_c * math.sin(angle)
            ax.plot(
                [x0, x1], [y0, y1], color=branch_color,
                lw=1.0, solid_capstyle="round", zorder=2,
            )

        # Draw arcs at internal nodes
        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal() or not clade.clades:
                continue
            cid = id(clade)
            if cid not in coords:
                continue
            child_angles = [
                coords[id(ch)]["angle"] for ch in clade.clades
                if id(ch) in coords
            ]
            if len(child_angles) < 2:
                continue
            min_a = min(child_angles)
            max_a = max(child_angles)
            r = coords[cid]["radius"]
            n_pts = max(20, int((max_a - min_a) * 50))
            angles = np.linspace(min_a, max_a, n_pts)
            xs = r * np.cos(angles)
            ys = r * np.sin(angles)
            ax.plot(xs, ys, color=branch_color, lw=0.8, zorder=2)

        # 95% HPD confidence interval arcs
        if hpd_intervals:
            for clade in tree.find_clades(order="preorder"):
                if clade.is_terminal() or clade == root:
                    continue
                cid = id(clade)
                if cid not in hpd_intervals or cid not in coords:
                    continue
                lo, hi = hpd_intervals[cid]
                angle = coords[cid]["angle"]
                # Convert age to radius
                r_lo = (self.root_age - hi) / self.root_age * max_radius * radius_scale
                r_hi = (self.root_age - lo) / self.root_age * max_radius * radius_scale
                if r_lo < 0:
                    r_lo = 0
                # Draw a thick translucent radial line for the CI
                x0 = r_lo * math.cos(angle)
                y0 = r_lo * math.sin(angle)
                x1 = r_hi * math.cos(angle)
                y1 = r_hi * math.sin(angle)
                ax.plot(
                    [x0, x1], [y0, y1],
                    color="#2b8cbe", alpha=0.3, lw=4,
                    solid_capstyle="round", zorder=1,
                )

        # Tip labels (italic)
        label_fs = config.ylabel_fontsize if config.ylabel_fontsize else 8
        if label_fs > 0:
            offset = max_radius * 0.03
            for tip in tips:
                cid = id(tip)
                if cid not in coords:
                    continue
                angle = coords[cid]["angle"]
                r = coords[cid]["radius"] + offset
                x = r * math.cos(angle)
                y = r * math.sin(angle)
                deg = math.degrees(angle)
                ha = "left"
                rotation = deg
                if 90 < deg < 270 or -270 < deg < -90:
                    rotation = deg + 180
                    ha = "right"
                ax.text(
                    x, y, tip.name, ha=ha, va="center",
                    fontsize=label_fs, fontstyle="italic",
                    rotation=rotation, rotation_mode="anchor", zorder=4,
                )

        # Color annotations
        if config.color_file:
            color_data = parse_color_file(config.color_file)
            for taxa_list, clr, lbl in color_data["ranges"]:
                mrca = resolve_mrca(tree, taxa_list)
                if mrca is not None:
                    draw_range_wedge(ax, tree, mrca, clr, coords)
            for taxon, lbl_color in color_data["labels"].items():
                for text_obj in ax.texts:
                    if text_obj.get_text() == taxon:
                        text_obj.set_color(lbl_color)
                        break

        if config.show_title:
            ax.set_title(
                config.title or "Chronogram",
                fontsize=config.title_fontsize,
            )

        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        if not self.json_output:
            print(f"Chronogram saved: {self.plot_output}")

    @staticmethod
    def _nice_tick_interval(root_age):
        """Choose a nice tick interval for age axis."""
        if root_age <= 10:
            return 2
        elif root_age <= 50:
            return 10
        elif root_age <= 150:
            return 20
        elif root_age <= 300:
            return 50
        else:
            return 100

    # ---- JSON output ----

    def _print_json(self, tree, parent_map, scale, root_to_tip, hpd_intervals=None):
        max_height = max(root_to_tip.values()) if root_to_tip else 1.0
        nodes = []
        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            cid = id(clade)
            dist = root_to_tip.get(cid, 0.0)
            age = self.root_age - (dist * scale)
            desc = sorted(t.name for t in clade.get_terminals())
            entry = {
                "descendants": desc,
                "age_ma": round(age, 4),
                "n_descendants": len(desc),
            }
            if hpd_intervals and cid in hpd_intervals:
                lo, hi = hpd_intervals[cid]
                entry["hpd_lower"] = round(lo, 4)
                entry["hpd_upper"] = round(hi, 4)
            nodes.append(entry)
        payload = {
            "root_age": self.root_age,
            "n_tips": tree.count_terminals(),
            "timescale": self.timescale,
            "node_ages": nodes,
            "plot_output": self.plot_output,
        }
        print_json(payload, sort_keys=False)
