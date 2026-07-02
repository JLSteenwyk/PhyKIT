from __future__ import annotations

import math
import sys

from .base import Tree
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


def _value_range(values):
    iterator = iter(values)
    try:
        first = next(iterator)
    except StopIteration:
        return None

    vmin = vmax = first
    for value in iterator:
        if value < vmin:
            vmin = value
        elif value > vmax:
            vmax = value
    return vmin, vmax


class ContMap(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.output_path = parsed["output_path"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        self.validate_tree(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="contMap visualization",
        )

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_values = self._parse_single_trait_data(
            self.trait_data_path, tree_tips
        )

        ordered_names = sorted(trait_values.keys())
        n = len(ordered_names)
        x = np.array([trait_values[name] for name in ordered_names])

        # Prune tree to shared taxa
        tips_to_prune = self._tips_to_prune_for_ordered_mapping(
            tree_tips, trait_values
        )
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

        # Label internal nodes
        node_labels = self._label_internal_nodes(tree_for_analysis)

        # Run fast ancestral reconstruction
        node_estimates, sigma2 = self._fast_anc(
            tree_for_analysis, x, ordered_names, node_labels
        )

        # Always produce contMap plot
        self._plot_contmap(
            tree_for_analysis, node_estimates, node_labels,
            trait_values, "trait", self.output_path
        )

        # Text output (always printed)
        self._print_text_output(n, sigma2)

        if self.json_output:
            result = {
                "n_tips": n,
                "sigma2": float(sigma2),
                "tip_values": {k: float(v) for k, v in sorted(trait_values.items())},
                "ancestral_estimates": {
                    label: float(node_estimates[label])
                    for label in sorted(node_estimates.keys())
                },
                "plot_output": self.output_path,
            }
            print_json(result)

    def _print_text_output(self, n: int, sigma2: float) -> None:
        print(
            "\n".join(
                [
                    "Continuous Trait Map (contMap)",
                    f"\nNumber of tips: {n}",
                    f"Sigma-squared (BM rate): {sigma2:.4f}",
                    f"Saved contMap plot: {self.output_path}",
                ]
            )
        )

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            output_path=args.output,
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def _parse_single_trait_data(
        self, path: str, tree_tips: list[str]
    ) -> dict[str, float]:
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

    def _label_internal_nodes(self, tree) -> dict:
        """Assign N1, N2, ... labels to unnamed internal nodes.

        Returns dict mapping id(clade) -> label for all internal nodes.
        """
        labels = {}
        counter = 1
        for clade in self._iter_preorder(tree.root):
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

    def _build_parent_map(self, tree) -> dict:
        parent_map = {}
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

    def _prepare_contmap_plot_data(
        self, tree, node_estimates, node_labels, trait_values
    ):
        parent_map = {}
        preorder_clades = []
        all_estimates = {}
        tips = []

        stack = [tree.root]
        while stack:
            clade = stack.pop()
            preorder_clades.append(clade)
            children = clade.clades
            if children:
                label = node_labels.get(id(clade))
                if label in node_estimates:
                    all_estimates[id(clade)] = node_estimates[label]
                for child in children:
                    parent_map[id(child)] = clade
                stack.extend(reversed(children))
            else:
                tips.append(clade)
                if clade.name in trait_values:
                    all_estimates[id(clade)] = trait_values[clade.name]

        return parent_map, preorder_clades, all_estimates, tips

    def _fast_anc(
        self, tree, x: np.ndarray, ordered_names: list[str],
        node_labels: dict,
    ) -> tuple[dict, float]:
        """Two-pass ML ancestral reconstruction (Felsenstein's algorithm).

        Pass 1 (postorder): compute subtree estimates and variances.
        Pass 2 (preorder): incorporate information from the rest of the
        tree so that every internal node gets the full-tree ML estimate.

        Returns (node_estimates, sigma2).
        """
        name_to_idx = {name: i for i, name in enumerate(ordered_names)}
        preorder_clades = list(self._iter_preorder(tree.root))
        postorder_clades = list(reversed(preorder_clades))
        parent_map = {}
        for clade in preorder_clades:
            for child in clade.clades:
                parent_map[id(child)] = clade

        # --- Pass 1: postorder (subtree estimates) ---
        est_down = {}
        var_down = {}

        for clade in postorder_clades:
            children = clade.clades
            clade_id = id(clade)
            if not children:
                idx = name_to_idx.get(clade.name)
                if idx is not None:
                    est_down[clade_id] = x[idx]
                    var_down[clade_id] = 0.0
            else:
                total_prec = 0.0
                weighted_sum = 0.0
                for child in children:
                    child_id = id(child)
                    child_estimate = est_down.get(child_id)
                    if child_estimate is None:
                        continue
                    v_i = child.branch_length if child.branch_length else 0.0
                    child_var = var_down[child_id]
                    denom = v_i + child_var
                    if denom == 0:
                        denom = 1e-10
                    prec = 1.0 / denom
                    total_prec += prec
                    weighted_sum += prec * child_estimate

                if total_prec:
                    est_down[clade_id] = weighted_sum / total_prec
                    var_down[clade_id] = 1.0 / total_prec

        # --- Pass 2: preorder (incorporate full-tree information) ---
        est_up = {}
        var_up = {}
        est_final = {}
        var_final = {}

        root = tree.root
        root_id = id(root)
        # Root has no information from above
        est_final[root_id] = est_down[root_id]
        var_final[root_id] = var_down[root_id]

        for clade in preorder_clades:
            clade_id = id(clade)
            if clade_id == root_id:
                continue
            parent = parent_map.get(clade_id)
            if parent is None:
                continue
            parent_id = id(parent)
            bl = clade.branch_length if clade.branch_length else 0.0

            # Gather precision from parent's "above" message
            prec_above = 0.0
            est_above_weighted = 0.0
            parent_estimate = est_up.get(parent_id)
            if parent_estimate is not None:
                v_up_parent = var_up[parent_id]
                if v_up_parent > 0:
                    p = 1.0 / v_up_parent
                    prec_above += p
                    est_above_weighted += p * parent_estimate
                elif v_up_parent == 0 and parent_id != root_id:
                    p = 1.0 / 1e-10
                    prec_above += p
                    est_above_weighted += p * parent_estimate

            # Gather precision from siblings
            for sibling in parent.clades:
                sibling_id = id(sibling)
                if sibling_id == clade_id:
                    continue
                sibling_estimate = est_down.get(sibling_id)
                if sibling_estimate is None:
                    continue
                sib_bl = sibling.branch_length if sibling.branch_length else 0.0
                sib_var = var_down[sibling_id] + sib_bl
                if sib_var == 0:
                    sib_var = 1e-10
                p = 1.0 / sib_var
                prec_above += p
                est_above_weighted += p * sibling_estimate

            # The "above" message for this node
            if prec_above > 0:
                up_estimate = est_above_weighted / prec_above
                est_up[clade_id] = up_estimate
                up_variance = bl + 1.0 / prec_above
                var_up[clade_id] = up_variance
            else:
                # No information from above (shouldn't happen in a valid tree)
                up_estimate = est_down.get(clade_id, 0.0)
                est_up[clade_id] = up_estimate
                up_variance = bl + 1e10
                var_up[clade_id] = up_variance

            # Combine down and up for final estimate
            down_estimate = est_down.get(clade_id)
            if down_estimate is not None:
                vd = var_down[clade_id]
                vu = up_variance
                if vd == 0:
                    vd = 1e-10
                if vu == 0:
                    vu = 1e-10
                prec_d = 1.0 / vd
                prec_u = 1.0 / vu
                total = prec_d + prec_u
                est_final[clade_id] = (
                    prec_d * down_estimate + prec_u * up_estimate
                ) / total
                var_final[clade_id] = 1.0 / total

        # Compute sigma^2 from phylogenetic independent contrasts
        sigma2 = self._compute_sigma2_from_contrasts(
            tree, x, est_down, var_down, ordered_names, postorder_clades
        )

        # Build output dicts keyed by node label
        node_estimates = {}
        for clade in preorder_clades:
            if clade.clades:
                clade_id = id(clade)
                label = node_labels.get(clade_id)
                final_estimate = est_final.get(clade_id)
                if label is not None and final_estimate is not None:
                    node_estimates[label] = final_estimate

        return node_estimates, sigma2

    def _compute_sigma2_from_contrasts(
        self, tree, x: np.ndarray, node_estimates: dict,
        node_variances: dict, ordered_names: list[str], postorder_clades=None,
    ) -> float:
        """Compute BM rate (sigma^2) from phylogenetic independent contrasts."""
        contrasts_sum = 0.0
        contrast_count = 0
        if postorder_clades is None:
            postorder_clades = self._iter_postorder(tree.root)
        for clade in postorder_clades:
            if not clade.clades:
                continue
            children = clade.clades
            first_child = None
            combined_est = 0.0
            combined_var = 0.0
            valid_count = 0

            for child in children:
                child_id = id(child)
                child_estimate = node_estimates.get(child_id)
                if child_estimate is None:
                    continue

                if valid_count == 0:
                    first_child = child
                    valid_count = 1
                    continue

                if valid_count == 1:
                    first_id = id(first_child)
                    v1 = (
                        first_child.branch_length or 0.0
                    ) + node_variances[first_id]
                    v2 = (child.branch_length or 0.0) + node_variances[child_id]
                    denom = v1 + v2
                    if denom > 0:
                        diff = node_estimates[first_id] - child_estimate
                        contrasts_sum += (diff * diff) / denom
                        contrast_count += 1
                        combined_est = (
                            v2 * node_estimates[first_id] + v1 * child_estimate
                        ) / denom
                        combined_var = (v1 * v2) / denom
                    else:
                        combined_est = node_estimates[first_id]
                        combined_var = 0.0
                    valid_count = 2
                    continue

                vk = (child.branch_length or 0.0) + node_variances[child_id]
                vc = combined_var
                denom_k = vk + vc
                if denom_k > 0:
                    diff = child_estimate - combined_est
                    contrasts_sum += (diff * diff) / denom_k
                    contrast_count += 1
                    combined_est = (
                        vc * child_estimate + vk * combined_est
                    ) / denom_k
                    combined_var = (vk * vc) / denom_k

        if contrast_count:
            return contrasts_sum / contrast_count
        return 0.0

    def _plot_contmap(
        self, tree, node_estimates, node_labels, trait_values,
        trait_name, output_path,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.colors import Normalize
        except ImportError:
            print(
                "matplotlib is required for contMap plotting. "
                "Install matplotlib and retry."
            )
            raise SystemExit(2)

        root = tree.root
        parent_map, preorder_clades, all_estimates, tips = (
            self._prepare_contmap_plot_data(
                tree, node_estimates, node_labels, trait_values
            )
        )
        from ...helpers.plot_config import compute_node_positions

        node_x, node_y = compute_node_positions(
            tree,
            parent_map,
            cladogram=self.plot_config.cladogram,
            preorder_clades=preorder_clades,
        )

        # Normalize trait values for colormap
        value_range = _value_range(all_estimates.values())
        if value_range is None:
            return
        vmin, vmax = value_range
        if vmin == vmax:
            vmax = vmin + 1.0
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap("coolwarm")

        config = self.plot_config
        config.resolve(n_rows=len(tips), n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        if self.plot_config.circular:
            from ...helpers.circular_layout import (
                compute_circular_coords,
                draw_circular_colored_arcs,
                draw_circular_gradient_branches,
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

            gradient_branches = []
            for clade in preorder_clades:
                if clade == root:
                    continue
                if id(clade) not in parent_map:
                    continue
                parent = parent_map[id(clade)]
                pid = id(parent)
                cid = id(clade)
                if pid not in coords or cid not in coords:
                    continue
                if pid not in all_estimates or cid not in all_estimates:
                    continue

                parent_val = all_estimates[pid]
                child_val = all_estimates[cid]

                gradient_branches.append(
                    (coords[pid], coords[cid], parent_val, child_val)
                )
            draw_circular_gradient_branches(
                ax, gradient_branches, cmap, vmin, vmax, lw=3,
            )

            # Arcs at internal nodes colored by trait value
            colored_arcs = []
            for clade in preorder_clades:
                if clade.is_terminal() or not clade.clades:
                    continue
                cid = id(clade)
                if cid not in coords or cid not in all_estimates:
                    continue
                child_angles = [coords[id(ch)]["angle"] for ch in clade.clades if id(ch) in coords]
                if len(child_angles) < 2:
                    continue
                node_val = all_estimates[cid]
                min_a = min(child_angles)
                max_a = max(child_angles)
                colored_arcs.append(
                    (0, 0, coords[cid]["radius"], min_a, max_a, cmap(norm(node_val)))
                )
            draw_circular_colored_arcs(ax, colored_arcs, lw=3)

            # Tip labels
            max_x = max(node_x.values()) if node_x else 1.0
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize is not None else 9
            if label_fontsize > 0:
                draw_circular_tip_labels(ax, tree, coords, fontsize=label_fontsize, offset=max_x * 0.03)

            # Apply color annotations (range + label only; branches are trait-colored)
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
                    ax.legend(handles=color_legend, loc="upper right", fontsize=8, frameon=True)

            # Colorbar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.15)
            cbar.set_label(trait_name)

            if config.show_title:
                ax.set_title(config.title or "Continuous Trait Map (contMap)", fontsize=config.title_fontsize)
        else:
            # --- Rectangular mode ---
            from matplotlib.collections import LineCollection

            n_seg = 50
            horizontal_segments = []
            horizontal_colors = []
            vertical_segments = []
            vertical_colors = []

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
                if pid not in all_estimates or cid not in all_estimates:
                    continue

                x0 = node_x[pid]
                x1 = node_x[cid]
                y0 = node_y[pid]
                y1 = node_y[cid]
                val0 = all_estimates[pid]
                val1 = all_estimates[cid]

                # Horizontal segment (at child's y): colored gradient
                for s in range(n_seg):
                    frac_s = s / n_seg
                    frac_e = (s + 1) / n_seg
                    xs = x0 + frac_s * (x1 - x0)
                    xe = x0 + frac_e * (x1 - x0)
                    v = val0 + ((frac_s + frac_e) / 2) * (val1 - val0)
                    horizontal_segments.append(((xs, y1), (xe, y1)))
                    horizontal_colors.append(cmap(norm(v)))

                # Vertical connector (at parent's x): parent's trait color
                vertical_segments.append(((x0, y0), (x0, y1)))
                vertical_colors.append(cmap(norm(val0)))

            if horizontal_segments:
                ax.add_collection(
                    LineCollection(
                        horizontal_segments,
                        colors=horizontal_colors,
                        linewidths=3,
                        capstyle="butt",
                    ),
                    autolim=True,
                )
            if vertical_segments:
                ax.add_collection(
                    LineCollection(
                        vertical_segments,
                        colors=vertical_colors,
                        linewidths=3,
                        capstyle="butt",
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

            # Apply color annotations (range + label only; branches are trait-colored)
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
                    ax.legend(handles=color_legend, loc="upper right", fontsize=8, frameon=True)

            # Colorbar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.15)
            cbar.set_label(trait_name)

            ax.set_xlabel("Branch length")
            ax.set_yticks([])
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            if config.show_title:
                ax.set_title(config.title or "Continuous Trait Map (contMap)", fontsize=config.title_fontsize)

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
