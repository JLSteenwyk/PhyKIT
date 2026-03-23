import sys
from typing import Dict, List

import numpy as np

from .base import Tree
from .stochastic_character_map import StochasticCharacterMap
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig, compute_node_x_cladogram
from ...helpers.circular_layout import (
    compute_circular_coords,
    draw_circular_tip_labels,
    draw_circular_colored_branch,
    draw_circular_colored_arc,
    draw_circular_multi_segment_branch,
)
from ...helpers.color_annotations import (
    parse_color_file,
    resolve_mrca,
    draw_range_rect,
    draw_range_wedge,
    build_color_legend_handles,
)
from ...errors import PhykitUserError


class DensityMap(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.trait_column = parsed["trait_column"]
        self.n_sim = parsed["n_sim"]
        self.seed = parsed["seed"]
        self.output_path = parsed["output_path"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            trait_column=args.trait,
            n_sim=getattr(args, "nsim", 100),
            seed=getattr(args, "seed", None),
            output_path=args.output,
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        tree = self.read_tree_file()

        # Use StochasticCharacterMap to do the heavy lifting
        scm = StochasticCharacterMap.__new__(StochasticCharacterMap)
        # Initialize the base Tree attributes manually
        scm.tree_file_path = self.tree_file_path
        scm.tree_format = "newick"

        scm._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)
        tip_states = scm._parse_discrete_trait_file(
            self.trait_data_path, self.trait_column, tree_tips
        )

        states = sorted(set(tip_states.values()))
        if len(states) < 2:
            raise PhykitUserError(
                ["At least 2 distinct character states are required."],
                code=2,
            )

        # Prune tree to shared taxa
        tips_in_tree = set(self.get_tip_names_from_tree(tree))
        tips_to_prune = [t for t in tips_in_tree if t not in tip_states]
        if tips_to_prune:
            for tip_name in tips_to_prune:
                tree.prune(tip_name)

        if self.plot_config.ladderize:
            tree.ladderize()

        # Fit Q matrix using ER model
        Q, loglik = scm._fit_q_matrix(tree, tip_states, states, "ER")
        k = len(states)
        pi = np.ones(k) / k

        # Compute conditional likelihoods once
        cond_liks, _ = scm._felsenstein_pruning(
            tree, tip_states, Q, pi, states
        )

        # Build parent lookup
        parent_map = scm._build_parent_map(tree)

        # Run stochastic mappings
        rng = np.random.default_rng(self.seed)
        mappings = []
        for _ in range(self.n_sim):
            mapping = scm._run_single_simulation(
                tree, tip_states, Q, pi, states, rng,
                cond_liks=cond_liks, parent_map=parent_map,
            )
            mappings.append(mapping)

        # Plot the density map
        self._plot_density_map(
            tree, mappings, states, self.output_path,
            parent_map=parent_map, scm=scm,
        )

        # Text output
        tips = list(tree.get_terminals())
        n_tips = len(tips)

        if self.json_output:
            result = {
                "n_tips": n_tips,
                "states": states,
                "n_sim": self.n_sim,
                "plot_output": self.output_path,
            }
            print_json(result)
        else:
            print("Density Map")
            print(
                f"\nStates: {len(states)} "
                f"({', '.join(states)})"
            )
            print(f"Number of tips: {n_tips}")
            print(f"Number of simulations: {self.n_sim}")
            print(f"Saved density map plot: {self.output_path}")

    @staticmethod
    def _get_state_at_position(history: List, t_query: float,
                               branch_length: float) -> int:
        """Find which state is active at position t_query along a branch.

        Walk through the (time, state) pairs until the cumulative position
        exceeds t_query.
        """
        current_state = history[0][1]
        for h_idx in range(1, len(history)):
            if history[h_idx][0] > t_query:
                break
            current_state = history[h_idx][1]
        return current_state

    def _compute_branch_posteriors(
        self, mappings: List[Dict], clade_id: int,
        branch_length: float, states: List[str],
        n_segments: int = 50,
    ) -> np.ndarray:
        """Compute posterior state probabilities at points along a branch.

        Returns an array of shape (n_segments, k) where each row gives the
        posterior probability of each state at the midpoint of that segment.
        """
        k = len(states)
        n_sim = len(mappings)
        posteriors = np.zeros((n_segments, k))

        if branch_length <= 0:
            # For zero-length branches, just count states at position 0
            for mapping in mappings:
                histories = mapping["branch_histories"]
                if clade_id in histories:
                    state = histories[clade_id][0][1]
                    posteriors[0, state] += 1
            if n_sim > 0:
                posteriors[0] /= n_sim
            # Fill all segments with the same value
            for seg in range(1, n_segments):
                posteriors[seg] = posteriors[0]
            return posteriors

        for seg in range(n_segments):
            frac_s = seg / n_segments
            frac_e = (seg + 1) / n_segments
            t_mid = (frac_s + frac_e) / 2.0 * branch_length

            for mapping in mappings:
                histories = mapping["branch_histories"]
                if clade_id not in histories:
                    continue
                history = histories[clade_id]
                state = self._get_state_at_position(
                    history, t_mid, branch_length
                )
                posteriors[seg, state] += 1

            if n_sim > 0:
                posteriors[seg] /= n_sim

        return posteriors

    def _plot_density_map(
        self, tree, mappings: List[Dict], states: List[str],
        output_path: str, parent_map: Dict = None,
        scm: StochasticCharacterMap = None,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.lines import Line2D
            from matplotlib.colors import LinearSegmentedColormap
            import matplotlib.colorbar as mcolorbar
        except ImportError:
            print(
                "matplotlib is required for density map plotting. "
                "Install matplotlib and retry."
            )
            raise SystemExit(2)

        if parent_map is None and scm is not None:
            parent_map = scm._build_parent_map(tree)

        k = len(states)
        n_segments = 50

        # Define state colors
        tab10 = plt.get_cmap("tab10")
        state_colors_rgb = []
        for i in range(k):
            state_colors_rgb.append(
                np.array(tab10(i / max(k - 1, 1))[:3])
            )

        # Compute node positions: x = distance from root, y = tip ordering
        node_x = {}
        node_y = {}

        tips = list(tree.get_terminals())
        for i, tip in enumerate(tips):
            node_y[id(tip)] = i

        root = tree.root
        if self.plot_config.cladogram:
            node_x = compute_node_x_cladogram(tree, parent_map)
        else:
            for clade in tree.find_clades(order="preorder"):
                if clade == root:
                    node_x[id(clade)] = 0.0
                else:
                    parent = scm._get_parent(tree, clade, parent_map)
                    if parent is not None:
                        t = clade.branch_length if clade.branch_length else 0.0
                        node_x[id(clade)] = node_x[id(parent)] + t

        for clade in tree.find_clades(order="postorder"):
            if not clade.is_terminal() and id(clade) not in node_y:
                child_ys = [
                    node_y[id(c)] for c in clade.clades if id(c) in node_y
                ]
                if child_ys:
                    node_y[id(clade)] = np.mean(child_ys)
                else:
                    node_y[id(clade)] = 0.0

        config = self.plot_config
        config.resolve(n_rows=len(tips), n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        if self.plot_config.circular:
            # --- Circular mode ---
            # Build parent_map keyed by id(child) -> parent clade
            circ_parent_map = {}
            for clade in tree.find_clades(order="preorder"):
                for child in clade.clades:
                    circ_parent_map[id(child)] = clade
            coords = compute_circular_coords(tree, node_x, circ_parent_map)
            ax.set_aspect("equal")
            ax.axis("off")

            for clade in tree.find_clades(order="preorder"):
                if clade == root:
                    continue
                parent = scm._get_parent(tree, clade, parent_map)
                if parent is None:
                    continue

                pid = id(parent)
                cid = id(clade)
                if pid not in coords or cid not in coords:
                    continue

                t = clade.branch_length if clade.branch_length else 0.0

                # Compute posterior probabilities along this branch
                posteriors = self._compute_branch_posteriors(
                    mappings, cid, t, states, n_segments=n_segments
                )

                # Build multi-segment data: each segment gets a blended color
                # We use draw_circular_colored_branch for each sub-segment
                if t > 0:
                    for seg in range(n_segments):
                        frac_s = seg / n_segments
                        frac_e = (seg + 1) / n_segments
                        # Compute weighted color
                        color = np.zeros(3)
                        for si in range(k):
                            color += posteriors[seg, si] * state_colors_rgb[si]
                        color = np.clip(color, 0, 1)

                        # Draw sub-segment along the radial branch
                        import math as _math
                        angle = coords[cid]["angle"]
                        r_p = coords[pid]["radius"]
                        r_c = coords[cid]["radius"]
                        r0 = r_p + (r_c - r_p) * frac_s
                        r1 = r_p + (r_c - r_p) * frac_e
                        x0 = r0 * _math.cos(angle)
                        y0 = r0 * _math.sin(angle)
                        x1 = r1 * _math.cos(angle)
                        y1 = r1 * _math.sin(angle)
                        ax.plot([x0, x1], [y0, y1], color=color,
                                linewidth=2.5, solid_capstyle="butt")
                else:
                    # Zero-length branch
                    color = np.zeros(3)
                    for si in range(k):
                        color += posteriors[0, si] * state_colors_rgb[si]
                    color = np.clip(color, 0, 1)
                    draw_circular_colored_branch(
                        ax, coords[pid], coords[cid], color=color, lw=2.5
                    )

            # Arcs at internal nodes colored by parent state (gray)
            for clade in tree.find_clades(order="preorder"):
                if clade.is_terminal() or not clade.clades:
                    continue
                cid = id(clade)
                if cid not in coords:
                    continue
                child_angles = [coords[id(ch)]["angle"] for ch in clade.clades if id(ch) in coords]
                if len(child_angles) < 2:
                    continue
                min_a = min(child_angles)
                max_a = max(child_angles)
                draw_circular_colored_arc(
                    ax, 0, 0, coords[cid]["radius"],
                    min_a, max_a, color="gray", lw=0.8,
                )

            # Tip labels
            max_x = max(node_x.values()) if node_x else 1.0
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize is not None else 9
            if label_fontsize > 0:
                draw_circular_tip_labels(ax, tree, coords, fontsize=label_fontsize, offset=max_x * 0.03)

            # Apply color annotations (range + label only; branches are trait-colored)
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

            # Add colorbar for 2-state case or legend for k-state
            if k == 2:
                cmap_colors = [state_colors_rgb[0], state_colors_rgb[1]]
                custom_cmap = LinearSegmentedColormap.from_list(
                    "density_grad", cmap_colors, N=256
                )
                sm = plt.cm.ScalarMappable(
                    cmap=custom_cmap,
                    norm=plt.Normalize(vmin=0, vmax=1),
                )
                sm.set_array([])
                cbar = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
                cbar.set_label(f"P({states[1]})", fontsize=9)
                cbar.set_ticks([0.0, 0.5, 1.0])
                cbar.set_ticklabels([
                    f"{states[0]} (1.0)", "0.5", f"{states[1]} (1.0)",
                ])

            handles = [
                Line2D([0], [0], color=state_colors_rgb[i],
                       linewidth=3, label=states[i])
                for i in range(k)
            ]
            if self.plot_config.color_file:
                color_legend = build_color_legend_handles(color_data)
                handles.extend(color_legend)
            ax.legend(handles=handles, title="States", loc="upper left",
                      fontsize=8, title_fontsize=9)

            if config.show_title:
                ax.set_title(config.title or "Density Map", fontsize=config.title_fontsize)
        else:
            # --- Rectangular mode ---
            for clade in tree.find_clades(order="preorder"):
                if clade == root:
                    continue
                parent = scm._get_parent(tree, clade, parent_map)
                if parent is None:
                    continue

                parent_x = node_x[id(parent)]
                parent_y = node_y[id(parent)]
                child_x = node_x[id(clade)]
                child_y = node_y[id(clade)]
                t = clade.branch_length if clade.branch_length else 0.0

                # Draw vertical connector at parent x
                ax.plot(
                    [parent_x, parent_x], [parent_y, child_y],
                    color="gray", linewidth=0.8, zorder=1
                )

                # Compute posterior probabilities along this branch
                clade_id = id(clade)
                posteriors = self._compute_branch_posteriors(
                    mappings, clade_id, t, states, n_segments=n_segments
                )

                # Draw horizontal branch as colored segments
                if t > 0:
                    for seg in range(n_segments):
                        frac_s = seg / n_segments
                        frac_e = (seg + 1) / n_segments
                        x0 = parent_x + frac_s * t
                        x1 = parent_x + frac_e * t

                        # Compute weighted color from state posteriors
                        color = np.zeros(3)
                        for si in range(k):
                            color += posteriors[seg, si] * state_colors_rgb[si]
                        color = np.clip(color, 0, 1)

                        ax.plot(
                            [x0, x1], [child_y, child_y],
                            color=color, linewidth=2.5,
                            solid_capstyle="butt", zorder=2,
                        )
                else:
                    # Zero-length branch: just draw a point
                    color = np.zeros(3)
                    for si in range(k):
                        color += posteriors[0, si] * state_colors_rgb[si]
                    color = np.clip(color, 0, 1)
                    ax.plot(
                        [parent_x, child_x], [child_y, child_y],
                        color=color, linewidth=2.5, zorder=2,
                    )

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

            # Add colorbar for 2-state case or legend for k-state
            if k == 2:
                # Create a gradient colorbar between the two state colors
                cmap_colors = [state_colors_rgb[0], state_colors_rgb[1]]
                custom_cmap = LinearSegmentedColormap.from_list(
                    "density_grad", cmap_colors, N=256
                )
                sm = plt.cm.ScalarMappable(
                    cmap=custom_cmap,
                    norm=plt.Normalize(vmin=0, vmax=1),
                )
                sm.set_array([])
                cbar = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
                cbar.set_label(
                    f"P({states[1]})",
                    fontsize=9,
                )
                cbar.set_ticks([0.0, 0.5, 1.0])
                cbar.set_ticklabels([
                    f"{states[0]} (1.0)",
                    "0.5",
                    f"{states[1]} (1.0)",
                ])

            # Always add a legend with state names and colors
            handles = [
                Line2D(
                    [0], [0], color=state_colors_rgb[i],
                    linewidth=3, label=states[i],
                )
                for i in range(k)
            ]
            if self.plot_config.color_file:
                color_legend = build_color_legend_handles(color_data)
                handles.extend(color_legend)
            ax.legend(
                handles=handles, title="States", loc="upper left",
                fontsize=8, title_fontsize=9,
            )

            ax.set_xlabel("Branch length")
            ax.set_yticks([])
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            if config.show_title:
                ax.set_title(config.title or "Density Map", fontsize=config.title_fontsize)

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
