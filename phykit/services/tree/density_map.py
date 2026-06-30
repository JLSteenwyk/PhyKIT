from __future__ import annotations

import sys

from .base import Tree
from .stochastic_character_map import StochasticCharacterMap
from ...errors import PhykitUserError


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


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

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

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
        tree = self.read_tree_file_unmodified()

        # Use StochasticCharacterMap to do the heavy lifting
        scm = StochasticCharacterMap.__new__(StochasticCharacterMap)
        # Initialize the base Tree attributes manually
        scm.tree_file_path = self.tree_file_path
        scm.tree_format = "newick"

        scm.validate_tree(tree, min_tips=3, require_branch_lengths=True, context="density map")

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

        copied_tree = False
        missing_tip_state_count = scm._count_missing_tip_states(tree, tip_states)
        if missing_tip_state_count is None or missing_tip_state_count:
            tree = self._fast_copy(tree)
            copied_tree = True
            scm._prune_tree_to_tip_states(tree, tip_states)

        if self.plot_config.ladderize:
            if not copied_tree:
                tree = self._fast_copy(tree)
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
        simulation_metadata = scm._build_simulation_metadata(
            tree, tip_states, states, parent_map
        )

        # Run stochastic mappings
        rng = np.random.default_rng(self.seed)
        mappings = []
        transition_cache = {}
        for _ in range(self.n_sim):
            mapping = scm._run_single_simulation(
                tree, tip_states, Q, pi, states, rng,
                cond_liks=cond_liks, parent_map=parent_map,
                transition_cache=transition_cache,
                simulation_metadata=simulation_metadata,
            )
            mappings.append(mapping)

        # Plot the density map
        self._plot_density_map(
            tree, mappings, states, self.output_path,
            parent_map=parent_map, scm=scm,
        )

        # Text output
        tips = self._terminal_clades(tree)
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
            self._print_text_output(states, n_tips)

    def _print_text_output(self, states, n_tips: int) -> None:
        print(
            "\n".join(
                [
                    "Density Map",
                    f"\nStates: {len(states)} ({', '.join(states)})",
                    f"Number of tips: {n_tips}",
                    f"Number of simulations: {self.n_sim}",
                    f"Saved density map plot: {self.output_path}",
                ]
            )
        )

    @staticmethod
    def _get_state_at_position(history: list, t_query: float,
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
        self, mappings: list[dict], clade_id: int,
        branch_length: float, states: list[str],
        n_segments: int = 50,
    ) -> np.ndarray:
        """Compute posterior state probabilities at points along a branch.

        Returns an array of shape (n_segments, k) where each row gives the
        posterior probability of each state at the midpoint of that segment.
        """
        k = len(states)
        n_sim = len(mappings)
        posteriors = np.zeros((n_segments, k))
        histories_for_branch = []
        for mapping in mappings:
            history = mapping["branch_histories"].get(clade_id)
            if history:
                histories_for_branch.append(history)

        if branch_length <= 0:
            # For zero-length branches, just count states at position 0
            for history in histories_for_branch:
                posteriors[0, history[0][1]] += 1
            if n_sim > 0:
                posteriors[0] /= n_sim
            # Fill all segments with the same value
            if n_segments > 1:
                posteriors[1:] = posteriors[0]
            return posteriors

        if histories_for_branch:
            midpoint_scale = branch_length / n_segments
            for history in histories_for_branch:
                history_idx = 1
                history_len = len(history)
                state = history[0][1]
                for seg in range(n_segments):
                    t_mid = (seg + 0.5) * midpoint_scale
                    while (
                        history_idx < history_len
                        and history[history_idx][0] <= t_mid
                    ):
                        state = history[history_idx][1]
                        history_idx += 1
                    posteriors[seg, state] += 1

        if n_sim > 0:
            posteriors /= n_sim

        return posteriors

    def _plot_density_map(
        self, tree, mappings: list[dict], states: list[str],
        output_path: str, parent_map: dict = None,
        scm: StochasticCharacterMap = None,
    ) -> None:
        from ...helpers.plot_config import compute_node_positions

        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.collections import LineCollection
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
        state_color_matrix = np.asarray(state_colors_rgb)

        root = tree.root
        preorder_clades = list(self._iter_preorder(root))
        tips = [clade for clade in preorder_clades if not clade.clades]
        node_x, node_y = compute_node_positions(
            tree,
            parent_map,
            cladogram=self.plot_config.cladogram,
            preorder_clades=preorder_clades,
        )

        config = self.plot_config
        config.resolve(n_rows=len(tips), n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        if self.plot_config.circular:
            # --- Circular mode ---
            from ...helpers.circular_layout import (
                _arc_fractions_array,
                compute_circular_coords,
                draw_circular_tip_labels,
            )
            import math as _math

            coords = compute_circular_coords(
                tree,
                node_x,
                parent_map,
                preorder_clades=preorder_clades,
                terminal_clades=tips,
            )
            ax.set_aspect("equal")
            ax.axis("off")

            branch_segments = []
            branch_colors = []
            frac_starts = np.arange(n_segments) / n_segments
            frac_ends = (np.arange(n_segments) + 1) / n_segments
            for clade in preorder_clades:
                if clade == root:
                    continue
                parent = parent_map.get(id(clade))
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

                color_values = np.clip(posteriors @ state_color_matrix, 0, 1)
                angle = coords[cid]["angle"]
                cos_angle = _math.cos(angle)
                sin_angle = _math.sin(angle)
                r_p = coords[pid]["radius"]
                r_c = coords[cid]["radius"]

                if t > 0:
                    radius_delta = r_c - r_p
                    r0s = r_p + radius_delta * frac_starts
                    r1s = r_p + radius_delta * frac_ends
                    for seg, color in enumerate(color_values):
                        branch_segments.append((
                            (r0s[seg] * cos_angle, r0s[seg] * sin_angle),
                            (r1s[seg] * cos_angle, r1s[seg] * sin_angle),
                        ))
                        branch_colors.append(color)
                else:
                    branch_segments.append((
                        (r_p * cos_angle, r_p * sin_angle),
                        (r_c * cos_angle, r_c * sin_angle),
                    ))
                    branch_colors.append(color_values[0])

            if branch_segments:
                ax.add_collection(
                    LineCollection(
                        branch_segments,
                        colors=branch_colors,
                        linewidths=2.5,
                        capstyle="butt",
                        zorder=2,
                    ),
                    autolim=True,
                )

            # Arcs at internal nodes colored by parent state (gray)
            arc_segments = []
            arc_fractions = _arc_fractions_array()
            tau = 2.0 * _math.pi
            for clade in preorder_clades:
                if not clade.clades:
                    continue
                cid = id(clade)
                if cid not in coords:
                    continue
                child_angles = [coords[id(ch)]["angle"] for ch in clade.clades if id(ch) in coords]
                if len(child_angles) < 2:
                    continue
                min_a = min(child_angles)
                max_a = max(child_angles)
                start = min_a % tau
                end = max_a % tau
                diff = (end - start) % tau
                if diff > _math.pi:
                    diff = diff - tau
                angles = start + diff * arc_fractions
                radius = coords[cid]["radius"]
                arc_segments.append(
                    np.column_stack((radius * np.cos(angles), radius * np.sin(angles)))
                )

            if arc_segments:
                ax.add_collection(
                    LineCollection(
                        arc_segments,
                        colors="gray",
                        linewidths=0.8,
                        capstyle="round",
                        zorder=1,
                    ),
                    autolim=True,
                )
            if branch_segments or arc_segments:
                ax.autoscale_view()

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
            vertical_segments = []
            branch_segments = []
            branch_colors = []
            frac_starts = np.arange(n_segments) / n_segments
            frac_ends = (np.arange(n_segments) + 1) / n_segments
            for clade in preorder_clades:
                if clade == root:
                    continue
                parent = parent_map.get(id(clade))
                if parent is None:
                    continue

                parent_x = node_x[id(parent)]
                parent_y = node_y[id(parent)]
                child_x = node_x[id(clade)]
                child_y = node_y[id(clade)]
                t = clade.branch_length if clade.branch_length else 0.0

                vertical_segments.append(((parent_x, parent_y), (parent_x, child_y)))

                # Compute posterior probabilities along this branch
                clade_id = id(clade)
                posteriors = self._compute_branch_posteriors(
                    mappings, clade_id, t, states, n_segments=n_segments
                )

                color_values = np.clip(posteriors @ state_color_matrix, 0, 1)
                if t > 0:
                    x0s = parent_x + frac_starts * t
                    x1s = parent_x + frac_ends * t
                    for seg, color in enumerate(color_values):
                        branch_segments.append((
                            (x0s[seg], child_y),
                            (x1s[seg], child_y),
                        ))
                        branch_colors.append(color)
                else:
                    branch_segments.append(((parent_x, child_y), (child_x, child_y)))
                    branch_colors.append(color_values[0])

            if vertical_segments:
                ax.add_collection(
                    LineCollection(
                        vertical_segments,
                        colors="gray",
                        linewidths=0.8,
                        zorder=1,
                    ),
                    autolim=True,
                )
            if branch_segments:
                ax.add_collection(
                    LineCollection(
                        branch_segments,
                        colors=branch_colors,
                        linewidths=2.5,
                        capstyle="butt",
                        zorder=2,
                    ),
                    autolim=True,
                )
            if vertical_segments or branch_segments:
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
    def _terminal_clades(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return list(tree.get_terminals())

        tips = []
        stack = [root]
        pop = stack.pop
        append = stack.append
        append_tip = tips.append
        while stack:
            clade = pop()
            children = clade.clades
            if children:
                child_count = len(children)
                if child_count == 2:
                    append(children[1])
                    append(children[0])
                else:
                    for index in range(child_count - 1, -1, -1):
                        append(children[index])
            else:
                append_tip(clade)
        return tips
