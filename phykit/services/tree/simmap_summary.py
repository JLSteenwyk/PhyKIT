"""
SIMMAP summary: run N stochastic character maps and provide per-branch
summaries of dwelling times, transition counts, and node posteriors.

Extends StochasticCharacterMap with detailed per-branch output matching
the output of phytools::describe.simmap() in R.
"""
from __future__ import annotations

import sys

from .stochastic_character_map import StochasticCharacterMap
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


class SimmapSummary(StochasticCharacterMap):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        from .base import Tree
        Tree.__init__(self, tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.trait_column = parsed["trait_column"]
        self.model = parsed["model"]
        self.nsim = parsed["nsim"]
        self.seed = parsed["seed"]
        self.plot_output = parsed["plot_output"]
        self.json_output = parsed["json_output"]
        self.csv_output = parsed["csv_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            trait_column=args.trait,
            model=getattr(args, "model", "ER"),
            nsim=getattr(args, "nsim", 100),
            seed=getattr(args, "seed", None),
            plot_output=getattr(args, "plot", None),
            json_output=getattr(args, "json", False),
            csv_output=getattr(args, "csv", None),
            plot_config=PlotConfig.from_args(args),
        )

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
    def _branch_clades_preorder(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            root = None

        if root is not None:
            try:
                return [
                    clade for clade in SimmapSummary._iter_preorder(root)
                    if clade != root
                ]
            except AttributeError:
                pass

        return [
            clade for clade in tree.find_clades(order="preorder")
            if clade != tree.root
        ]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        self.validate_tree(tree, min_tips=3, require_branch_lengths=True, context="SIMMAP summary")

        tree_tips = self.get_tip_names_from_tree(tree)
        tip_states = self._parse_discrete_trait_file(
            self.trait_data_path, self.trait_column, tree_tips
        )

        states = sorted(set(tip_states.values()))
        if len(states) < 2:
            raise PhykitUserError(
                ["At least 2 distinct character states are required."],
                code=2,
            )

        copied_tree = False
        missing_tip_state_count = self._count_missing_tip_states(tree, tip_states)
        if missing_tip_state_count is None or missing_tip_state_count:
            tree = self._fast_copy(tree)
            copied_tree = True
            self._prune_tree_to_tip_states(tree, tip_states)

        if self.plot_config.ladderize:
            if not copied_tree:
                tree = self._fast_copy(tree)
            tree.ladderize()

        # Fit Q matrix
        Q, loglik = self._fit_q_matrix(tree, tip_states, states, self.model)
        k = len(states)
        pi = np.ones(k) / k

        # Compute conditional likelihoods once
        cond_liks, _ = self._felsenstein_pruning(
            tree, tip_states, Q, pi, states
        )

        parent_map = self._build_parent_map(tree)
        simulation_metadata = self._build_simulation_metadata(
            tree, tip_states, states, parent_map
        )

        # Run stochastic mappings
        rng = np.random.default_rng(self.seed)
        mappings = []
        transition_cache = {}
        for _ in range(self.nsim):
            mapping = self._run_single_simulation(
                tree, tip_states, Q, pi, states, rng,
                cond_liks=cond_liks, parent_map=parent_map,
                transition_cache=transition_cache,
                simulation_metadata=simulation_metadata,
            )
            mappings.append(mapping)

        # Per-branch summary
        branch_summary = self._summarize_per_branch(
            mappings, states, tree, parent_map
        )
        node_posteriors = self._summarize_node_posteriors(
            mappings, states, tree
        )
        tree_summary = self._summarize_simulations(mappings, states, tree)

        # Output
        if self.csv_output:
            self._write_csv(
                branch_summary, node_posteriors, states, tree,
                parent_map, self.csv_output,
            )

        if self.plot_output:
            self._plot_posterior_pie(
                tree, node_posteriors, states, branch_summary,
                parent_map, self.plot_output,
            )

        if self.json_output:
            self._print_json_output(
                Q, loglik, states, branch_summary,
                node_posteriors, tree_summary, tree, parent_map,
            )
        else:
            self._print_text(
                Q, loglik, states, branch_summary,
                node_posteriors, tree_summary, tree, parent_map,
            )

    def _summarize_per_branch(
        self, mappings: list[dict], states: list[str],
        tree, parent_map: dict,
    ) -> dict[int, dict]:
        """Per-branch: proportion of time in each state, expected transitions."""
        k = len(states)
        nsim = len(mappings)

        branch_clades = self._branch_clades_preorder(tree)
        zeros = np.zeros
        branch_data = {}
        for clade in branch_clades:
            cid = id(clade)
            branch_data[cid] = [
                clade.branch_length if clade.branch_length else 1e-8,
                zeros(k),
                0.0,
            ]

        for mapping in mappings:
            histories = mapping["branch_histories"]
            for cid, history in histories.items():
                data = branch_data.get(cid)
                if data is None:
                    continue
                t = data[0]
                dwelling = data[1]

                if len(history) == 1:
                    dwelling[history[0][1]] += t
                    continue

                transitions = 0
                prev_time, prev_state = history[0]
                for h_idx in range(1, len(history)):
                    time_val, state = history[h_idx]
                    dwelling[prev_state] += time_val - prev_time
                    if prev_state != state:
                        transitions += 1
                    prev_time = time_val
                    prev_state = state
                dwelling[prev_state] += t - prev_time
                data[2] += transitions

        # Compute means and proportions
        result = {}
        for clade in branch_clades:
            cid = id(clade)
            t, dwelling, transitions = branch_data[cid]
            mean_dwell = dwelling / nsim
            total_dwell = mean_dwell.sum()
            proportions = (
                mean_dwell / total_dwell if total_dwell > 0
                else np.zeros(k)
            )
            mean_trans = float(transitions / nsim)

            result[cid] = {
                "mean_dwelling": mean_dwell,
                "proportions": proportions,
                "mean_transitions": mean_trans,
                "branch_length": t,
            }

        return result

    def _summarize_node_posteriors(
        self, mappings: list[dict], states: list[str], tree,
    ) -> dict[int, np.ndarray]:
        """Posterior probability of each state at each internal node."""
        k = len(states)
        nsim = len(mappings)
        posteriors = {}

        try:
            root = tree.root
            root.clades
        except AttributeError:
            root = None

        if root is not None:
            try:
                stack = [root]
                pop = stack.pop
                append = stack.append
                zeros = np.zeros
                id_ = id
                while stack:
                    clade = pop()
                    children = clade.clades
                    if children:
                        posteriors[id_(clade)] = zeros(k)
                        child_count = len(children)
                        if child_count == 2:
                            append(children[1])
                            append(children[0])
                        else:
                            for index in range(child_count - 1, -1, -1):
                                append(children[index])
            except AttributeError:
                posteriors = {}

        if root is None or not posteriors:
            for clade in tree.find_clades(order="preorder"):
                if clade.is_terminal():
                    continue
                cid = id(clade)
                posteriors[cid] = np.zeros(k)

        posterior_items = tuple(posteriors.items())
        for mapping in mappings:
            node_states = mapping["node_states"]
            get_state = node_states.get
            for cid, posterior in posterior_items:
                state = get_state(cid)
                if state is not None:
                    posterior[state] += 1

        inv_nsim = 1.0 / nsim
        for posterior in posteriors.values():
            posterior *= inv_nsim

        return posteriors

    @staticmethod
    def _collect_clade_tip_names(tree) -> dict[int, tuple]:
        direct_tip_names = SimmapSummary._collect_clade_tip_names_direct(tree)
        if direct_tip_names is not None:
            return direct_tip_names

        clade_tip_names: dict[int, tuple] = {}
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                clade_tip_names[id(clade)] = (clade.name,)
            else:
                tips = []
                for child in clade.clades:
                    tips.extend(clade_tip_names.get(id(child), ()))
                clade_tip_names[id(clade)] = tuple(sorted(tips))
        return clade_tip_names

    @staticmethod
    def _collect_clade_tip_names_direct(tree):
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

        clade_tip_names: dict[int, tuple] = {}
        for clade in reversed(clades):
            children = clade.clades
            if not children:
                clade_tip_names[id(clade)] = (clade.name,)
            else:
                tips = []
                for child in children:
                    tips.extend(clade_tip_names.get(id(child), ()))
                clade_tip_names[id(clade)] = tuple(sorted(tips))
        return clade_tip_names

    @staticmethod
    def _branch_label(clade, tree, parent_map, clade_tip_names=None) -> str:
        """Create a human-readable label for a branch."""
        if clade.is_terminal():
            return clade.name
        if clade_tip_names is None:
            tips = sorted(t.name for t in clade.get_terminals())
        else:
            tips = clade_tip_names.get(id(clade), ())
        if len(tips) <= 3:
            return f"({', '.join(tips)})"
        return f"({tips[0]}, ..., {tips[-1]}; {len(tips)} taxa)"

    def _print_text(
        self, Q, loglik, states, branch_summary,
        node_posteriors, tree_summary, tree, parent_map,
    ) -> None:
        clade_tip_names = self._collect_clade_tip_names(tree)
        try:
            preorder = list(self._iter_preorder(tree.root))
        except AttributeError:
            preorder = list(tree.find_clades(order="preorder"))
        lines = [
            f"SIMMAP Summary ({self.nsim} stochastic maps)",
            f"Model: {self.model}",
            f"Log-likelihood: {loglik:.4f}",
            f"States: {', '.join(states)}",
        ]
        try:
            # Q matrix
            lines.append("\nQ matrix:")
            header = "              " + "".join(f"{s:>12s}" for s in states)
            lines.append(header)
            q_row_format = "%-14s" + "%12.4f" * len(states)
            for s, q_row in zip(states, Q):
                lines.append(q_row_format % ((s,) + tuple(q_row)))

            # Tree-level summary
            mean_trans = tree_summary["mean_transitions"]
            total_trans = 0.0
            lines.append("\nExpected transitions (tree-wide):")
            for si, trans_row in zip(states, mean_trans):
                for sj, value in zip(states, trans_row):
                    if si != sj and value > 0.001:
                        lines.append("  %s -> %s:  %.2f" % (si, sj, value))
                        total_trans += value
            lines.append("  %-28s%.2f" % ("Total:", total_trans))

            # Per-branch summary
            lines.append("\nPer-branch dwelling time proportions:")
            col_header = (
                f"{'Branch':<35s} {'Length':>8s}"
                + "".join(f" {s:>10s}" for s in states)
                + f" {'E[trans]':>10s}"
            )
            lines.append(col_header)
            lines.append("-" * len(col_header))

            branch_rows = []
            node_rows = []
            branch_row_format = (
                "%-35s %8.4f"
                + " %10.4f" * len(states)
                + " %10.2f"
            )
            node_row_format = "%-35s" + " %10.4f" * len(states)
            for clade in preorder:
                cid = id(clade)

                if clade != tree.root and cid in branch_summary:
                    bs = branch_summary[cid]
                    label = self._branch_label(
                        clade, tree, parent_map, clade_tip_names
                    )
                    if len(label) > 33:
                        label = label[:30] + "..."
                    row = branch_row_format % (
                        (label, bs["branch_length"])
                        + tuple(bs["proportions"])
                        + (bs["mean_transitions"],)
                    )
                    branch_rows.append(row)

                if clade.clades and cid in node_posteriors:
                    label = self._branch_label(
                        clade, tree, parent_map, clade_tip_names
                    )
                    if len(label) > 33:
                        label = label[:30] + "..."
                    row = node_row_format % (
                        (label,) + tuple(node_posteriors[cid])
                    )
                    node_rows.append(row)

            lines.extend(branch_rows)

            # Node posteriors
            lines.append("\nNode posterior probabilities:")
            node_header = (
                f"{'Node':<35s}"
                + "".join(f" {s:>10s}" for s in states)
            )
            lines.append(node_header)
            lines.append("-" * len(node_header))
            lines.extend(node_rows)

            print("\n".join(lines))

        except BrokenPipeError:
            return

    def _print_json_output(
        self, Q, loglik, states, branch_summary,
        node_posteriors, tree_summary, tree, parent_map,
    ) -> None:
        k = len(states)
        clade_tip_names = self._collect_clade_tip_names(tree)

        # Q matrix
        q_dict = {}
        for si, q_row in zip(states, Q):
            q_dict[si] = {
                sj: round(float(value), 6) for sj, value in zip(states, q_row)
            }

        # Tree-level transitions
        mean_trans = tree_summary["mean_transitions"]
        trans_dict = {}
        total_trans = 0.0
        for i in range(k):
            for j in range(k):
                if i != j:
                    key = f"{states[i]} -> {states[j]}"
                    trans_dict[key] = round(float(mean_trans[i, j]), 4)
                    total_trans += float(mean_trans[i, j])

        # Per-branch
        branches = []
        nodes = []
        try:
            preorder_clades = list(self._iter_preorder(tree.root))
        except AttributeError:
            preorder_clades = list(tree.find_clades(order="preorder"))

        for clade in preorder_clades:
            cid = id(clade)

            if clade != tree.root and cid in branch_summary:
                bs = branch_summary[cid]
                label = self._branch_label(clade, tree, parent_map, clade_tip_names)
                branches.append(
                    {
                        "branch": label,
                        "branch_length": round(float(bs["branch_length"]), 6),
                        "dwelling_proportions": {
                            states[i]: round(float(bs["proportions"][i]), 4)
                            for i in range(k)
                        },
                        "expected_transitions": round(bs["mean_transitions"], 4),
                    }
                )

            if clade.clades and cid in node_posteriors:
                label = self._branch_label(clade, tree, parent_map, clade_tip_names)
                nodes.append(
                    {
                        "node": label,
                        "posteriors": {
                            states[i]: round(float(node_posteriors[cid][i]), 4)
                            for i in range(k)
                        },
                    }
                )

        payload = {
            "model": self.model,
            "nsim": self.nsim,
            "log_likelihood": round(float(loglik), 4),
            "states": states,
            "q_matrix": q_dict,
            "expected_transitions": trans_dict,
            "total_expected_transitions": round(total_trans, 4),
            "branches": branches,
            "node_posteriors": nodes,
        }
        print_json(payload, sort_keys=False)

    def _write_csv(
        self, branch_summary, node_posteriors, states, tree,
        parent_map, csv_path,
    ) -> None:
        k = len(states)
        clade_tip_names = self._collect_clade_tip_names(tree)
        try:
            preorder_clades = list(self._iter_preorder(tree.root))
        except AttributeError:
            preorder_clades = list(tree.find_clades(order="preorder"))

        with open(csv_path, "w") as f:
            # Branch dwelling proportions
            header = "branch,branch_length," + ",".join(
                f"prop_{s}" for s in states
            ) + ",expected_transitions\n"
            f.write(header)

            node_rows = []
            for clade in preorder_clades:
                cid = id(clade)

                if clade != tree.root and cid in branch_summary:
                    bs = branch_summary[cid]
                    label = self._branch_label(
                        clade, tree, parent_map, clade_tip_names
                    )
                    label = label.replace(",", ";")
                    proportions = ",".join(
                        f"{bs['proportions'][i]:.4f}" for i in range(k)
                    )
                    f.write(
                        f"{label},{bs['branch_length']:.6f},"
                        f"{proportions},{bs['mean_transitions']:.4f}\n"
                    )

                if clade.clades and cid in node_posteriors:
                    label = self._branch_label(
                        clade, tree, parent_map, clade_tip_names
                    )
                    label = label.replace(",", ";")
                    posteriors = ",".join(
                        f"{node_posteriors[cid][i]:.4f}" for i in range(k)
                    )
                    node_rows.append(f"{label},{posteriors}\n")

            # Node posteriors
            f.write("\n")
            f.write("node," + ",".join(
                f"posterior_{s}" for s in states
            ) + "\n")
            f.writelines(node_rows)

    def _plot_posterior_pie(
        self, tree, node_posteriors, states, branch_summary,
        parent_map, output_path,
    ) -> None:
        """Plot tree with pie charts at internal nodes showing posteriors."""
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.collections import LineCollection
            from matplotlib.lines import Line2D
        except ImportError:
            print("matplotlib is required for plotting.")
            raise SystemExit(2)

        from ...helpers.plot_config import (
            cleanup_tree_axes,
            compute_node_positions,
            draw_tip_labels,
        )

        k = len(states)
        cmap = plt.get_cmap("tab10")
        state_colors = [cmap(i / max(k - 1, 1)) for i in range(k)]

        preorder_clades = list(self._iter_preorder(tree.root))
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

        # Draw branches colored by dominant state proportion
        root = tree.root
        horizontal_segments = []
        horizontal_colors = []
        vertical_segments = []
        for clade in preorder_clades:
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

            # Color branch by dominant state
            if cid in branch_summary:
                props = branch_summary[cid]["proportions"]
                dom = int(np.argmax(props))
                color = state_colors[dom]
                alpha = float(max(props[dom], 0.4))
            else:
                color = "gray"
                alpha = 1.0

            horizontal_segments.append(((x0, y1), (x1, y1)))
            horizontal_colors.append(color if alpha == 1.0 else (*color[:3], alpha))
            vertical_segments.append(((x0, y0), (x0, y1)))

        if horizontal_segments:
            ax.add_collection(
                LineCollection(
                    horizontal_segments,
                    colors=horizontal_colors,
                    linewidths=2,
                    zorder=2,
                ),
                autolim=True,
            )
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
        if horizontal_segments or vertical_segments:
            ax.autoscale_view()

        # Tip labels
        label_fs = (
            config.ylabel_fontsize if config.ylabel_fontsize else 9
        )
        draw_tip_labels(ax, tree, node_x, node_y, fontsize=label_fs)

        # Pie charts at internal nodes
        n_tips = len(tips)
        pie_size = min(0.06, 0.8 / max(n_tips, 1))

        fig.subplots_adjust(left=0.05, right=0.85, top=0.92, bottom=0.12)
        fig.canvas.draw()

        for clade in preorder_clades:
            if not clade.clades or clade == root:
                continue
            cid = id(clade)
            if cid not in node_posteriors:
                continue

            probs = node_posteriors[cid]
            cx = node_x.get(cid, 0)
            cy = node_y.get(cid, 0)

            disp = ax.transData.transform((cx, cy))
            fig_coord = fig.transFigure.inverted().transform(disp)
            fx, fy = fig_coord

            inset = fig.add_axes(
                [fx - pie_size / 2, fy - pie_size / 2,
                 pie_size, pie_size],
                zorder=10,
            )
            wedge_vals = []
            wedge_colors = []
            for i in range(k):
                if probs[i] > 1e-6:
                    wedge_vals.append(probs[i])
                    wedge_colors.append(state_colors[i])
            if wedge_vals:
                inset.pie(
                    wedge_vals, colors=wedge_colors, startangle=90,
                    wedgeprops={"edgecolor": "black", "linewidth": 0.5},
                )
            inset.set_aspect("equal")
            inset.axis("off")

        # Legend
        handles = [
            Line2D([0], [0], color=state_colors[i], linewidth=3,
                   label=states[i])
            for i in range(k)
        ]
        legend_loc = config.legend_position or "upper left"
        if legend_loc != "none":
            ax.legend(
                handles=handles, title="States", loc=legend_loc,
                fontsize=8, title_fontsize=9,
            )

        cleanup_tree_axes(ax)
        if config.show_title:
            ax.set_title(
                config.title or "SIMMAP Summary — Node Posteriors",
                fontsize=config.title_fontsize,
            )

        fig.savefig(output_path, dpi=config.dpi)
        plt.close(fig)
        print(f"Plot saved: {output_path}")
