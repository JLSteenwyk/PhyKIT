"""
SIMMAP summary: run N stochastic character maps and provide per-branch
summaries of dwelling times, transition counts, and node posteriors.

Extends StochasticCharacterMap with detailed per-branch output matching
the output of phytools::describe.simmap() in R.
"""
import sys
from typing import Dict, List, Tuple

import numpy as np

from .stochastic_character_map import StochasticCharacterMap
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig, compute_node_x_cladogram
from ...errors import PhykitUserError


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

    def process_args(self, args) -> Dict:
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

    def run(self) -> None:
        tree = self.read_tree_file()
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

        # Prune tree to shared taxa
        tips_in_tree = set(self.get_tip_names_from_tree(tree))
        tips_to_prune = [t for t in tips_in_tree if t not in tip_states]
        if tips_to_prune:
            for tip_name in tips_to_prune:
                tree.prune(tip_name)

        if self.plot_config.ladderize:
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

        # Run stochastic mappings
        rng = np.random.default_rng(self.seed)
        mappings = []
        for _ in range(self.nsim):
            mapping = self._run_single_simulation(
                tree, tip_states, Q, pi, states, rng,
                cond_liks=cond_liks, parent_map=parent_map,
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
        self, mappings: List[Dict], states: List[str],
        tree, parent_map: Dict,
    ) -> Dict[int, Dict]:
        """Per-branch: proportion of time in each state, expected transitions."""
        k = len(states)
        nsim = len(mappings)

        # Initialize accumulators
        branch_dwelling = {}   # clade_id -> (nsim, k) array
        branch_transitions = {}  # clade_id -> (nsim,) array

        for clade in tree.find_clades(order="preorder"):
            if clade == tree.root:
                continue
            cid = id(clade)
            branch_dwelling[cid] = np.zeros((nsim, k))
            branch_transitions[cid] = np.zeros(nsim)

        for sim_idx, mapping in enumerate(mappings):
            histories = mapping["branch_histories"]
            for clade in tree.find_clades(order="preorder"):
                if clade == tree.root:
                    continue
                cid = id(clade)
                if cid not in histories:
                    continue

                history = histories[cid]
                t = clade.branch_length if clade.branch_length else 1e-8

                # Dwelling times
                for h_idx in range(len(history)):
                    state = history[h_idx][1]
                    start_t = history[h_idx][0]
                    end_t = (
                        history[h_idx + 1][0]
                        if h_idx + 1 < len(history)
                        else t
                    )
                    branch_dwelling[cid][sim_idx, state] += end_t - start_t

                # Count transitions on this branch
                n_trans = sum(
                    1 for h_idx in range(1, len(history))
                    if history[h_idx][1] != history[h_idx - 1][1]
                )
                branch_transitions[cid][sim_idx] = n_trans

        # Compute means and proportions
        result = {}
        for clade in tree.find_clades(order="preorder"):
            if clade == tree.root:
                continue
            cid = id(clade)
            t = clade.branch_length if clade.branch_length else 1e-8
            mean_dwell = branch_dwelling[cid].mean(axis=0)
            total_dwell = mean_dwell.sum()
            proportions = (
                mean_dwell / total_dwell if total_dwell > 0
                else np.zeros(k)
            )
            mean_trans = float(branch_transitions[cid].mean())

            result[cid] = {
                "mean_dwelling": mean_dwell,
                "proportions": proportions,
                "mean_transitions": mean_trans,
                "branch_length": t,
            }

        return result

    def _summarize_node_posteriors(
        self, mappings: List[Dict], states: List[str], tree,
    ) -> Dict[int, np.ndarray]:
        """Posterior probability of each state at each internal node."""
        k = len(states)
        nsim = len(mappings)
        posteriors = {}

        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            cid = id(clade)
            posteriors[cid] = np.zeros(k)

        for mapping in mappings:
            node_states = mapping["node_states"]
            for cid in posteriors:
                if cid in node_states:
                    posteriors[cid][node_states[cid]] += 1

        for cid in posteriors:
            posteriors[cid] /= nsim

        return posteriors

    @staticmethod
    def _branch_label(clade, tree, parent_map) -> str:
        """Create a human-readable label for a branch."""
        if clade.is_terminal():
            return clade.name
        tips = sorted(t.name for t in clade.get_terminals())
        if len(tips) <= 3:
            return f"({', '.join(tips)})"
        return f"({tips[0]}, ..., {tips[-1]}; {len(tips)} taxa)"

    def _print_text(
        self, Q, loglik, states, branch_summary,
        node_posteriors, tree_summary, tree, parent_map,
    ) -> None:
        k = len(states)
        try:
            print(f"SIMMAP Summary ({self.nsim} stochastic maps)")
            print(f"Model: {self.model}")
            print(f"Log-likelihood: {loglik:.4f}")
            print(f"States: {', '.join(states)}")

            # Q matrix
            print("\nQ matrix:")
            header = "              " + "".join(f"{s:>12s}" for s in states)
            print(header)
            for i, s in enumerate(states):
                row = f"{s:14s}" + "".join(
                    f"{Q[i, j]:12.4f}" for j in range(k)
                )
                print(row)

            # Tree-level summary
            mean_trans = tree_summary["mean_transitions"]
            total_trans = 0.0
            print("\nExpected transitions (tree-wide):")
            for i in range(k):
                for j in range(k):
                    if i != j and mean_trans[i, j] > 0.001:
                        print(
                            f"  {states[i]} -> {states[j]}:"
                            f"  {mean_trans[i, j]:.2f}"
                        )
                        total_trans += mean_trans[i, j]
            print(f"  {'Total:':28s}{total_trans:.2f}")

            # Per-branch summary
            print(f"\nPer-branch dwelling time proportions:")
            col_header = (
                f"{'Branch':<35s} {'Length':>8s}"
                + "".join(f" {s:>10s}" for s in states)
                + f" {'E[trans]':>10s}"
            )
            print(col_header)
            print("-" * len(col_header))

            for clade in tree.find_clades(order="preorder"):
                if clade == tree.root:
                    continue
                cid = id(clade)
                if cid not in branch_summary:
                    continue
                bs = branch_summary[cid]
                label = self._branch_label(clade, tree, parent_map)
                if len(label) > 33:
                    label = label[:30] + "..."
                row = f"{label:<35s} {bs['branch_length']:>8.4f}"
                for i in range(k):
                    row += f" {bs['proportions'][i]:>10.4f}"
                row += f" {bs['mean_transitions']:>10.2f}"
                print(row)

            # Node posteriors
            print(f"\nNode posterior probabilities:")
            node_header = (
                f"{'Node':<35s}"
                + "".join(f" {s:>10s}" for s in states)
            )
            print(node_header)
            print("-" * len(node_header))

            for clade in tree.find_clades(order="preorder"):
                if clade.is_terminal():
                    continue
                cid = id(clade)
                if cid not in node_posteriors:
                    continue
                label = self._branch_label(clade, tree, parent_map)
                if len(label) > 33:
                    label = label[:30] + "..."
                row = f"{label:<35s}"
                for i in range(k):
                    row += f" {node_posteriors[cid][i]:>10.4f}"
                print(row)

        except BrokenPipeError:
            return

    def _print_json_output(
        self, Q, loglik, states, branch_summary,
        node_posteriors, tree_summary, tree, parent_map,
    ) -> None:
        k = len(states)

        # Q matrix
        q_dict = {}
        for i, si in enumerate(states):
            q_dict[si] = {
                states[j]: round(float(Q[i, j]), 6) for j in range(k)
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
        for clade in tree.find_clades(order="preorder"):
            if clade == tree.root:
                continue
            cid = id(clade)
            if cid not in branch_summary:
                continue
            bs = branch_summary[cid]
            label = self._branch_label(clade, tree, parent_map)
            entry = {
                "branch": label,
                "branch_length": round(float(bs["branch_length"]), 6),
                "dwelling_proportions": {
                    states[i]: round(float(bs["proportions"][i]), 4)
                    for i in range(k)
                },
                "expected_transitions": round(bs["mean_transitions"], 4),
            }
            branches.append(entry)

        # Node posteriors
        nodes = []
        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            cid = id(clade)
            if cid not in node_posteriors:
                continue
            label = self._branch_label(clade, tree, parent_map)
            entry = {
                "node": label,
                "posteriors": {
                    states[i]: round(float(node_posteriors[cid][i]), 4)
                    for i in range(k)
                },
            }
            nodes.append(entry)

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
        with open(csv_path, "w") as f:
            # Branch dwelling proportions
            header = "branch,branch_length," + ",".join(
                f"prop_{s}" for s in states
            ) + ",expected_transitions\n"
            f.write(header)
            for clade in tree.find_clades(order="preorder"):
                if clade == tree.root:
                    continue
                cid = id(clade)
                if cid not in branch_summary:
                    continue
                bs = branch_summary[cid]
                label = self._branch_label(clade, tree, parent_map)
                label = label.replace(",", ";")
                row = f"{label},{bs['branch_length']:.6f}"
                for i in range(k):
                    row += f",{bs['proportions'][i]:.4f}"
                row += f",{bs['mean_transitions']:.4f}\n"
                f.write(row)

            # Node posteriors
            f.write("\n")
            f.write("node," + ",".join(
                f"posterior_{s}" for s in states
            ) + "\n")
            for clade in tree.find_clades(order="preorder"):
                if clade.is_terminal():
                    continue
                cid = id(clade)
                if cid not in node_posteriors:
                    continue
                label = self._branch_label(clade, tree, parent_map)
                label = label.replace(",", ";")
                row = label
                for i in range(k):
                    row += f",{node_posteriors[cid][i]:.4f}"
                row += "\n"
                f.write(row)

    def _plot_posterior_pie(
        self, tree, node_posteriors, states, branch_summary,
        parent_map, output_path,
    ) -> None:
        """Plot tree with pie charts at internal nodes showing posteriors."""
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.lines import Line2D
        except ImportError:
            print("matplotlib is required for plotting.")
            raise SystemExit(2)

        k = len(states)
        cmap = plt.get_cmap("tab10")
        state_colors = [cmap(i / max(k - 1, 1)) for i in range(k)]

        tips = list(tree.get_terminals())
        node_x = {}
        node_y = {}

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
                    parent = parent_map.get(id(clade))
                    if parent is not None:
                        t = clade.branch_length if clade.branch_length else 0.0
                        node_x[id(clade)] = node_x[id(parent)] + t

        for clade in tree.find_clades(order="postorder"):
            if not clade.is_terminal() and id(clade) not in node_y:
                child_ys = [
                    node_y[id(c)] for c in clade.clades
                    if id(c) in node_y
                ]
                node_y[id(clade)] = (
                    float(np.mean(child_ys)) if child_ys else 0.0
                )

        config = self.plot_config
        config.resolve(n_rows=len(tips), n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        # Draw branches colored by dominant state proportion
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

            # Color branch by dominant state
            if cid in branch_summary:
                props = branch_summary[cid]["proportions"]
                dom = int(np.argmax(props))
                color = state_colors[dom]
                alpha = float(max(props[dom], 0.4))
            else:
                color = "gray"
                alpha = 1.0

            ax.plot([x0, x1], [y1, y1], color=color, lw=2, alpha=alpha)
            ax.plot([x0, x0], [y0, y1], color="gray", lw=0.8)

        # Tip labels
        max_x = max(node_x.values()) if node_x else 1.0
        offset = max_x * 0.02
        label_fs = (
            config.ylabel_fontsize if config.ylabel_fontsize else 9
        )
        if label_fs > 0:
            for tip in tips:
                ax.text(
                    node_x[id(tip)] + offset, node_y[id(tip)],
                    tip.name, va="center", fontsize=label_fs,
                )

        # Pie charts at internal nodes
        n_tips = len(tips)
        pie_size = min(0.06, 0.8 / max(n_tips, 1))

        fig.subplots_adjust(left=0.05, right=0.85, top=0.92, bottom=0.12)
        fig.canvas.draw()

        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal() or clade == root:
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

        ax.set_xlabel("Branch length")
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        if config.show_title:
            ax.set_title(
                config.title or "SIMMAP Summary — Node Posteriors",
                fontsize=config.title_fontsize,
            )

        fig.savefig(output_path, dpi=config.dpi)
        plt.close(fig)
        print(f"Plot saved: {output_path}")
