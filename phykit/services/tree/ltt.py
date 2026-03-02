import math
from typing import Dict, List, Tuple

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class LTT(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.verbose = parsed["verbose"]
        self.json_output = parsed["json_output"]
        self.plot_output = parsed["plot_output"]

    def run(self) -> None:
        tree = self.read_tree_file()
        self._validate_tree(tree)

        gamma, p_value, bt, g = self._compute_gamma(tree)
        ltt_data = self._compute_ltt(tree)

        if self.json_output:
            self._output_json(gamma, p_value, ltt_data, bt, g)
            return

        self._output_text(gamma, p_value, ltt_data, bt, g)

        if self.plot_output:
            self._plot_ltt(ltt_data, self.plot_output, gamma=gamma, p_value=p_value)

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            verbose=getattr(args, "verbose", False),
            json_output=getattr(args, "json", False),
            plot_output=getattr(args, "plot_output", None),
        )

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for gamma statistic."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                raise PhykitUserError(
                    ["All branches in the tree must have lengths."],
                    code=2,
                )

    @staticmethod
    def _compute_gamma(tree):
        """Pybus & Harvey (2000) gamma statistic, matching ape::gammaStat().

        Replicates the exact algorithm from R's ape source (gammaStat.R):
            N <- length(phy$tip.label)
            bt <- sort(branching.times(phy))
            g <- rev(c(bt[1], diff(bt)))
            ST <- sum((2:N) * g)
            stat <- sum(cumsum((2:(N-1)) * g[-(N-1)])) / (N-2)
            m <- ST / 2
            s <- ST * sqrt(1 / (12 * (N - 2)))
            (stat - m) / s
        """
        tips = list(tree.get_terminals())
        N = len(tips)

        if N < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for gamma statistic."],
                code=2,
            )

        # Get branching times (node ages = distance from present/tips)
        root = tree.root
        max_height = max(tree.distance(root, tip) for tip in tips)

        bt = []
        for clade in tree.find_clades(order="level"):
            if clade.is_terminal():
                continue
            node_dist_from_root = tree.distance(root, clade)
            node_age = max_height - node_dist_from_root
            bt.append(node_age)

        bt.sort()  # ascending: most recent nodes first

        # Internode intervals reversed (ape: g <- rev(c(bt[1], diff(bt))))
        # This orders intervals from past to present (root interval first)
        g_unreversed = [bt[0]] + [bt[i] - bt[i - 1] for i in range(1, len(bt))]
        g = list(reversed(g_unreversed))

        # ST = sum((2:N) * g)
        ST = sum((k + 2) * g[k] for k in range(N - 1))

        # stat = sum(cumsum((2:(N-1)) * g[-(N-1)])) / (N-2)
        # g[-(N-1)] in R removes the last element
        g_partial = g[:-1]
        partial = [(k + 2) * g_partial[k] for k in range(N - 2)]
        cumsum_partial = []
        running = 0.0
        for val in partial:
            running += val
            cumsum_partial.append(running)

        stat = sum(cumsum_partial) / (N - 2)

        m = ST / 2
        s = ST * math.sqrt(1.0 / (12 * (N - 2)))
        gamma = (stat - m) / s

        # Two-tailed p-value from standard normal
        from scipy.stats import norm

        p_value = 2 * norm.sf(abs(gamma))

        return gamma, p_value, bt, g

    @staticmethod
    def _compute_ltt(tree):
        """Compute lineage-through-time data.

        Returns list of (time_from_root, n_lineages) tuples.
        Time runs from 0 (root) to tree_height (present).
        """
        # Get branching times as distances from root
        root = tree.root
        max_height = max(
            tree.distance(root, tip) for tip in tree.get_terminals()
        )

        branching_times_from_root = []
        for clade in tree.find_clades(order="level"):
            if clade.is_terminal():
                continue
            if clade == root:
                continue
            branching_times_from_root.append(tree.distance(root, clade))

        branching_times_from_root.sort()

        # LTT: start with 2 lineages at root (time=0)
        ltt = [(0.0, 2)]
        n_lineages = 2
        for bt in branching_times_from_root:
            n_lineages += 1
            ltt.append((bt, n_lineages))
        # End at present
        ltt.append((max_height, n_lineages))

        return ltt

    @staticmethod
    def _plot_ltt(ltt_data, output_path, gamma=None, p_value=None):
        """Plot lineage-through-time as a step function."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        times = [pt[0] for pt in ltt_data]
        lineages = [pt[1] for pt in ltt_data]

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.step(times, lineages, where="post", linewidth=2, color="black")
        ax.set_xlabel("Time from root", fontsize=12)
        ax.set_ylabel("Number of lineages", fontsize=12)
        ax.set_yscale("log")

        if gamma is not None:
            label = f"\u03b3 = {gamma:.4f}"
            if p_value is not None:
                label += f" (p = {p_value:.4e})"
            ax.text(
                0.05,
                0.95,
                label,
                transform=ax.transAxes,
                fontsize=11,
                verticalalignment="top",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
            )

        fig.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

    def _output_text(self, gamma, p_value, ltt_data, bt, g):
        try:
            print(f"{round(gamma, 4)}\t{round(p_value, 4)}")
            if self.verbose:
                print("\nBranching times (node ages):")
                for i, t in enumerate(bt):
                    print(f"  {i + 1}\t{t:.6f}")
                print("\nLineage-through-time:")
                print("  time_from_root\tn_lineages")
                for time_val, n_lin in ltt_data:
                    print(f"  {time_val:.6f}\t{n_lin}")
        except BrokenPipeError:
            pass

    def _output_json(self, gamma, p_value, ltt_data, bt, g):
        result = dict(
            gamma=float(gamma),
            p_value=float(p_value),
            branching_times=[float(t) for t in bt],
            internode_intervals=[float(v) for v in g],
            ltt=[
                dict(time_from_root=float(t), n_lineages=int(n))
                for t, n in ltt_data
            ],
        )
        print_json(result)
