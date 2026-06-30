from __future__ import annotations

import math

from .base import Tree
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def _max_tip_height(tips, depths, root_depth):
    iterator = iter(tips)
    first = next(iterator)
    local_depths = depths
    max_height = local_depths[first] - root_depth
    for tip in iterator:
        height = local_depths[tip] - root_depth
        if height > max_height:
            max_height = height
    return max_height


def _ltt_from_internal_depths(internal_depths, max_height):
    ltt = [(0.0, 2)]
    n_lineages = 2
    append = ltt.append
    iterator = iter(internal_depths)
    next(iterator, None)
    for bt_from_root in iterator:
        n_lineages += 1
        append((bt_from_root, n_lineages))
    append((max_height, n_lineages))
    return ltt


class LTT(Tree):
    _DEPTH_DATA_UNSET = object()

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.verbose = parsed["verbose"]
        self.json_output = parsed["json_output"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        self.validate_tree(tree, min_tips=3, context="gamma statistic")

        tips = self._terminal_clades(tree)
        depth_data = self._depths_from_root(tree)
        gamma, p_value, bt, g, ltt_data = self._compute_gamma_and_ltt(
            tree, tips=tips, depth_data=depth_data
        )

        if self.json_output:
            self._output_json(gamma, p_value, ltt_data, bt, g)
            return

        self._output_text(gamma, p_value, ltt_data, bt, g)

        if self.plot_output:
            self._plot_ltt(ltt_data, self.plot_output, gamma=gamma, p_value=p_value)

    def process_args(self, args) -> dict[str, str]:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree,
            verbose=getattr(args, "verbose", False),
            json_output=getattr(args, "json", False),
            plot_output=getattr(args, "plot_output", None),
            plot_config=PlotConfig.from_args(args),
        )

    @staticmethod
    def _terminal_clades(tree):
        terminals = LTT._terminal_clades_fast(tree)
        if terminals is not None:
            return terminals
        return list(tree.get_terminals())

    @staticmethod
    def _terminal_clades_fast(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        terminals = []
        stack = [root]
        pop = stack.pop
        append = stack.append
        append_terminal = terminals.append
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    if len(children) == 2:
                        left, right = children
                        append(right)
                        append(left)
                    else:
                        for index in range(len(children) - 1, -1, -1):
                            append(children[index])
                else:
                    append_terminal(clade)
        except AttributeError:
            return None
        return terminals

    @staticmethod
    def _depths_from_root(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            root = None
        else:
            depths = {root: 0.0}
            stack = [root]
            try:
                while stack:
                    clade = stack.pop()
                    parent_depth = depths[clade]
                    for child in reversed(clade.clades):
                        depths[child] = parent_depth + (child.branch_length or 0.0)
                        stack.append(child)
            except (AttributeError, TypeError):
                pass
            else:
                return root, depths, 0.0

        try:
            depths = tree.depths()
            root = tree.root
            root_depth = depths[root]
        except (AttributeError, KeyError):
            return None
        return root, depths, root_depth

    @staticmethod
    def _internal_depths_from_root(tree, depth_data, include_root):
        if depth_data is None:
            return None
        try:
            root, depths, root_depth = depth_data
            root.clades
        except (AttributeError, TypeError, ValueError):
            return None

        internal_depths = []
        stack = [root]
        pop = stack.pop
        append = stack.append
        append_depth = internal_depths.append
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    if include_root or clade != root:
                        append_depth(depths[clade] - root_depth)
                    child_count = len(children)
                    if child_count == 2:
                        append(children[1])
                        append(children[0])
                    else:
                        for index in range(child_count - 1, -1, -1):
                            append(children[index])
        except (AttributeError, KeyError):
            return None

        return internal_depths

    @staticmethod
    def _compute_gamma(tree, tips=None, depth_data=_DEPTH_DATA_UNSET):
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
        if tips is None:
            tips = LTT._terminal_clades(tree)
        N = len(tips)

        if N < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for gamma statistic."],
                code=2,
            )

        # Get branching times (node ages = distance from present/tips)
        if depth_data is LTT._DEPTH_DATA_UNSET:
            depth_data = LTT._depths_from_root(tree)
        if depth_data is None:
            root = tree.root
            max_height = max(tree.distance(root, tip) for tip in tips)
        else:
            root, depths, root_depth = depth_data
            max_height = _max_tip_height(tips, depths, root_depth)

        internal_depths = LTT._internal_depths_from_root(
            tree, depth_data, include_root=True
        )
        if internal_depths is None:
            bt = []
            for clade in tree.find_clades(order="level"):
                if clade.is_terminal():
                    continue
                if depth_data is None:
                    node_dist_from_root = tree.distance(root, clade)
                else:
                    node_dist_from_root = depths[clade] - root_depth
                node_age = max_height - node_dist_from_root
                bt.append(node_age)
        else:
            bt = [
                max_height - node_dist_from_root
                for node_dist_from_root in internal_depths
            ]

        bt.sort()  # ascending: most recent nodes first

        # Internode intervals reversed (ape: g <- rev(c(bt[1], diff(bt))))
        # This orders intervals from past to present (root interval first)
        g_unreversed = [bt[0]] + [bt[i] - bt[i - 1] for i in range(1, len(bt))]
        g = list(reversed(g_unreversed))

        # ST = sum((2:N) * g)
        ST = sum((k + 2) * g[k] for k in range(N - 1))

        # stat = sum(cumsum((2:(N-1)) * g[-(N-1)])) / (N-2)
        # g[-(N-1)] in R removes the last element
        running = 0.0
        stat_sum = 0.0
        for k in range(N - 2):
            running += (k + 2) * g[k]
            stat_sum += running

        stat = stat_sum / (N - 2)

        m = ST / 2
        s = ST * math.sqrt(1.0 / (12 * (N - 2)))
        gamma = (stat - m) / s

        # Two-tailed p-value from standard normal:
        # 2 * norm.sf(abs(gamma)) == erfc(abs(gamma) / sqrt(2)).
        p_value = math.erfc(abs(gamma) / math.sqrt(2.0))

        return gamma, p_value, bt, g

    @staticmethod
    def _compute_gamma_and_ltt(tree, tips=None, depth_data=_DEPTH_DATA_UNSET):
        if tips is None:
            tips = LTT._terminal_clades(tree)
        N = len(tips)

        if N < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for gamma statistic."],
                code=2,
            )

        depth_data = (
            LTT._depths_from_root(tree)
            if depth_data is LTT._DEPTH_DATA_UNSET
            else depth_data
        )
        if depth_data is None:
            gamma, p_value, bt, g = LTT._compute_gamma(
                tree, tips=tips, depth_data=None
            )
            ltt_data = LTT._compute_ltt(tree, tips=tips, depth_data=None)
            return gamma, p_value, bt, g, ltt_data

        root, depths, root_depth = depth_data
        max_height = _max_tip_height(tips, depths, root_depth)

        internal_depths = LTT._internal_depths_from_root(
            tree, depth_data, include_root=True
        )
        if internal_depths is None:
            gamma, p_value, bt, g = LTT._compute_gamma(
                tree, tips=tips, depth_data=depth_data
            )
            ltt_data = LTT._compute_ltt(tree, tips=tips, depth_data=depth_data)
            return gamma, p_value, bt, g, ltt_data

        internal_depths.sort()
        bt = [max_height - depth for depth in reversed(internal_depths)]

        # Internode intervals reversed (ape: g <- rev(c(bt[1], diff(bt)))).
        g_unreversed = [bt[0]] + [bt[i] - bt[i - 1] for i in range(1, len(bt))]
        g = list(reversed(g_unreversed))

        ST = sum((k + 2) * g[k] for k in range(N - 1))

        running = 0.0
        stat_sum = 0.0
        for k in range(N - 2):
            running += (k + 2) * g[k]
            stat_sum += running

        stat = stat_sum / (N - 2)
        m = ST / 2
        s = ST * math.sqrt(1.0 / (12 * (N - 2)))
        gamma = (stat - m) / s
        p_value = math.erfc(abs(gamma) / math.sqrt(2.0))

        ltt = _ltt_from_internal_depths(internal_depths, max_height)

        return gamma, p_value, bt, g, ltt

    @staticmethod
    def _compute_ltt(tree, tips=None, depth_data=_DEPTH_DATA_UNSET):
        """Compute lineage-through-time data.

        Returns list of (time_from_root, n_lineages) tuples.
        Time runs from 0 (root) to tree_height (present).
        """
        # Get branching times as distances from root
        if tips is None:
            tips = LTT._terminal_clades(tree)

        depth_data = (
            LTT._depths_from_root(tree)
            if depth_data is LTT._DEPTH_DATA_UNSET
            else depth_data
        )
        if depth_data is None:
            root = tree.root
            max_height = max(
                tree.distance(root, tip) for tip in tips
            )
        else:
            root, depths, root_depth = depth_data
            max_height = _max_tip_height(tips, depths, root_depth)

        branching_times_from_root = LTT._internal_depths_from_root(
            tree, depth_data, include_root=False
        )
        if branching_times_from_root is None:
            branching_times_from_root = []
            for clade in tree.find_clades(order="level"):
                if clade.is_terminal():
                    continue
                if clade == root:
                    continue
                if depth_data is None:
                    branching_times_from_root.append(tree.distance(root, clade))
                else:
                    branching_times_from_root.append(depths[clade] - root_depth)

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

    def _plot_ltt(self, ltt_data, output_path, gamma=None, p_value=None):
        """Plot lineage-through-time as a step function."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        config = self.plot_config
        config.resolve(n_rows=len(ltt_data), n_cols=None)
        default_colors = ["black"]
        colors = config.merge_colors(default_colors)

        times = [pt[0] for pt in ltt_data]
        lineages = [pt[1] for pt in ltt_data]

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        ax.step(times, lineages, where="post", linewidth=2, color=colors[0])
        ax.set_xlabel("Time from root", fontsize=config.axis_fontsize or 12)
        ax.set_ylabel("Number of lineages", fontsize=config.axis_fontsize or 12)
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

        if config.show_title and config.title:
            ax.set_title(config.title, fontsize=config.title_fontsize)

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def _output_text(self, gamma, p_value, ltt_data, bt, g):
        try:
            lines = [f"{round(gamma, 4)}\t{round(p_value, 4)}"]
            if self.verbose:
                lines.append("\nBranching times (node ages):")
                lines.extend(f"  {i + 1}\t{t:.6f}" for i, t in enumerate(bt))
                lines.append("\nLineage-through-time:")
                lines.append("  time_from_root\tn_lineages")
                lines.extend(f"  {time_val:.6f}\t{n_lin}" for time_val, n_lin in ltt_data)
            print("\n".join(lines))
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
