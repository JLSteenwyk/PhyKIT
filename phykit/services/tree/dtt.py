"""
Disparity through time (DTT).

Calculates and plots disparity-through-time for a phylogenetic tree
and phenotypic data, following Harmon et al. (2003). Computes the
Morphological Disparity Index (MDI) comparing observed DTT to a
Brownian motion null expectation.
"""
import sys
from typing import Dict, List, Tuple

import numpy as np

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig
from ...helpers.trait_parsing import parse_multi_trait_file
from ...errors import PhykitUserError


class Dtt(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.trait_column = parsed["trait_column"]
        self.index = parsed["index"]
        self.nsim = parsed["nsim"]
        self.seed = parsed["seed"]
        self.plot_output = parsed["plot_output"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.traits,
            trait_column=getattr(args, "trait", None),
            index=getattr(args, "index", "avg_sq"),
            nsim=getattr(args, "nsim", 0),
            seed=getattr(args, "seed", None),
            plot_output=getattr(args, "plot_output", None),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        tree = self.read_tree_file()
        self.validate_tree(
            tree, min_tips=3, require_branch_lengths=True,
            context="disparity through time",
        )

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(
            self.trait_data_path, tree_tips
        )

        # Get trait data
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)

        if self.trait_column:
            if self.trait_column not in trait_names:
                raise PhykitUserError(
                    [
                        f"Trait column '{self.trait_column}' not found.",
                        f"Available: {', '.join(trait_names)}",
                    ],
                    code=2,
                )
            col_idx = trait_names.index(self.trait_column)
            data = np.array([traits[name][col_idx] for name in ordered_names])
            data = data.reshape(-1, 1)
        else:
            p = len(trait_names)
            data = np.array([
                [traits[name][j] for j in range(p)]
                for name in ordered_names
            ])

        # Prune tree to trait taxa
        import pickle
        tree_copy = pickle.loads(pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL))
        tips_in_tree = [t.name for t in tree_copy.get_terminals()]
        tips_to_prune = [t for t in tips_in_tree if t not in traits]
        if tips_to_prune:
            tree_copy = self.prune_tree_using_taxa_list(tree_copy, tips_to_prune)

        # Compute DTT
        times, dtt_values = self._compute_dtt(tree_copy, data, ordered_names)

        # Compute MDI and null envelope via simulation
        mdi = None
        mdi_p = None
        sim_dtt = None
        if self.nsim > 0:
            sim_dtt, mdi, mdi_p = self._simulate_null(
                tree_copy, data, ordered_names, times
            )

        # Output
        if self.plot_output:
            self._plot_dtt(times, dtt_values, sim_dtt, mdi)

        if self.json_output:
            self._print_json(times, dtt_values, mdi, mdi_p, sim_dtt)
        else:
            self._print_text(times, dtt_values, mdi, mdi_p)

    def _compute_disparity(self, data: np.ndarray) -> float:
        """Compute disparity (average squared Euclidean distance among all pairs)."""
        n = data.shape[0]
        if n < 2:
            return 0.0

        if self.index == "avg_manhattan":
            total = 0.0
            count = 0
            for i in range(n):
                for j in range(i + 1, n):
                    total += np.sum(np.abs(data[i] - data[j]))
                    count += 1
            return total / count if count > 0 else 0.0
        else:
            # avg_sq: average squared Euclidean distance
            total = 0.0
            count = 0
            for i in range(n):
                for j in range(i + 1, n):
                    total += np.sum((data[i] - data[j]) ** 2)
                    count += 1
            return total / count if count > 0 else 0.0

    def _compute_dtt(
        self, tree, data: np.ndarray, ordered_names: List[str]
    ) -> Tuple[List[float], List[float]]:
        """Compute DTT curve following Harmon et al. (2003).

        At each branching time, compute the mean relative disparity of
        subclades with >= 2 species. Times are relative (0 = root, 1 = tips).
        """
        root = tree.root
        tip_heights = [tree.distance(root, t) for t in tree.get_terminals()]
        tree_height = max(tip_heights)
        if tree_height == 0:
            return [0.0, 1.0], [1.0, 0.0]

        total_disp = self._compute_disparity(data)
        if total_disp == 0:
            return [0.0, 1.0], [1.0, 0.0]

        name_to_idx = {name: i for i, name in enumerate(ordered_names)}

        # Precompute disparity for every clade (internal node)
        clade_disp = {}
        for clade in tree.find_clades(order="preorder"):
            tips = [t.name for t in clade.get_terminals()
                    if t.name in name_to_idx]
            if len(tips) >= 2:
                indices = [name_to_idx[n] for n in tips]
                clade_disp[id(clade)] = self._compute_disparity(data[indices])

        # Precompute node times (distance from root, relative)
        node_time = {}
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            node_time[id(clade)] = tree.distance(root, clade) / tree_height
            for child in clade.clades:
                parent_map[id(child)] = clade

        # Branching times (all internal nodes including root, sorted)
        internal_nodes = [
            c for c in tree.find_clades(order="preorder")
            if not c.is_terminal()
        ]
        internal_nodes.sort(key=lambda c: node_time[id(c)])

        times = [0.0]
        dtt_values = [1.0]

        for node in internal_nodes:
            t = node_time[id(node)]

            # At time t (just after this node branches), find all
            # lineages: branches where parent_time <= t < child_tip_time.
            # Each lineage corresponds to the clade below the child end
            # of that branch.
            lineage_disparities = []

            for clade in tree.find_clades(order="preorder"):
                if clade == root:
                    continue
                parent = parent_map.get(id(clade))
                if parent is None:
                    continue

                parent_t = node_time[id(parent)]
                # For a terminal, its "node time" = tree height = 1.0
                if clade.is_terminal():
                    child_t = 1.0
                else:
                    child_t = node_time[id(clade)]

                # This branch spans time t if parent_t <= t < child_t
                if parent_t <= t + 1e-10 and child_t > t + 1e-10:
                    # The subclade below this branch
                    tips = [tip.name for tip in clade.get_terminals()
                            if tip.name in name_to_idx]
                    if len(tips) >= 2:
                        lineage_disparities.append(
                            clade_disp.get(id(clade), 0.0)
                        )

            if lineage_disparities:
                rel_disp = float(np.mean(lineage_disparities)) / total_disp
            else:
                rel_disp = 0.0

            times.append(t)
            dtt_values.append(rel_disp)

        # Terminal point (only if last entry isn't already at/near 0)
        if dtt_values and dtt_values[-1] > 1e-10:
            times.append(1.0)
            dtt_values.append(0.0)

        return times, dtt_values

    def _compute_mdi(
        self, times: List[float], dtt_values: List[float],
        null_median: np.ndarray,
    ) -> float:
        """Compute MDI: area between observed DTT and null median."""
        times_arr = np.array(times)
        dtt_arr = np.array(dtt_values)

        # Interpolate null median to match observed times
        null_times = np.linspace(0, 1, len(null_median))
        null_interp = np.interp(times_arr, null_times, null_median)

        # MDI = area between observed and null (trapezoidal integration)
        diff = dtt_arr - null_interp
        mdi = float(np.trapezoid(diff, times_arr))
        return mdi

    def _simulate_null(
        self, tree, data, ordered_names, obs_times,
    ) -> Tuple[np.ndarray, float, float]:
        """Simulate BM data and compute null DTT distribution."""
        from .vcv_utils import build_vcv_matrix

        rng = np.random.default_rng(self.seed)
        n = len(ordered_names)
        p = data.shape[1] if data.ndim > 1 else 1

        vcv = build_vcv_matrix(tree, ordered_names)

        # Cholesky for simulation
        try:
            L = np.linalg.cholesky(vcv)
        except np.linalg.LinAlgError:
            # Add small diagonal for PD
            vcv += np.eye(n) * 1e-8
            L = np.linalg.cholesky(vcv)

        # Compute total disparity and mean
        total_disp = self._compute_disparity(data)

        # Simulate nsim datasets
        n_time_points = len(obs_times)
        sim_dtt_matrix = np.zeros((self.nsim, n_time_points))

        for s in range(self.nsim):
            # Simulate BM data with same VCV structure
            if p == 1:
                z = rng.standard_normal(n)
                sim_data = (L @ z).reshape(-1, 1)
            else:
                Z = rng.standard_normal((n, p))
                sim_data = L @ Z

            sim_times, sim_values = self._compute_dtt(
                tree, sim_data, ordered_names
            )

            # Interpolate to observed time points
            sim_interp = np.interp(obs_times, sim_times, sim_values)
            sim_dtt_matrix[s, :] = sim_interp

        # Null median
        null_median = np.median(sim_dtt_matrix, axis=0)

        # MDI
        obs_arr = np.array(
            self._compute_dtt(tree, data, ordered_names)[1]
        )
        # Interpolate observed to consistent times
        obs_times_arr = np.array(obs_times)
        mdi = float(np.trapezoid(
            np.interp(obs_times_arr, obs_times, obs_arr) - null_median,
            obs_times_arr,
        ))

        # MDI p-value: proportion of simulated MDIs >= observed
        sim_mdis = np.zeros(self.nsim)
        for s in range(self.nsim):
            sim_mdi = float(np.trapezoid(
                sim_dtt_matrix[s, :] - null_median, obs_times_arr
            ))
            sim_mdis[s] = sim_mdi

        mdi_p = float(np.mean(np.abs(sim_mdis) >= np.abs(mdi)))

        return sim_dtt_matrix, mdi, mdi_p

    def _plot_dtt(self, times, dtt_values, sim_dtt, mdi):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for DTT plotting.")
            return

        config = self.plot_config
        config.resolve(n_rows=10, n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        # Null envelope
        if sim_dtt is not None and sim_dtt.shape[0] > 0:
            null_median = np.median(sim_dtt, axis=0)
            null_lo = np.percentile(sim_dtt, 2.5, axis=0)
            null_hi = np.percentile(sim_dtt, 97.5, axis=0)

            ax.fill_between(
                times, null_lo, null_hi,
                color="#cccccc", alpha=0.4, label="95% BM null envelope",
                zorder=1,
            )
            ax.plot(
                times, null_median,
                color="#888888", lw=1, linestyle="--",
                label="BM null median", zorder=2,
            )

        # Observed DTT
        ax.plot(
            times, dtt_values,
            color="#2b8cbe", lw=2, solid_capstyle="round",
            label="Observed", zorder=3,
        )

        # BM expectation line (linear decline from 1 to 0)
        ax.plot(
            [0, 1], [1, 0],
            color="#d62728", lw=1, linestyle=":",
            alpha=0.5, label="BM expectation", zorder=1,
        )

        ax.set_xlabel("Relative time", fontsize=10)
        ax.set_ylabel("Relative disparity", fontsize=10)
        ax.set_xlim(0, 1)
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=8, loc="upper right")

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        title = config.title or "Disparity Through Time"
        if mdi is not None:
            title += f" (MDI = {mdi:.4f})"
        if config.show_title:
            ax.set_title(title, fontsize=config.title_fontsize)

        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"DTT plot saved: {self.plot_output}")

    def _print_text(self, times, dtt_values, mdi, mdi_p):
        try:
            print("Disparity Through Time (DTT)")
            print(f"Index: {self.index}")
            print(f"N time points: {len(times)}")
            if mdi is not None:
                print(f"MDI: {mdi:.6f}")
            if mdi_p is not None:
                print(f"MDI p-value: {mdi_p:.4f}")
            print()
            print(f"{'Time':>10s} {'Relative disparity':>20s}")
            print("-" * 32)
            for t, d in zip(times, dtt_values):
                print(f"{t:>10.6f} {d:>20.6f}")
        except BrokenPipeError:
            return

    def _print_json(self, times, dtt_values, mdi, mdi_p, sim_dtt):
        payload = {
            "index": self.index,
            "n_time_points": len(times),
            "times": [round(t, 6) for t in times],
            "dtt": [round(d, 6) for d in dtt_values],
        }
        if mdi is not None:
            payload["mdi"] = round(mdi, 6)
        if mdi_p is not None:
            payload["mdi_p_value"] = round(mdi_p, 4)
        if self.plot_output:
            payload["plot_output"] = self.plot_output
        print_json(payload, sort_keys=False)
