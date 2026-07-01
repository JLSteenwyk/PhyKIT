"""
Disparity through time (DTT).

Calculates and plots disparity-through-time for a phylogenetic tree
and phenotypic data, following Harmon et al. (2003). Computes the
Morphological Disparity Index (MDI) comparing observed DTT to a
Brownian motion null expectation.
"""
from __future__ import annotations

import sys

from .base import Tree
from ...errors import PhykitUserError


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


def _permutation_p_value_ge(permutations: np.ndarray, observed: float) -> float:
    if permutations.size == 0:
        return float("nan")
    return float(np.count_nonzero(permutations >= observed) / permutations.size)


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def parse_multi_trait_file(*args, **kwargs):
    from ...helpers.trait_parsing import (
        parse_multi_trait_file as _parse_multi_trait_file,
    )

    return _parse_multi_trait_file(*args, **kwargs)


def trait_column_from_rows(*args, **kwargs):
    from ...helpers.trait_parsing import (
        trait_column_from_rows as _trait_column_from_rows,
    )

    return _trait_column_from_rows(*args, **kwargs)


def trait_matrix_from_rows(*args, **kwargs):
    from ...helpers.trait_parsing import (
        trait_matrix_from_rows as _trait_matrix_from_rows,
    )

    return _trait_matrix_from_rows(*args, **kwargs)


def _trapezoid(y, x):
    return np.trapezoid(y, x) if hasattr(np, "trapezoid") else np.trapz(y, x)


def _batch_row_sum_squares(values: np.ndarray) -> np.ndarray:
    rows = values.reshape(values.shape[0], -1)
    return np.einsum("ij,ij->i", rows, rows)


def _mean_selected_columns(values: np.ndarray, positions: list[int]) -> np.ndarray:
    if len(positions) == 1:
        return values[:, positions[0]]
    if len(positions) <= 8:
        selected_sum = values[:, positions[0]].copy()
        for position in positions[1:]:
            selected_sum += values[:, position]
        selected_sum /= float(len(positions))
        return selected_sum
    return np.mean(values[:, positions], axis=1)


def _max_terminal_depth(terminals, depths, root_depth):
    iterator = iter(terminals)
    first = next(iterator)
    local_depths = depths
    max_height = local_depths[first] - root_depth
    for terminal in iterator:
        height = local_depths[terminal] - root_depth
        if height > max_height:
            max_height = height
    return max_height


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

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

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
        tree = self.read_tree_file_unmodified()
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
            data = trait_column_from_rows(traits, ordered_names, col_idx).reshape(
                -1, 1
            )
        else:
            data = trait_matrix_from_rows(traits, ordered_names)

        # Prune tree to trait taxa
        tips_to_prune = [t for t in tree_tips if t not in traits]
        tree_for_analysis = tree
        if tips_to_prune:
            tree_for_analysis = self._fast_copy(tree)
            tree_for_analysis = self.prune_tree_using_taxa_list(
                tree_for_analysis,
                tips_to_prune,
            )

        # Compute DTT
        times, dtt_values = self._compute_dtt(
            tree_for_analysis,
            data,
            ordered_names,
        )

        # Compute MDI and null envelope via simulation
        mdi = None
        mdi_p = None
        sim_dtt = None
        if self.nsim > 0:
            sim_dtt, mdi, mdi_p = self._simulate_null(
                tree_for_analysis,
                data,
                ordered_names,
                times,
                dtt_values,
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

        count = n * (n - 1) / 2
        if self.index == "avg_manhattan":
            sorted_data = np.sort(data, axis=0)
            coeffs = 2 * np.arange(n) - n + 1
            total = float(np.sum(sorted_data * coeffs[:, None]))
            return total / count
        else:
            # avg_sq: average squared Euclidean distance
            flattened = data.ravel()
            sum_sq = float(np.dot(flattened, flattened))
            sums = data.sum(axis=0)
            total = n * sum_sq - float(np.dot(sums, sums))
            return total / count

    @staticmethod
    def _root_relative_depths(tree):
        try:
            depths = tree.depths()
            root = tree.root
            root_depth = depths[root]
        except (AttributeError, KeyError):
            return None
        return root, depths, root_depth

    @staticmethod
    def _standard_dtt_traversal_context(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        all_clades = []
        terminals = []
        depths = {}
        parent_map = {}
        stack = [(root, 0.0)]
        pop = stack.pop
        append_clade = all_clades.append
        append_terminal = terminals.append
        try:
            while stack:
                clade, depth = pop()
                append_clade(clade)
                depths[clade] = depth
                children = clade.clades
                if children:
                    for child in reversed(children):
                        parent_map[id(child)] = clade
                        stack.append(
                            (child, depth + (child.branch_length or 0.0))
                        )
                else:
                    append_terminal(clade)
        except (AttributeError, TypeError):
            return None

        return root, all_clades, terminals, depths, parent_map

    def _prepare_dtt_context(self, tree, ordered_names: list[str]) -> dict:
        standard_context = self._standard_dtt_traversal_context(tree)
        if standard_context is None:
            depth_data = self._root_relative_depths(tree)
            all_clades = list(tree.find_clades(order="preorder"))
            terminals = [clade for clade in all_clades if clade.is_terminal()]
            if depth_data is None:
                root = tree.root
                tip_heights = [tree.distance(root, t) for t in terminals]
                tree_height = max(tip_heights)
            else:
                root, depths, root_depth = depth_data
                tree_height = _max_terminal_depth(terminals, depths, root_depth)
            parent_map = {}
        else:
            root, all_clades, terminals, depths, parent_map = standard_context
            root_depth = 0.0
            tree_height = _max_terminal_depth(terminals, depths, root_depth)
        if tree_height == 0:
            return {"zero_height": True}

        name_to_idx = {name: i for i, name in enumerate(ordered_names)}
        clade_tip_indices = {}
        for clade in reversed(all_clades):
            children = clade.clades
            if children:
                clade_tip_indices[id(clade)] = tuple(
                    idx
                    for child in children
                    for idx in clade_tip_indices.get(id(child), ())
                )
            elif clade.name in name_to_idx:
                clade_tip_indices[id(clade)] = (name_to_idx[clade.name],)
            else:
                clade_tip_indices[id(clade)] = ()

        node_time = {}
        for clade in all_clades:
            if standard_context is None and depth_data is None:
                node_time[id(clade)] = tree.distance(root, clade) / tree_height
            else:
                node_time[id(clade)] = (depths[clade] - root_depth) / tree_height
            for child in clade.clades:
                parent_map[id(child)] = clade

        internal_nodes = [
            c for c in all_clades
            if c.clades
        ]
        internal_nodes.sort(key=lambda c: node_time[id(c)])

        lineage_starts = []
        lineage_ends = []
        for clade in all_clades:
            if clade == root:
                continue

            clade_id = id(clade)
            if len(clade_tip_indices[clade_id]) < 2:
                continue

            parent = parent_map.get(clade_id)
            if parent is None:
                continue

            parent_t = node_time[id(parent)]
            child_t = 1.0 if not clade.clades else node_time[clade_id]
            lineage_starts.append((parent_t, clade_id))
            lineage_ends.append((child_t, clade_id))

        lineage_starts.sort()
        lineage_ends.sort()

        times = [0.0]
        lineage_clade_ids = []
        active_lineages = set()
        start_idx = 0
        end_idx = 0
        eps = 1e-10

        for node in internal_nodes:
            t = node_time[id(node)]
            time_limit = t + eps
            while (
                start_idx < len(lineage_starts)
                and lineage_starts[start_idx][0] <= time_limit
            ):
                active_lineages.add(lineage_starts[start_idx][1])
                start_idx += 1

            while (
                end_idx < len(lineage_ends)
                and lineage_ends[end_idx][0] <= time_limit
            ):
                active_lineages.discard(lineage_ends[end_idx][1])
                end_idx += 1

            times.append(t)
            lineage_clade_ids.append(list(active_lineages))

        clade_index_arrays = {
            clade_id: np.asarray(indices, dtype=np.intp)
            for clade_id, indices in clade_tip_indices.items()
            if len(indices) >= 2
        }
        terminal_index_by_clade_id = {
            id(clade): name_to_idx[clade.name]
            for clade in terminals
            if clade.name in name_to_idx
        }
        child_clade_ids = {
            id(clade): tuple(id(child) for child in clade.clades)
            for clade in all_clades
            if clade.clades
        }

        return {
            "zero_height": False,
            "times": times,
            "lineage_clade_ids": lineage_clade_ids,
            "clade_index_arrays": clade_index_arrays,
            "postorder_clade_ids": [id(clade) for clade in reversed(all_clades)],
            "terminal_index_by_clade_id": terminal_index_by_clade_id,
            "child_clade_ids": child_clade_ids,
        }

    def _compute_dtt(
        self, tree, data: np.ndarray, ordered_names: list[str], context=None
    ) -> tuple[list[float], list[float]]:
        """Compute DTT curve following Harmon et al. (2003).

        At each branching time, compute the mean relative disparity of
        subclades with >= 2 species. Times are relative (0 = root, 1 = tips).
        """
        if context is None:
            context = self._prepare_dtt_context(tree, ordered_names)

        if context["zero_height"]:
            return [0.0, 1.0], [1.0, 0.0]

        total_disp = self._compute_disparity(data)
        if total_disp == 0:
            return [0.0, 1.0], [1.0, 0.0]

        clade_disp = {
            clade_id: self._compute_disparity(data[indices])
            for clade_id, indices in context["clade_index_arrays"].items()
        }

        times = list(context["times"])
        dtt_values = [1.0]

        for lineage_ids in context["lineage_clade_ids"]:
            lineage_disparities = [
                clade_disp.get(clade_id, 0.0)
                for clade_id in lineage_ids
            ]
            if lineage_disparities:
                rel_disp = (
                    sum(lineage_disparities)
                    / len(lineage_disparities)
                    / total_disp
                )
            else:
                rel_disp = 0.0
            dtt_values.append(rel_disp)

        # Terminal point (only if last entry isn't already at/near 0)
        if dtt_values and dtt_values[-1] > 1e-10:
            times.append(1.0)
            dtt_values.append(0.0)

        return times, dtt_values

    def _compute_mdi(
        self, times: list[float], dtt_values: list[float],
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
        mdi = float(_trapezoid(diff, times_arr))
        return mdi

    def _simulate_null(
        self, tree, data, ordered_names, obs_times, obs_values=None,
    ) -> tuple[np.ndarray, float, float]:
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

        # Simulate nsim datasets
        n_time_points = len(obs_times)
        sim_dtt_matrix = np.zeros((self.nsim, n_time_points))
        dtt_context = self._prepare_dtt_context(tree, ordered_names)

        if self.index == "avg_sq" and not dtt_context["zero_height"]:
            return self._simulate_null_avg_sq_batch(
                L,
                data,
                obs_times,
                dtt_context,
                obs_values=obs_values,
            )

        for s in range(self.nsim):
            # Simulate BM data with same VCV structure
            if p == 1:
                z = rng.standard_normal(n)
                sim_data = (L @ z).reshape(-1, 1)
            else:
                Z = rng.standard_normal((n, p))
                sim_data = L @ Z

            sim_times, sim_values = self._compute_dtt(
                tree, sim_data, ordered_names, context=dtt_context
            )

            # Interpolate to observed time points
            sim_interp = np.interp(obs_times, sim_times, sim_values)
            sim_dtt_matrix[s, :] = sim_interp

        # Null median
        null_median = np.median(sim_dtt_matrix, axis=0)

        # MDI
        obs_arr = np.array(
            self._compute_dtt(tree, data, ordered_names, context=dtt_context)[1]
        )
        # Interpolate observed to consistent times
        obs_times_arr = np.array(obs_times)
        mdi = float(_trapezoid(
            np.interp(obs_times_arr, obs_times, obs_arr) - null_median,
            obs_times_arr,
        ))

        # MDI p-value: proportion of simulated MDIs >= observed
        sim_mdis = np.zeros(self.nsim)
        for s in range(self.nsim):
            sim_mdi = float(_trapezoid(
                sim_dtt_matrix[s, :] - null_median, obs_times_arr
            ))
            sim_mdis[s] = sim_mdi

        mdi_p = _permutation_p_value_ge(np.abs(sim_mdis), np.abs(mdi))

        return sim_dtt_matrix, mdi, mdi_p

    def _simulate_null_avg_sq_batch(
        self,
        cholesky_factor: np.ndarray,
        data: np.ndarray,
        obs_times: list[float],
        dtt_context: dict,
        obs_values=None,
    ) -> tuple[np.ndarray, float, float]:
        rng = np.random.default_rng(self.seed)
        n = cholesky_factor.shape[0]
        p = data.shape[1] if data.ndim > 1 else 1

        if p == 1:
            z = rng.standard_normal((self.nsim, n))
            sim_data = (z @ cholesky_factor.T)[:, :, None]
        else:
            z = rng.standard_normal((self.nsim, n, p))
            sim_data = np.matmul(cholesky_factor, z)

        pair_count = n * (n - 1) / 2
        sim_sums = np.sum(sim_data, axis=1)
        total_disp = (
            n * _batch_row_sum_squares(sim_data)
            - _batch_row_sum_squares(sim_sums)
        ) / pair_count

        clade_ids = list(dtt_context["clade_index_arrays"].keys())
        clade_id_to_pos = {
            clade_id: idx for idx, clade_id in enumerate(clade_ids)
        }
        clade_disp = self._batch_clade_disparities_avg_sq(
            sim_data,
            dtt_context,
            clade_id_to_pos,
        )

        sim_times = np.asarray(dtt_context["times"], dtype=float)
        sim_values = np.empty(
            (self.nsim, 1 + len(dtt_context["lineage_clade_ids"])),
            dtype=float,
        )
        sim_values[:, 0] = 1.0
        valid_total = total_disp != 0.0
        for time_idx, lineage_ids in enumerate(
            dtt_context["lineage_clade_ids"],
            start=1,
        ):
            positions = [
                clade_id_to_pos[clade_id]
                for clade_id in lineage_ids
                if clade_id in clade_id_to_pos
            ]
            if positions:
                mean_disp = _mean_selected_columns(clade_disp, positions)
                rel_disp = np.zeros(self.nsim)
                np.divide(
                    mean_disp,
                    total_disp,
                    out=rel_disp,
                    where=valid_total,
                )
                sim_values[:, time_idx] = rel_disp
            else:
                sim_values[:, time_idx] = 0.0

        obs_times_arr = np.asarray(obs_times, dtype=float)
        if (
            len(obs_times_arr) == len(sim_times) + 1
            and abs(obs_times_arr[-1] - 1.0) < 1e-10
        ):
            extended_times = np.empty(len(sim_times) + 1, dtype=float)
            extended_times[:-1] = sim_times
            extended_times[-1] = 1.0

            extended_values = np.empty(
                (self.nsim, sim_values.shape[1] + 1),
                dtype=float,
            )
            extended_values[:, :-1] = sim_values
            last_values = sim_values[:, -1]
            final_values = extended_values[:, -1]
            final_values[:] = last_values
            final_values[last_values > 1e-10] = 0.0
            sim_times = extended_times
            sim_values = extended_values

        if (
            len(obs_times_arr) == len(sim_times)
            and np.array_equal(obs_times_arr, sim_times)
            and (np.diff(sim_times) > 0.0).all()
        ):
            sim_dtt_matrix = sim_values.copy()
            if not valid_total.all():
                sim_dtt_matrix[~valid_total, :] = 1.0 - obs_times_arr
        else:
            sim_dtt_matrix = np.empty((self.nsim, len(obs_times_arr)), dtype=float)
            for sim_idx in range(self.nsim):
                if not valid_total[sim_idx]:
                    sim_dtt_matrix[sim_idx, :] = np.interp(
                        obs_times_arr,
                        [0.0, 1.0],
                        [1.0, 0.0],
                    )
                else:
                    sim_dtt_matrix[sim_idx, :] = np.interp(
                        obs_times_arr,
                        sim_times,
                        sim_values[sim_idx],
                    )

        null_median = np.median(sim_dtt_matrix, axis=0)
        if obs_values is None:
            obs_values = self._compute_dtt(
                None,
                data,
                [],
                context=dtt_context,
            )[1]
        obs_values = np.asarray(obs_values, dtype=float)
        mdi = float(_trapezoid(
            obs_values - null_median,
            obs_times_arr,
        ))

        sim_mdis = np.asarray(_trapezoid(sim_dtt_matrix - null_median, obs_times_arr))
        mdi_p = _permutation_p_value_ge(np.abs(sim_mdis), np.abs(mdi))

        return sim_dtt_matrix, mdi, mdi_p

    def _batch_clade_disparities_avg_sq(
        self,
        sim_data: np.ndarray,
        dtt_context: dict,
        clade_id_to_pos: dict[int, int],
    ) -> np.ndarray:
        """Compute per-simulation avg-squared disparities for prepared clades."""
        sim_count = sim_data.shape[0]
        clade_disp = np.zeros((sim_count, len(clade_id_to_pos)))
        postorder_ids = dtt_context.get("postorder_clade_ids")
        terminal_indices = dtt_context.get("terminal_index_by_clade_id")
        child_clade_ids = dtt_context.get("child_clade_ids")
        if not postorder_ids or terminal_indices is None or child_clade_ids is None:
            for clade_id, clade_pos in clade_id_to_pos.items():
                indices = dtt_context["clade_index_arrays"][clade_id]
                m = len(indices)
                if m < 2:
                    continue
                subset = sim_data[:, indices, :]
                subset_pair_count = m * (m - 1) / 2
                subset_sums = np.sum(subset, axis=1)
                clade_disp[:, clade_pos] = (
                    m * _batch_row_sum_squares(subset)
                    - _batch_row_sum_squares(subset_sums)
                ) / subset_pair_count
            return clade_disp

        counts = {}
        sums = {}
        sums_sq = {}
        for clade_id in postorder_ids:
            terminal_idx = terminal_indices.get(clade_id)
            if terminal_idx is not None:
                values = sim_data[:, terminal_idx, :]
                counts[clade_id] = 1
                sums[clade_id] = values
                sums_sq[clade_id] = _batch_row_sum_squares(values)
                continue

            children = [
                child_id
                for child_id in child_clade_ids.get(clade_id, ())
                if child_id in counts
            ]
            if not children:
                continue

            count = 0
            total_sum = None
            total_sum_sq = None
            for child_id in children:
                count += counts[child_id]
                if total_sum is None:
                    total_sum = sums[child_id].copy()
                    total_sum_sq = sums_sq[child_id].copy()
                else:
                    total_sum += sums[child_id]
                    total_sum_sq += sums_sq[child_id]

            counts[clade_id] = count
            sums[clade_id] = total_sum
            sums_sq[clade_id] = total_sum_sq

            clade_pos = clade_id_to_pos.get(clade_id)
            if clade_pos is None or count < 2:
                continue

            pair_count = count * (count - 1) / 2
            clade_disp[:, clade_pos] = (
                count * total_sum_sq
                - _batch_row_sum_squares(total_sum)
            ) / pair_count

        return clade_disp

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
            lines = [
                "Disparity Through Time (DTT)",
                f"Index: {self.index}",
                f"N time points: {len(times)}",
            ]
            if mdi is not None:
                lines.append(f"MDI: {mdi:.6f}")
            if mdi_p is not None:
                lines.append(f"MDI p-value: {mdi_p:.4f}")
            lines.extend(
                [
                    "",
                    f"{'Time':>10s} {'Relative disparity':>20s}",
                    "-" * 32,
                ]
            )
            lines.extend(f"{t:>10.6f} {d:>20.6f}" for t, d in zip(times, dtt_values))
            print("\n".join(lines))
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
