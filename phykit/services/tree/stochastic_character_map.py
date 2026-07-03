from __future__ import annotations

import sys

from .base import Tree
from ...errors import PhykitUserError


class _LazyNumpy:
    _module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = _np
            self._module = module

        value = getattr(module, name)
        setattr(self, name, value)
        return value


np = _LazyNumpy()


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class StochasticCharacterMap(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.trait_column = parsed["trait_column"]
        self.model = parsed["model"]
        self.nsim = parsed["nsim"]
        self.seed = parsed["seed"]
        self.plot_output = parsed["plot_output"]
        self.json_output = parsed["json_output"]
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
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        self.validate_tree(tree, min_tips=3, require_branch_lengths=True, context="stochastic character mapping")

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

        # Compute conditional likelihoods once (they don't change between sims)
        cond_liks, _ = self._felsenstein_pruning(
            tree, tip_states, Q, pi, states
        )

        # Build parent lookup dict once
        parent_map = self._build_parent_map(tree)
        simulation_metadata = self._build_simulation_metadata(
            tree, tip_states, states, parent_map
        )
        sampling_nodes = self._prepare_ancestral_sampling_nodes(
            simulation_metadata[0],
            Q,
            cond_liks,
            len(states),
            transition_cache={},
        )
        branch_history_context = self._prepare_branch_history_context(
            Q,
            len(states),
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
                sampling_nodes=sampling_nodes,
                branch_history_context=branch_history_context,
            )
            mappings.append(mapping)

        # Summarize
        summary = self._summarize_simulations(mappings, states, tree)

        # Plot
        if self.plot_output:
            self._plot_stochastic_map(
                tree, mappings[0], states, self.plot_output,
                parent_map=parent_map,
            )

        if self.json_output:
            result = self._format_result(
                Q, loglik, states, summary, self.plot_output
            )
            print_json(result)
        else:
            self._print_text_output(Q, loglik, states, summary)

    def _parse_discrete_trait_file(
        self, path: str, column: str, tree_tips: list[str]
    ) -> dict[str, str]:
        from ...helpers.discrete_models import parse_discrete_traits

        return parse_discrete_traits(path, tree_tips, trait_column=column)

    @staticmethod
    def _count_missing_tip_states(tree, tip_states: dict[str, str]) -> int | None:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        missing_tip_state_count = 0
        stack = [root]
        pop = stack.pop
        extend = stack.extend
        contains = tip_states.__contains__
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    extend(children)
                elif not contains(clade.name):
                    missing_tip_state_count += 1
        except AttributeError:
            return None

        return missing_tip_state_count

    @staticmethod
    def _prune_tree_to_tip_states(tree, tip_states: dict[str, str]) -> int:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            root = None

        if root is not None:
            tips_to_prune = []
            target_ids = set()
            terminal_count = 0
            stack = [root]
            pop = stack.pop
            append = stack.append
            append_prune = tips_to_prune.append
            add_target = target_ids.add
            contains = tip_states.__contains__
            try:
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
                        terminal_count += 1
                        if not contains(clade.name):
                            append_prune(clade)
                            add_target(id(clade))
            except AttributeError:
                tips_to_prune = None

            if tips_to_prune is not None:
                if (
                    tips_to_prune
                    and len(tips_to_prune) < terminal_count
                    and Tree._prune_terminal_objects_batch_standard_tree(
                        tree,
                        target_ids,
                    )
                ):
                    return len(tips_to_prune)
                for tip in tips_to_prune:
                    tree.prune(tip)
                return len(tips_to_prune)

        tips_to_prune = [
            tip for tip in tree.get_terminals() if tip.name not in tip_states
        ]
        for tip in tips_to_prune:
            tree.prune(tip)
        return len(tips_to_prune)

    def _build_q_matrix(self, params, k, model):
        from ...helpers.discrete_models import build_q_matrix

        return build_q_matrix(params, k, model)

    def _matrix_exp(self, Q, t):
        from ...helpers.discrete_models import matrix_exp

        return matrix_exp(Q, t)

    def _felsenstein_pruning(self, tree, tip_states, Q, pi, states):
        from ...helpers.discrete_models import felsenstein_pruning

        return felsenstein_pruning(tree, tip_states, Q, pi, states)

    def _fit_q_matrix(self, tree, tip_states, states, model):
        from ...helpers.discrete_models import fit_q_matrix

        return fit_q_matrix(tree, tip_states, states, model)

    @staticmethod
    def _sample_categorical(probs: np.ndarray, rng) -> int:
        return StochasticCharacterMap._sample_categorical_cdf(
            np.cumsum(probs),
            rng,
        )

    @staticmethod
    def _sample_categorical_cdf(cdf: np.ndarray, rng) -> int:
        value = rng.random()
        length = len(cdf)
        if length == 2:
            if value < cdf[0]:
                return 0
            return 1
        if length == 3:
            if value < cdf[0]:
                return 0
            if value < cdf[1]:
                return 1
            return 2
        if length == 4:
            if value < cdf[0]:
                return 0
            if value < cdf[1]:
                return 1
            if value < cdf[2]:
                return 2
            return 3
        if length == 5:
            if value < cdf[2]:
                if value < cdf[0]:
                    return 0
                if value < cdf[1]:
                    return 1
                return 2
            if value < cdf[3]:
                return 3
            return 4
        if length == 6:
            if value < cdf[2]:
                if value < cdf[0]:
                    return 0
                if value < cdf[1]:
                    return 1
                return 2
            if value < cdf[4]:
                if value < cdf[3]:
                    return 3
                return 4
            return 5
        if length == 7:
            if value < cdf[3]:
                if value < cdf[1]:
                    if value < cdf[0]:
                        return 0
                    return 1
                if value < cdf[2]:
                    return 2
                return 3
            if value < cdf[5]:
                if value < cdf[4]:
                    return 4
                return 5
            return 6
        if length == 8:
            if value < cdf[3]:
                if value < cdf[1]:
                    if value < cdf[0]:
                        return 0
                    return 1
                if value < cdf[2]:
                    return 2
                return 3
            if value < cdf[5]:
                if value < cdf[4]:
                    return 4
                return 5
            if value < cdf[6]:
                return 6
            return 7

        idx = int(np.searchsorted(cdf, value, side="right"))
        if idx == length:
            idx -= 1
        return idx

    def _sample_ancestral_states(
        self, tree, tip_states: dict[str, str],
        Q: np.ndarray, pi: np.ndarray,
        cond_liks: dict, states: list[str], rng,
        parent_map: dict = None, transition_cache: dict = None,
        simulation_nodes=None, sampling_nodes=None,
    ) -> dict:
        k = len(states)
        node_states = {}

        # Sample root
        root_lik = cond_liks[id(tree.root)]
        root_probs = pi * root_lik
        root_probs /= root_probs.sum()
        root_state_idx = self._sample_categorical(root_probs, rng)
        node_states[id(tree.root)] = root_state_idx

        if simulation_nodes is None:
            if parent_map is None:
                parent_map = self._build_parent_map(tree)
            simulation_nodes, _ = self._build_simulation_metadata(
                tree, tip_states, states, parent_map
            )
        if sampling_nodes is None:
            sampling_nodes = self._prepare_ancestral_sampling_nodes(
                simulation_nodes,
                Q,
                cond_liks,
                k,
                transition_cache,
            )

        for clade_id, parent_id, tip_state_idx, cdfs_by_parent in sampling_nodes:
            parent_state = node_states[parent_id]

            if tip_state_idx is not None:
                node_states[clade_id] = tip_state_idx
            else:
                node_states[clade_id] = self._sample_categorical_cdf(
                    cdfs_by_parent[parent_state],
                    rng,
                )

        return node_states

    def _prepare_ancestral_sampling_nodes(
        self,
        simulation_nodes,
        Q: np.ndarray,
        cond_liks: dict,
        k: int,
        transition_cache: dict = None,
    ):
        sampling_nodes = []
        uniform = np.ones(k) / k
        uniform_cdf = np.cumsum(uniform)
        for clade_id, parent_id, branch_length, tip_state_idx in simulation_nodes:
            if tip_state_idx is not None:
                sampling_nodes.append((clade_id, parent_id, tip_state_idx, None))
                continue

            if transition_cache is None:
                P = self._matrix_exp(Q, branch_length)
            else:
                P = transition_cache.get(branch_length)
                if P is None:
                    P = self._matrix_exp(Q, branch_length)
                    transition_cache[branch_length] = P

            child_lik = cond_liks[clade_id]
            if k < 5:
                cdfs_by_parent = np.empty((k, k), dtype=float)
                for parent_state in range(k):
                    probs = P[parent_state, :] * child_lik
                    total = probs.sum()
                    if total > 0:
                        cdfs_by_parent[parent_state, :] = np.cumsum(probs / total)
                    else:
                        cdfs_by_parent[parent_state, :] = uniform_cdf
            else:
                weighted = P * child_lik
                totals = weighted.sum(axis=1)
                if (totals > 0).all():
                    cdfs_by_parent = np.cumsum(
                        weighted / totals[:, None],
                        axis=1,
                    )
                else:
                    positive = totals > 0
                    cdfs_by_parent = np.empty((k, k), dtype=float)
                    if positive.any():
                        cdfs_by_parent[positive] = np.cumsum(
                            weighted[positive] / totals[positive, None],
                            axis=1,
                        )
                    cdfs_by_parent[~positive] = uniform_cdf
            sampling_nodes.append(
                (clade_id, parent_id, tip_state_idx, cdfs_by_parent)
            )
        return sampling_nodes

    def _build_parent_map(self, tree) -> dict:
        parent_map = {}
        try:
            root = tree.root
            root.clades
        except AttributeError:
            root = None

        if root is not None:
            stack = [root]
            try:
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
            except AttributeError:
                parent_map = {}

        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

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

    def _get_parent(self, tree, target, parent_map=None):
        if parent_map is not None:
            return parent_map.get(id(target))
        for clade in tree.find_clades(order="preorder"):
            if target in clade.clades:
                return clade
        return None

    def _build_simulation_metadata(self, tree, tip_states, states, parent_map):
        state_idx = {state: i for i, state in enumerate(states)}
        simulation_nodes = []
        branch_nodes = []
        try:
            root = tree.root
            root.clades
        except AttributeError:
            root = None

        if root is not None:
            try:
                for clade in self._iter_preorder(root):
                    children = clade.clades
                    if clade is not root:
                        parent = parent_map.get(id(clade))
                        if parent is not None:
                            clade_id = id(clade)
                            parent_id = id(parent)
                            branch_length = (
                                clade.branch_length
                                if clade.branch_length
                                else 1e-8
                            )
                            tip_state_idx = (
                                state_idx[tip_states[clade.name]]
                                if not children and clade.name in tip_states
                                else None
                            )
                            simulation_nodes.append(
                                (
                                    clade_id,
                                    parent_id,
                                    branch_length,
                                    tip_state_idx,
                                )
                            )
                            branch_nodes.append(
                                (clade_id, parent_id, branch_length)
                            )
                return simulation_nodes, branch_nodes
            except AttributeError:
                simulation_nodes = []
                branch_nodes = []

        for clade in tree.find_clades(order="preorder"):
            if clade == tree.root:
                continue
            parent = parent_map.get(id(clade))
            if parent is None:
                continue
            clade_id = id(clade)
            parent_id = id(parent)
            branch_length = clade.branch_length if clade.branch_length else 1e-8
            tip_state_idx = (
                state_idx[tip_states[clade.name]]
                if clade.is_terminal() and clade.name in tip_states
                else None
            )
            simulation_nodes.append(
                (clade_id, parent_id, branch_length, tip_state_idx)
            )
            branch_nodes.append((clade_id, parent_id, branch_length))
        return simulation_nodes, branch_nodes

    def _simulate_branch_history(
        self, Q: np.ndarray, start_state: int, end_state: int,
        branch_length: float, k: int, rng, max_attempts: int = 1000,
        rates: np.ndarray = None, transition_probs_by_state: list[np.ndarray] = None,
        transition_cdfs_by_state: list[np.ndarray] = None,
        deterministic_next_states: list[int] = None,
    ) -> list[tuple[float, int]]:
        if branch_length <= 0:
            return [(0.0, start_state)]

        if rates is None:
            rates = -Q.diagonal()
        if transition_cdfs_by_state is None:
            if transition_probs_by_state is None:
                transition_probs_by_state = self._transition_probs_by_state(Q, k)
            transition_cdfs_by_state = self._transition_cdfs_by_state(
                transition_probs_by_state
            )

        for _ in range(max_attempts):
            history = [(0.0, start_state)]
            current_state = start_state
            current_time = 0.0

            while current_time < branch_length:
                rate = rates[current_state]
                if rate <= 0:
                    break
                wait_time = rng.exponential(1.0 / rate)
                if current_time + wait_time >= branch_length:
                    break
                current_time += wait_time
                # Choose next state
                trans_cdf = transition_cdfs_by_state[current_state]
                if trans_cdf is None:
                    break
                if deterministic_next_states is not None:
                    new_state = deterministic_next_states[current_state]
                    if new_state >= 0:
                        rng.random()
                    else:
                        new_state = self._sample_categorical_cdf(trans_cdf, rng)
                else:
                    new_state = self._sample_categorical_cdf(trans_cdf, rng)
                current_state = new_state
                history.append((current_time, current_state))

            if current_state == end_state:
                return history

        # Fallback: direct assignment
        if start_state == end_state:
            return [(0.0, start_state)]
        else:
            midpoint = branch_length / 2.0
            return [(0.0, start_state), (midpoint, end_state)]

    @staticmethod
    def _transition_probs_by_state(Q: np.ndarray, k: int) -> list[np.ndarray]:
        transition_probs = []
        for state in range(k):
            probs = Q[state, :].copy()
            probs[state] = 0.0
            total = probs.sum()
            transition_probs.append((probs / total) if total > 0 else None)
        return transition_probs

    @staticmethod
    def _transition_cdfs_by_state(
        transition_probs_by_state: list[np.ndarray],
    ) -> list[np.ndarray]:
        return [
            np.cumsum(probs) if probs is not None else None
            for probs in transition_probs_by_state
        ]

    def _prepare_branch_history_context(self, Q: np.ndarray, k: int):
        transition_probs_by_state = self._transition_probs_by_state(Q, k)
        return (
            -Q.diagonal(),
            self._transition_cdfs_by_state(transition_probs_by_state),
            self._deterministic_transition_targets(transition_probs_by_state),
        )

    @staticmethod
    def _deterministic_transition_targets(
        transition_probs_by_state: list[np.ndarray],
    ) -> list[int]:
        targets = []
        for probs in transition_probs_by_state:
            if probs is None:
                targets.append(-1)
                continue
            nonzero = np.flatnonzero(probs > 0.0)
            if len(nonzero) == 1:
                targets.append(int(nonzero[0]))
            else:
                targets.append(-1)
        return targets

    def _run_single_simulation(
        self, tree, tip_states: dict[str, str],
        Q: np.ndarray, pi: np.ndarray,
        states: list[str], rng,
        cond_liks: dict = None, parent_map: dict = None,
        transition_cache: dict = None,
        simulation_metadata=None, sampling_nodes=None,
        branch_history_context=None,
    ) -> dict:
        k = len(states)

        if cond_liks is None:
            cond_liks, _ = self._felsenstein_pruning(
                tree, tip_states, Q, pi, states
            )
        if parent_map is None:
            parent_map = self._build_parent_map(tree)
        if simulation_metadata is None:
            simulation_metadata = self._build_simulation_metadata(
                tree, tip_states, states, parent_map
            )
        simulation_nodes, branch_nodes = simulation_metadata
        node_states = self._sample_ancestral_states(
            tree, tip_states, Q, pi, cond_liks, states, rng,
            parent_map=parent_map, transition_cache=transition_cache,
            simulation_nodes=simulation_nodes,
            sampling_nodes=sampling_nodes,
        )

        branch_histories = {}
        if branch_history_context is None:
            branch_history_context = self._prepare_branch_history_context(
                Q,
                k,
            )
        if len(branch_history_context) == 2:
            rates, transition_cdfs_by_state = branch_history_context
            deterministic_next_states = None
        else:
            rates, transition_cdfs_by_state, deterministic_next_states = (
                branch_history_context
            )
        for clade_id, parent_id, t in branch_nodes:
            start_state = node_states[parent_id]
            end_state = node_states[clade_id]

            history = self._simulate_branch_history(
                Q, start_state, end_state, t, k, rng,
                rates=rates,
                transition_cdfs_by_state=transition_cdfs_by_state,
                deterministic_next_states=deterministic_next_states,
            )
            branch_histories[clade_id] = history

        return {
            "node_states": node_states,
            "branch_histories": branch_histories,
        }

    def _summarize_simulations(
        self, mappings: list[dict], states: list[str], tree
    ) -> dict:
        k = len(states)
        nsim = len(mappings)

        total_dwelling = np.zeros(k)
        total_transitions = np.zeros((k, k))
        node_posteriors = {}

        branch_lengths = {}
        internal_nodes = []
        try:
            root = tree.root
            root.clades
        except AttributeError:
            root = None

        if root is not None:
            stack = [root]
            try:
                while stack:
                    clade = stack.pop()
                    children = clade.clades
                    cid = id(clade)
                    if clade is not root:
                        branch_lengths[cid] = (
                            clade.branch_length
                            if clade.branch_length
                            else 1e-8
                        )
                    if children:
                        internal_nodes.append(cid)
                        node_posteriors[cid] = np.zeros(k)
                        stack.extend(reversed(children))
            except AttributeError:
                branch_lengths = {}
                internal_nodes = []
                node_posteriors = {}

        if root is None or not internal_nodes:
            for clade in tree.find_clades(order="preorder"):
                cid = id(clade)
                if clade != tree.root:
                    branch_lengths[cid] = (
                        clade.branch_length if clade.branch_length else 1e-8
                    )
                if clade.is_terminal():
                    continue
                internal_nodes.append(cid)
                node_posteriors[cid] = np.zeros(k)

        for mapping in mappings:
            histories = mapping["branch_histories"]
            node_states = mapping["node_states"]

            for clade_id, history in histories.items():
                t = branch_lengths.get(clade_id)
                if t is None:
                    continue

                if len(history) == 1:
                    total_dwelling[history[0][1]] += t
                    continue

                prev_time, prev_state = history[0]
                for h_idx in range(1, len(history)):
                    time_val, state = history[h_idx]
                    total_dwelling[prev_state] += time_val - prev_time
                    if prev_state != state:
                        total_transitions[prev_state, state] += 1
                    prev_time = time_val
                    prev_state = state
                total_dwelling[prev_state] += t - prev_time

            # Node posteriors
            for node_id in internal_nodes:
                if node_id in node_states:
                    state = node_states[node_id]
                    node_posteriors[node_id][state] += 1

        # Average across simulations
        mean_dwelling = total_dwelling / nsim
        mean_transitions = total_transitions / nsim

        # Normalize node posteriors
        for node_id in node_posteriors:
            node_posteriors[node_id] /= nsim

        return {
            "mean_dwelling_times": mean_dwelling,
            "mean_transitions": mean_transitions,
            "node_posteriors": node_posteriors,
        }

    def _print_text_output(
        self, Q: np.ndarray, loglik: float,
        states: list[str], summary: dict
    ) -> None:
        lines = [
            "Stochastic Character Mapping (SIMMAP)",
            f"\nModel: {self.model}",
            f"Number of simulations: {self.nsim}",
        ]
        append = lines.append

        append("\nFitted Q matrix:")
        # Header
        header = "              " + "".join(f"{s:>12s}" for s in states)
        append(header)
        q_row_format = "%-14s" + "%12.4f" * len(states)
        for s, q_row in zip(states, Q):
            append(q_row_format % ((s,) + tuple(q_row)))

        append(f"\nLog-likelihood: {loglik:.4f}")

        mean_dwelling = summary["mean_dwelling_times"]
        total_dwelling = mean_dwelling.sum()
        append("\nMean dwelling times:")
        for s, value in zip(states, mean_dwelling):
            pct = 100.0 * value / total_dwelling if total_dwelling > 0 else 0
            append(f"  {s:14s} {value:8.2f} ({pct:.1f}%)")

        mean_trans = summary["mean_transitions"]
        append("\nMean transitions:")
        total_trans = 0.0
        for si, trans_row in zip(states, mean_trans):
            for sj, value in zip(states, trans_row):
                if si != sj and value > 0:
                    append(
                        f"  {si} -> {sj}:"
                        f"  {value:.2f}"
                    )
                    total_trans += value
        append(f"  {'Total:':28s}{total_trans:.2f}")
        print("\n".join(lines))

    def _format_result(
        self, Q: np.ndarray, loglik: float,
        states: list[str], summary: dict,
        plot_output: str
    ) -> dict:
        q_dict = {}
        for si, q_row in zip(states, Q):
            q_dict[si] = {
                sj: float(value) for sj, value in zip(states, q_row)
            }

        mean_dwelling = summary["mean_dwelling_times"]
        total_dwelling = float(mean_dwelling.sum())

        dwelling_dict = {
            state: float(value) for state, value in zip(states, mean_dwelling)
        }
        proportion_dict = {
            state: float(value / total_dwelling)
            if total_dwelling > 0 else 0.0
            for state, value in zip(states, mean_dwelling)
        }

        mean_trans = summary["mean_transitions"]
        trans_dict = {}
        total_trans = 0.0
        for si, trans_row in zip(states, mean_trans):
            for sj, value in zip(states, trans_row):
                if si != sj:
                    scalar = float(value)
                    key = f"{si} -> {sj}"
                    trans_dict[key] = scalar
                    total_trans += scalar

        result = {
            "model": self.model,
            "nsim": self.nsim,
            "q_matrix": q_dict,
            "log_likelihood": float(loglik),
            "states": states,
            "mean_dwelling_times": dwelling_dict,
            "mean_dwelling_proportions": proportion_dict,
            "mean_transitions": trans_dict,
            "mean_total_transitions": total_trans,
        }

        if plot_output:
            result["plot_output"] = plot_output

        return result

    def _plot_stochastic_map(
        self, tree, mapping: dict, states: list[str],
        output_path: str, parent_map: dict = None,
    ) -> None:
        from ...helpers.plot_config import compute_node_positions

        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.collections import LineCollection
            from matplotlib.lines import Line2D
        except ImportError:
            print(
                "matplotlib is required for stochastic character map plotting. "
                "Install matplotlib and retry."
            )
            raise SystemExit(2)

        if parent_map is None:
            parent_map = self._build_parent_map(tree)

        k = len(states)
        cmap = plt.get_cmap("tab10")
        state_colors = {
            i: cmap(i / max(k - 1, 1)) for i in range(k)
        }

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

        histories = mapping["branch_histories"]

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
                angle = coords[cid]["angle"]
                r_p = coords[pid]["radius"]
                r_c = coords[cid]["radius"]
                radius_delta = r_c - r_p
                cos_angle = _math.cos(angle)
                sin_angle = _math.sin(angle)

                if cid in histories:
                    history = histories[cid]
                    for h_idx in range(len(history)):
                        state = history[h_idx][1]
                        seg_start = history[h_idx][0]
                        if h_idx + 1 < len(history):
                            seg_end = history[h_idx + 1][0]
                        else:
                            seg_end = t
                        if t > 0:
                            start_frac = seg_start / t
                            end_frac = seg_end / t
                        else:
                            start_frac = 0.0
                            end_frac = 1.0
                        r0 = r_p + radius_delta * start_frac
                        r1 = r_p + radius_delta * end_frac
                        branch_segments.append((
                            (r0 * cos_angle, r0 * sin_angle),
                            (r1 * cos_angle, r1 * sin_angle),
                        ))
                        branch_colors.append(state_colors[state])
                else:
                    branch_segments.append((
                        (r_p * cos_angle, r_p * sin_angle),
                        (r_c * cos_angle, r_c * sin_angle),
                    ))
                    branch_colors.append("gray")

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

            # Arcs at internal nodes
            node_states = mapping.get("node_states", {})
            arc_segments = []
            arc_colors = []
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
                # Use the node's sampled state for color if available
                arc_color = "gray"
                if cid in node_states:
                    arc_color = state_colors[node_states[cid]]

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
                arc_colors.append(arc_color)

            if arc_segments:
                ax.add_collection(
                    LineCollection(
                        arc_segments,
                        colors=arc_colors,
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

            # Legend
            handles = [
                Line2D([0], [0], color=state_colors[i], linewidth=3, label=states[i])
                for i in range(k)
            ]
            if self.plot_config.color_file:
                color_legend = build_color_legend_handles(color_data)
                handles.extend(color_legend)
            ax.legend(handles=handles, title="States", loc="upper left",
                      fontsize=8, title_fontsize=9)

            if config.show_title:
                ax.set_title(config.title or "Stochastic Character Map", fontsize=config.title_fontsize)
        else:
            # --- Rectangular mode ---
            vertical_segments = []
            branch_segments = []
            branch_colors = []
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

                # Draw horizontal branch segments colored by history
                clade_id = id(clade)
                if clade_id in histories:
                    history = histories[clade_id]
                    for h_idx in range(len(history)):
                        state = history[h_idx][1]
                        seg_start = history[h_idx][0]
                        if h_idx + 1 < len(history):
                            seg_end = history[h_idx + 1][0]
                        else:
                            seg_end = t
                        x0 = parent_x + seg_start
                        x1 = parent_x + seg_end
                        color = state_colors[state]
                        branch_segments.append(((x0, child_y), (x1, child_y)))
                        branch_colors.append(color)
                else:
                    branch_segments.append(((parent_x, child_y), (child_x, child_y)))
                    branch_colors.append("gray")

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
                        tip.name, va="center", fontsize=label_fontsize
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

            # Legend
            handles = [
                Line2D([0], [0], color=state_colors[i], linewidth=3, label=states[i])
                for i in range(k)
            ]
            if self.plot_config.color_file:
                color_legend = build_color_legend_handles(color_data)
                handles.extend(color_legend)
            ax.legend(
                handles=handles, title="States", loc="upper left",
                fontsize=8, title_fontsize=9
            )

            ax.set_xlabel("Branch length")
            ax.set_yticks([])
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            if config.show_title:
                ax.set_title(config.title or "Stochastic Character Map", fontsize=config.title_fontsize)

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved stochastic character map plot: {output_path}")
