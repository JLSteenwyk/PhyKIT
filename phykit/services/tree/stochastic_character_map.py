import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


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
        )

    def run(self) -> None:
        tree = self.read_tree_file()
        self._validate_tree(tree)

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

        # Run stochastic mappings
        rng = np.random.default_rng(self.seed)
        mappings = []
        for _ in range(self.nsim):
            mapping = self._run_single_simulation(
                tree, tip_states, Q, pi, states, rng,
                cond_liks=cond_liks, parent_map=parent_map,
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

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for stochastic character mapping."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                raise PhykitUserError(
                    ["All branches in the tree must have lengths."],
                    code=2,
                )

    def _parse_discrete_trait_file(
        self, path: str, column: str, tree_tips: List[str]
    ) -> Dict[str, str]:
        try:
            with open(path) as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        data_lines = []
        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            data_lines.append(stripped)

        if len(data_lines) < 2:
            raise PhykitUserError(
                ["Trait file must have a header row and at least one data row."],
                code=2,
            )

        header_parts = data_lines[0].split("\t")
        if len(header_parts) < 2:
            raise PhykitUserError(
                ["Header must have at least 2 columns (taxon + at least 1 trait)."],
                code=2,
            )

        trait_names = header_parts[1:]
        if column not in trait_names:
            raise PhykitUserError(
                [
                    f"Column '{column}' not found in trait file.",
                    f"Available columns: {', '.join(trait_names)}",
                ],
                code=2,
            )

        col_idx = trait_names.index(column)
        n_cols = len(header_parts)

        traits = {}
        for line_idx, line in enumerate(data_lines[1:], 2):
            parts = line.split("\t")
            if len(parts) != n_cols:
                raise PhykitUserError(
                    [
                        f"Line {line_idx} has {len(parts)} columns; expected {n_cols}.",
                    ],
                    code=2,
                )
            taxon = parts[0]
            value = parts[1 + col_idx].strip()
            if not value:
                raise PhykitUserError(
                    [f"Missing trait value for taxon '{taxon}' on line {line_idx}."],
                    code=2,
                )
            traits[taxon] = value

        tree_tip_set = set(tree_tips)
        trait_taxa_set = set(traits.keys())
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

    def _build_q_matrix(
        self, params: np.ndarray, k: int, model: str
    ) -> np.ndarray:
        Q = np.zeros((k, k))
        if model == "ER":
            rate = params[0]
            Q[:] = rate
            np.fill_diagonal(Q, 0.0)
        elif model == "SYM":
            idx = 0
            for i in range(k):
                for j in range(i + 1, k):
                    Q[i, j] = params[idx]
                    Q[j, i] = params[idx]
                    idx += 1
        elif model == "ARD":
            idx = 0
            for i in range(k):
                for j in range(k):
                    if i != j:
                        Q[i, j] = params[idx]
                        idx += 1
        # Set diagonal
        for i in range(k):
            Q[i, i] = -np.sum(Q[i, :])
        return Q

    def _matrix_exp(self, Q: np.ndarray, t: float) -> np.ndarray:
        return expm(Q * t)

    def _felsenstein_pruning(
        self, tree, tip_states: Dict[str, str], Q: np.ndarray,
        pi: np.ndarray, states: List[str]
    ) -> Tuple[Dict, float]:
        k = len(states)
        state_idx = {s: i for i, s in enumerate(states)}
        cond_liks = {}

        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                lik = np.zeros(k)
                if clade.name in tip_states:
                    lik[state_idx[tip_states[clade.name]]] = 1.0
                cond_liks[id(clade)] = lik
            else:
                lik = np.ones(k)
                for child in clade.clades:
                    t = child.branch_length if child.branch_length else 1e-8
                    P = self._matrix_exp(Q, t)
                    child_lik = cond_liks[id(child)]
                    lik *= P @ child_lik
                cond_liks[id(clade)] = lik

        root_lik = cond_liks[id(tree.root)]
        total_lik = np.sum(pi * root_lik)
        if total_lik <= 0:
            loglik = -1e20
        else:
            loglik = np.log(total_lik)

        return cond_liks, loglik

    def _fit_q_matrix(
        self, tree, tip_states: Dict[str, str],
        states: List[str], model: str
    ) -> Tuple[np.ndarray, float]:
        k = len(states)

        if model == "ER":
            n_params = 1
        elif model == "SYM":
            n_params = k * (k - 1) // 2
        elif model == "ARD":
            n_params = k * (k - 1)
        else:
            raise PhykitUserError(
                [f"Unknown model '{model}'. Use ER, SYM, or ARD."],
                code=2,
            )

        pi = np.ones(k) / k

        def neg_loglik(params):
            Q = self._build_q_matrix(np.abs(params), k, model)
            _, ll = self._felsenstein_pruning(tree, tip_states, Q, pi, states)
            return -ll

        bounds = [(1e-8, 100.0)] * n_params

        # Multi-start optimization for robustness
        starting_values = [0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0]
        best_negll = np.inf
        best_params = np.ones(n_params) * 0.1

        for sv in starting_values:
            x0 = np.ones(n_params) * sv
            for method in ["L-BFGS-B", "Nelder-Mead"]:
                try:
                    kwargs = {"method": method}
                    if method == "L-BFGS-B":
                        kwargs["bounds"] = bounds
                    result = minimize(neg_loglik, x0, **kwargs)
                    if result.fun < best_negll:
                        best_negll = result.fun
                        best_params = np.abs(result.x)
                except (ValueError, np.linalg.LinAlgError):
                    continue

        # Refine best result with Nelder-Mead
        try:
            result = minimize(
                neg_loglik, best_params, method="Nelder-Mead",
                options={"maxiter": 10000, "xatol": 1e-10, "fatol": 1e-10},
            )
            if result.fun < best_negll:
                best_params = np.abs(result.x)
        except (ValueError, np.linalg.LinAlgError):
            pass

        Q = self._build_q_matrix(best_params, k, model)
        _, loglik = self._felsenstein_pruning(tree, tip_states, Q, pi, states)

        return Q, loglik

    def _sample_ancestral_states(
        self, tree, tip_states: Dict[str, str],
        Q: np.ndarray, pi: np.ndarray,
        cond_liks: Dict, states: List[str], rng,
        parent_map: Dict = None,
    ) -> Dict:
        k = len(states)
        state_idx = {s: i for i, s in enumerate(states)}
        node_states = {}

        # Sample root
        root_lik = cond_liks[id(tree.root)]
        root_probs = pi * root_lik
        root_probs /= root_probs.sum()
        root_state_idx = rng.choice(k, p=root_probs)
        node_states[id(tree.root)] = root_state_idx

        # Pre-order traversal
        for clade in tree.find_clades(order="preorder"):
            if clade == tree.root:
                continue
            parent = self._get_parent(tree, clade, parent_map)
            if parent is None:
                continue
            parent_state = node_states[id(parent)]

            t = clade.branch_length if clade.branch_length else 1e-8

            if clade.is_terminal() and clade.name in tip_states:
                node_states[id(clade)] = state_idx[tip_states[clade.name]]
            else:
                P = self._matrix_exp(Q, t)
                child_lik = cond_liks[id(clade)]
                probs = P[parent_state, :] * child_lik
                total = probs.sum()
                if total > 0:
                    probs /= total
                else:
                    probs = np.ones(k) / k
                node_states[id(clade)] = rng.choice(k, p=probs)

        return node_states

    def _build_parent_map(self, tree) -> Dict:
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    def _get_parent(self, tree, target, parent_map=None):
        if parent_map is not None:
            return parent_map.get(id(target))
        for clade in tree.find_clades(order="preorder"):
            if target in clade.clades:
                return clade
        return None

    def _simulate_branch_history(
        self, Q: np.ndarray, start_state: int, end_state: int,
        branch_length: float, k: int, rng, max_attempts: int = 1000
    ) -> List[Tuple[float, int]]:
        if branch_length <= 0:
            return [(0.0, start_state)]

        rates = -np.diag(Q)

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
                trans_probs = Q[current_state, :].copy()
                trans_probs[current_state] = 0.0
                total = trans_probs.sum()
                if total <= 0:
                    break
                trans_probs /= total
                new_state = rng.choice(k, p=trans_probs)
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

    def _run_single_simulation(
        self, tree, tip_states: Dict[str, str],
        Q: np.ndarray, pi: np.ndarray,
        states: List[str], rng,
        cond_liks: Dict = None, parent_map: Dict = None,
    ) -> Dict:
        k = len(states)
        state_idx = {s: i for i, s in enumerate(states)}

        if cond_liks is None:
            cond_liks, _ = self._felsenstein_pruning(
                tree, tip_states, Q, pi, states
            )
        if parent_map is None:
            parent_map = self._build_parent_map(tree)
        node_states = self._sample_ancestral_states(
            tree, tip_states, Q, pi, cond_liks, states, rng,
            parent_map=parent_map,
        )

        branch_histories = {}
        for clade in tree.find_clades(order="preorder"):
            if clade == tree.root:
                continue
            parent = self._get_parent(tree, clade, parent_map)
            if parent is None:
                continue

            t = clade.branch_length if clade.branch_length else 1e-8
            start_state = node_states[id(parent)]
            end_state = node_states[id(clade)]

            history = self._simulate_branch_history(
                Q, start_state, end_state, t, k, rng
            )
            branch_histories[id(clade)] = history

        return {
            "node_states": node_states,
            "branch_histories": branch_histories,
        }

    def _summarize_simulations(
        self, mappings: List[Dict], states: List[str], tree
    ) -> Dict:
        k = len(states)
        nsim = len(mappings)

        total_dwelling = np.zeros((nsim, k))
        total_transitions = np.zeros((nsim, k, k))
        node_posteriors = {}

        # Collect node ids for internal nodes
        internal_nodes = []
        for clade in tree.find_clades(order="preorder"):
            if not clade.is_terminal():
                internal_nodes.append(id(clade))

        for sim_idx, mapping in enumerate(mappings):
            histories = mapping["branch_histories"]
            node_states = mapping["node_states"]

            for clade in tree.find_clades(order="preorder"):
                if clade == tree.root:
                    continue
                clade_id = id(clade)
                if clade_id not in histories:
                    continue

                history = histories[clade_id]
                t = clade.branch_length if clade.branch_length else 1e-8

                # Compute dwelling times for this branch
                for h_idx in range(len(history)):
                    state = history[h_idx][1]
                    start_t = history[h_idx][0]
                    if h_idx + 1 < len(history):
                        end_t = history[h_idx + 1][0]
                    else:
                        end_t = t
                    total_dwelling[sim_idx, state] += end_t - start_t

                # Count transitions
                for h_idx in range(1, len(history)):
                    from_state = history[h_idx - 1][1]
                    to_state = history[h_idx][1]
                    if from_state != to_state:
                        total_transitions[sim_idx, from_state, to_state] += 1

            # Node posteriors
            for node_id in internal_nodes:
                if node_id not in node_posteriors:
                    node_posteriors[node_id] = np.zeros(k)
                if node_id in node_states:
                    state = node_states[node_id]
                    node_posteriors[node_id][state] += 1

        # Average across simulations
        mean_dwelling = total_dwelling.mean(axis=0)
        mean_transitions = total_transitions.mean(axis=0)

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
        states: List[str], summary: Dict
    ) -> None:
        k = len(states)
        print("Stochastic Character Mapping (SIMMAP)")
        print(f"\nModel: {self.model}")
        print(f"Number of simulations: {self.nsim}")

        print("\nFitted Q matrix:")
        # Header
        header = "              " + "".join(f"{s:>12s}" for s in states)
        print(header)
        for i, s in enumerate(states):
            row = f"{s:14s}" + "".join(f"{Q[i, j]:12.4f}" for j in range(k))
            print(row)

        print(f"\nLog-likelihood: {loglik:.4f}")

        mean_dwelling = summary["mean_dwelling_times"]
        total_dwelling = mean_dwelling.sum()
        print("\nMean dwelling times:")
        for i, s in enumerate(states):
            pct = 100.0 * mean_dwelling[i] / total_dwelling if total_dwelling > 0 else 0
            print(f"  {s:14s} {mean_dwelling[i]:8.2f} ({pct:.1f}%)")

        mean_trans = summary["mean_transitions"]
        print("\nMean transitions:")
        total_trans = 0.0
        for i in range(k):
            for j in range(k):
                if i != j and mean_trans[i, j] > 0:
                    print(
                        f"  {states[i]} -> {states[j]}:"
                        f"  {mean_trans[i, j]:.2f}"
                    )
                    total_trans += mean_trans[i, j]
        print(f"  {'Total:':28s}{total_trans:.2f}")

    def _format_result(
        self, Q: np.ndarray, loglik: float,
        states: List[str], summary: Dict,
        plot_output: str
    ) -> Dict:
        k = len(states)
        q_dict = {}
        for i, si in enumerate(states):
            q_dict[si] = {}
            for j, sj in enumerate(states):
                q_dict[si][sj] = float(Q[i, j])

        mean_dwelling = summary["mean_dwelling_times"]
        total_dwelling = float(mean_dwelling.sum())

        dwelling_dict = {
            states[i]: float(mean_dwelling[i]) for i in range(k)
        }
        proportion_dict = {
            states[i]: float(mean_dwelling[i] / total_dwelling)
            if total_dwelling > 0 else 0.0
            for i in range(k)
        }

        mean_trans = summary["mean_transitions"]
        trans_dict = {}
        total_trans = 0.0
        for i in range(k):
            for j in range(k):
                if i != j:
                    key = f"{states[i]} -> {states[j]}"
                    trans_dict[key] = float(mean_trans[i, j])
                    total_trans += float(mean_trans[i, j])

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
        self, tree, mapping: Dict, states: List[str],
        output_path: str, parent_map: Dict = None,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
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

        # Compute node positions: x = distance from root, y = assigned by tips
        node_x = {}
        node_y = {}

        tips = list(tree.get_terminals())
        for i, tip in enumerate(tips):
            node_y[id(tip)] = i

        # Compute x positions (distance from root)
        root = tree.root
        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                node_x[id(clade)] = 0.0
            else:
                parent = self._get_parent(tree, clade, parent_map)
                if parent is not None:
                    t = clade.branch_length if clade.branch_length else 0.0
                    node_x[id(clade)] = node_x[id(parent)] + t

        # Compute y positions for internal nodes (average of children)
        for clade in tree.find_clades(order="postorder"):
            if not clade.is_terminal() and id(clade) not in node_y:
                child_ys = [
                    node_y[id(c)] for c in clade.clades if id(c) in node_y
                ]
                if child_ys:
                    node_y[id(clade)] = np.mean(child_ys)
                else:
                    node_y[id(clade)] = 0.0

        fig, ax = plt.subplots(figsize=(10, max(4, len(tips) * 0.4)))

        histories = mapping["branch_histories"]

        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            parent = self._get_parent(tree, clade, parent_map)
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
                    ax.plot(
                        [x0, x1], [child_y, child_y],
                        color=color, linewidth=2.5, solid_capstyle="butt",
                        zorder=2
                    )
            else:
                ax.plot(
                    [parent_x, child_x], [child_y, child_y],
                    color="gray", linewidth=2.5, zorder=2
                )

        # Tip labels
        max_x = max(node_x.values()) if node_x else 0
        offset = max_x * 0.02
        for tip in tips:
            ax.text(
                node_x[id(tip)] + offset, node_y[id(tip)],
                tip.name, va="center", fontsize=9
            )

        # Legend
        handles = [
            Line2D([0], [0], color=state_colors[i], linewidth=3, label=states[i])
            for i in range(k)
        ]
        ax.legend(
            handles=handles, title="States", loc="upper left",
            fontsize=8, title_fontsize=9
        )

        ax.set_xlabel("Branch length")
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.set_title("Stochastic Character Map")
        fig.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved stochastic character map plot: {output_path}")
