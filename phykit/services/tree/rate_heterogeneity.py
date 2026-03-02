import copy
import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class RateHeterogeneity(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.regime_data_path = parsed["regime_data_path"]
        self.n_sim = parsed["n_sim"]
        self.seed = parsed["seed"]
        self.plot_output = parsed["plot_output"]
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file()
        self._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_values = self._parse_trait_file(self.trait_data_path, tree_tips)
        regime_assignments = self._parse_regime_file(
            self.regime_data_path, tree_tips
        )

        # Use intersection of all three sets
        shared = set(trait_values.keys()) & set(regime_assignments.keys())
        if len(shared) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa among tree, trait, and regime files.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        trait_values = {k: trait_values[k] for k in shared}
        regime_assignments = {k: regime_assignments[k] for k in shared}

        # Prune tree to shared taxa
        tree_copy = copy.deepcopy(tree)
        tip_names_in_tree = [t.name for t in tree_copy.get_terminals()]
        tips_to_prune = [t for t in tip_names_in_tree if t not in shared]
        if tips_to_prune:
            tree_copy = self.prune_tree_using_taxa_list(tree_copy, tips_to_prune)

        ordered_names = sorted(trait_values.keys())
        n = len(ordered_names)
        y = np.array([trait_values[name] for name in ordered_names])

        regimes = sorted(set(regime_assignments.values()))
        k = len(regimes)

        if k < 2:
            raise PhykitUserError(
                ["At least 2 distinct regimes are required."],
                code=2,
            )

        # Assign regimes to branches via Fitch parsimony
        parent_map = self._build_parent_map(tree_copy)
        branch_regimes = self._assign_branch_regimes(
            tree_copy, regime_assignments, parent_map
        )

        # Build per-regime VCV matrices
        per_regime_vcv = self._build_per_regime_vcv(
            tree_copy, ordered_names, regimes, branch_regimes, parent_map
        )

        # Build total VCV
        C_total = sum(per_regime_vcv.values())

        # Fit single-rate model
        sigma2_single, anc_single, ll_single = self._fit_single_rate(y, C_total)
        aic_single = -2.0 * ll_single + 2.0 * 2  # 2 params: sigma2, anc

        # Fit multi-rate model
        sigma2_multi, anc_multi, ll_multi = self._fit_multi_rate(
            y, per_regime_vcv, regimes
        )
        aic_multi = -2.0 * ll_multi + 2.0 * (k + 1)  # k sigma2 + anc

        # LRT
        lrt_stat = 2.0 * (ll_multi - ll_single)
        lrt_stat = max(lrt_stat, 0.0)
        df = k - 1
        chi2_p = float(chi2.sf(lrt_stat, df=df))

        # Parametric bootstrap
        sim_p = None
        n_sim_done = None
        if self.n_sim > 0:
            sim_p, n_sim_done = self._parametric_bootstrap(
                y, C_total, per_regime_vcv, regimes, lrt_stat, self.n_sim, self.seed
            )

        # Plot
        if self.plot_output:
            self._plot_regime_tree(
                tree_copy, branch_regimes, regimes, parent_map, self.plot_output
            )

        # Build sigma2_multi dict
        sigma2_multi_dict = {regimes[i]: float(sigma2_multi[i]) for i in range(k)}

        # Compute R²_regime: 1 - (σ²_multi_weighted / σ²_single)
        regime_tip_counts = {}
        for r_name in regimes:
            regime_tip_counts[r_name] = sum(
                1 for t in regime_assignments.values() if t == r_name
            )
        sig2_weighted = sum(
            (regime_tip_counts[r_name] / n) * sigma2_multi_dict[r_name]
            for r_name in regimes
        )
        if sigma2_single > 0:
            r2_regime = 1.0 - sig2_weighted / sigma2_single
        else:
            r2_regime = float("nan")

        if self.json_output:
            result = self._format_result(
                n_tips=n,
                regimes=regimes,
                sigma2_single=sigma2_single,
                anc_single=anc_single,
                ll_single=ll_single,
                aic_single=aic_single,
                sigma2_multi_dict=sigma2_multi_dict,
                anc_multi=anc_multi,
                ll_multi=ll_multi,
                aic_multi=aic_multi,
                lrt_stat=lrt_stat,
                df=df,
                chi2_p=chi2_p,
                sim_p=sim_p,
                n_sim_done=n_sim_done,
                r2_regime=r2_regime,
            )
            if self.plot_output:
                result["plot_output"] = self.plot_output
            print_json(result)
        else:
            self._print_text_output(
                n_tips=n,
                regimes=regimes,
                sigma2_single=sigma2_single,
                anc_single=anc_single,
                ll_single=ll_single,
                aic_single=aic_single,
                sigma2_multi_dict=sigma2_multi_dict,
                anc_multi=anc_multi,
                ll_multi=ll_multi,
                aic_multi=aic_multi,
                lrt_stat=lrt_stat,
                df=df,
                chi2_p=chi2_p,
                sim_p=sim_p,
                n_sim_done=n_sim_done,
                r2_regime=r2_regime,
            )

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            regime_data_path=args.regime_data,
            n_sim=getattr(args, "nsim", 0),
            seed=getattr(args, "seed", None),
            plot_output=getattr(args, "plot", None),
            json_output=getattr(args, "json", False),
        )

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for rate heterogeneity test."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                raise PhykitUserError(
                    ["All branches in the tree must have lengths."],
                    code=2,
                )

    def _parse_trait_file(
        self, path: str, tree_tips: List[str]
    ) -> Dict[str, float]:
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

        traits = {}
        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise PhykitUserError(
                    [
                        f"Line {line_num} in trait file has {len(parts)} columns; expected 2.",
                        "Each line should be: taxon_name<tab>trait_value",
                    ],
                    code=2,
                )
            taxon, value_str = parts
            try:
                traits[taxon] = float(value_str)
            except ValueError:
                raise PhykitUserError(
                    [
                        f"Non-numeric trait value '{value_str}' for taxon '{taxon}' on line {line_num}.",
                    ],
                    code=2,
                )

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

    def _parse_regime_file(
        self, path: str, tree_tips: List[str]
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

        regimes = {}
        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise PhykitUserError(
                    [
                        f"Line {line_num} in regime file has {len(parts)} columns; expected 2.",
                        "Each line should be: taxon_name<tab>regime_label",
                    ],
                    code=2,
                )
            taxon, regime = parts
            regimes[taxon] = regime

        tree_tip_set = set(tree_tips)
        regime_taxa_set = set(regimes.keys())
        shared = tree_tip_set & regime_taxa_set

        tree_only = tree_tip_set - regime_taxa_set
        regime_only = regime_taxa_set - tree_tip_set

        if tree_only:
            print(
                f"Warning: {len(tree_only)} taxa in tree but not in regime file: "
                f"{', '.join(sorted(tree_only))}",
                file=sys.stderr,
            )
        if regime_only:
            print(
                f"Warning: {len(regime_only)} taxa in regime file but not in tree: "
                f"{', '.join(sorted(regime_only))}",
                file=sys.stderr,
            )

        if len(shared) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa between tree and regime file.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        return {taxon: regimes[taxon] for taxon in shared}

    def _build_parent_map(self, tree) -> Dict:
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    def _assign_branch_regimes(
        self, tree, tip_regimes: Dict[str, str], parent_map: Dict
    ) -> Dict:
        """Assign regime labels to all branches using Fitch parsimony.

        Returns dict mapping id(clade) -> regime_label for every non-root clade.
        """
        # Postorder pass: build state sets
        state_sets = {}
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                if clade.name in tip_regimes:
                    state_sets[id(clade)] = {tip_regimes[clade.name]}
                else:
                    state_sets[id(clade)] = set()
            else:
                child_sets = [
                    state_sets[id(c)] for c in clade.clades
                    if id(c) in state_sets and state_sets[id(c)]
                ]
                if not child_sets:
                    state_sets[id(clade)] = set()
                else:
                    intersection = child_sets[0]
                    for cs in child_sets[1:]:
                        intersection = intersection & cs
                    if intersection:
                        state_sets[id(clade)] = intersection
                    else:
                        union = set()
                        for cs in child_sets:
                            union = union | cs
                        state_sets[id(clade)] = union

        # Preorder pass: resolve ambiguities
        node_regimes = {}
        root = tree.root
        root_set = state_sets.get(id(root), set())
        if root_set:
            node_regimes[id(root)] = sorted(root_set)[0]
        else:
            node_regimes[id(root)] = sorted(set(tip_regimes.values()))[0]

        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            clade_set = state_sets.get(id(clade), set())
            if not clade_set:
                # Fallback: use parent regime
                parent = parent_map.get(id(clade))
                if parent and id(parent) in node_regimes:
                    node_regimes[id(clade)] = node_regimes[id(parent)]
                continue

            parent = parent_map.get(id(clade))
            if parent and id(parent) in node_regimes:
                parent_regime = node_regimes[id(parent)]
                if parent_regime in clade_set:
                    node_regimes[id(clade)] = parent_regime
                else:
                    node_regimes[id(clade)] = sorted(clade_set)[0]
            else:
                node_regimes[id(clade)] = sorted(clade_set)[0]

        # Branch regime = regime of the child node
        branch_regimes = {}
        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            branch_regimes[id(clade)] = node_regimes.get(id(clade), "unknown")

        return branch_regimes

    def _build_root_to_tip_paths(
        self, tree, ordered_names: List[str], parent_map: Dict
    ) -> Dict[str, List]:
        """Build root-to-tip paths for each tip.

        Returns dict mapping tip_name -> list of (clade_id, branch_length)
        tuples from root to tip.
        """
        # Map tip names to clade objects
        tip_map = {}
        for tip in tree.get_terminals():
            if tip.name in ordered_names:
                tip_map[tip.name] = tip

        paths = {}
        for name in ordered_names:
            tip = tip_map[name]
            path = []
            current = tip
            while id(current) in parent_map:
                bl = current.branch_length if current.branch_length else 0.0
                path.append((id(current), bl))
                current = parent_map[id(current)]
            path.reverse()
            paths[name] = path

        return paths

    def _build_per_regime_vcv(
        self, tree, ordered_names: List[str], regimes: List[str],
        branch_regimes: Dict, parent_map: Dict,
    ) -> Dict[str, np.ndarray]:
        """Build separate VCV matrix for each regime.

        C_r[i,j] = portion of shared path in regime r between tips i and j.
        """
        n = len(ordered_names)
        paths = self._build_root_to_tip_paths(tree, ordered_names, parent_map)

        per_regime_vcv = {r: np.zeros((n, n)) for r in regimes}

        for i in range(n):
            path_i = paths[ordered_names[i]]
            # Diagonal: sum branch lengths where regime == r
            for clade_id, bl in path_i:
                r = branch_regimes.get(clade_id, regimes[0])
                per_regime_vcv[r][i, i] += bl

            for j in range(i + 1, n):
                path_j = paths[ordered_names[j]]
                # Find shared prefix (matching clade_ids from root)
                min_len = min(len(path_i), len(path_j))
                for s in range(min_len):
                    if path_i[s][0] == path_j[s][0]:
                        clade_id = path_i[s][0]
                        bl = path_i[s][1]
                        r = branch_regimes.get(clade_id, regimes[0])
                        per_regime_vcv[r][i, j] += bl
                        per_regime_vcv[r][j, i] += bl
                    else:
                        break

        return per_regime_vcv

    def _fit_single_rate(
        self, y: np.ndarray, C_total: np.ndarray
    ) -> Tuple[float, float, float]:
        """Fit single-rate BM model.

        Returns (sigma2, ancestral_state, log_likelihood).
        """
        n = len(y)
        try:
            C_inv = np.linalg.inv(C_total)
        except np.linalg.LinAlgError:
            raise PhykitUserError(
                [
                    "Singular VCV matrix: cannot invert.",
                    "Check that the tree has valid branch lengths.",
                ],
                code=2,
            )

        ones = np.ones(n)
        a_hat = float(ones @ C_inv @ y) / float(ones @ C_inv @ ones)
        e = y - a_hat
        sigma2 = float(e @ C_inv @ e) / n

        sign, logdet = np.linalg.slogdet(C_total)
        if sign <= 0 or sigma2 <= 0:
            return sigma2, a_hat, float("-inf")

        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sigma2) + logdet + n)
        return float(sigma2), float(a_hat), float(ll)

    def _fit_multi_rate(
        self, y: np.ndarray, per_regime_vcv: Dict[str, np.ndarray],
        regimes: List[str],
    ) -> Tuple[np.ndarray, float, float]:
        """Fit multi-rate BM model by optimizing per-regime sigma2 values.

        Returns (sigma2_array, ancestral_state, log_likelihood).
        """
        n = len(y)
        k = len(regimes)
        regime_matrices = [per_regime_vcv[r] for r in regimes]

        def neg_log_likelihood(log_sigma2s):
            sigma2s = np.exp(log_sigma2s)
            V = np.zeros((n, n))
            for i_r in range(k):
                V += sigma2s[i_r] * regime_matrices[i_r]

            try:
                sign, logdet = np.linalg.slogdet(V)
                if sign <= 0:
                    return 1e20
                V_inv = np.linalg.inv(V)
            except np.linalg.LinAlgError:
                return 1e20

            ones = np.ones(n)
            denom = float(ones @ V_inv @ ones)
            if denom == 0:
                return 1e20
            a_hat = float(ones @ V_inv @ y) / denom
            e = y - a_hat
            quad = float(e @ V_inv @ e)

            ll = -0.5 * (n * np.log(2 * np.pi) + logdet + quad)
            return -ll

        # Multi-start optimization
        starting_scales = [0.001, 0.01, 0.1, 1.0, 10.0]
        best_negll = np.inf
        best_params = np.zeros(k)

        for sv in starting_scales:
            x0 = np.log(np.ones(k) * sv)
            try:
                result = minimize(
                    neg_log_likelihood, x0, method="L-BFGS-B",
                    options={"maxiter": 5000},
                )
                if result.fun < best_negll:
                    best_negll = result.fun
                    best_params = result.x
            except (ValueError, np.linalg.LinAlgError):
                continue

        # Refine with Nelder-Mead
        try:
            result = minimize(
                neg_log_likelihood, best_params, method="Nelder-Mead",
                options={"maxiter": 10000, "xatol": 1e-10, "fatol": 1e-10},
            )
            if result.fun < best_negll:
                best_negll = result.fun
                best_params = result.x
        except (ValueError, np.linalg.LinAlgError):
            pass

        sigma2s = np.exp(best_params)

        # Compute final ancestral state and log-likelihood
        V = np.zeros((n, n))
        for i_r in range(k):
            V += sigma2s[i_r] * regime_matrices[i_r]

        try:
            V_inv = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            return sigma2s, 0.0, float("-inf")

        ones = np.ones(n)
        a_hat = float(ones @ V_inv @ y) / float(ones @ V_inv @ ones)
        e = y - a_hat
        sign, logdet = np.linalg.slogdet(V)
        quad = float(e @ V_inv @ e)
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet + quad)

        return sigma2s, float(a_hat), float(ll)

    def _parametric_bootstrap(
        self, y: np.ndarray, C_total: np.ndarray,
        per_regime_vcv: Dict[str, np.ndarray], regimes: List[str],
        observed_lrt: float, n_sim: int, seed=None,
    ) -> Tuple[float, int]:
        """Parametric bootstrap p-value for the LRT.

        Simulate data under the single-rate null, refit both models,
        compute LRT, and count exceedances.
        """
        n = len(y)
        rng = np.random.default_rng(seed)

        # Fit null model to get parameters
        sigma2_null, anc_null, _ = self._fit_single_rate(y, C_total)
        if sigma2_null <= 0:
            return 1.0, n_sim

        # Cholesky for simulation
        C_scaled = sigma2_null * C_total
        try:
            L = np.linalg.cholesky(C_scaled)
        except np.linalg.LinAlgError:
            # Add small jitter
            C_scaled += np.eye(n) * 1e-8
            L = np.linalg.cholesky(C_scaled)

        count = 0
        for _ in range(n_sim):
            z = rng.standard_normal(n)
            y_sim = anc_null + L @ z

            _, _, ll_single_sim = self._fit_single_rate(y_sim, C_total)
            _, _, ll_multi_sim = self._fit_multi_rate(
                y_sim, per_regime_vcv, regimes
            )
            lrt_sim = 2.0 * (ll_multi_sim - ll_single_sim)
            lrt_sim = max(lrt_sim, 0.0)
            if lrt_sim >= observed_lrt:
                count += 1

        p_value = (count + 1) / (n_sim + 1)
        return float(p_value), n_sim

    def _plot_regime_tree(
        self, tree, branch_regimes: Dict, regimes: List[str],
        parent_map: Dict, output_path: str,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.lines import Line2D
        except ImportError:
            print(
                "matplotlib is required for regime tree plotting. "
                "Install matplotlib and retry."
            )
            raise SystemExit(2)

        k = len(regimes)
        cmap = plt.get_cmap("tab10")
        regime_colors = {
            regimes[i]: cmap(i / max(k - 1, 1)) for i in range(k)
        }

        node_x = {}
        node_y = {}
        tips = list(tree.get_terminals())

        for i, tip in enumerate(tips):
            node_y[id(tip)] = i

        root = tree.root
        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                node_x[id(clade)] = 0.0
            elif id(clade) in parent_map:
                parent = parent_map[id(clade)]
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

        fig, ax = plt.subplots(figsize=(10, max(4, len(tips) * 0.4)))

        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            if id(clade) not in parent_map:
                continue

            parent = parent_map[id(clade)]
            parent_x = node_x[id(parent)]
            parent_y = node_y[id(parent)]
            child_y = node_y[id(clade)]
            child_x = node_x[id(clade)]

            regime = branch_regimes.get(id(clade), regimes[0])
            color = regime_colors[regime]

            # Vertical connector
            ax.plot(
                [parent_x, parent_x], [parent_y, child_y],
                color="gray", linewidth=0.8, zorder=1,
            )

            # Horizontal branch
            ax.plot(
                [parent_x, child_x], [child_y, child_y],
                color=color, linewidth=2.5, solid_capstyle="butt", zorder=2,
            )

        max_x = max(node_x.values()) if node_x else 0
        offset = max_x * 0.02
        for tip in tips:
            ax.text(
                node_x[id(tip)] + offset, node_y[id(tip)],
                tip.name, va="center", fontsize=9,
            )

        handles = [
            Line2D([0], [0], color=regime_colors[r], linewidth=3, label=r)
            for r in regimes
        ]
        ax.legend(
            handles=handles, title="Regimes", loc="upper left",
            fontsize=8, title_fontsize=9,
        )

        ax.set_xlabel("Branch length")
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.set_title("Regime Tree (Rate Heterogeneity)")
        fig.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved regime tree plot: {output_path}")

    def _print_text_output(
        self, *, n_tips, regimes, sigma2_single, anc_single, ll_single,
        aic_single, sigma2_multi_dict, anc_multi, ll_multi, aic_multi,
        lrt_stat, df, chi2_p, sim_p, n_sim_done, r2_regime,
    ) -> None:
        print("Rate Heterogeneity Test (Multi-rate Brownian Motion)")

        print(f"\nRegimes: {len(regimes)} ({', '.join(regimes)})")
        print(f"Number of tips: {n_tips}")

        print("\nSingle-rate model (H0):")
        print(f"  Sigma-squared:      {sigma2_single:.4f}")
        print(f"  Ancestral state:    {anc_single:.4f}")
        print(f"  Log-likelihood:     {ll_single:.2f}")
        print(f"  AIC:                {aic_single:.2f}")

        print("\nMulti-rate model (H1):")
        print(f"  {'Regime':<20s}{'Sigma-squared':>14s}")
        for r in regimes:
            print(f"  {r:<20s}{sigma2_multi_dict[r]:>14.4f}")
        print(f"  Ancestral state:    {anc_multi:.4f}")
        print(f"  Log-likelihood:     {ll_multi:.2f}")
        print(f"  AIC:                {aic_multi:.2f}")

        print("\nLikelihood ratio test:")
        print(f"  LRT statistic:      {lrt_stat:.4f}")
        print(f"  Degrees of freedom: {df}")
        print(f"  Chi-squared p-value: {chi2_p:.4f}")
        if sim_p is not None:
            print(f"  Simulated p-value:  {sim_p:.4f} ({n_sim_done} simulations)")

        print(f"\nEffect size:")
        print(f"  R2_regime: {r2_regime:.4f}")

    def _format_result(
        self, *, n_tips, regimes, sigma2_single, anc_single, ll_single,
        aic_single, sigma2_multi_dict, anc_multi, ll_multi, aic_multi,
        lrt_stat, df, chi2_p, sim_p, n_sim_done, r2_regime,
    ) -> Dict:
        return {
            "n_tips": n_tips,
            "regimes": regimes,
            "single_rate": {
                "sigma2": float(sigma2_single),
                "ancestral_state": float(anc_single),
                "log_likelihood": float(ll_single),
                "aic": float(aic_single),
            },
            "multi_rate": {
                "sigma2": sigma2_multi_dict,
                "ancestral_state": float(anc_multi),
                "log_likelihood": float(ll_multi),
                "aic": float(aic_multi),
            },
            "lrt": {
                "statistic": float(lrt_stat),
                "df": df,
                "chi2_p_value": float(chi2_p),
                "sim_p_value": float(sim_p) if sim_p is not None else None,
                "n_sim": n_sim_done,
            },
            "r_squared_regime": float(r2_regime),
        }
