from __future__ import annotations

import math
import sys
from collections import Counter

from .base import Tree
from ...errors import PhykitUserError


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()
_CHO_FACTOR = None
_CHO_SOLVE = None
_MINIMIZE = None


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def cho_factor(*args, **kwargs):
    global _CHO_FACTOR

    if _CHO_FACTOR is None:
        from scipy.linalg import cho_factor as _cho_factor

        _CHO_FACTOR = _cho_factor

    return _CHO_FACTOR(*args, **kwargs)


def cho_solve(*args, **kwargs):
    global _CHO_SOLVE

    if _CHO_SOLVE is None:
        from scipy.linalg import cho_solve as _cho_solve

        _CHO_SOLVE = _cho_solve

    return _CHO_SOLVE(*args, **kwargs)


def minimize(*args, **kwargs):
    global _MINIMIZE

    if _MINIMIZE is None:
        from scipy.optimize import minimize as _minimize

        _MINIMIZE = _minimize

    return _MINIMIZE(*args, **kwargs)


def _chi2_sf(lrt_stat: float, df: int) -> float:
    if df == 1:
        return math.erfc(math.sqrt(lrt_stat / 2.0))
    if df == 2:
        return math.exp(-lrt_stat / 2.0)

    from scipy.stats import chi2
    return float(chi2.sf(lrt_stat, df=df))


def _merge_nonempty_child_state_sets(state_sets, children):
    intersection = None
    union = None
    for child in children:
        child_set = state_sets.get(id(child))
        if not child_set:
            continue
        if intersection is None:
            intersection = child_set
            union = child_set
        else:
            intersection = intersection & child_set
            union = union | child_set

    if intersection is None:
        return set()
    return intersection if intersection else union


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
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        self.validate_tree(tree, min_tips=3, require_branch_lengths=True, context="rate heterogeneity test")

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_values = self._parse_trait_file(self.trait_data_path, tree_tips)
        regime_assignments = self._parse_regime_file(
            self.regime_data_path, tree_tips
        )

        # Use intersection of all three sets
        shared = set(trait_values) & set(regime_assignments)
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

        tips_to_prune = [t for t in tree_tips if t not in shared]
        needs_working_copy = bool(tips_to_prune) or bool(
            self.plot_output and self.plot_config.ladderize
        )
        tree_for_analysis = self._fast_copy(tree) if needs_working_copy else tree
        if tips_to_prune:
            tree_for_analysis = self.prune_tree_using_taxa_list(
                tree_for_analysis, tips_to_prune
            )

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

        preorder_clades = list(self._iter_preorder(tree_for_analysis.root))
        postorder_clades = list(reversed(preorder_clades))

        # Assign regimes to branches via Fitch parsimony
        parent_map = self._build_parent_map(tree_for_analysis, preorder_clades)
        branch_regimes = self._assign_branch_regimes(
            tree_for_analysis, regime_assignments, parent_map,
            preorder_clades, postorder_clades,
        )

        # Build per-regime VCV matrices
        per_regime_vcv = self._build_per_regime_vcv(
            tree_for_analysis, ordered_names, regimes, branch_regimes, parent_map,
            preorder_clades,
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
        chi2_p = _chi2_sf(lrt_stat, df=df)

        # Parametric bootstrap
        sim_p = None
        n_sim_done = None
        if self.n_sim > 0:
            sim_p, n_sim_done = self._parametric_bootstrap(
                y, C_total, per_regime_vcv, regimes, lrt_stat, self.n_sim, self.seed
            )

        # Plot
        if self.plot_output:
            if self.plot_config.ladderize:
                tree_for_analysis.ladderize()
                parent_map = self._build_parent_map(tree_for_analysis)
            self._plot_regime_tree(
                tree_for_analysis, branch_regimes, regimes, parent_map, self.plot_output
            )

        # Build sigma2_multi dict
        sigma2_multi_dict = {regimes[i]: float(sigma2_multi[i]) for i in range(k)}

        # Compute R²_regime: 1 - (σ²_multi_weighted / σ²_single)
        regime_tip_counts = self._count_regime_tips(regime_assignments, regimes)
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

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            regime_data_path=args.regime_data,
            n_sim=getattr(args, "nsim", 0),
            seed=getattr(args, "seed", None),
            plot_output=getattr(args, "plot", None),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    @staticmethod
    def _count_regime_tips(
        regime_assignments: dict[str, str],
        regimes: list[str],
    ) -> dict[str, int]:
        counts = Counter(regime_assignments.values())
        return {regime: counts[regime] for regime in regimes}

    def _parse_trait_file(
        self, path: str, tree_tips: list[str]
    ) -> dict[str, float]:
        try:
            traits = {}
            with open(path) as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t", 2)
                    if len(parts) != 2:
                        column_count = line.count("\t") + 1
                        raise PhykitUserError(
                            [
                                f"Line {line_num} in trait file has {column_count} columns; expected 2.",
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
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        tree_tip_set = set(tree_tips)
        if (
            len(tree_tip_set) >= 3
            and len(tree_tip_set) == len(traits)
            and tree_tip_set == traits.keys()
        ):
            return traits

        trait_taxa_set = set(traits)
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
        self, path: str, tree_tips: list[str]
    ) -> dict[str, str]:
        try:
            regimes = {}
            with open(path, "rb") as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line[0] == 35:
                        continue
                    parts = line.split(b"\t", 2)
                    if len(parts) != 2:
                        column_count = line.count(b"\t") + 1
                        raise PhykitUserError(
                            [
                                f"Line {line_num} in regime file has {column_count} columns; expected 2.",
                                "Each line should be: taxon_name<tab>regime_label",
                            ],
                            code=2,
                        )
                    taxon, regime = parts
                    taxon = taxon.decode()
                    regime = regime.decode()
                    regimes[taxon] = regime
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        tree_tip_set = set(tree_tips)
        if (
            len(tree_tip_set) >= 3
            and len(tree_tip_set) == len(regimes)
            and tree_tip_set == regimes.keys()
        ):
            return dict(regimes)

        regime_taxa_set = set(regimes)
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

    @staticmethod
    def _iter_postorder(root):
        clades = []
        stack = [root]
        pop = stack.pop
        append_clade = clades.append
        extend = stack.extend
        while stack:
            clade = pop()
            append_clade(clade)
            children = clade.clades
            if children:
                extend(children)
        yield from reversed(clades)

    def _build_parent_map(self, tree, preorder_clades=None) -> dict:
        parent_map = {}
        if preorder_clades is not None:
            for clade in preorder_clades:
                for child in clade.clades:
                    parent_map[id(child)] = clade
            return parent_map

        stack = [tree.root]
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

    def _assign_branch_regimes(
        self, tree, tip_regimes: dict[str, str], parent_map: dict,
        preorder_clades=None, postorder_clades=None,
    ) -> dict:
        """Assign regime labels to all branches using Fitch parsimony.

        Returns dict mapping id(clade) -> regime_label for every non-root clade.
        """
        if preorder_clades is None:
            preorder_clades = list(self._iter_preorder(tree.root))
        if postorder_clades is None:
            postorder_clades = reversed(preorder_clades)

        # Postorder pass: build state sets
        state_sets = {}
        for clade in postorder_clades:
            if clade.is_terminal():
                if clade.name in tip_regimes:
                    state_sets[id(clade)] = {tip_regimes[clade.name]}
                else:
                    state_sets[id(clade)] = set()
            else:
                children = clade.clades
                if len(children) == 2:
                    left = state_sets.get(id(children[0]))
                    right = state_sets.get(id(children[1]))
                    if left and right:
                        intersection = left & right
                        state_sets[id(clade)] = (
                            intersection if intersection else left | right
                        )
                    elif left:
                        state_sets[id(clade)] = left
                    elif right:
                        state_sets[id(clade)] = right
                    else:
                        state_sets[id(clade)] = set()
                else:
                    state_sets[id(clade)] = _merge_nonempty_child_state_sets(
                        state_sets,
                        children,
                    )

        # Preorder pass: resolve ambiguities
        node_regimes = {}
        root = tree.root
        root_set = state_sets.get(id(root), set())
        if root_set:
            node_regimes[id(root)] = min(root_set)
        else:
            node_regimes[id(root)] = sorted(set(tip_regimes.values()))[0]

        for clade in preorder_clades:
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
                    node_regimes[id(clade)] = min(clade_set)
            else:
                node_regimes[id(clade)] = min(clade_set)

        # Branch regime = regime of the child node
        branch_regimes = {}
        for clade in preorder_clades:
            if clade == root:
                continue
            branch_regimes[id(clade)] = node_regimes.get(id(clade), "unknown")

        return branch_regimes

    def _build_root_to_tip_paths(
        self, tree, ordered_names: list[str], parent_map: dict,
        preorder_clades=None,
    ) -> dict[str, list]:
        """Build root-to-tip paths for each tip.

        Returns dict mapping tip_name -> list of (clade_id, branch_length)
        tuples from root to tip.
        """
        # Map tip names to clade objects
        ordered_set = set(ordered_names)
        tip_map = {}
        if preorder_clades is None:
            preorder_clades = self._iter_preorder(tree.root)
        for clade in preorder_clades:
            if not clade.clades and clade.name in ordered_set:
                tip_map[clade.name] = clade

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
        self, tree, ordered_names: list[str], regimes: list[str],
        branch_regimes: dict, parent_map: dict, preorder_clades=None,
    ) -> dict[str, np.ndarray]:
        """Build separate VCV matrix for each regime.

        C_r[i,j] = portion of shared path in regime r between tips i and j.
        """
        n = len(ordered_names)
        paths = self._build_root_to_tip_paths(
            tree, ordered_names, parent_map, preorder_clades
        )

        per_regime_vcv = {r: np.zeros((n, n)) for r in regimes}

        clade_indices = {}
        clade_lengths = {}
        for idx, name in enumerate(ordered_names):
            for clade_id, bl in paths[name]:
                clade_indices.setdefault(clade_id, []).append(idx)
                clade_lengths.setdefault(clade_id, bl)

        for clade_id, indices in clade_indices.items():
            r = branch_regimes.get(clade_id, regimes[0])
            branch_length = clade_lengths[clade_id]
            if len(indices) == 1:
                per_regime_vcv[r][indices[0], indices[0]] += branch_length
                continue
            idx = np.asarray(indices, dtype=np.intp)
            per_regime_vcv[r][np.ix_(idx, idx)] += branch_length

        return per_regime_vcv

    def _fit_single_rate(
        self, y: np.ndarray, C_total: np.ndarray
    ) -> tuple[float, float, float]:
        """Fit single-rate BM model.

        Returns (sigma2, ancestral_state, log_likelihood).
        """
        try:
            return self._fit_single_rate_cholesky(y, C_total)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._fit_single_rate_inverse(y, C_total)

    def _fit_single_rate_cholesky(
        self, y: np.ndarray, C_total: np.ndarray
    ) -> tuple[float, float, float]:
        """Cholesky-backed single-rate BM fit for positive-definite VCVs."""
        n = len(y)
        factor = cho_factor(C_total, lower=True, check_finite=False)
        ones = np.ones(n)
        solve_rhs = np.empty((n, 2), dtype=np.result_type(y, C_total))
        solve_rhs[:, 0] = 1.0
        solve_rhs[:, 1] = y
        solved = cho_solve(factor, solve_rhs, check_finite=False)
        C_inv_ones = solved[:, 0]
        C_inv_y = solved[:, 1]

        a_hat = float(ones @ C_inv_y) / float(ones @ C_inv_ones)
        e = y - a_hat
        C_inv_e = C_inv_y - C_inv_ones * a_hat
        sigma2 = float(e @ C_inv_e) / n

        logdet = 2.0 * float(np.log(np.diag(factor[0])).sum())
        if sigma2 <= 0:
            return sigma2, a_hat, float("-inf")

        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sigma2) + logdet + n)
        return float(sigma2), float(a_hat), float(ll)

    def _fit_single_rate_inverse(
        self, y: np.ndarray, C_total: np.ndarray
    ) -> tuple[float, float, float]:
        """Inverse-based single-rate BM fit retained as a fallback."""
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
        self, y: np.ndarray, per_regime_vcv: dict[str, np.ndarray],
        regimes: list[str],
    ) -> tuple[np.ndarray, float, float]:
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

            _, ll = self._multi_rate_log_likelihood(y, V, ones)
            if not np.isfinite(ll):
                return 1e20
            return -ll

        ones = np.ones(n)

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

        a_hat, ll = self._multi_rate_log_likelihood(y, V, ones)
        if not np.isfinite(ll):
            return sigma2s, 0.0, float("-inf")

        return sigma2s, float(a_hat), float(ll)

    def _multi_rate_log_likelihood(
        self, y: np.ndarray, V: np.ndarray, ones: np.ndarray
    ) -> tuple[float, float]:
        """GLS ancestral state and log-likelihood for a multi-rate covariance."""
        try:
            return self._multi_rate_log_likelihood_cholesky(y, V, ones)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._multi_rate_log_likelihood_inverse(y, V, ones)

    def _multi_rate_log_likelihood_cholesky(
        self, y: np.ndarray, V: np.ndarray, ones: np.ndarray
    ) -> tuple[float, float]:
        n = len(y)
        factor = cho_factor(V, lower=True, check_finite=False)
        solve_rhs = np.empty((n, 2), dtype=np.result_type(y, V))
        solve_rhs[:, 0] = ones
        solve_rhs[:, 1] = y
        solved = cho_solve(factor, solve_rhs, check_finite=False)
        V_inv_ones = solved[:, 0]
        V_inv_y = solved[:, 1]

        denom = float(ones @ V_inv_ones)
        if denom == 0:
            return 0.0, float("-inf")

        a_hat = float(ones @ V_inv_y) / denom
        e = y - a_hat
        V_inv_e = V_inv_y - V_inv_ones * a_hat
        quad = float(e @ V_inv_e)
        logdet = 2.0 * float(np.log(np.diag(factor[0])).sum())
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet + quad)
        return a_hat, ll

    def _multi_rate_log_likelihood_inverse(
        self, y: np.ndarray, V: np.ndarray, ones: np.ndarray
    ) -> tuple[float, float]:
        n = len(y)
        try:
            V_inv = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            return 0.0, float("-inf")

        denom = float(ones @ V_inv @ ones)
        if denom == 0:
            return 0.0, float("-inf")

        a_hat = float(ones @ V_inv @ y) / denom
        e = y - a_hat
        sign, logdet = np.linalg.slogdet(V)
        if sign <= 0:
            return 0.0, float("-inf")
        quad = float(e @ V_inv @ e)
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet + quad)
        return a_hat, ll

    def _parametric_bootstrap(
        self, y: np.ndarray, C_total: np.ndarray,
        per_regime_vcv: dict[str, np.ndarray], regimes: list[str],
        observed_lrt: float, n_sim: int, seed=None,
    ) -> tuple[float, int]:
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
        self, tree, branch_regimes: dict, regimes: list[str],
        parent_map: dict, output_path: str,
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
                "matplotlib is required for regime tree plotting. "
                "Install matplotlib and retry."
            )
            raise SystemExit(2)

        k = len(regimes)
        cmap = plt.get_cmap("tab10")
        regime_colors = {
            regimes[i]: cmap(i / max(k - 1, 1)) for i in range(k)
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

        if self.plot_config.circular:
            # --- Circular mode ---
            import math
            from ...helpers.circular_layout import (
                compute_circular_coords,
                draw_circular_tip_labels,
            )

            coords = compute_circular_coords(
                tree,
                node_x,
                parent_map,
                preorder_clades=preorder_clades,
                terminal_clades=tips,
            )
            ax.set_aspect("equal")
            ax.axis("off")

            # Draw each branch individually with its regime color
            radial_segments = []
            radial_colors = []
            for clade in preorder_clades:
                if clade == root:
                    continue
                cid = id(clade)
                if cid not in parent_map:
                    continue
                parent = parent_map[cid]
                if id(parent) not in coords or cid not in coords:
                    continue

                regime = branch_regimes.get(cid, regimes[0])
                color = regime_colors[regime]

                angle = coords[cid]["angle"]
                r_p = coords[id(parent)]["radius"]
                r_c = coords[cid]["radius"]
                x0 = r_p * math.cos(angle)
                y0 = r_p * math.sin(angle)
                x1 = r_c * math.cos(angle)
                y1 = r_c * math.sin(angle)
                radial_segments.append(((x0, y0), (x1, y1)))
                radial_colors.append(color)

            if radial_segments:
                ax.add_collection(
                    LineCollection(
                        radial_segments,
                        colors=radial_colors,
                        linewidths=2.5,
                        capstyle="round",
                        zorder=2,
                    ),
                    autolim=True,
                )

            # Draw arcs at internal nodes colored by regime
            arc_segments = []
            arc_colors = []
            arc_fractions = [idx / 60 for idx in range(61)]
            for clade in preorder_clades:
                if not clade.clades:
                    continue
                cid = id(clade)
                pc = coords[cid]
                child_angles = [coords[id(ch)]["angle"] for ch in clade.clades]
                if len(child_angles) < 2:
                    continue
                start_a = min(child_angles)
                end_a = max(child_angles)
                span = (end_a - start_a) % (2.0 * math.pi)
                if span > math.pi:
                    start_a, end_a = end_a, start_a

                arc_color = "gray"
                if cid in branch_regimes:
                    arc_color = regime_colors.get(branch_regimes[cid], "gray")

                start = start_a % (2.0 * math.pi)
                end = end_a % (2.0 * math.pi)
                diff = (end - start) % (2.0 * math.pi)
                if diff > math.pi:
                    diff = diff - 2.0 * math.pi
                radius = pc["radius"]
                arc_segments.append([
                    (
                        radius * math.cos(start + diff * fraction),
                        radius * math.sin(start + diff * fraction),
                    )
                    for fraction in arc_fractions
                ])
                arc_colors.append(arc_color)

            if arc_segments:
                ax.add_collection(
                    LineCollection(
                        arc_segments,
                        colors=arc_colors,
                        linewidths=1.5,
                        capstyle="round",
                        zorder=1,
                    ),
                    autolim=True,
                )
            if radial_segments or arc_segments:
                ax.autoscale_view()

            # Tip labels
            max_x = max(node_x.values()) if node_x else 1.0
            label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize is not None else 9
            if label_fontsize > 0:
                draw_circular_tip_labels(ax, tree, coords, fontsize=label_fontsize, offset=max_x * 0.02)

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
                Line2D([0], [0], color=regime_colors[r], linewidth=3, label=r)
                for r in regimes
            ]
            if self.plot_config.color_file:
                color_legend = build_color_legend_handles(color_data)
                handles.extend(color_legend)
            ax.legend(
                handles=handles, title="Regimes", loc="upper left",
                fontsize=8, title_fontsize=9,
            )

            if config.show_title:
                ax.set_title(config.title or "Regime Tree (Rate Heterogeneity)", fontsize=config.title_fontsize)
            fig.tight_layout()
            fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
            plt.close(fig)
            print(f"Saved regime tree plot: {output_path}")

        else:
            # --- Rectangular mode ---
            from matplotlib.collections import LineCollection

            vertical_segments = []
            horizontal_segments = []
            horizontal_colors = []
            for clade in preorder_clades:
                if clade == root:
                    continue
                cid = id(clade)
                if cid not in parent_map:
                    continue

                parent = parent_map[cid]
                pid = id(parent)
                if pid not in node_x or cid not in node_x:
                    continue

                parent_x = node_x[pid]
                parent_y = node_y[pid]
                child_y = node_y[cid]
                child_x = node_x[cid]

                regime = branch_regimes.get(cid, regimes[0])
                color = regime_colors[regime]

                vertical_segments.append(((parent_x, parent_y), (parent_x, child_y)))
                horizontal_segments.append(((parent_x, child_y), (child_x, child_y)))
                horizontal_colors.append(color)

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
            if horizontal_segments:
                ax.add_collection(
                    LineCollection(
                        horizontal_segments,
                        colors=horizontal_colors,
                        linewidths=2.5,
                        capstyle="butt",
                        zorder=2,
                    ),
                    autolim=True,
                )
            ax.autoscale_view()

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

            handles = [
                Line2D([0], [0], color=regime_colors[r], linewidth=3, label=r)
                for r in regimes
            ]
            if self.plot_config.color_file:
                color_legend = build_color_legend_handles(color_data)
                handles.extend(color_legend)
            ax.legend(
                handles=handles, title="Regimes", loc="upper left",
                fontsize=8, title_fontsize=9,
            )

            ax.set_xlabel("Branch length")
            ax.set_yticks([])
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            if config.show_title:
                ax.set_title(config.title or "Regime Tree (Rate Heterogeneity)", fontsize=config.title_fontsize)
            fig.tight_layout()
            fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
            plt.close(fig)
            print(f"Saved regime tree plot: {output_path}")

    def _print_text_output(
        self, *, n_tips, regimes, sigma2_single, anc_single, ll_single,
        aic_single, sigma2_multi_dict, anc_multi, ll_multi, aic_multi,
        lrt_stat, df, chi2_p, sim_p, n_sim_done, r2_regime,
    ) -> None:
        lines = [
            "Rate Heterogeneity Test (Multi-rate Brownian Motion)",
            f"\nRegimes: {len(regimes)} ({', '.join(regimes)})",
            f"Number of tips: {n_tips}",
            "\nSingle-rate model (H0):",
            f"  Sigma-squared:      {sigma2_single:.4f}",
            f"  Ancestral state:    {anc_single:.4f}",
            f"  Log-likelihood:     {ll_single:.2f}",
            f"  AIC:                {aic_single:.2f}",
            "\nMulti-rate model (H1):",
            f"  {'Regime':<20s}{'Sigma-squared':>14s}",
        ]
        lines.extend(f"  {r:<20s}{sigma2_multi_dict[r]:>14.4f}" for r in regimes)
        lines.extend(
            [
                f"  Ancestral state:    {anc_multi:.4f}",
                f"  Log-likelihood:     {ll_multi:.2f}",
                f"  AIC:                {aic_multi:.2f}",
                "\nLikelihood ratio test:",
                f"  LRT statistic:      {lrt_stat:.4f}",
                f"  Degrees of freedom: {df}",
                f"  Chi-squared p-value: {chi2_p:.4f}",
            ]
        )
        if sim_p is not None:
            lines.append(f"  Simulated p-value:  {sim_p:.4f} ({n_sim_done} simulations)")

        lines.extend([f"\nEffect size:", f"  R2_regime: {r2_regime:.4f}"])
        print("\n".join(lines))

    def _format_result(
        self, *, n_tips, regimes, sigma2_single, anc_single, ll_single,
        aic_single, sigma2_multi_dict, anc_multi, ll_multi, aic_multi,
        lrt_stat, df, chi2_p, sim_p, n_sim_done, r2_regime,
    ) -> dict:
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
