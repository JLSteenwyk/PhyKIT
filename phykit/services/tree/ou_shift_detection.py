import copy
import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize_scalar
from sklearn.linear_model import lars_path

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class OUShiftDetection(Tree):
    """Automatic OU shift detection using LASSO (l1ou approach).

    Discovers where on the phylogeny the adaptive optimum changed,
    using the LASSO-based approach from Khabbazian et al. (2016).
    No regime file is needed — only a tree and trait data.

    Key algorithmic choices matching R's l1ou package:
    - Two-pass LASSO: pass 1 with alpha≈0, pass 2 with re-estimated alpha
    - Per-candidate alpha optimization (like R's phylolm)
    - R-style pBIC with Fisher information correction
    - R-style single-pass backward elimination per LASSO candidate
    """

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.criterion = parsed["criterion"]
        self.max_shifts = parsed["max_shifts"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> Dict:
        criterion = getattr(args, "criterion", "pBIC")
        if criterion not in ("pBIC", "BIC", "AICc"):
            raise PhykitUserError(
                [
                    f"Unknown criterion '{criterion}'.",
                    "Available criteria: pBIC, BIC, AICc",
                ],
                code=2,
            )
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            criterion=criterion,
            max_shifts=getattr(args, "max_shifts", None),
            json_output=getattr(args, "json", False),
        )

    def run(self) -> None:
        tree = self.read_tree_file()
        self._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)
        traits = self._parse_trait_file(self.trait_data_path, tree_tips)

        shared = set(traits.keys())

        # Prune tree to shared taxa
        tree_copy = copy.deepcopy(tree)
        tip_names_in_tree = [t.name for t in tree_copy.get_terminals()]
        to_prune = [t for t in tip_names_in_tree if t not in shared]
        if to_prune:
            for t in to_prune:
                tree_copy.prune(t)

        # Store working data as instance attributes for internal methods
        self._ordered_names = sorted(shared)
        self._x = np.array([traits[name] for name in self._ordered_names])
        self._n = len(self._ordered_names)
        self._parent_map = self._build_parent_map(tree_copy)
        self._lineage_info = self._build_lineage_info_no_regime(
            tree_copy, self._ordered_names, self._parent_map
        )
        self._tree_height = max(
            sum(bl for _, bl, _, _ in self._lineage_info[name])
            for name in self._ordered_names
        )
        self._edges = self._enumerate_edges(tree_copy, self._parent_map)
        self._n_edges = len(self._edges)

        # Precompute shared path lengths for vectorized VCV construction
        self._S, self._tip_heights_arr = self._precompute_shared_path_lengths()

        max_shifts = self.max_shifts
        if max_shifts is None:
            max_shifts = self._n // 2
        max_shifts = min(max_shifts, self._n - 2)

        if self._n_edges == 0:
            null = self._fit_and_score_config([])
            if null is None:
                null = self._empty_result()
            result = self._build_final_result(null, tree_copy)
            self._output(result)
            return

        # ── Two-pass LASSO (Khabbazian et al. 2016) ──────────────
        # Pass 1: BM approximation (alpha ≈ 0)
        best = self._lasso_pass_with_selection(
            alpha_for_lasso=1e-6, max_shifts=max_shifts
        )

        if best is None:
            best = self._fit_and_score_config([])
        if best is None:
            best = self._empty_result()

        alpha_pass1 = best.get("alpha", 0.0)

        # Pass 2: OU with estimated alpha (if alpha changed meaningfully)
        if alpha_pass1 > 1e-4:
            best2 = self._lasso_pass_with_selection(
                alpha_for_lasso=alpha_pass1, max_shifts=max_shifts
            )
            if best2 is not None and best2[self.criterion] < best[self.criterion]:
                best = best2

        result = self._build_final_result(best, tree_copy)
        self._output(result)

    # ── Tree & data parsing (adapted from OUwie) ─────────────────────

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for OU shift detection."],
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

    # ── Tree structure helpers ───────────────────────────────────────

    def _build_parent_map(self, tree) -> Dict:
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    def _build_lineage_info_no_regime(
        self, tree, ordered_names: List[str], parent_map: Dict,
    ) -> Dict[str, List[Tuple]]:
        """Build root-to-tip lineage info without regime labels.

        Returns dict mapping tip_name -> list of tuples:
            (clade_id, branch_length, dist_from_root_start, dist_from_root_end)
        """
        tip_map = {}
        for tip in tree.get_terminals():
            if tip.name in ordered_names:
                tip_map[tip.name] = tip

        lineage_info = {}
        for name in ordered_names:
            tip = tip_map[name]
            path = []
            current = tip
            while id(current) in parent_map:
                bl = current.branch_length if current.branch_length else 0.0
                path.append((id(current), bl))
                current = parent_map[id(current)]
            path.reverse()

            enriched = []
            cum = 0.0
            for clade_id, bl in path:
                d_start = cum
                d_end = cum + bl
                enriched.append((clade_id, bl, d_start, d_end))
                cum += bl

            lineage_info[name] = enriched

        return lineage_info

    # ── Edge enumeration ─────────────────────────────────────────────

    def _enumerate_edges(self, tree, parent_map: Dict) -> List[Tuple]:
        """List all non-root edges as candidate shift locations.

        Returns list of (clade_id, clade_object) tuples for every
        non-root clade in the tree.
        """
        edges = []
        root = tree.root
        for clade in tree.find_clades(order="preorder"):
            if clade == root:
                continue
            edges.append((id(clade), clade))
        return edges

    def _count_descendants(self, tree, edges: List[Tuple]) -> Dict:
        """Count the number of tip descendants for each edge."""
        counts = {}
        for clade_id, clade in edges:
            tips = list(clade.get_terminals())
            counts[clade_id] = len(tips)
        return counts

    # ── Precomputed phylogenetic data ────────────────────────────────

    def _precompute_shared_path_lengths(self) -> Tuple[np.ndarray, np.ndarray]:
        """Precompute shared path length matrix S and tip height vector.

        S[i,j] = total branch length shared between root-to-tip_i and
                  root-to-tip_j paths (= distance from root to their LCA).
        """
        n = self._n
        ordered_names = self._ordered_names
        lineage_info = self._lineage_info

        S = np.zeros((n, n))
        tip_heights = np.zeros(n)

        for i in range(n):
            path_i = lineage_info[ordered_names[i]]
            tip_heights[i] = sum(bl for _, bl, _, _ in path_i)
            S[i, i] = tip_heights[i]
            for j in range(i + 1, n):
                path_j = lineage_info[ordered_names[j]]
                s = 0.0
                for idx in range(min(len(path_i), len(path_j))):
                    if path_i[idx][0] == path_j[idx][0]:
                        s += path_i[idx][1]
                    else:
                        break
                S[i, j] = s
                S[j, i] = s

        return S, tip_heights

    def _build_ou_vcv_fast(self, alpha: float, sigma2: float = 1.0) -> np.ndarray:
        """Vectorized OU VCV matrix using precomputed shared path lengths.

        V[i,j] = (sigma2/(2*alpha)) * exp(-alpha*(T_i+T_j-2*s_ij))
                                     * (1 - exp(-2*alpha*s_ij))
        """
        S = self._S
        T = self._tip_heights_arr

        if alpha < 1e-10:
            return sigma2 * S

        D = T[:, None] + T[None, :] - 2.0 * S
        return (sigma2 / (2.0 * alpha)) * np.exp(-alpha * D) * (
            1.0 - np.exp(-2.0 * alpha * S)
        )

    # ── OU1 alpha estimation ─────────────────────────────────────────

    def _fit_ou1_for_alpha(
        self, x: np.ndarray = None, ordered_names: List[str] = None,
        lineage_info: Dict = None, tree_height: float = None,
    ) -> float:
        """Estimate alpha from a single-regime OU1 fit."""
        x = x if x is not None else self._x
        n = len(x)
        th = tree_height if tree_height is not None else self._tree_height
        upper = 100.0 / th if th > 0 else 100.0

        use_fast = hasattr(self, "_S")

        def neg_ll(alpha):
            if use_fast:
                V = self._build_ou_vcv_fast(alpha, 1.0)
            else:
                V = self._build_ou_vcv_single_alpha_no_regime(
                    ordered_names, lineage_info, alpha, 1.0
                )
            try:
                V_inv = np.linalg.inv(V)
            except np.linalg.LinAlgError:
                return 1e20

            ones = np.ones(n)
            denom = float(ones @ V_inv @ ones)
            if abs(denom) < 1e-300:
                return 1e20
            theta = float(ones @ V_inv @ x) / denom
            e = x - theta
            sig2 = float(e @ V_inv @ e) / n
            if sig2 <= 0:
                return 1e20

            sign, logdet = np.linalg.slogdet(V)
            if sign <= 0:
                return 1e20

            ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)
            return -ll

        alpha_hat, _ = self._optimize_parameter(neg_ll, (1e-8, upper))
        return alpha_hat

    # ── Two-pass LASSO with per-candidate model selection ────────────

    def _lasso_pass_with_selection(
        self, alpha_for_lasso: float, max_shifts: int,
    ) -> Dict:
        """Run one LASSO pass: build design, extract configs, evaluate each.

        For each candidate configuration from the LASSO path:
        1. Fit full model with per-candidate alpha optimization
        2. Apply R-style backward elimination (single pass, ≥3 shifts)
        3. Score with selected criterion
        Returns the best-scoring model.
        """
        edges = self._edges
        n = self._n

        # Build design matrix and VCV at the LASSO alpha
        W = self._build_shift_weight_matrix(
            self._ordered_names, self._lineage_info, edges,
            alpha_for_lasso, self._tree_height,
        )
        V = self._build_ou_vcv_fast(alpha_for_lasso)
        X_star, y_star = self._transform_to_independent(self._x, W, V)

        # Extract candidate configurations from LASSO path
        configs = self._extract_lasso_configs(X_star, y_star, max_shifts)

        # Evaluate null model
        best = self._fit_and_score_config([])

        # Evaluate each LASSO candidate
        for config in configs:
            fitted = self._fit_and_score_config(config)
            if fitted is None:
                continue

            # R-style backward elimination (single pass, only for ≥3 shifts)
            if fitted["n_shifts"] >= 3:
                fitted = self._backward_eliminate_single_pass(fitted)

            if best is None or fitted[self.criterion] < best[self.criterion]:
                best = fitted

        return best

    def _extract_lasso_configs(
        self, X_star: np.ndarray, y_star: np.ndarray, max_shifts: int,
    ) -> List[List[int]]:
        """Extract unique shift configurations from LASSO path.

        The root optimum (column 0) is unpenalized via residualization.
        Returns list of shift-edge-index lists.
        """
        X_intercept = X_star[:, 0:1]
        X_shifts = X_star[:, 1:]
        p = X_shifts.shape[1]

        if p == 0:
            return []

        # Residualize against intercept (project out root optimum)
        XtX_inv = 1.0 / (X_intercept.T @ X_intercept)[0, 0]
        beta_int_y = XtX_inv * (X_intercept.T @ y_star)[0]
        y_resid = y_star - X_intercept.flatten() * beta_int_y

        X_resid = np.zeros_like(X_shifts)
        for j in range(p):
            beta_int_j = XtX_inv * (X_intercept.T @ X_shifts[:, j:j + 1])[0, 0]
            X_resid[:, j] = X_shifts[:, j] - X_intercept.flatten() * beta_int_j

        col_norms = np.linalg.norm(X_resid, axis=0)
        col_norms[col_norms < 1e-12] = 1.0
        X_resid_norm = X_resid / col_norms

        max_iter = min(max_shifts, p, self._n - 1)
        if max_iter < 1:
            return []

        try:
            _, _, coefs_path = lars_path(
                X_resid_norm, y_resid, method="lasso", max_iter=max_iter,
            )
        except Exception:
            return []

        configs = []
        seen = {frozenset()}

        for step in range(coefs_path.shape[1]):
            coefs = coefs_path[:, step] / col_norms
            nonzero_idx = np.where(np.abs(coefs) > 1e-10)[0]
            config = frozenset(nonzero_idx.tolist())
            if config in seen or len(nonzero_idx) == 0 or len(nonzero_idx) > max_shifts:
                continue
            seen.add(config)
            configs.append(nonzero_idx.tolist())

        return configs

    # ── Per-candidate model fitting ──────────────────────────────────

    def _fit_and_score_config(self, shift_edge_indices: List[int]) -> Dict:
        """Fit OU model with given shifts, jointly optimizing alpha.

        Uses R's l1ou approach: indicator (simpX) design matrix for model
        fitting and scoring (phylolm style), not OU weight matrix. The OU
        weight matrix is only used during the LASSO step.

        For each candidate shift configuration:
        1. Build indicator design matrix (0/1 for shift edge on path)
        2. Optimize alpha by profiled concentrated likelihood
        3. GLS estimates of optima (intercept = root, shifts = deltas)
        4. Compute IC scores using indicator matrix for pBIC vcov

        Returns dict with model parameters and scores, or None on failure.
        """
        k = len(shift_edge_indices)
        n = self._n
        x = self._x
        upper = 100.0 / self._tree_height if self._tree_height > 0 else 100.0

        if k == 0:
            return self._fit_null_model()

        # Build indicator design matrix (R's simpX approach)
        W_ind = self._build_indicator_design_matrix(shift_edge_indices)

        def neg_ll(alpha):
            V = self._build_ou_vcv_fast(alpha)
            try:
                V_inv = np.linalg.inv(V)
            except np.linalg.LinAlgError:
                return 1e20

            theta = self._gls_theta_hat(x, W_ind, V_inv)
            e = x - W_ind @ theta
            sig2 = float(e @ V_inv @ e) / n
            if sig2 <= 0:
                return 1e20

            sign, logdet = np.linalg.slogdet(V)
            if sign <= 0:
                return 1e20

            return 0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)

        alpha_hat = self._optimize_alpha(neg_ll, (1e-8, upper))

        # Final fit at optimal alpha with indicator matrix
        V = self._build_ou_vcv_fast(alpha_hat)
        try:
            V_inv = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            return None

        theta = self._gls_theta_hat(x, W_ind, V_inv)
        e = x - W_ind @ theta
        sig2 = float(e @ V_inv @ e) / n
        if sig2 <= 0:
            return None

        sign, logdet = np.linalg.slogdet(V)
        if sign <= 0:
            return None

        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)

        # theta[0] = root optimum (intercept)
        # theta[j+1] = shift delta for edge j (added to root)
        # Actual optimum for shifted regime = theta[0] + theta[j+1]
        optima = {}
        for idx_in_reduced, j in enumerate(shift_edge_indices):
            optima[j] = float(theta[0] + theta[idx_in_reduced + 1])

        return dict(
            n_shifts=k,
            shift_edge_indices=list(shift_edge_indices),
            alpha=float(alpha_hat),
            sigma2=float(sig2),
            theta_root=float(theta[0]),
            optima=optima,
            log_likelihood=float(ll),
            pBIC=float(self._compute_pbic(n, ll, k, W_ind, V_inv, sig2, x)),
            BIC=float(self._compute_bic(n, ll, k)),
            AICc=float(self._compute_aicc(n, ll, k)),
        )

    def _fit_null_model(self) -> Dict:
        """Fit the null (no-shift) OU1 model with alpha optimization."""
        n = self._n
        x = self._x

        alpha_hat = self._fit_ou1_for_alpha()

        V = self._build_ou_vcv_fast(alpha_hat)
        try:
            V_inv = np.linalg.inv(V)
        except np.linalg.LinAlgError:
            return None

        W = np.ones((n, 1))
        theta = self._gls_theta_hat(x, W, V_inv)
        e = x - W @ theta
        sig2 = float(e @ V_inv @ e) / n
        if sig2 <= 0:
            return None

        sign, logdet = np.linalg.slogdet(V)
        if sign <= 0:
            return None

        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)

        return dict(
            n_shifts=0,
            shift_edge_indices=[],
            alpha=float(alpha_hat),
            sigma2=float(sig2),
            theta_root=float(theta[0]),
            optima={},
            log_likelihood=float(ll),
            pBIC=float(self._compute_pbic(n, ll, 0, W, V_inv, sig2, x)),
            BIC=float(self._compute_bic(n, ll, 0)),
            AICc=float(self._compute_aicc(n, ll, 0)),
        )

    # ── Backward elimination (R-style) ──────────────────────────────

    def _backward_eliminate_single_pass(self, fitted: Dict) -> Dict:
        """R-style backward elimination: single greedy pass.

        For each shift in the current configuration, tests removing it
        (with full model re-fit including alpha optimization). If removing
        it improves the criterion score, keeps the removal. One pass only.
        """
        current = fitted
        shift_indices = list(current["shift_edge_indices"])
        current_score = current[self.criterion]

        for idx_to_remove in list(shift_indices):
            if idx_to_remove not in shift_indices:
                continue  # already removed earlier in this pass

            reduced = [j for j in shift_indices if j != idx_to_remove]
            candidate = self._fit_and_score_config(reduced)
            if candidate is None:
                continue

            if candidate[self.criterion] < current_score:
                shift_indices = list(reduced)
                current = candidate
                current_score = candidate[self.criterion]

        return current

    # ── Alpha optimization ──────────────────────────────────────────

    def _optimize_alpha(self, neg_ll_func, bounds) -> float:
        """Fast 1-D bounded optimization for alpha."""
        lo, hi = bounds
        try:
            res = minimize_scalar(
                neg_ll_func, bounds=(lo, hi), method="bounded",
                options={"xatol": 1e-4},
            )
            return res.x
        except (ValueError, RuntimeError):
            return (lo + hi) / 2.0

    def _optimize_parameter(self, neg_ll_func, bounds, niter=10):
        """Multi-bracket 1-D optimization (for OU1 alpha estimation)."""
        lo_bound, hi_bound = bounds
        bounds_lo = np.linspace(
            lo_bound, hi_bound - (hi_bound - lo_bound) / niter, niter
        )
        bounds_hi = np.linspace(
            lo_bound + (hi_bound - lo_bound) / niter, hi_bound, niter
        )

        best_ll = -np.inf
        best_param = (lo_bound + hi_bound) / 2.0

        for lo, hi in zip(bounds_lo, bounds_hi):
            if lo >= hi:
                continue
            try:
                res = minimize_scalar(
                    neg_ll_func, bounds=(lo, hi), method="bounded"
                )
                ll_val = -res.fun
                if ll_val > best_ll:
                    best_ll = ll_val
                    best_param = res.x
            except (ValueError, RuntimeError):
                continue

        return best_param, best_ll

    # ── OU VCV (slow path for backward compatibility) ────────────────

    def _build_ou_vcv_single_alpha_no_regime(
        self, ordered_names: List[str], lineage_info: Dict,
        alpha: float, sigma2: float,
    ) -> np.ndarray:
        """Build OU VCV matrix without regime labels (loop-based)."""
        n = len(ordered_names)
        V = np.zeros((n, n))

        tip_heights = {}
        for name in ordered_names:
            tip_heights[name] = sum(bl for _, bl, _, _ in lineage_info[name])

        for i in range(n):
            for j in range(i, n):
                T_i = tip_heights[ordered_names[i]]
                T_j = tip_heights[ordered_names[j]]

                if i == j:
                    s_ij = T_i
                else:
                    path_i = lineage_info[ordered_names[i]]
                    path_j = lineage_info[ordered_names[j]]
                    s_ij = 0.0
                    min_len = min(len(path_i), len(path_j))
                    for s in range(min_len):
                        if path_i[s][0] == path_j[s][0]:
                            s_ij += path_i[s][1]
                        else:
                            break

                if alpha < 1e-10:
                    val = sigma2 * s_ij
                else:
                    d_i = T_i - s_ij
                    d_j = T_j - s_ij
                    val = (sigma2 / (2.0 * alpha)) * np.exp(
                        -alpha * (d_i + d_j)
                    ) * (1.0 - np.exp(-2.0 * alpha * s_ij))

                V[i, j] = val
                V[j, i] = val

        return V

    # ── Design matrices ────────────────────────────────────────────

    def _build_shift_weight_matrix(
        self, ordered_names: List[str], lineage_info: Dict,
        edges: List[Tuple], alpha: float, tree_height: float,
    ) -> np.ndarray:
        """Build n x (p+1) OU-weight design matrix for LASSO shift detection.

        Column 0: root/background optimum contribution
        Columns 1..p: one per candidate shift edge
        Column 0 = 1 - sum(columns 1..p), ensuring rows always sum to 1.
        """
        n = len(ordered_names)
        p = len(edges)
        W = np.zeros((n, p + 1))

        edge_col = {}
        for j, (clade_id, _) in enumerate(edges):
            edge_col[clade_id] = j + 1

        for i, name in enumerate(ordered_names):
            path = lineage_info[name]
            T_i = sum(bl for _, bl, _, _ in path)

            for clade_id, bl, d_start, d_end in path:
                if clade_id in edge_col:
                    col = edge_col[clade_id]
                    if alpha < 1e-10:
                        w = bl / T_i if T_i > 0 else 0.0
                    else:
                        w = (1.0 - np.exp(-alpha * bl)) * np.exp(
                            -alpha * (T_i - d_end)
                        )
                    W[i, col] = w

            W[i, 0] = 1.0 - np.sum(W[i, 1:])

        return W

    def _build_indicator_design_matrix(
        self, shift_edge_indices: List[int],
    ) -> np.ndarray:
        """Build n x (k+1) indicator (simpX) design matrix for scoring.

        This matches R's l1ou ``generate_design_matrix(tree, "simpX")``:
        - Column 0: intercept (all 1s)
        - Column j+1: 1 if shift edge j is on the root-to-tip_i path, else 0

        R uses this indicator matrix (not OU weights) for phylolm fitting
        and pBIC computation when scoring candidate shift configurations.
        """
        n = self._n
        k = len(shift_edge_indices)
        W = np.zeros((n, k + 1))
        W[:, 0] = 1.0  # intercept

        if k == 0:
            return W

        # Map selected shift edge indices to their clade IDs
        shift_clade_ids = {}
        for col_idx, edge_idx in enumerate(shift_edge_indices):
            clade_id = self._edges[edge_idx][0]
            shift_clade_ids[clade_id] = col_idx + 1

        # For each tip, check which shift edges are on its root-to-tip path
        for i, name in enumerate(self._ordered_names):
            path = self._lineage_info[name]
            for clade_id, bl, d_start, d_end in path:
                if clade_id in shift_clade_ids:
                    W[i, shift_clade_ids[clade_id]] = 1.0

        return W

    # ── Cholesky transform ───────────────────────────────────────────

    def _transform_to_independent(
        self, x: np.ndarray, W: np.ndarray, V: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Cholesky-transform to remove phylogenetic correlation."""
        try:
            L = np.linalg.cholesky(V)
            L_inv = np.linalg.inv(L)
        except np.linalg.LinAlgError:
            V_reg = V + np.eye(len(x)) * 1e-8 * np.trace(V) / len(x)
            L = np.linalg.cholesky(V_reg)
            L_inv = np.linalg.inv(L)

        y_star = L_inv @ x
        X_star = L_inv @ W

        return X_star, y_star

    # ── GLS ──────────────────────────────────────────────────────────

    def _gls_theta_hat(
        self, x: np.ndarray, W: np.ndarray, V_inv: np.ndarray
    ) -> np.ndarray:
        """GLS estimate of theta: (W'V^{-1}W)^{-1} W'V^{-1}x"""
        WtVi = W.T @ V_inv
        try:
            theta = np.linalg.solve(WtVi @ W, WtVi @ x)
        except np.linalg.LinAlgError:
            theta = np.linalg.lstsq(WtVi @ W, WtVi @ x, rcond=None)[0]
        return theta

    # ── Information criteria (R-style) ──────────────────────────────

    def _compute_pbic(
        self, n: int, log_likelihood: float, k: int,
        W: np.ndarray, V_inv: np.ndarray, sigma2: float,
        x: np.ndarray,
    ) -> float:
        """Phylogenetic BIC matching R's l1ou implementation.

        pBIC = 2*k*log(nEdges-1) - 2*LL + 2*log(n) - log(det(scaled_vcov))

        where scaled_vcov = sigma2 * (W'V^{-1}W)^{-1} * (n-d)/(var(y)*n)
        and d = 1 + k (intercept + shift magnitudes).
        """
        n_edges = self._n_edges if hasattr(self, "_n_edges") else 0

        # Combinatorial penalty for choosing k edges from (nEdges-1) candidates
        df1 = 2.0 * k * np.log(max(n_edges - 1, 1)) if k > 0 else 0.0

        # Fisher information correction
        d = 1 + k  # intercept + shift magnitudes
        WtVi = W.T @ V_inv
        info = WtVi @ W

        try:
            vcov = sigma2 * np.linalg.inv(info)
        except np.linalg.LinAlgError:
            return float("inf")

        var_y = np.var(x, ddof=1)
        if var_y < 1e-300:
            return float("inf")

        if n - d <= 0:
            return float("inf")

        scaled = vcov * (n - d) / (var_y * n)
        sign, ld = np.linalg.slogdet(scaled)
        if sign <= 0:
            return float("inf")

        return df1 - 2.0 * log_likelihood + 2.0 * np.log(n) - ld

    def _compute_bic(self, n: int, log_likelihood: float, k: int) -> float:
        """BIC matching R's l1ou: -2*LL + (2k+3)*log(n).

        Penalty includes k for shift selection + (k+3) for regression
        parameters (k shift magnitudes + intercept + alpha + sigma2).
        """
        return -2.0 * log_likelihood + (2 * k + 3) * np.log(n)

    def _compute_aicc(self, n: int, log_likelihood: float, k: int) -> float:
        """AICc matching R's l1ou with p = 2k + 3 effective parameters."""
        p = 2 * k + 3
        if p >= n - 1:
            return float("inf")
        return -2.0 * log_likelihood + 2.0 * p + 2.0 * p * (p + 1) / (n - p - 1)

    # ── Edge description ─────────────────────────────────────────────

    def _describe_edge(self, clade, tree) -> str:
        """Human-readable edge name."""
        tips = [t.name for t in clade.get_terminals()]
        if len(tips) == 1:
            return f"terminal branch to {tips[0]}"
        elif len(tips) <= 3:
            return f"stem of ({', '.join(sorted(tips))})"
        else:
            return f"stem of ({', '.join(sorted(tips[:2]))}, ... +{len(tips)-2} more)"

    def _identify_shift_edges(
        self, best: Dict, edges: List[Tuple], tree
    ) -> List[Dict]:
        """Map selected shift edges to human-readable descriptions."""
        shifts = []
        for j in best["shift_edge_indices"]:
            clade_id, clade = edges[j]
            desc = self._describe_edge(clade, tree)
            optimum = best["optima"].get(j, best["theta_root"])
            shifts.append(dict(
                edge_index=j,
                description=desc,
                optimum=float(optimum),
            ))
        return shifts

    # ── Final result assembly ────────────────────────────────────────

    def _build_final_result(self, best: Dict, tree) -> Dict:
        """Assemble final result dictionary for output."""
        shifts = self._identify_shift_edges(best, self._edges, tree)

        return dict(
            n_tips=self._n,
            n_shifts=best["n_shifts"],
            criterion=self.criterion,
            alpha=best["alpha"],
            sigma2=best["sigma2"],
            theta_root=best["theta_root"],
            log_likelihood=best["log_likelihood"],
            pBIC=best["pBIC"],
            BIC=best["BIC"],
            AICc=best["AICc"],
            shifts=shifts,
        )

    def _empty_result(self) -> Dict:
        """Fallback result when no model could be fit."""
        return dict(
            n_shifts=0,
            shift_edge_indices=[],
            alpha=0.0,
            sigma2=0.0,
            theta_root=0.0,
            optima={},
            log_likelihood=float("-inf"),
            pBIC=float("inf"),
            BIC=float("inf"),
            AICc=float("inf"),
        )

    # ── Output methods ───────────────────────────────────────────────

    def _output(self, result: Dict) -> None:
        if self.json_output:
            self._print_json_output(result)
        else:
            self._print_text_output(result)

    def _print_text_output(self, result: Dict) -> None:
        try:
            print("=" * 60)
            print("OU Shift Detection (l1ou)")
            print("=" * 60)
            print(f"Number of tips:       {result['n_tips']}")
            print(f"Number of shifts:     {result['n_shifts']}")
            print(f"Selection criterion:  {result['criterion']}")
            print(f"Alpha (OU strength):  {result['alpha']:.6f}")
            print(f"Sigma² (BM rate):     {result['sigma2']:.6f}")
            print(f"Root optimum (θ₀):    {result['theta_root']:.6f}")
            print(f"Log-likelihood:       {result['log_likelihood']:.4f}")
            print(f"pBIC:                 {result['pBIC']:.4f}")
            print(f"BIC:                  {result['BIC']:.4f}")
            print(f"AICc:                 {result['AICc']:.4f}")

            if result["shifts"]:
                print()
                print("Detected shifts:")
                print("-" * 60)
                for i, shift in enumerate(result["shifts"], 1):
                    print(f"  Shift {i}: {shift['description']}")
                    print(f"           New optimum: {shift['optimum']:.6f}")
            else:
                print()
                print("No shifts detected — single-regime OU is best.")
            print("=" * 60)
        except BrokenPipeError:
            pass

    def _print_json_output(self, result: Dict) -> None:
        print_json(result, sort_keys=False)
