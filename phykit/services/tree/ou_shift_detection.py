from __future__ import annotations

import sys

from .base import Tree
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyPickle:
    _module = None

    def _load(self):
        module = self._module
        if module is None:
            import pickle as module

            self._module = module
        return module

    def dumps(self, *args, **kwargs):
        module = self._load()
        dumps = module.dumps
        self.dumps = dumps
        if "loads" not in self.__dict__:
            self.loads = module.loads

        return dumps(*args, **kwargs)

    def loads(self, *args, **kwargs):
        module = self._load()
        loads = module.loads
        self.loads = loads
        if "dumps" not in self.__dict__:
            self.dumps = module.dumps

        return loads(*args, **kwargs)

    def __getattr__(self, name):
        value = getattr(self._load(), name)
        setattr(self, name, value)
        return value



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


pickle = _LazyPickle()
np = _LazyNumpy()
_CHO_FACTOR = None
_CHO_SOLVE = None
_SOLVE_TRIANGULAR = None
_MINIMIZE_SCALAR = None


def cho_factor(*args, **kwargs):
    global _CHO_FACTOR
    if _CHO_FACTOR is None:
        from scipy.linalg import cho_factor as _CHO_FACTOR

    return _CHO_FACTOR(*args, **kwargs)


def cho_solve(*args, **kwargs):
    global _CHO_SOLVE
    if _CHO_SOLVE is None:
        from scipy.linalg import cho_solve as _CHO_SOLVE

    return _CHO_SOLVE(*args, **kwargs)


def solve_triangular(*args, **kwargs):
    global _SOLVE_TRIANGULAR
    if _SOLVE_TRIANGULAR is None:
        from scipy.linalg import solve_triangular as _SOLVE_TRIANGULAR

    return _SOLVE_TRIANGULAR(*args, **kwargs)


def minimize_scalar(*args, **kwargs):
    global _MINIMIZE_SCALAR
    if _MINIMIZE_SCALAR is None:
        from scipy.optimize import minimize_scalar as _MINIMIZE_SCALAR

    return _MINIMIZE_SCALAR(*args, **kwargs)


def _lars_path(*args, **kwargs):
    from sklearn.linear_model import lars_path

    return lars_path(*args, **kwargs)


def _column_l2_norms(matrix):
    return np.sqrt(np.einsum("ij,ij->j", matrix, matrix))


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

    def process_args(self, args) -> dict:
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
        tree = self.read_tree_file_unmodified()
        self.validate_tree(tree, min_tips=3, require_branch_lengths=True, context="OU shift detection")

        tree_tips = self.get_tip_names_from_tree(tree)
        traits = self._parse_trait_file(self.trait_data_path, tree_tips)

        to_prune, ordered_names = self._prepare_shared_trait_data(
            tree_tips,
            traits,
        )
        tree_for_analysis = self._fast_copy(tree) if to_prune else tree
        if to_prune:
            for t in to_prune:
                tree_for_analysis.prune(t)

        # Store working data as instance attributes for internal methods
        self._ordered_names = ordered_names
        self._x = np.array([traits[name] for name in self._ordered_names])
        self._n = len(self._ordered_names)
        preorder_clades = list(self._iter_preorder(tree_for_analysis.root))
        self._parent_map = self._build_parent_map(tree_for_analysis, preorder_clades)
        self._lineage_info = self._build_lineage_info_no_regime(
            tree_for_analysis, self._ordered_names, self._parent_map, preorder_clades
        )
        self._lineage_rows_by_clade_id = None
        self._tree_height = max(
            sum(bl for _, bl, _, _ in self._lineage_info[name])
            for name in self._ordered_names
        )
        self._edges = self._enumerate_edges(
            tree_for_analysis, self._parent_map, preorder_clades
        )
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
            result = self._build_final_result(null, tree_for_analysis)
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

        result = self._build_final_result(best, tree_for_analysis)
        self._output(result)

    # ── Tree & data parsing (adapted from OUwie) ─────────────────────

    @staticmethod
    def _prepare_shared_trait_data(tree_tips, traits):
        if len(traits) == len(tree_tips) and all(
            name in traits for name in tree_tips
        ):
            return [], sorted(traits)

        shared = set(traits.keys())
        return [name for name in tree_tips if name not in shared], sorted(shared)

    def _parse_trait_file(
        self, path: str, tree_tips: list[str]
    ) -> dict[str, float]:
        try:
            traits = {}
            with open(path) as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line[0] == "#":
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

        if (
            len(tree_tips) >= 3
            and len(tree_tips) == len(traits)
            and next(iter(traits)) == tree_tips[0]
            and next(reversed(traits)) == tree_tips[-1]
            and list(traits) == tree_tips
        ):
            return traits

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

    # ── Tree structure helpers ───────────────────────────────────────

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
        while stack:
            clade = stack.pop()
            clades.append(clade)
            stack.extend(clade.clades)
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

    def _build_lineage_info_no_regime(
        self, tree, ordered_names: list[str], parent_map: dict,
        preorder_clades=None,
    ) -> dict[str, list[tuple]]:
        """Build root-to-tip lineage info without regime labels.

        Returns dict mapping tip_name -> list of tuples:
            (clade_id, branch_length, dist_from_root_start, dist_from_root_end)
        """
        ordered_set = set(ordered_names)
        tip_map = {}
        if preorder_clades is None:
            preorder_clades = self._iter_preorder(tree.root)
        for clade in preorder_clades:
            if not clade.clades and clade.name in ordered_set:
                tip_map[clade.name] = clade

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

    def _enumerate_edges(
        self, tree, parent_map: dict, preorder_clades=None
    ) -> list[tuple]:
        """List all non-root edges as candidate shift locations.

        Returns list of (clade_id, clade_object) tuples for every
        non-root clade in the tree.
        """
        edges = []
        root = tree.root
        if preorder_clades is None:
            preorder_clades = self._iter_preorder(root)
        for clade in preorder_clades:
            if clade == root:
                continue
            edges.append((id(clade), clade))
        return edges

    def _count_descendants(self, tree, edges: list[tuple]) -> dict:
        """Count the number of tip descendants for each edge."""
        descendant_counts = {}
        for clade in self._iter_postorder(tree.root):
            children = clade.clades
            if not children:
                descendant_counts[id(clade)] = 1
            elif len(children) == 2:
                descendant_counts[id(clade)] = (
                    descendant_counts[id(children[0])]
                    + descendant_counts[id(children[1])]
                )
            else:
                total = 0
                for child in children:
                    total += descendant_counts[id(child)]
                descendant_counts[id(clade)] = total
        counts = {}
        for clade_id, clade in edges:
            counts[clade_id] = descendant_counts.get(clade_id, 0)
        return counts

    # ── Precomputed phylogenetic data ────────────────────────────────

    def _precompute_shared_path_lengths(self) -> tuple[np.ndarray, np.ndarray]:
        """Precompute shared path length matrix S and tip height vector.

        S[i,j] = total branch length shared between root-to-tip_i and
                  root-to-tip_j paths (= distance from root to their LCA).
        """
        S, tip_heights, rows_by_clade_id = self._shared_path_lengths_from_lineage_info(
            self._ordered_names,
            self._lineage_info,
        )
        self._lineage_rows_by_clade_id = rows_by_clade_id
        return S, tip_heights

    @staticmethod
    def _shared_path_lengths_from_lineage_info(
        ordered_names: list[str],
        lineage_info: dict,
    ) -> tuple[np.ndarray, np.ndarray, dict[int, np.ndarray]]:
        n = len(ordered_names)
        S = np.zeros((n, n))
        tip_heights = np.zeros(n)
        row_lists = {}
        branch_lengths = {}

        for row_idx, name in enumerate(ordered_names):
            tip_height = 0.0
            for clade_id, bl, *_ in lineage_info[name]:
                tip_height += bl
                row_lists.setdefault(clade_id, []).append(row_idx)
                branch_lengths[clade_id] = bl
            tip_heights[row_idx] = tip_height

        rows_by_clade_id = {}
        for clade_id, rows in row_lists.items():
            row_idx = np.asarray(rows, dtype=np.intp)
            rows_by_clade_id[clade_id] = row_idx
            if len(row_idx) == n:
                S += branch_lengths[clade_id]
            elif len(row_idx) == 1:
                S[row_idx[0], row_idx[0]] += branch_lengths[clade_id]
            else:
                S[np.ix_(row_idx, row_idx)] += branch_lengths[clade_id]

        return S, tip_heights, rows_by_clade_id

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
        self, x: np.ndarray = None, ordered_names: list[str] = None,
        lineage_info: dict = None, tree_height: float = None,
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
            ones = np.ones(n)
            _, sig2, ll = self._gls_profile_likelihood(x, ones[:, None], V)
            if sig2 <= 0 or not np.isfinite(ll):
                return 1e20
            return -ll

        alpha_hat, _ = self._optimize_parameter(neg_ll, (1e-8, upper))
        return alpha_hat

    # ── Two-pass LASSO with per-candidate model selection ────────────

    def _lasso_pass_with_selection(
        self, alpha_for_lasso: float, max_shifts: int,
    ) -> dict:
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
    ) -> list[list[int]]:
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

        col_norms = _column_l2_norms(X_resid)
        col_norms[col_norms < 1e-12] = 1.0
        X_resid_norm = X_resid / col_norms

        max_iter = min(max_shifts, p, self._n - 1)
        if max_iter < 1:
            return []

        try:
            _, _, coefs_path = _lars_path(
                X_resid_norm, y_resid, method="lasso", max_iter=max_iter,
            )
        except Exception:
            return []

        configs = []
        seen = {frozenset()}

        for step in range(coefs_path.shape[1]):
            coefs = coefs_path[:, step] / col_norms
            nonzero_idx = np.flatnonzero(np.abs(coefs) > 1e-10)
            config = frozenset(nonzero_idx.tolist())
            if config in seen or len(nonzero_idx) == 0 or len(nonzero_idx) > max_shifts:
                continue
            seen.add(config)
            configs.append(nonzero_idx.tolist())

        return configs

    # ── Per-candidate model fitting ──────────────────────────────────

    def _fit_and_score_config(self, shift_edge_indices: list[int]) -> dict:
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
            _, sig2, ll = self._gls_profile_likelihood(x, W_ind, V)
            if sig2 <= 0 or not np.isfinite(ll):
                return 1e20
            return -ll

        alpha_hat = self._optimize_alpha(neg_ll, (1e-8, upper))

        # Final fit at optimal alpha with indicator matrix
        V = self._build_ou_vcv_fast(alpha_hat)
        theta, sig2, ll = self._gls_profile_likelihood(x, W_ind, V)
        if sig2 <= 0 or not np.isfinite(ll):
            return None
        pbic = self._compute_pbic_from_vcv(n, ll, k, W_ind, V, sig2, x)
        if pbic is None:
            return None

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
            pBIC=float(pbic),
            BIC=float(self._compute_bic(n, ll, k)),
            AICc=float(self._compute_aicc(n, ll, k)),
        )

    def _fit_null_model(self) -> dict:
        """Fit the null (no-shift) OU1 model with alpha optimization."""
        n = self._n
        x = self._x

        alpha_hat = self._fit_ou1_for_alpha()

        V = self._build_ou_vcv_fast(alpha_hat)
        W = np.ones((n, 1))
        theta, sig2, ll = self._gls_profile_likelihood(x, W, V)
        if sig2 <= 0 or not np.isfinite(ll):
            return None
        pbic = self._compute_pbic_from_vcv(n, ll, 0, W, V, sig2, x)
        if pbic is None:
            return None

        return dict(
            n_shifts=0,
            shift_edge_indices=[],
            alpha=float(alpha_hat),
            sigma2=float(sig2),
            theta_root=float(theta[0]),
            optima={},
            log_likelihood=float(ll),
            pBIC=float(pbic),
            BIC=float(self._compute_bic(n, ll, 0)),
            AICc=float(self._compute_aicc(n, ll, 0)),
        )

    # ── Backward elimination (R-style) ──────────────────────────────

    def _backward_eliminate_single_pass(self, fitted: dict) -> dict:
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
        self, ordered_names: list[str], lineage_info: dict,
        alpha: float, sigma2: float,
    ) -> np.ndarray:
        """Build OU VCV matrix without regime labels."""
        S, tip_heights, _ = self._shared_path_lengths_from_lineage_info(
            ordered_names,
            lineage_info,
        )
        if alpha < 1e-10:
            return sigma2 * S

        D = tip_heights[:, None] + tip_heights[None, :] - 2.0 * S
        return (sigma2 / (2.0 * alpha)) * np.exp(-alpha * D) * (
            1.0 - np.exp(-2.0 * alpha * S)
        )

    # ── Design matrices ────────────────────────────────────────────

    def _get_lineage_rows_by_clade_id(self) -> dict[int, np.ndarray]:
        rows_by_clade_id = getattr(self, "_lineage_rows_by_clade_id", None)
        if rows_by_clade_id is not None:
            return rows_by_clade_id

        row_lists = {}
        for row_idx, name in enumerate(self._ordered_names):
            for clade_id, *_ in self._lineage_info[name]:
                row_lists.setdefault(clade_id, []).append(row_idx)

        rows_by_clade_id = {
            clade_id: np.asarray(rows, dtype=np.intp)
            for clade_id, rows in row_lists.items()
        }
        self._lineage_rows_by_clade_id = rows_by_clade_id
        return rows_by_clade_id

    def _build_shift_weight_matrix(
        self, ordered_names: list[str], lineage_info: dict,
        edges: list[tuple], alpha: float, tree_height: float,
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

            W[i, 0] = 1.0 - W[i, 1:].sum()

        return W

    def _build_indicator_design_matrix(
        self, shift_edge_indices: list[int],
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

        # Map selected shift edge indices to their clade IDs.
        shift_clade_ids = {}
        for col_idx, edge_idx in enumerate(shift_edge_indices):
            clade_id = self._edges[edge_idx][0]
            shift_clade_ids[clade_id] = col_idx + 1

        rows_by_clade_id = self._get_lineage_rows_by_clade_id()
        for clade_id, col_idx in shift_clade_ids.items():
            rows = rows_by_clade_id.get(clade_id)
            if rows is not None:
                W[rows, col_idx] = 1.0

        return W

    # ── Cholesky transform ───────────────────────────────────────────

    def _transform_to_independent(
        self, x: np.ndarray, W: np.ndarray, V: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Cholesky-transform to remove phylogenetic correlation."""
        try:
            L = np.linalg.cholesky(V)
            y_star = solve_triangular(L, x, lower=True, check_finite=False)
            X_star = solve_triangular(L, W, lower=True, check_finite=False)
        except np.linalg.LinAlgError:
            V_reg = V + np.eye(len(x)) * 1e-8 * np.trace(V) / len(x)
            L = np.linalg.cholesky(V_reg)
            y_star = solve_triangular(L, x, lower=True, check_finite=False)
            X_star = solve_triangular(L, W, lower=True, check_finite=False)

        return X_star, y_star

    # ── GLS ──────────────────────────────────────────────────────────

    def _gls_profile_likelihood(
        self, x: np.ndarray, W: np.ndarray, V: np.ndarray
    ) -> tuple[np.ndarray, float, float]:
        """Profiled GLS likelihood with Cholesky solve and inverse fallback."""
        try:
            return self._gls_profile_likelihood_cholesky(x, W, V)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._gls_profile_likelihood_inverse(x, W, V)

    def _gls_profile_likelihood_cholesky(
        self, x: np.ndarray, W: np.ndarray, V: np.ndarray
    ) -> tuple[np.ndarray, float, float]:
        n = len(x)
        factor = cho_factor(V, lower=True, check_finite=False)
        n_predictors = W.shape[1]
        solve_rhs = np.empty(
            (n, n_predictors + 1),
            dtype=np.result_type(W, x, V),
        )
        solve_rhs[:, :n_predictors] = W
        solve_rhs[:, n_predictors] = x
        solved = cho_solve(factor, solve_rhs, check_finite=False)
        V_inv_W = solved[:, :n_predictors]
        V_inv_x = solved[:, n_predictors]

        info = W.T @ V_inv_W
        rhs = W.T @ V_inv_x
        try:
            theta = np.linalg.solve(info, rhs)
        except np.linalg.LinAlgError:
            theta = np.linalg.lstsq(info, rhs, rcond=None)[0]

        e = x - W @ theta
        V_inv_e = V_inv_x - V_inv_W @ theta
        sig2 = float(e @ V_inv_e) / n
        if sig2 <= 0:
            return theta, sig2, float("-inf")

        logdet = 2.0 * float(np.log(np.diagonal(factor[0])).sum())
        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)
        return theta, sig2, float(ll)

    def _gls_profile_likelihood_inverse(
        self, x: np.ndarray, W: np.ndarray, V: np.ndarray
    ) -> tuple[np.ndarray, float, float]:
        n = len(x)
        V_inv = np.linalg.inv(V)
        theta = self._gls_theta_hat(x, W, V_inv)
        e = x - W @ theta
        sig2 = float(e @ V_inv @ e) / n
        if sig2 <= 0:
            return theta, sig2, float("-inf")

        sign, logdet = np.linalg.slogdet(V)
        if sign <= 0:
            return theta, sig2, float("-inf")

        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet + n)
        return theta, sig2, float(ll)

    def _invert_vcv(self, V: np.ndarray):
        try:
            factor = cho_factor(V, lower=True, check_finite=False)
            return cho_solve(factor, np.eye(V.shape[0]), check_finite=False)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            try:
                return np.linalg.inv(V)
            except np.linalg.LinAlgError:
                return None

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
        info = W.T @ V_inv @ W
        return self._compute_pbic_from_info(
            n, log_likelihood, k, info, sigma2, x
        )

    def _compute_pbic_from_vcv(
        self, n: int, log_likelihood: float, k: int,
        W: np.ndarray, V: np.ndarray, sigma2: float,
        x: np.ndarray,
    ):
        """Compute pBIC without forming the full covariance inverse."""
        try:
            factor = cho_factor(V, lower=True, check_finite=False)
            V_inv_W = cho_solve(factor, W, check_finite=False)
            info = W.T @ V_inv_W
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            try:
                V_inv = np.linalg.inv(V)
            except np.linalg.LinAlgError:
                return None
            info = W.T @ V_inv @ W

        return self._compute_pbic_from_info(
            n, log_likelihood, k, info, sigma2, x
        )

    def _compute_pbic_from_info(
        self, n: int, log_likelihood: float, k: int,
        info: np.ndarray, sigma2: float, x: np.ndarray,
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

        var_y = np.var(x, ddof=1)
        if var_y < 1e-300:
            return float("inf")

        if n - d <= 0:
            return float("inf")

        scale = sigma2 * (n - d) / (var_y * n)
        if scale == 0.0:
            return float("inf")

        sign_info, logdet_info = np.linalg.slogdet(info)
        sign_scale = 1.0 if scale > 0.0 or d % 2 == 0 else -1.0
        if sign_info * sign_scale <= 0:
            return float("inf")

        ld = d * np.log(abs(scale)) - logdet_info
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
        tips = Tree.calculate_terminal_names_fast(clade)
        if tips is None:
            tips = [t.name for t in clade.get_terminals()]
        if len(tips) == 1:
            return f"terminal branch to {tips[0]}"
        elif len(tips) <= 3:
            return f"stem of ({', '.join(sorted(tips))})"
        else:
            return f"stem of ({', '.join(sorted(tips[:2]))}, ... +{len(tips)-2} more)"

    def _identify_shift_edges(
        self, best: dict, edges: list[tuple], tree
    ) -> list[dict]:
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

    def _build_final_result(self, best: dict, tree) -> dict:
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

    def _empty_result(self) -> dict:
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

    def _output(self, result: dict) -> None:
        if self.json_output:
            self._print_json_output(result)
        else:
            self._print_text_output(result)

    def _print_text_output(self, result: dict) -> None:
        try:
            lines = [
                "=" * 60,
                "OU Shift Detection (l1ou)",
                "=" * 60,
                f"Number of tips:       {result['n_tips']}",
                f"Number of shifts:     {result['n_shifts']}",
                f"Selection criterion:  {result['criterion']}",
                f"Alpha (OU strength):  {result['alpha']:.6f}",
                f"Sigma² (BM rate):     {result['sigma2']:.6f}",
                f"Root optimum (θ₀):    {result['theta_root']:.6f}",
                f"Log-likelihood:       {result['log_likelihood']:.4f}",
                f"pBIC:                 {result['pBIC']:.4f}",
                f"BIC:                  {result['BIC']:.4f}",
                f"AICc:                 {result['AICc']:.4f}",
            ]

            shifts = result["shifts"]
            if shifts:
                lines.extend(["", "Detected shifts:", "-" * 60])
                append = lines.append
                for i, shift in enumerate(shifts, 1):
                    append(
                        f"  Shift {i}: {shift['description']}\n"
                        f"           New optimum: {shift['optimum']:.6f}"
                    )
            else:
                lines.extend(["", "No shifts detected — single-regime OU is best."])
            lines.append("=" * 60)
            print("\n".join(lines))
        except BrokenPipeError:
            pass

    def _print_json_output(self, result: dict) -> None:
        print_json(result, sort_keys=False)
