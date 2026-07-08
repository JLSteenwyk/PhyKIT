"""
Phylogenetic ANOVA / MANOVA using RRPP.

Implements the Residual Randomization Permutation Procedure (RRPP)
of Adams & Collyer (2018) for testing group differences while
accounting for phylogenetic non-independence.

Auto-detects univariate (ANOVA) vs multivariate (MANOVA) based on
the number of response trait columns, with optional user override.
"""
from __future__ import annotations

import sys
from operator import itemgetter

from .base import Tree
from ...errors import PhykitUserError


class _LazyNumpy:
    def __init__(self):
        self._module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = _np
            self._module = module

        attr = getattr(module, name)
        setattr(self, name, attr)
        return attr


np = _LazyNumpy()
_SOLVE_TRIANGULAR = None


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def solve_triangular(*args, **kwargs):
    global _SOLVE_TRIANGULAR
    if _SOLVE_TRIANGULAR is None:
        from scipy.linalg import solve_triangular as _SOLVE_TRIANGULAR

    return _SOLVE_TRIANGULAR(*args, **kwargs)


def _permutation_p_value_and_z(observed: float, permutations: np.ndarray) -> tuple[float, float]:
    p_value = float(np.count_nonzero(permutations >= observed) / permutations.size)
    perm_mean = float(permutations.mean())
    perm_std = float(permutations.std())
    z_score = (observed - perm_mean) / perm_std if perm_std > 0 else 0.0
    return p_value, float(z_score)


def _pillai_trace_from_eigenvalues(eigenvalues) -> float:
    return float((eigenvalues / (1.0 + eigenvalues)).sum())


def _singular_value_total_variance(singular_values) -> float:
    return float(np.dot(singular_values, singular_values))


class PhyloAnova(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.group_column = parsed["group_column"]
        self.method = parsed["method"]
        self.permutations = parsed["permutations"]
        self.pairwise = parsed["pairwise"]
        self.plot_output = parsed["plot_output"]
        self.plot_type = parsed["plot_type"]
        self.json_output = parsed["json_output"]
        self.seed = parsed["seed"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.traits,
            group_column=getattr(args, "group_column", None),
            method=getattr(args, "method", "auto"),
            permutations=getattr(args, "permutations", 1000),
            pairwise=getattr(args, "pairwise", False),
            plot_output=getattr(args, "plot_output", None),
            plot_type=getattr(args, "plot_type", "auto"),
            json_output=getattr(args, "json", False),
            seed=getattr(args, "seed", None),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        from .vcv_utils import build_vcv_matrix

        tree = self.read_tree_file_unmodified()
        if self._needs_default_branch_lengths(tree):
            tree = self._fast_copy(tree)
        self.validate_tree(tree, min_tips=3, assign_default_branch_length=1e-8, context="phylogenetic ANOVA")

        tree_tips = self.get_tip_names_from_tree(tree)
        header, traits = self._parse_trait_file(self.trait_data_path, tree_tips)

        # Identify group column
        group_col = self.group_column or header[0]
        if group_col not in header:
            raise PhykitUserError(
                [
                    f"Group column '{group_col}' not found in trait file.",
                    f"Available columns: {', '.join(header)}",
                ],
                code=2,
            )
        group_idx = header.index(group_col)

        # Separate response columns (numeric) from group column
        response_names = [h for i, h in enumerate(header) if i != group_idx]
        if not response_names:
            raise PhykitUserError(
                ["No response trait columns found (only the group column)."],
                code=2,
            )

        # Determine method
        n_responses = len(response_names)
        if self.method == "auto":
            method = "manova" if n_responses > 1 else "anova"
        else:
            method = self.method
            if method == "anova" and n_responses > 1:
                raise PhykitUserError(
                    [
                        f"ANOVA requires 1 response trait but found {n_responses}.",
                        "Use --method manova for multiple response traits.",
                    ],
                    code=2,
                )

        # Build ordered data
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)

        groups, Y = self._groups_and_response_matrix(
            traits, ordered_names, group_idx, len(header)
        )
        unique_groups = sorted(set(groups))
        n_groups = len(unique_groups)

        if n_groups < 2:
            raise PhykitUserError(
                ["At least 2 groups are required for ANOVA/MANOVA."],
                code=2,
            )

        # Build design matrix (intercept + dummy variables)
        X_full = self._build_design_matrix(groups, unique_groups, n)
        X_reduced = np.ones((n, 1))

        # Build phylogenetic VCV
        vcv = build_vcv_matrix(tree, ordered_names)

        # Cholesky transform
        try:
            Y_star, X_full_star, X_red_star = self._cholesky_transform(
                vcv, Y, X_full, X_reduced
            )
        except np.linalg.LinAlgError:
            raise PhykitUserError(
                [
                    "VCV matrix is not positive definite.",
                    "Check that the tree has valid branch lengths.",
                ],
                code=2,
            )

        # Fit models
        if method == "anova":
            result = self._run_anova(
                Y_star, X_full_star, X_red_star, n, n_groups,
            )
        else:
            result = self._run_manova(
                Y_star, X_full_star, X_red_star, n, n_groups,
            )

        # Pairwise comparisons
        pairwise_results = None
        if self.pairwise:
            pairwise_results = self._run_pairwise(
                Y_star, groups, unique_groups, X_red_star,
            )

        # Output
        self._print_results(
            method, result, unique_groups, response_names,
            pairwise_results, n, n_groups,
        )

        if self.json_output:
            self._print_json(
                method, result, unique_groups, response_names,
                pairwise_results, n, n_groups,
            )

        # Plot
        if self.plot_output:
            plot_type = self.plot_type
            if plot_type == "auto":
                plot_type = "boxplot" if method == "anova" else "phylomorphospace"
            if plot_type == "boxplot":
                self._plot_boxplot(
                    tree, ordered_names, Y, groups, unique_groups,
                    response_names, self.plot_output,
                )
            else:
                self._plot_phylomorphospace(
                    tree, ordered_names, Y, groups, unique_groups,
                    response_names, self.plot_output,
                )

    @staticmethod
    def _needs_default_branch_lengths(tree) -> bool:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return True

        stack = [root]
        pop = stack.pop
        extend = stack.extend
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if clade is not root and clade.branch_length is None:
                    return True
                if children:
                    extend(children)
        except AttributeError:
            return True

        return False

    @staticmethod
    def _groups_and_response_matrix(
        traits: dict[str, list],
        ordered_names: list[str],
        group_idx: int,
        n_cols: int,
    ) -> tuple[list[str], np.ndarray]:
        response_indices = tuple(idx for idx in range(n_cols) if idx != group_idx)
        n_taxa = len(ordered_names)
        groups = []
        groups_append = groups.append
        local_traits = traits

        if len(response_indices) == 1:
            response_idx = response_indices[0]
            Y = np.empty((n_taxa, 1), dtype=float)
            for row_idx, name in enumerate(ordered_names):
                row = local_traits[name]
                groups_append(row[group_idx])
                Y[row_idx, 0] = row[response_idx]
            return groups, Y

        get_response_values = itemgetter(*response_indices)
        Y = np.empty((n_taxa, len(response_indices)), dtype=float)
        for row_idx, name in enumerate(ordered_names):
            row = local_traits[name]
            groups_append(row[group_idx])
            Y[row_idx] = get_response_values(row)

        return groups, Y

    def _parse_trait_file(
        self, path: str, tree_tips: list[str]
    ) -> tuple[list[str], dict[str, list]]:
        """Parse TSV with header. Group column kept as string; others as float."""
        try:
            with open(path) as f:
                header_parts = None
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue
                    header_parts = stripped.split("\t")
                    break

                if header_parts is None:
                    raise PhykitUserError(
                        ["Trait file must have a header row and at least one data row."],
                        code=2,
                    )

                n_cols = len(header_parts)
                if n_cols < 3:
                    raise PhykitUserError(
                        ["Header must have at least 3 columns (taxon + group + trait)."],
                        code=2,
                    )
                header = header_parts[1:]

                group_col_name = self.group_column or header[0]
                group_idx = (
                    header.index(group_col_name)
                    if group_col_name in header
                    else None
                )
                numeric_cols = [
                    (i, col_name)
                    for i, col_name in enumerate(header)
                    if i != group_idx
                ]
                float_ = float

                traits = {}
                saw_data = False
                line_idx = 2
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue
                    saw_data = True
                    parts = stripped.split("\t")
                    if len(parts) != n_cols:
                        raise PhykitUserError(
                            [f"Line {line_idx} has {len(parts)} columns; expected {n_cols}."],
                            code=2,
                        )
                    taxon = parts[0]
                    values = [None] * len(header)
                    if group_idx is not None:
                        values[group_idx] = parts[group_idx + 1].strip()
                    for i, col_name in numeric_cols:
                        val_str = parts[i + 1]
                        try:
                            values[i] = float_(val_str)
                        except ValueError:
                            raise PhykitUserError(
                                [
                                    f"Non-numeric trait value '{val_str}' for taxon "
                                    f"'{taxon}' (trait '{col_name}') on line {line_idx}.",
                                ],
                                code=2,
                            )
                    traits[taxon] = values
                    line_idx += 1

                if not saw_data:
                    raise PhykitUserError(
                        ["Trait file must have a header row and at least one data row."],
                        code=2,
                    )
        except FileNotFoundError:
            raise PhykitUserError(
                [f"{path} corresponds to no such file or directory."],
                code=2,
            )

        tree_tip_set = set(tree_tips)
        if (
            len(tree_tip_set) >= 3
            and len(tree_tip_set) == len(traits)
            and tree_tip_set == traits.keys()
        ):
            return header, traits

        trait_taxa_set = set(traits)
        shared = tree_tip_set & trait_taxa_set

        if tree_tip_set - trait_taxa_set:
            print(
                f"Warning: {len(tree_tip_set - trait_taxa_set)} taxa in tree "
                f"but not in trait file",
                file=sys.stderr,
            )
        if trait_taxa_set - tree_tip_set:
            print(
                f"Warning: {len(trait_taxa_set - tree_tip_set)} taxa in trait "
                f"file but not in tree",
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

        filtered = {taxon: traits[taxon] for taxon in shared}
        return header, filtered

    @staticmethod
    def _build_design_matrix(
        groups: list[str], unique_groups: list[str], n: int
    ) -> np.ndarray:
        """Build design matrix: intercept + (g-1) dummy columns."""
        # Use treatment coding: intercept + dummies for groups 1..g-1
        X = np.zeros((n, len(unique_groups)))
        X[:, 0] = 1.0  # intercept
        group_to_idx = {
            group: idx for idx, group in enumerate(unique_groups)
        }
        for i, g in enumerate(groups):
            idx = group_to_idx[g]
            if idx > 0:
                X[i, idx] = 1.0
        return X

    @staticmethod
    def _cholesky_transform(
        vcv: np.ndarray, Y: np.ndarray, X_full: np.ndarray, X_reduced: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Whiten response/design matrices without explicitly inverting L."""
        L = np.linalg.cholesky(vcv)
        return (
            solve_triangular(L, Y, lower=True, check_finite=False),
            solve_triangular(L, X_full, lower=True, check_finite=False),
            solve_triangular(L, X_reduced, lower=True, check_finite=False),
        )

    @staticmethod
    def _projection_basis(X: np.ndarray) -> np.ndarray:
        """Return an orthonormal basis for repeated projections onto X."""
        Q, _ = np.linalg.qr(X, mode="reduced")
        return Q

    @staticmethod
    def _project(Q: np.ndarray, Y: np.ndarray) -> np.ndarray:
        return Q @ (Q.T @ Y)

    def _run_anova(
        self,
        Y_star: np.ndarray,
        X_full_star: np.ndarray,
        X_red_star: np.ndarray,
        n: int,
        n_groups: int,
    ) -> dict:
        """Univariate phylogenetic ANOVA with RRPP."""
        y = Y_star.ravel() if Y_star.ndim > 1 else Y_star
        if Y_star.ndim > 1:
            y = Y_star[:, 0]

        Q_red = self._projection_basis(X_red_star)
        Q_full = self._projection_basis(X_full_star)

        # Fit models
        hat_red = self._project(Q_red, y)
        resid_red = y - hat_red

        hat_full = self._project(Q_full, y)
        resid_full = y - hat_full

        ss_total = float(resid_red @ resid_red)
        ss_resid = float(resid_full @ resid_full)
        ss_model = ss_total - ss_resid

        df_model = n_groups - 1
        df_resid = n - n_groups
        ms_model = ss_model / df_model if df_model > 0 else 0.0
        ms_resid = ss_resid / df_resid if df_resid > 0 else 0.0
        f_obs = ms_model / ms_resid if ms_resid > 0 else 0.0

        # RRPP permutations
        rng = np.random.RandomState(self.seed)
        f_perms = np.zeros(self.permutations)
        for p in range(self.permutations):
            perm_idx = rng.permutation(n)
            resid_perm = resid_red[perm_idx]
            y_perm = hat_red + resid_perm

            projected_coeffs = Q_full.T @ y_perm
            ss_resid_p = float(
                y_perm @ y_perm - projected_coeffs @ projected_coeffs
            )
            ss_model_p = ss_total - ss_resid_p
            ms_model_p = ss_model_p / df_model if df_model > 0 else 0.0
            ms_resid_p = ss_resid_p / df_resid if df_resid > 0 else 0.0
            f_perms[p] = ms_model_p / ms_resid_p if ms_resid_p > 0 else 0.0

        p_value, z_score = _permutation_p_value_and_z(f_obs, f_perms)

        return dict(
            ss_model=ss_model,
            ss_resid=ss_resid,
            ss_total=ss_total,
            df_model=df_model,
            df_resid=df_resid,
            df_total=n - 1,
            ms_model=ms_model,
            ms_resid=ms_resid,
            f_stat=f_obs,
            z_score=z_score,
            p_value=p_value,
            test_statistic_name="F",
        )

    def _run_manova(
        self,
        Y_star: np.ndarray,
        X_full_star: np.ndarray,
        X_red_star: np.ndarray,
        n: int,
        n_groups: int,
    ) -> dict:
        """Multivariate phylogenetic MANOVA with RRPP using Pillai's trace."""
        Q_red = self._projection_basis(X_red_star)
        Q_full = self._projection_basis(X_full_star)

        # Fit models
        hat_red = self._project(Q_red, Y_star)
        resid_red = Y_star - hat_red

        hat_full = self._project(Q_full, Y_star)
        resid_full = Y_star - hat_full

        # SSCP matrices
        SS_total = resid_red.T @ resid_red
        SS_resid = resid_full.T @ resid_full
        SS_model = SS_total - SS_resid

        df_model = n_groups - 1
        df_resid = n - n_groups

        # Pillai's trace
        try:
            SS_resid_inv = np.linalg.inv(SS_resid)
            H_E_inv = SS_model @ SS_resid_inv
            eigenvalues = np.real(np.linalg.eigvals(H_E_inv))
            pillai = _pillai_trace_from_eigenvalues(eigenvalues)
        except np.linalg.LinAlgError:
            pillai = 0.0
            eigenvalues = np.array([0.0])

        # Approximate F for Pillai's trace
        p = Y_star.shape[1]  # number of response variables
        s = min(df_model, p)
        m = (abs(df_model - p) - 1) / 2.0
        nn = (df_resid - p - 1) / 2.0
        if s > 0 and (s + 2 * nn + 2) > 0:
            f_approx = (
                (pillai / s)
                * ((2 * nn + s + 1) / (2 * m + s + 1))
            )
        else:
            f_approx = 0.0

        # RRPP permutations
        rng = np.random.RandomState(self.seed)
        pillai_perms = np.zeros(self.permutations)
        for perm in range(self.permutations):
            perm_idx = rng.permutation(n)
            resid_perm = resid_red[perm_idx]
            Y_perm = hat_red + resid_perm

            projected_coeffs = Q_full.T @ Y_perm
            SS_resid_p = (
                Y_perm.T @ Y_perm - projected_coeffs.T @ projected_coeffs
            )
            SS_model_p = SS_total - SS_resid_p

            try:
                SS_resid_inv_p = np.linalg.inv(SS_resid_p)
                H_E_inv_p = SS_model_p @ SS_resid_inv_p
                eig_p = np.real(np.linalg.eigvals(H_E_inv_p))
                pillai_perms[perm] = _pillai_trace_from_eigenvalues(eig_p)
            except np.linalg.LinAlgError:
                pillai_perms[perm] = 0.0

        p_value, z_score = _permutation_p_value_and_z(pillai, pillai_perms)

        # Scalar SS for display (trace of SSCP)
        ss_model = float(np.trace(SS_model))
        ss_resid = float(np.trace(SS_resid))
        ss_total = float(np.trace(SS_total))

        return dict(
            ss_model=ss_model,
            ss_resid=ss_resid,
            ss_total=ss_total,
            df_model=df_model,
            df_resid=df_resid,
            df_total=n - 1,
            ms_model=ss_model / df_model if df_model > 0 else 0.0,
            ms_resid=ss_resid / df_resid if df_resid > 0 else 0.0,
            f_stat=f_approx,
            pillai_trace=pillai,
            z_score=z_score,
            p_value=p_value,
            test_statistic_name="Pillai's trace",
        )

    def _run_pairwise(
        self,
        Y_star: np.ndarray,
        groups: list[str],
        unique_groups: list[str],
        X_red_star: np.ndarray,
    ) -> list[dict]:
        """Pairwise group comparisons via RRPP."""
        rng = np.random.RandomState(self.seed)
        results = []

        # Reduced model residuals
        beta_red = np.linalg.lstsq(X_red_star, Y_star, rcond=None)[0]
        hat_red = X_red_star @ beta_red
        resid_red = Y_star - hat_red

        groups_arr = np.array(groups)

        for i in range(len(unique_groups)):
            for j in range(i + 1, len(unique_groups)):
                g1, g2 = unique_groups[i], unique_groups[j]
                mask = (groups_arr == g1) | (groups_arr == g2)
                idx = np.flatnonzero(mask)

                Y_sub = Y_star[idx]
                groups_sub = groups_arr[idx]
                pos1 = np.flatnonzero(groups_sub == g1)
                pos2 = np.flatnonzero(groups_sub == g2)

                # Observed distance between group means
                mean1 = Y_sub[pos1].mean(axis=0)
                mean2 = Y_sub[pos2].mean(axis=0)
                mean_diff = mean1 - mean2
                d_obs = float(np.sqrt(mean_diff @ mean_diff))

                # Permute residuals for pairwise test
                resid_sub = resid_red[idx]
                hat_sub = hat_red[idx]
                hat_diff = (
                    hat_sub[pos1].mean(axis=0) - hat_sub[pos2].mean(axis=0)
                )

                perm_indices = np.empty(
                    (self.permutations, len(idx)), dtype=np.intp
                )
                for p in range(self.permutations):
                    perm_indices[p] = rng.permutation(len(idx))

                contrast = np.zeros(len(idx), dtype=resid_sub.dtype)
                contrast[pos1] = 1.0 / len(pos1)
                contrast[pos2] = -1.0 / len(pos2)
                weights = np.empty(
                    (self.permutations, len(idx)), dtype=resid_sub.dtype
                )
                weights[np.arange(self.permutations)[:, None], perm_indices] = (
                    contrast
                )
                resid_diffs = weights @ resid_sub
                perm_diffs = hat_diff + resid_diffs
                d_perms = np.sqrt(np.einsum("ij,ij->i", perm_diffs, perm_diffs))

                p_val, z_val = _permutation_p_value_and_z(d_obs, d_perms)

                results.append(dict(
                    group1=g1, group2=g2,
                    distance=d_obs, z_score=z_val, p_value=p_val,
                ))

        return results

    def _print_results(
        self, method, result, unique_groups, response_names,
        pairwise_results, n, n_groups,
    ) -> None:
        method_label = "Phylogenetic ANOVA" if method == "anova" else "Phylogenetic MANOVA"
        try:
            lines = [f"{method_label} (RRPP, {self.permutations} permutations)"]
            append = lines.append
            if method == "anova":
                append(f"Response: {response_names[0]}")
            else:
                append(f"Response traits: {', '.join(response_names)}")
            append(f"Groups ({n_groups}): {', '.join(unique_groups)}")
            append(f"N taxa: {n}")
            append("")
            append(
                f"{'Source':<15} {'Df':>5} {'SS':>12} {'MS':>12} "
                f"{'F':>10} {'Z':>10} {'p-value':>10}"
            )
            append("-" * 76)
            append(
                f"{'group':<15} {result['df_model']:>5d} "
                f"{result['ss_model']:>12.4f} {result['ms_model']:>12.4f} "
                f"{result['f_stat']:>10.4f} {result['z_score']:>10.4f} "
                f"{result['p_value']:>10.4f}"
            )
            append(
                f"{'residuals':<15} {result['df_resid']:>5d} "
                f"{result['ss_resid']:>12.4f} {result['ms_resid']:>12.4f}"
            )
            append(
                f"{'total':<15} {result['df_total']:>5d} "
                f"{result['ss_total']:>12.4f}"
            )

            if "pillai_trace" in result:
                append("")
                append(f"Pillai's trace: {result['pillai_trace']:.4f}")

            if pairwise_results:
                append("")
                append("Pairwise comparisons:")
                append(
                    f"  {'Comparison':<30} {'d':>10} {'Z':>10} {'p-value':>10}"
                )
                append("  " + "-" * 62)
                for pw in pairwise_results:
                    label = f"{pw['group1']} vs {pw['group2']}"
                    append(
                        f"  {label:<30} {pw['distance']:>10.4f} "
                        f"{pw['z_score']:>10.4f} {pw['p_value']:>10.4f}"
                    )
            print("\n".join(lines))
        except BrokenPipeError:
            return

    def _print_json(
        self, method, result, unique_groups, response_names,
        pairwise_results, n, n_groups,
    ) -> None:
        payload = {
            "method": method,
            "permutations": self.permutations,
            "n_taxa": n,
            "n_groups": n_groups,
            "groups": unique_groups,
            "response_traits": response_names,
            "anova_table": {
                "df_model": result["df_model"],
                "df_resid": result["df_resid"],
                "df_total": result["df_total"],
                "ss_model": round(result["ss_model"], 6),
                "ss_resid": round(result["ss_resid"], 6),
                "ss_total": round(result["ss_total"], 6),
                "ms_model": round(result["ms_model"], 6),
                "ms_resid": round(result["ms_resid"], 6),
                "f_stat": round(result["f_stat"], 6),
                "z_score": round(result["z_score"], 4),
                "p_value": round(result["p_value"], 4),
            },
        }
        if "pillai_trace" in result:
            payload["anova_table"]["pillai_trace"] = round(
                result["pillai_trace"], 6
            )
        if pairwise_results:
            payload["pairwise"] = [
                {
                    "group1": pw["group1"],
                    "group2": pw["group2"],
                    "distance": round(pw["distance"], 6),
                    "z_score": round(pw["z_score"], 4),
                    "p_value": round(pw["p_value"], 4),
                }
                for pw in pairwise_results
            ]
        print_json(payload, sort_keys=False)

    def _plot_boxplot(
        self, tree, ordered_names, Y, groups, unique_groups,
        response_names, output_path,
    ) -> None:
        """Violin + boxplot for univariate ANOVA."""
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.collections import LineCollection
        except ImportError:
            print("matplotlib is required for plotting.")
            return

        config = self.plot_config
        config.resolve(n_rows=len(unique_groups), n_cols=None)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        y = Y[:, 0] if Y.ndim > 1 else Y
        default_colors = ["#2b8cbe", "#d62728", "#969696", "#31a354",
                          "#756bb1", "#e6550d", "#636363", "#de2d26"]
        colors = config.merge_colors(default_colors)

        data_by_group = []
        group_labels = []
        group_colors = []
        groups_arr = np.asarray(groups)
        for i, g in enumerate(unique_groups):
            vals = y[groups_arr == g]
            data_by_group.append(vals)
            group_labels.append(g)
            group_colors.append(colors[i % len(colors)])

        positions = list(range(1, len(unique_groups) + 1))

        # Violin
        vp = ax.violinplot(
            data_by_group, positions=positions,
            showmeans=False, showmedians=False, showextrema=False,
        )
        for i, body in enumerate(vp["bodies"]):
            body.set_facecolor(group_colors[i])
            body.set_alpha(0.3)
            body.set_edgecolor(group_colors[i])

        # Boxplot
        bp = ax.boxplot(
            data_by_group, positions=positions,
            widths=0.15, patch_artist=True,
            showfliers=True, zorder=3,
        )
        for i, box in enumerate(bp["boxes"]):
            box.set_facecolor(group_colors[i])
            box.set_alpha(0.7)
        for element in ["whiskers", "caps", "medians"]:
            for line in bp[element]:
                line.set_color("black")

        ax.set_xticks(positions)
        ax.set_xticklabels(group_labels)
        label_fs = config.xlabel_fontsize or 10
        ax.tick_params(axis="x", labelsize=label_fs)

        ylabel = response_names[0] if response_names else "Trait value"
        ylabel_fs = config.ylabel_fontsize if config.ylabel_fontsize else 10
        ax.set_ylabel(ylabel, fontsize=ylabel_fs)

        if config.show_title:
            ax.set_title(
                config.title or "Phylogenetic ANOVA",
                fontsize=config.title_fontsize,
            )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.subplots_adjust(
            left=0.10,
            right=0.97,
            top=0.90 if config.show_title else 0.96,
            bottom=0.14,
        )
        fig.savefig(output_path, dpi=config.dpi)
        plt.close(fig)
        print(f"Plot saved: {output_path}")

    def _plot_phylomorphospace(
        self, tree, ordered_names, Y, groups, unique_groups,
        response_names, output_path,
    ) -> None:
        """Phylomorphospace colored by group for multivariate MANOVA."""
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.collections import LineCollection
        except ImportError:
            print("matplotlib is required for plotting.")
            return

        config = self.plot_config
        config.resolve(n_rows=len(ordered_names), n_cols=None)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        default_colors = ["#2b8cbe", "#d62728", "#969696", "#31a354",
                          "#756bb1", "#e6550d", "#636363", "#de2d26"]
        colors = config.merge_colors(default_colors)

        # Use first 2 traits as axes (or PCA if more)
        if Y.shape[1] == 2:
            pc = Y
            xlabel = response_names[0]
            ylabel = response_names[1]
        else:
            # PCA
            Y_centered = Y - Y.mean(axis=0)
            U, S, Vt = np.linalg.svd(Y_centered, full_matrices=False)
            pc = Y_centered @ Vt[:2].T
            total_var = _singular_value_total_variance(S)
            var_explained = S[:2] ** 2 / total_var * 100
            xlabel = f"PC1 ({var_explained[0]:.1f}%)"
            ylabel = f"PC2 ({var_explained[1]:.1f}%)"

        parent_map, node_coords, preorder_clades = (
            self._prepare_phylomorphospace_overlay(tree, ordered_names, pc)
        )

        # Draw phylogeny branches
        branch_segments = []
        for clade in preorder_clades:
            if clade == tree.root:
                continue
            pid = parent_map.get(id(clade))
            if pid is None or id(pid) not in node_coords or id(clade) not in node_coords:
                continue
            p_coord = node_coords[id(pid)]
            c_coord = node_coords[id(clade)]
            branch_segments.append([(p_coord[0], p_coord[1]), (c_coord[0], c_coord[1])])
        if branch_segments:
            ax.add_collection(
                LineCollection(
                    branch_segments,
                    colors="gray",
                    linewidths=0.5,
                    alpha=0.5,
                    zorder=1,
                ),
                autolim=True,
        )

        # Plot points colored by group
        groups_arr = np.asarray(groups)
        for i, g in enumerate(unique_groups):
            mask = groups_arr == g
            color = colors[i % len(colors)]
            ax.scatter(
                pc[mask, 0], pc[mask, 1],
                c=color, label=g, s=50, edgecolors="black",
                linewidths=0.5, zorder=3,
            )

        ax.set_xlabel(xlabel, fontsize=config.axis_fontsize or 10)
        ax.set_ylabel(ylabel, fontsize=config.axis_fontsize or 10)
        ax.legend(fontsize=8, frameon=True)

        if config.show_title:
            ax.set_title(
                config.title or "Phylogenetic MANOVA — Phylomorphospace",
                fontsize=config.title_fontsize,
            )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.subplots_adjust(
            left=0.10,
            right=0.97,
            top=0.90 if config.show_title else 0.96,
            bottom=0.13,
        )
        fig.savefig(output_path, dpi=config.dpi)
        plt.close(fig)
        print(f"Plot saved: {output_path}")

    @staticmethod
    def _prepare_phylomorphospace_overlay(tree, ordered_names, pc):
        direct_result = PhyloAnova._prepare_phylomorphospace_overlay_direct(
            tree, ordered_names, pc
        )
        if direct_result is not None:
            return direct_result

        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade

        # Map tip names to indices
        name_to_idx = {name: i for i, name in enumerate(ordered_names)}

        # Simple ancestral reconstruction (weighted mean of descendants)
        node_coords = {}
        for tip in tree.get_terminals():
            if tip.name in name_to_idx:
                idx = name_to_idx[tip.name]
                node_coords[id(tip)] = pc[idx]

        for clade in tree.find_clades(order="postorder"):
            if id(clade) in node_coords:
                continue
            child_coords = [
                node_coords[id(c)] for c in clade.clades
                if id(c) in node_coords
            ]
            if child_coords:
                node_coords[id(clade)] = np.mean(child_coords, axis=0)

        return parent_map, node_coords, list(tree.find_clades(order="preorder"))

    @staticmethod
    def _prepare_phylomorphospace_overlay_direct(tree, ordered_names, pc):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        preorder_clades = []
        parent_map = {}
        stack = [root]
        try:
            while stack:
                clade = stack.pop()
                preorder_clades.append(clade)
                children = clade.clades
                if children:
                    for child in children:
                        parent_map[id(child)] = clade
                    child_count = len(children)
                    if child_count == 2:
                        stack.append(children[1])
                        stack.append(children[0])
                    else:
                        for idx in range(child_count - 1, -1, -1):
                            stack.append(children[idx])
        except AttributeError:
            return None

        name_to_idx = {name: i for i, name in enumerate(ordered_names)}
        node_coords = {}
        node_coords_get = node_coords.get
        id_ = id
        for clade in reversed(preorder_clades):
            children = clade.clades
            if not children:
                idx = name_to_idx.get(clade.name)
                if idx is not None:
                    node_coords[id_(clade)] = pc[idx]
                continue

            child_count = len(children)
            if child_count == 2:
                left = node_coords_get(id_(children[0]))
                right = node_coords_get(id_(children[1]))
                if left is not None:
                    if right is not None:
                        node_coords[id_(clade)] = (left + right) * 0.5
                    else:
                        node_coords[id_(clade)] = left.copy()
                elif right is not None:
                    node_coords[id_(clade)] = right.copy()
                continue
            if child_count == 1:
                coord = node_coords_get(id_(children[0]))
                if coord is not None:
                    node_coords[id_(clade)] = coord.copy()
                continue

            coord_sum = None
            count = 0
            for child in children:
                coord = node_coords_get(id_(child))
                if coord is None:
                    continue
                if coord_sum is None:
                    coord_sum = coord.copy()
                else:
                    coord_sum += coord
                count += 1
            if count:
                node_coords[id_(clade)] = coord_sum / count

        return parent_map, node_coords, preorder_clades
