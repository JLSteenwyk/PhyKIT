from __future__ import annotations

import math

from .base import Tree
from ...errors import PhykitUserError


_STDTR = None
_FDTRC = None
_CHO_FACTOR = None
_CHO_SOLVE = None


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def parse_multi_trait_file(*args, **kwargs):
    from ...helpers.trait_parsing import (
        parse_multi_trait_file as _parse_multi_trait_file,
    )

    return _parse_multi_trait_file(*args, **kwargs)


def response_predictor_arrays(*args, **kwargs):
    from ...helpers.trait_parsing import (
        response_predictor_arrays as _response_predictor_arrays,
    )

    return _response_predictor_arrays(*args, **kwargs)


def subset_traits_to_ordered_shared_taxa(*args, **kwargs):
    from ...helpers.trait_parsing import (
        subset_traits_to_ordered_shared_taxa as _subset_traits_to_ordered_shared_taxa,
    )

    return _subset_traits_to_ordered_shared_taxa(*args, **kwargs)


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


def _stdtr(df: int, values):
    global _STDTR

    if _STDTR is None:
        from scipy.special import stdtr as _stdtr

        _STDTR = _stdtr

    return _STDTR(df, values)


def _t_two_tailed_p_values(t_stats: np.ndarray, df: int) -> np.ndarray:
    return 2.0 * _stdtr(df, -np.abs(t_stats))


def _t_two_tailed_p_value(t_stat: float, df: int) -> float:
    return float(2.0 * _stdtr(df, -abs(t_stat)))


def _f_sf(f_stat: float, dfn: int, dfd: int) -> float:
    if dfn == 1 and f_stat >= 0:
        return _t_two_tailed_p_value(math.sqrt(f_stat), dfd)

    global _FDTRC

    if _FDTRC is None:
        from scipy.special import fdtrc as _fdtrc

        _FDTRC = _fdtrc

    return float(_FDTRC(dfn, dfd, f_stat))


class PhylogeneticRegression(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.response = parsed["response"]
        self.predictors = parsed["predictors"]
        self.method = parsed["method"]
        self.json_output = parsed["json_output"]
        self.gene_trees_path = parsed["gene_trees_path"]

    def run(self) -> None:
        from .vcv_utils import (
            build_discordance_vcv,
            build_vcv_matrix,
            parse_gene_trees,
        )

        tree = self.read_tree_file_unmodified()
        self.validate_tree(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="phylogenetic regression",
        )

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(
            self.trait_data_path, tree_tips
        )

        trait_name_to_idx = {}
        for idx, name in enumerate(trait_names):
            if name not in trait_name_to_idx:
                trait_name_to_idx[name] = idx

        # Validate column names
        all_columns = [self.response] + list(self.predictors)
        for col in all_columns:
            if col not in trait_name_to_idx:
                raise PhykitUserError(
                    [
                        f"Column '{col}' not found in trait file.",
                        f"Available columns: {', '.join(trait_names)}",
                    ],
                    code=2,
                )

        if self.response in self.predictors:
            raise PhykitUserError(
                [
                    f"Response variable '{self.response}' must not also be a predictor.",
                ],
                code=2,
            )

        ordered_names = sorted(traits.keys())

        if self.gene_trees_path:
            gene_trees = parse_gene_trees(self.gene_trees_path)
            vcv, vcv_meta = build_discordance_vcv(tree, gene_trees, ordered_names)
            shared = vcv_meta["shared_taxa"]
            traits, ordered_names = subset_traits_to_ordered_shared_taxa(
                traits, ordered_names, shared
            )
        else:
            vcv = build_vcv_matrix(tree, ordered_names)
            vcv_meta = None

        n = len(ordered_names)
        k = len(self.predictors)

        if n <= k + 1:
            raise PhykitUserError(
                [
                    f"Insufficient observations: n={n} but need at least {k + 2} "
                    f"(more than k+1={k + 1} parameters).",
                ],
                code=2,
            )

        # Build response vector y and design matrix X (with intercept)
        resp_idx = trait_name_to_idx[self.response]
        pred_indices = [trait_name_to_idx[p] for p in self.predictors]

        y, X = response_predictor_arrays(
            traits, ordered_names, resp_idx, pred_indices
        )

        lambda_val = None
        log_likelihood_lambda = None

        if self.method == "lambda":
            from ...helpers.pgls_utils import (
                apply_lambda,
                estimate_lambda,
                max_lambda as compute_max_lambda,
            )

            max_lam = compute_max_lambda(tree)
            lambda_val, log_likelihood_lambda = estimate_lambda(
                y, X, vcv, max_lam
            )
            vcv = apply_lambda(vcv, lambda_val)

        (
            beta_hat,
            residuals,
            sigma2,
            var_beta,
            r_squared,
            adj_r_squared,
            f_stat,
            f_p_value,
            r2_total,
            r2_phylo,
            ll,
            aic,
        ) = self._fit_model(y, X, vcv, k, n)

        # Standard errors, t-stats, p-values
        se = np.sqrt(np.diag(var_beta))
        df_resid = n - k - 1
        t_stats = beta_hat / se
        p_values = _t_two_tailed_p_values(t_stats, df_resid)

        # Fitted values
        fitted = X @ beta_hat

        # Build result
        coef_names = ["(Intercept)"] + list(self.predictors)
        formula = f"{self.response} ~ {' + '.join(self.predictors)}"

        result = self._format_result(
            coef_names=coef_names,
            beta_hat=beta_hat,
            se=se,
            t_stats=t_stats,
            p_values=p_values,
            sigma2=sigma2,
            df_resid=df_resid,
            r_squared=r_squared,
            adj_r_squared=adj_r_squared,
            f_stat=f_stat,
            f_p_value=f_p_value,
            ll=ll,
            aic=aic,
            formula=formula,
            k=k,
            n=n,
            residuals=residuals,
            fitted=fitted,
            ordered_names=ordered_names,
            lambda_val=lambda_val,
            r2_total=r2_total,
            r2_phylo=r2_phylo,
        )

        if vcv_meta is not None:
            result["vcv_metadata"] = vcv_meta

        if self.json_output:
            print_json(result)
        else:
            self._print_text_output(
                coef_names=coef_names,
                beta_hat=beta_hat,
                se=se,
                t_stats=t_stats,
                p_values=p_values,
                sigma2=sigma2,
                df_resid=df_resid,
                r_squared=r_squared,
                adj_r_squared=adj_r_squared,
                f_stat=f_stat,
                f_p_value=f_p_value,
                ll=ll,
                aic=aic,
                formula=formula,
                k=k,
                n=n,
                lambda_val=lambda_val,
                r2_total=r2_total,
                r2_phylo=r2_phylo,
            )

    def process_args(self, args) -> dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            response=args.response,
            predictors=args.predictors,
            method=getattr(args, "method", "BM"),
            json_output=getattr(args, "json", False),
            gene_trees_path=getattr(args, "gene_trees", None),
        )

    def _build_vcv_matrix(
        self, tree, ordered_names: list[str]
    ) -> np.ndarray:
        from .vcv_utils import build_vcv_matrix
        return build_vcv_matrix(tree, ordered_names)

    def _fit_model(
        self,
        y: np.ndarray,
        X: np.ndarray,
        vcv: np.ndarray,
        k: int,
        n: int,
    ):
        try:
            factor = cho_factor(vcv, lower=True, check_finite=False)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._fit_model_inverse(y, X, vcv, k, n)

        n_predictors = X.shape[1]
        solve_rhs = np.empty(
            (n, n_predictors + 2),
            dtype=np.result_type(X, y, vcv),
        )
        solve_rhs[:, :n_predictors] = X
        solve_rhs[:, n_predictors] = y
        solve_rhs[:, n_predictors + 1] = 1.0
        solved = cho_solve(factor, solve_rhs, check_finite=False)
        C_inv_X = solved[:, :n_predictors]
        C_inv_y = solved[:, n_predictors]
        C_inv_ones = solved[:, n_predictors + 1]

        XtCiX = X.T @ C_inv_X
        try:
            normal_rhs = np.empty(
                (n_predictors, n_predictors + 1),
                dtype=XtCiX.dtype,
            )
            normal_rhs[:, 0] = X.T @ C_inv_y
            normal_rhs[:, 1:] = np.eye(n_predictors, dtype=XtCiX.dtype)
            normal_solved = np.linalg.solve(XtCiX, normal_rhs)
        except np.linalg.LinAlgError:
            raise PhykitUserError(
                [
                    "Singular design matrix: cannot estimate coefficients.",
                    "Check that predictors are not collinear.",
                ],
                code=2,
            )

        beta_hat = normal_solved[:, 0]
        XtCiX_inv = normal_solved[:, 1:]
        residuals = y - X @ beta_hat
        C_inv_residuals = C_inv_y - C_inv_X @ beta_hat

        df_resid = n - X.shape[1]
        sigma2 = float(residuals @ C_inv_residuals) / max(df_resid, 1)
        var_beta = sigma2 * XtCiX_inv

        ones = np.ones(n)
        y_bar_gls = float(ones @ C_inv_y) / float(ones @ C_inv_ones)
        y_centered = y - y_bar_gls
        C_inv_y_centered = C_inv_y - y_bar_gls * C_inv_ones

        ss_tot = float(y_centered @ C_inv_y_centered)
        ss_res = float(residuals @ C_inv_residuals)
        stats = self._compute_model_stats_from_sums(
            ss_tot,
            ss_res,
            y,
            residuals,
            C_inv_residuals,
            k,
            n,
        )

        sigma2_ml = float(residuals @ C_inv_residuals) / n
        logdet_C = 2.0 * float(np.log(np.diagonal(factor[0])).sum())
        ll = -0.5 * (
            n * np.log(2 * np.pi) + n * np.log(sigma2_ml) + logdet_C + n
        )
        aic = -2.0 * ll + 2.0 * (k + 2)  # k predictors + intercept + sigma2

        return beta_hat, residuals, sigma2, var_beta, *stats, ll, aic

    def _fit_model_inverse(
        self,
        y: np.ndarray,
        X: np.ndarray,
        vcv: np.ndarray,
        k: int,
        n: int,
    ):
        try:
            C_inv = np.linalg.inv(vcv)
        except np.linalg.LinAlgError:
            raise PhykitUserError(
                [
                    "Singular VCV matrix: cannot invert.",
                    "Check that the tree has valid branch lengths.",
                ],
                code=2,
            )

        from ...helpers.pgls_utils import fit_gls

        beta_hat, residuals, sigma2, var_beta = fit_gls(y, X, C_inv)
        fitted = X @ beta_hat
        stats = self._compute_model_stats(y, fitted, residuals, C_inv, k, n)

        sign, logdet_C = np.linalg.slogdet(vcv)
        sigma2_ml = float(residuals @ C_inv @ residuals) / n
        ll = -0.5 * (
            n * np.log(2 * np.pi) + n * np.log(sigma2_ml) + logdet_C + n
        )
        aic = -2.0 * ll + 2.0 * (k + 2)  # k predictors + intercept + sigma2

        return beta_hat, residuals, sigma2, var_beta, *stats, ll, aic

    def _compute_model_stats(
        self,
        y: np.ndarray,
        fitted: np.ndarray,
        residuals: np.ndarray,
        C_inv: np.ndarray,
        k: int,
        n: int,
    ) -> tuple[float, float, float, float, float, float]:
        """Compute R², adjusted R², F-statistic, F p-value, R²_total, R²_phylo.

        Uses GLS-weighted sums of squares:
          SS_res = e' C_inv e
          SS_tot = (y - y_bar_gls)' C_inv (y - y_bar_gls)
        where y_bar_gls = (1' C_inv y) / (1' C_inv 1)

        Three-way R² variance decomposition:
          R²_pred  = standard GLS R² (predictor contribution given phylogeny)
          R²_total = 1 - σ²_gls_ml / σ²_ols_ml (phylogeny + predictor)
          R²_phylo = R²_total - R²_pred (phylogeny's unique contribution)
        """
        ones = np.ones(n)
        y_bar_gls = float(ones @ C_inv @ y) / float(ones @ C_inv @ ones)
        y_centered = y - y_bar_gls

        ss_tot = float(y_centered @ C_inv @ y_centered)
        ss_res = float(residuals @ C_inv @ residuals)

        return self._compute_model_stats_from_sums(
            ss_tot,
            ss_res,
            y,
            residuals,
            C_inv @ residuals,
            k,
            n,
        )

    def _compute_model_stats_from_sums(
        self,
        ss_tot: float,
        ss_res: float,
        y: np.ndarray,
        residuals: np.ndarray,
        C_inv_residuals: np.ndarray,
        k: int,
        n: int,
    ) -> tuple[float, float, float, float, float, float]:
        if ss_tot == 0:
            r_squared = 0.0
        else:
            r_squared = 1.0 - ss_res / ss_tot

        df_resid = n - k - 1
        if n - 1 == 0:
            adj_r_squared = r_squared
        else:
            adj_r_squared = 1.0 - (1.0 - r_squared) * (n - 1) / df_resid

        # F-statistic: (SS_tot - SS_res)/k / (SS_res / df_resid)
        ss_reg = ss_tot - ss_res
        if k == 0 or ss_res == 0:
            f_stat = 0.0
            f_p_value = 1.0
        else:
            ms_reg = ss_reg / k
            ms_res = ss_res / df_resid
            f_stat = ms_reg / ms_res
            f_p_value = _f_sf(f_stat, k, df_resid)

        # Three-way R² variance decomposition
        sig2_gls_ml = float(residuals @ C_inv_residuals) / n
        sig2_ols_ml = float(np.var(y))
        if sig2_ols_ml == 0:
            r2_total = 0.0
            r2_phylo = 0.0
        else:
            r2_total = 1.0 - sig2_gls_ml / sig2_ols_ml
            r2_phylo = r2_total - r_squared

        return r_squared, adj_r_squared, f_stat, f_p_value, r2_total, r2_phylo

    def _signif_code(self, p: float) -> str:
        if p < 0.001:
            return "***"
        elif p < 0.01:
            return "**"
        elif p < 0.05:
            return "*"
        elif p < 0.1:
            return "."
        else:
            return " "

    def _print_text_output(
        self, *, coef_names, beta_hat, se, t_stats, p_values,
        sigma2, df_resid, r_squared, adj_r_squared,
        f_stat, f_p_value, ll, aic, formula, k, n, lambda_val,
        r2_total, r2_phylo,
    ) -> None:
        lines = [
            "Phylogenetic Generalized Least Squares (PGLS)",
            f"\nFormula: {formula}",
        ]

        if lambda_val is not None:
            lines.append(f"\nEstimated lambda: {lambda_val:.4f}")

        # Header
        lines.extend(
            [
                "\nCoefficients:",
                (
                    f"{'':20s}{'Estimate':>12s}{'Std.Error':>12s}"
                    f"{'t-value':>12s}{'p-value':>12s}"
                ),
            ]
        )

        for name, beta, std_error, t_stat, p_value in zip(
            coef_names, beta_hat, se, t_stats, p_values
        ):
            sig = self._signif_code(p_value)
            lines.append(
                f"{name:20s}{beta:12.4f}{std_error:12.4f}"
                f"{t_stat:12.4f}{p_value:12.6f}    {sig}"
            )

        rse = np.sqrt(sigma2)
        lines.extend(
            [
                "---",
                "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1",
                (
                    f"\nResidual standard error: {rse:.4f} on {df_resid} "
                    "degrees of freedom"
                ),
                (
                    f"Multiple R-squared: {r_squared:.4f}    "
                    f"Adjusted R-squared: {adj_r_squared:.4f}"
                ),
                f"R-squared (total):   {r2_total:.4f}   (phylo + predictor)",
                f"R-squared (phylo):   {r2_phylo:.4f}   (phylogeny contribution)",
                (
                    f"F-statistic: {f_stat:.2f} on {k} and {df_resid} DF    "
                    f"p-value: {f_p_value:.6f}"
                ),
                f"Log-likelihood: {ll:.4f}    AIC: {aic:.4f}",
            ]
        )
        print("\n".join(lines))

    def _format_result(
        self, *, coef_names, beta_hat, se, t_stats, p_values,
        sigma2, df_resid, r_squared, adj_r_squared,
        f_stat, f_p_value, ll, aic, formula, k, n,
        residuals, fitted, ordered_names, lambda_val,
        r2_total, r2_phylo,
    ) -> dict:
        coefficients = {
            name: {
                "estimate": float(beta),
                "std_error": float(std_error),
                "t_value": float(t_stat),
                "p_value": float(p_value),
            }
            for name, beta, std_error, t_stat, p_value in zip(
                coef_names, beta_hat, se, t_stats, p_values
            )
        }

        result = {
            "formula": formula,
            "coefficients": coefficients,
            "residual_standard_error": float(np.sqrt(sigma2)),
            "df_residual": df_resid,
            "r_squared": float(r_squared),
            "adj_r_squared": float(adj_r_squared),
            "r_squared_total": float(r2_total),
            "r_squared_phylo": float(r2_phylo),
            "f_statistic": float(f_stat),
            "f_p_value": float(f_p_value),
            "log_likelihood": float(ll),
            "aic": float(aic),
            "n_observations": n,
            "n_predictors": k,
            "residuals": dict(zip(ordered_names, residuals.tolist())),
            "fitted_values": dict(zip(ordered_names, fitted.tolist())),
        }

        if lambda_val is not None:
            result["lambda"] = float(lambda_val)

        return result
