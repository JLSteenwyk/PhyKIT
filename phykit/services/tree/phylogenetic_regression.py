import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import t as t_dist
from scipy.stats import f as f_dist

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


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
        from .vcv_utils import build_vcv_matrix, build_discordance_vcv, parse_gene_trees

        tree = self.read_tree_file()
        self._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_names, traits = self._parse_multi_trait_file(
            self.trait_data_path, tree_tips
        )

        # Validate column names
        all_columns = [self.response] + list(self.predictors)
        for col in all_columns:
            if col not in trait_names:
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
            if set(shared) != set(ordered_names):
                traits = {k: traits[k] for k in shared}
                ordered_names = shared
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
        resp_idx = trait_names.index(self.response)
        pred_indices = [trait_names.index(p) for p in self.predictors]

        y = np.array([traits[name][resp_idx] for name in ordered_names])
        X_pred = np.array(
            [[traits[name][j] for j in pred_indices] for name in ordered_names]
        )
        # Design matrix: intercept + predictors
        X = np.column_stack([np.ones(n), X_pred])

        lambda_val = None
        log_likelihood_lambda = None

        if self.method == "lambda":
            max_lam = self._max_lambda(tree)
            lambda_val, log_likelihood_lambda = self._estimate_lambda(
                y, X, vcv, max_lam
            )
            # Transform VCV
            diag_vals = np.diag(vcv).copy()
            vcv = vcv * lambda_val
            np.fill_diagonal(vcv, diag_vals)

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

        # GLS estimation
        beta_hat, residuals, sigma2, var_beta = self._fit_gls(y, X, C_inv)

        # Standard errors, t-stats, p-values
        se = np.sqrt(np.diag(var_beta))
        df_resid = n - k - 1
        t_stats = beta_hat / se
        p_values = 2.0 * t_dist.sf(np.abs(t_stats), df=df_resid)

        # Fitted values
        fitted = X @ beta_hat

        # Model statistics
        r_squared, adj_r_squared, f_stat, f_p_value, r2_total, r2_phylo = self._compute_model_stats(
            y, fitted, residuals, C_inv, k, n
        )

        # Log-likelihood and AIC
        sign, logdet_C = np.linalg.slogdet(vcv)
        sigma2_ml = float(residuals @ C_inv @ residuals) / n
        ll = -0.5 * (
            n * np.log(2 * np.pi) + n * np.log(sigma2_ml) + logdet_C + n
        )
        aic = -2.0 * ll + 2.0 * (k + 2)  # k predictors + intercept + sigma2

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

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            response=args.response,
            predictors=args.predictors,
            method=getattr(args, "method", "BM"),
            json_output=getattr(args, "json", False),
            gene_trees_path=getattr(args, "gene_trees", None),
        )

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for phylogenetic regression."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                raise PhykitUserError(
                    ["All branches in the tree must have lengths."],
                    code=2,
                )

    def _parse_multi_trait_file(
        self, path: str, tree_tips: List[str]
    ) -> Tuple[List[str], Dict[str, List[float]]]:
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

        # Filter out comments and blank lines
        data_lines = []
        for line in lines:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            data_lines.append(stripped)

        if len(data_lines) < 2:
            raise PhykitUserError(
                [
                    "Multi-trait file must have a header row and at least one data row.",
                ],
                code=2,
            )

        # First data line is the header
        header_parts = data_lines[0].split("\t")
        n_cols = len(header_parts)
        if n_cols < 2:
            raise PhykitUserError(
                [
                    "Header must have at least 2 columns (taxon + at least 1 trait).",
                ],
                code=2,
            )
        trait_names = header_parts[1:]

        traits = {}
        for line_idx, line in enumerate(data_lines[1:], 2):
            parts = line.split("\t")
            if len(parts) != n_cols:
                raise PhykitUserError(
                    [
                        f"Line {line_idx} has {len(parts)} columns; expected {n_cols}.",
                        f"Each line should have: taxon_name<tab>{'<tab>'.join(['trait'] * len(trait_names))}",
                    ],
                    code=2,
                )
            taxon = parts[0]
            values = []
            for i, val_str in enumerate(parts[1:]):
                try:
                    values.append(float(val_str))
                except ValueError:
                    raise PhykitUserError(
                        [
                            f"Non-numeric trait value '{val_str}' for taxon '{taxon}' "
                            f"(trait '{trait_names[i]}') on line {line_idx}.",
                        ],
                        code=2,
                    )
            traits[taxon] = values

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

        filtered = {taxon: traits[taxon] for taxon in shared}
        return trait_names, filtered

    def _build_vcv_matrix(
        self, tree, ordered_names: List[str]
    ) -> np.ndarray:
        from .vcv_utils import build_vcv_matrix
        return build_vcv_matrix(tree, ordered_names)

    def _max_lambda(self, tree) -> float:
        tips = tree.get_terminals()
        root = tree.root
        tip_heights = [tree.distance(root, tip) for tip in tips]
        max_tip_height = max(tip_heights)
        min_tip_height = min(tip_heights)

        is_ultrametric = (max_tip_height - min_tip_height) / max_tip_height < 1e-6

        if not is_ultrametric:
            return 1.0

        max_parent_height = 0.0
        for clade in tree.find_clades(order="level"):
            if clade == root:
                continue
            node_height = tree.distance(root, clade)
            parent_height = node_height - (clade.branch_length or 0.0)
            if parent_height > max_parent_height:
                max_parent_height = parent_height

        if max_parent_height == 0.0:
            return 1.0

        return max_tip_height / max_parent_height

    def _fit_gls(
        self, y: np.ndarray, X: np.ndarray, C_inv: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, float, np.ndarray]:
        """GLS estimation: beta_hat = (X' C_inv X)^{-1} X' C_inv y

        Returns (beta_hat, residuals, sigma2_reml, var_beta).
        """
        n, k_plus_1 = X.shape
        XtCiX = X.T @ C_inv @ X
        try:
            XtCiX_inv = np.linalg.inv(XtCiX)
        except np.linalg.LinAlgError:
            raise PhykitUserError(
                [
                    "Singular design matrix: cannot estimate coefficients.",
                    "Check that predictors are not collinear.",
                ],
                code=2,
            )

        beta_hat = XtCiX_inv @ X.T @ C_inv @ y
        residuals = y - X @ beta_hat

        # REML-style sigma^2 for inference
        df_resid = n - k_plus_1
        sigma2 = float(residuals @ C_inv @ residuals) / df_resid

        var_beta = sigma2 * XtCiX_inv

        return beta_hat, residuals, sigma2, var_beta

    def _pgls_log_likelihood(
        self, y: np.ndarray, X: np.ndarray, C: np.ndarray
    ) -> float:
        """Concentrated log-likelihood with beta and sigma^2 profiled out."""
        n = len(y)
        C_inv = np.linalg.inv(C)

        XtCiX = X.T @ C_inv @ X
        XtCiX_inv = np.linalg.inv(XtCiX)
        beta_hat = XtCiX_inv @ X.T @ C_inv @ y
        e = y - X @ beta_hat

        sigma2_ml = float(e @ C_inv @ e) / n

        sign, logdet_C = np.linalg.slogdet(C)
        if sign <= 0 or sigma2_ml <= 0:
            return -1e20

        ll = -0.5 * (
            n * np.log(2 * np.pi) + n * np.log(sigma2_ml) + logdet_C + n
        )
        return float(ll)

    def _estimate_lambda(
        self,
        y: np.ndarray,
        X: np.ndarray,
        vcv: np.ndarray,
        max_lambda: float,
    ) -> Tuple[float, float]:
        """Optimize Pagel's lambda via ML using multi-interval bounded search."""
        diag_vals = np.diag(vcv).copy()
        niter = 10

        def neg_ll(lam):
            C_lam = vcv * lam
            np.fill_diagonal(C_lam, diag_vals)
            try:
                ll = self._pgls_log_likelihood(y, X, C_lam)
                return -ll
            except (np.linalg.LinAlgError, FloatingPointError, ValueError):
                return 1e10

        bounds_lo = np.linspace(0, max_lambda - max_lambda / niter, niter)
        bounds_hi = np.linspace(max_lambda / niter, max_lambda, niter)

        best_ll = -np.inf
        lambda_hat = 0.0
        for lo, hi in zip(bounds_lo, bounds_hi):
            res = minimize_scalar(neg_ll, bounds=(lo, hi), method="bounded")
            ll_val = -res.fun
            if ll_val > best_ll:
                best_ll = ll_val
                lambda_hat = res.x

        # Compute log-likelihood at fitted lambda
        C_fitted = vcv * lambda_hat
        np.fill_diagonal(C_fitted, diag_vals)
        ll_fitted = self._pgls_log_likelihood(y, X, C_fitted)

        return float(lambda_hat), float(ll_fitted)

    def _compute_model_stats(
        self,
        y: np.ndarray,
        fitted: np.ndarray,
        residuals: np.ndarray,
        C_inv: np.ndarray,
        k: int,
        n: int,
    ) -> Tuple[float, float, float, float, float, float]:
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
            f_p_value = float(f_dist.sf(f_stat, dfn=k, dfd=df_resid))

        # Three-way R² variance decomposition
        sig2_gls_ml = float(residuals @ C_inv @ residuals) / n
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
        print("Phylogenetic Generalized Least Squares (PGLS)")
        print(f"\nFormula: {formula}")

        if lambda_val is not None:
            print(f"\nEstimated lambda: {lambda_val:.4f}")

        print("\nCoefficients:")
        # Header
        print(
            f"{'':20s}{'Estimate':>12s}{'Std.Error':>12s}"
            f"{'t-value':>12s}{'p-value':>12s}"
        )

        for i, name in enumerate(coef_names):
            sig = self._signif_code(p_values[i])
            print(
                f"{name:20s}{beta_hat[i]:12.4f}{se[i]:12.4f}"
                f"{t_stats[i]:12.4f}{p_values[i]:12.6f}    {sig}"
            )

        print("---")
        print("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1")

        rse = np.sqrt(sigma2)
        print(
            f"\nResidual standard error: {rse:.4f} on {df_resid} degrees of freedom"
        )
        print(
            f"Multiple R-squared: {r_squared:.4f}    "
            f"Adjusted R-squared: {adj_r_squared:.4f}"
        )
        print(f"R-squared (total):   {r2_total:.4f}   (phylo + predictor)")
        print(f"R-squared (phylo):   {r2_phylo:.4f}   (phylogeny contribution)")
        print(
            f"F-statistic: {f_stat:.2f} on {k} and {df_resid} DF    "
            f"p-value: {f_p_value:.6f}"
        )
        print(f"Log-likelihood: {ll:.4f}    AIC: {aic:.4f}")

    def _format_result(
        self, *, coef_names, beta_hat, se, t_stats, p_values,
        sigma2, df_resid, r_squared, adj_r_squared,
        f_stat, f_p_value, ll, aic, formula, k, n,
        residuals, fitted, ordered_names, lambda_val,
        r2_total, r2_phylo,
    ) -> Dict:
        coefficients = {}
        for i, name in enumerate(coef_names):
            coefficients[name] = {
                "estimate": float(beta_hat[i]),
                "std_error": float(se[i]),
                "t_value": float(t_stats[i]),
                "p_value": float(p_values[i]),
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
            "residuals": {
                ordered_names[i]: float(residuals[i])
                for i in range(len(ordered_names))
            },
            "fitted_values": {
                ordered_names[i]: float(fitted[i])
                for i in range(len(ordered_names))
            },
        }

        if lambda_val is not None:
            result["lambda"] = float(lambda_val)

        return result
