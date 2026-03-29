import copy
import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm as norm_dist

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class PhyloLogistic(Tree):
    """Phylogenetic logistic regression (Ives & Garland 2010).

    Fits a logistic regression for binary (0/1) response data on a
    phylogenetic tree using maximum penalized likelihood estimation
    (logistic_MPLE) with Firth's bias-correction penalty.

    The phylogenetic correlation is modelled through the OU-transformed
    variance-covariance matrix and a jointly estimated alpha parameter
    that controls phylogenetic signal.
    """

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.response = parsed["response"]
        self.predictors = parsed["predictors"]
        self.method = parsed["method"]
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        from .vcv_utils import build_vcv_matrix

        tree = self.read_tree_file()
        self.validate_tree(tree, min_tips=3, require_branch_lengths=True, context="phylogenetic logistic regression")

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
        X = np.column_stack([np.ones(n), X_pred])

        # Validate binary response
        unique_vals = set(y)
        if not unique_vals.issubset({0.0, 1.0}):
            raise PhykitUserError(
                [
                    f"Response variable '{self.response}' must contain only "
                    f"0 and 1 values for logistic regression.",
                    f"Found values: {sorted(unique_vals)}",
                ],
                code=2,
            )

        if len(unique_vals) < 2:
            raise PhykitUserError(
                [
                    f"Response variable '{self.response}' must contain both "
                    f"0 and 1 values.",
                    f"Found only: {sorted(unique_vals)}",
                ],
                code=2,
            )

        # Fit the model
        result = self._fit(tree, y, X, ordered_names, build_vcv_matrix)

        if self.json_output:
            print_json(result)
        else:
            self._print_text_output(result)

    def process_args(self, args) -> Dict[str, str]:
        predictors_raw = getattr(args, "predictor", "")
        if isinstance(predictors_raw, list):
            predictors = predictors_raw
        else:
            predictors = [p.strip() for p in predictors_raw.split(",") if p.strip()]
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            response=args.response,
            predictors=predictors,
            method=getattr(args, "method", "logistic_MPLE"),
            json_output=getattr(args, "json", False),
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
                        f"Each line should have: taxon_name<tab>"
                        f"{'<tab>'.join(['trait'] * len(trait_names))}",
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

    # ----------------------------------------------------------------
    # VCV construction with OU-transformed branch lengths
    # ----------------------------------------------------------------

    def _build_logistic_vcv(
        self, tree, alpha: float, ordered_names: List[str], build_vcv_fn
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Build the OU-transformed VCV + diagonal correction.

        For a given alpha, transform the phylogenetic tree:
        - For each branch with length bl, the transformed length is:
          (1 - exp(-2*alpha*bl)) / (2*alpha)
        - Build the VCV matrix from these transformed branch lengths
        - Compute diagonal correction: exp(-2*alpha*dist_root_to_tip_i)
        - Add diagonal correction to VCV diagonal

        Returns (vcv, diag_correction).
        """
        if alpha < 1e-10:
            # BM limit
            vcv = build_vcv_fn(tree, ordered_names)
            diag_corr = np.ones(len(ordered_names))
            return vcv, diag_corr

        # Transform branch lengths
        tree_t = copy.deepcopy(tree)
        for clade in tree_t.find_clades():
            if clade.branch_length is not None and clade.branch_length > 0:
                bl = clade.branch_length
                clade.branch_length = (1 - np.exp(-2 * alpha * bl)) / (2 * alpha)

        vcv = build_vcv_fn(tree_t, ordered_names)

        # Diagonal correction
        diag_corr = np.zeros(len(ordered_names))
        for i, name in enumerate(ordered_names):
            dist = tree.distance(tree.root, name)  # original tree distances
            diag_corr[i] = np.exp(-2 * alpha * dist)

        # Add diagonal correction to VCV
        np.fill_diagonal(vcv, np.diag(vcv) + diag_corr)

        return vcv, diag_corr

    def _mean_root_tip_distance(self, tree, ordered_names: List[str]) -> float:
        """Mean root-to-tip distance across all tips."""
        dists = [tree.distance(tree.root, name) for name in ordered_names]
        return float(np.mean(dists))

    # ----------------------------------------------------------------
    # Starting values
    # ----------------------------------------------------------------

    def _logistic_starting_values(
        self, y: np.ndarray, X: np.ndarray, btol: float = 10
    ) -> np.ndarray:
        """Standard logistic regression via IRLS as starting values."""
        n, p = X.shape
        beta = np.zeros(p)

        for _ in range(50):
            eta = X @ beta
            eta = np.clip(eta, -btol, btol)
            mu = 1.0 / (1.0 + np.exp(-eta))
            mu = np.clip(mu, 1e-6, 1 - 1e-6)
            W = np.diag(mu * (1 - mu))
            z = eta + (y - mu) / (mu * (1 - mu))
            try:
                beta_new = np.linalg.solve(X.T @ W @ X, X.T @ W @ z)
            except np.linalg.LinAlgError:
                break
            if np.max(np.abs(beta_new - beta)) < 1e-8:
                beta = beta_new
                break
            beta = beta_new

        return beta

    # ----------------------------------------------------------------
    # Penalized log-likelihood (MPLE)
    # ----------------------------------------------------------------

    def _compute_info_matrix(
        self, X: np.ndarray, mu: np.ndarray, vcv: np.ndarray
    ) -> np.ndarray:
        """Compute the Fisher information matrix for phylogenetic logistic regression.

        I(beta) = X^T @ W^{1/2} @ V^{-1} @ W^{1/2} @ X
        where W = diag(mu*(1-mu)) and V is the OU-transformed VCV.
        """
        W_sqrt = np.diag(np.sqrt(mu * (1 - mu)))
        try:
            V_inv = np.linalg.inv(vcv)
        except np.linalg.LinAlgError:
            return np.eye(X.shape[1]) * 1e-10

        return X.T @ W_sqrt @ V_inv @ W_sqrt @ X

    def _neg_pen_loglik(
        self,
        params: np.ndarray,
        X: np.ndarray,
        y: np.ndarray,
        tree,
        ordered_names: List[str],
        build_vcv_fn,
        btol: float = 10,
    ) -> float:
        """Negative penalized log-likelihood for logistic MPLE.

        params = [beta_0, beta_1, ..., beta_k, log(1/alpha)]

        The penalized log-likelihood is:
            PLL = LL + 0.5 * log|I(beta)|
        where LL is the standard Bernoulli log-likelihood and I is the
        Fisher information matrix accounting for phylogenetic correlation.
        """
        n, dk = X.shape
        beta = params[:dk]
        log_inv_alpha = params[dk]
        alpha = 1.0 / np.exp(log_inv_alpha)

        # Linear predictor
        eta = X @ beta
        if np.any(np.abs(eta) >= btol):
            return 1e10

        mu = 1.0 / (1.0 + np.exp(-np.clip(eta, -btol, btol)))
        mu = np.clip(mu, 1e-16, 1 - 1e-16)

        # Build VCV
        try:
            vcv, diag_corr = self._build_logistic_vcv(
                tree, alpha, ordered_names, build_vcv_fn
            )
        except Exception:
            return 1e10

        # Standard Bernoulli log-likelihood
        ll = float(np.sum(y * np.log(mu) + (1 - y) * np.log(1 - mu)))
        if not np.isfinite(ll):
            return 1e10

        # Information matrix for Firth correction
        # I(beta) = X^T @ W^{1/2} @ V^{-1} @ W^{1/2} @ X
        info_matrix = self._compute_info_matrix(X, mu, vcv)

        # Firth penalty: 0.5 * log|I(beta)|
        try:
            sign, logdet = np.linalg.slogdet(info_matrix)
            if sign <= 0:
                return 1e10
        except np.linalg.LinAlgError:
            return 1e10

        pen_ll = ll + 0.5 * logdet
        return -pen_ll

    # ----------------------------------------------------------------
    # Fit
    # ----------------------------------------------------------------

    def _fit(
        self,
        tree,
        y: np.ndarray,
        X: np.ndarray,
        ordered_names: List[str],
        build_vcv_fn,
    ) -> Dict:
        """Run the optimisation and compute results."""
        n, p = X.shape
        k = p - 1  # number of predictors (excluding intercept)
        btol = 10

        # Starting values
        beta0 = self._logistic_starting_values(y, X, btol)
        eta0 = X @ beta0
        if np.any(np.abs(eta0) >= btol):
            beta0 = np.zeros(p)
            n1 = np.sum(y == 1)
            n0 = np.sum(y == 0)
            if n1 > 0 and n0 > 0:
                beta0[0] = np.log(n1 / n0)
                if np.any(np.abs(X @ beta0) >= btol):
                    beta0[0] = 0.0

        # Alpha bounds
        mean_rt = self._mean_root_tip_distance(tree, ordered_names)
        log_alpha_bound = 4.0
        lL_init = np.log(max(mean_rt, 1e-10))
        lL_lower = lL_init - log_alpha_bound
        lL_upper = lL_init + log_alpha_bound

        x0 = np.concatenate([beta0, [lL_init]])

        # Bounds: beta unbounded, log(1/alpha) bounded
        bounds = [(None, None)] * p + [(lL_lower, lL_upper)]

        result = minimize(
            self._neg_pen_loglik,
            x0,
            args=(X, y, tree, ordered_names, build_vcv_fn, btol),
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 500, "ftol": 2.2e-4, "gtol": 1e-5},
        )

        beta_hat = result.x[:p]
        alpha_hat = 1.0 / np.exp(result.x[p])
        convergence = result.status

        # Final quantities
        eta = X @ beta_hat
        mu = 1.0 / (1.0 + np.exp(-np.clip(eta, -btol, btol)))
        mu = np.clip(mu, 1e-16, 1 - 1e-16)

        # Unpenalized log-likelihood
        ll = float(np.sum(y * np.log(mu) + (1 - y) * np.log(1 - mu)))

        # Penalized log-likelihood and standard errors
        vcv, diag_corr = self._build_logistic_vcv(
            tree, alpha_hat, ordered_names, build_vcv_fn
        )
        info_matrix = self._compute_info_matrix(X, mu, vcv)

        try:
            sign, logdet = np.linalg.slogdet(info_matrix)
            if sign > 0:
                pen_ll = ll + 0.5 * logdet
            else:
                pen_ll = ll
        except np.linalg.LinAlgError:
            pen_ll = ll

        # Standard errors from inverse of info matrix
        try:
            I_inv = np.linalg.inv(info_matrix)
            se = np.sqrt(np.abs(np.diag(I_inv)))
        except np.linalg.LinAlgError:
            se = np.full(p, np.nan)

        z_stats = beta_hat / se
        p_values = 2.0 * norm_dist.sf(np.abs(z_stats))

        # AIC: -2*ll + 2*(p+1) where p coefficients + alpha
        aic = -2.0 * ll + 2.0 * (p + 1)

        coef_names = ["(Intercept)"] + list(self.predictors)
        formula = f"{self.response} ~ {' + '.join(self.predictors)}"

        return self._format_result(
            method=self.method,
            coef_names=coef_names,
            beta_hat=beta_hat,
            se=se,
            z_stats=z_stats,
            p_values=p_values,
            ll=ll,
            pen_ll=pen_ll,
            aic=aic,
            formula=formula,
            n=n,
            k=k,
            ordered_names=ordered_names,
            alpha=alpha_hat,
            convergence=convergence,
        )

    # ----------------------------------------------------------------
    # Output formatting
    # ----------------------------------------------------------------

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

    def _print_text_output(self, result: Dict) -> None:
        print("Phylogenetic Logistic Regression (Ives & Garland 2010)")
        print("=" * 54)
        print(f"Response: {self.response}")
        print(f"Predictor(s): {', '.join(self.predictors)}")
        print(f"Formula: {result['formula']}")
        print(f"Method: {result['method']}")
        print(f"Taxa: {result['n_taxa']}")
        print(f"Alpha: {result['alpha']:.6f}")

        print("\nCoefficients:")
        print(
            f"{'':20s}{'Estimate':>12s}{'Std.Error':>12s}"
            f"{'z-value':>12s}{'Pr(>|z|)':>12s}"
        )

        coefficients = result["coefficients"]
        for name in coefficients:
            coef = coefficients[name]
            sig = self._signif_code(coef["p_value"])
            print(
                f"{name:20s}{coef['estimate']:12.4f}{coef['std_error']:12.4f}"
                f"{coef['z_value']:12.4f}{coef['p_value']:12.6f}    {sig}"
            )

        print("---")
        print("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1")

        print(f"\nLog-likelihood: {result['log_likelihood']:.4f}")
        print(f"Penalized log-likelihood: {result['penalized_log_likelihood']:.4f}")
        print(f"AIC: {result['aic']:.4f}")

    def _format_result(
        self,
        *,
        method: str,
        coef_names: List[str],
        beta_hat: np.ndarray,
        se: np.ndarray,
        z_stats: np.ndarray,
        p_values: np.ndarray,
        ll: float,
        pen_ll: float,
        aic: float,
        formula: str,
        n: int,
        k: int,
        ordered_names: List[str],
        alpha: float,
        convergence: int,
    ) -> Dict:
        coefficients = {}
        for i, name in enumerate(coef_names):
            coefficients[name] = {
                "estimate": float(beta_hat[i]),
                "std_error": float(se[i]),
                "z_value": float(z_stats[i]),
                "p_value": float(p_values[i]),
            }

        result = {
            "method": method,
            "response": self.response,
            "predictors": list(self.predictors),
            "n_taxa": n,
            "alpha": float(alpha),
            "coefficients": coefficients,
            "log_likelihood": float(ll),
            "penalized_log_likelihood": float(pen_ll),
            "aic": float(aic),
            "convergence": int(convergence),
            "formula": formula,
        }

        return result
