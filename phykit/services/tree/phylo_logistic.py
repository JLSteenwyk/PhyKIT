from __future__ import annotations

from math import exp, sqrt

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
_SPECIAL_ERFC = None


def _binary_response_class_counts(y: np.ndarray) -> tuple[int, int]:
    n1 = int(np.count_nonzero(y))
    return n1, int(y.size - n1)


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


def special_erfc(*args, **kwargs):
    global _SPECIAL_ERFC

    if _SPECIAL_ERFC is None:
        from scipy.special import erfc as _erfc

        _SPECIAL_ERFC = _erfc

    return _SPECIAL_ERFC(*args, **kwargs)


def _normal_two_tailed_p_values(z_stats: np.ndarray) -> np.ndarray:
    z_values = np.asarray(z_stats, dtype=float)
    p_values = special_erfc(
        np.abs(z_values) / sqrt(2.0),
    )
    return np.asarray(p_values, dtype=np.float64).reshape(z_values.shape)


def _standard_errors_from_info_matrix(info_matrix: np.ndarray) -> np.ndarray:
    try:
        factor = cho_factor(info_matrix, lower=True, check_finite=False)
        info_inv_diag = np.diag(
            cho_solve(
                factor,
                np.eye(info_matrix.shape[0]),
                check_finite=False,
            )
        )
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        try:
            info_inv_diag = np.diag(np.linalg.inv(info_matrix))
        except np.linalg.LinAlgError:
            return np.full(info_matrix.shape[0], np.nan)

    return np.sqrt(np.abs(info_inv_diag))


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

        tree = self.read_tree_file_unmodified()
        self.validate_tree(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="phylogenetic logistic regression",
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

    def process_args(self, args) -> dict[str, str]:
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

    # ----------------------------------------------------------------
    # VCV construction with OU-transformed branch lengths
    # ----------------------------------------------------------------

    @staticmethod
    def _root_tip_distances(tree, ordered_names: list[str]) -> np.ndarray:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            root = None

        if root is not None:
            try:
                distance_by_name = {}
                stack = [(root, 0.0)]
                while stack:
                    clade, distance = stack.pop()
                    children = clade.clades
                    if children:
                        for child in reversed(children):
                            stack.append(
                                (
                                    child,
                                    distance + (child.branch_length or 0.0),
                                )
                            )
                    else:
                        distance_by_name[clade.name] = distance
                return np.array([distance_by_name[name] for name in ordered_names])
            except (AttributeError, KeyError, TypeError):
                pass

        try:
            depths = tree.depths()
            root_depth = depths[tree.root]
            terminals = tree.get_terminals()
            terminal_by_name = {terminal.name: terminal for terminal in terminals}
            return np.array([
                depths[terminal_by_name[name]] - root_depth for name in ordered_names
            ])
        except (AttributeError, KeyError, TypeError):
            return np.array([
                tree.distance(tree.root, name) for name in ordered_names
            ])

    def _build_logistic_vcv(
        self,
        tree,
        alpha: float,
        ordered_names: list[str],
        build_vcv_fn,
        root_tip_distances: np.ndarray | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
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

        def transform_branch_length(branch_length):
            if branch_length <= 0:
                return 0.0
            return (1 - exp(-2 * alpha * branch_length)) / (2 * alpha)

        from .vcv_utils import build_transformed_vcv_matrix, build_vcv_matrix

        if build_vcv_fn is build_vcv_matrix:
            vcv = build_transformed_vcv_matrix(
                tree, ordered_names, transform_branch_length
            )
        else:
            import pickle

            tree_t = pickle.loads(pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL))
            for clade in tree_t.find_clades():
                if clade.branch_length is not None and clade.branch_length > 0:
                    clade.branch_length = transform_branch_length(clade.branch_length)
            vcv = build_vcv_fn(tree_t, ordered_names)

        if root_tip_distances is None:
            root_tip_distances = self._root_tip_distances(tree, ordered_names)
        diag_corr = np.exp(-2 * alpha * root_tip_distances)

        # Add diagonal correction to VCV
        np.fill_diagonal(vcv, np.diag(vcv) + diag_corr)

        return vcv, diag_corr

    def _mean_root_tip_distance(self, tree, ordered_names: list[str]) -> float:
        """Mean root-to-tip distance across all tips."""
        return float(np.mean(self._root_tip_distances(tree, ordered_names)))

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
            weights = mu * (1 - mu)
            z = eta + (y - mu) / weights
            try:
                beta_new = np.linalg.solve(
                    X.T @ (weights[:, None] * X),
                    X.T @ (weights * z),
                )
            except np.linalg.LinAlgError:
                break
            if np.abs(beta_new - beta).max() < 1e-8:
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
        try:
            return self._compute_info_matrix_cholesky(X, mu, vcv)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._compute_info_matrix_inverse(X, mu, vcv)

    def _compute_info_matrix_cholesky(
        self, X: np.ndarray, mu: np.ndarray, vcv: np.ndarray
    ) -> np.ndarray:
        weights = np.sqrt(mu * (1 - mu))
        X_weighted = X * weights[:, None]
        factor = cho_factor(vcv, lower=True, check_finite=False)
        V_inv_X_weighted = cho_solve(factor, X_weighted, check_finite=False)
        return X_weighted.T @ V_inv_X_weighted

    def _compute_info_matrix_inverse(
        self, X: np.ndarray, mu: np.ndarray, vcv: np.ndarray
    ) -> np.ndarray:
        weights = np.sqrt(mu * (1 - mu))
        X_weighted = X * weights[:, None]
        try:
            V_inv = np.linalg.inv(vcv)
        except np.linalg.LinAlgError:
            return np.eye(X.shape[1]) * 1e-10

        return X_weighted.T @ V_inv @ X_weighted

    def _neg_pen_loglik(
        self,
        params: np.ndarray,
        X: np.ndarray,
        y: np.ndarray,
        tree,
        ordered_names: list[str],
        build_vcv_fn,
        root_tip_distances: np.ndarray,
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
        if (np.abs(eta) >= btol).any():
            return 1e10

        mu = 1.0 / (1.0 + np.exp(-np.clip(eta, -btol, btol)))
        mu = np.clip(mu, 1e-16, 1 - 1e-16)

        # Build VCV
        try:
            vcv, diag_corr = self._build_logistic_vcv(
                tree, alpha, ordered_names, build_vcv_fn, root_tip_distances
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
        ordered_names: list[str],
        build_vcv_fn,
    ) -> dict:
        """Run the optimisation and compute results."""
        n, p = X.shape
        k = p - 1  # number of predictors (excluding intercept)
        btol = 10

        # Starting values
        beta0 = self._logistic_starting_values(y, X, btol)
        eta0 = X @ beta0
        if (np.abs(eta0) >= btol).any():
            beta0 = np.zeros(p)
            n1, n0 = _binary_response_class_counts(y)
            if n1 > 0 and n0 > 0:
                beta0[0] = np.log(n1 / n0)
                if (np.abs(X @ beta0) >= btol).any():
                    beta0[0] = 0.0

        root_tip_distances = self._root_tip_distances(tree, ordered_names)
        mean_rt = float(np.mean(root_tip_distances))
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
            args=(X, y, tree, ordered_names, build_vcv_fn, root_tip_distances, btol),
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
            tree, alpha_hat, ordered_names, build_vcv_fn, root_tip_distances
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

        se = _standard_errors_from_info_matrix(info_matrix)

        z_stats = beta_hat / se
        p_values = _normal_two_tailed_p_values(z_stats)

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

    def _print_text_output(self, result: dict) -> None:
        lines = [
            "Phylogenetic Logistic Regression (Ives & Garland 2010)",
            "=" * 54,
            f"Response: {self.response}",
            f"Predictor(s): {', '.join(self.predictors)}",
            f"Formula: {result['formula']}",
            f"Method: {result['method']}",
            f"Taxa: {result['n_taxa']}",
            f"Alpha: {result['alpha']:.6f}",
            "\nCoefficients:",
            f"{'':20s}{'Estimate':>12s}{'Std.Error':>12s}"
            f"{'z-value':>12s}{'Pr(>|z|)':>12s}",
        ]

        append = lines.append
        for name, coef in result["coefficients"].items():
            p_value = coef["p_value"]
            if p_value < 0.001:
                sig = "***"
            elif p_value < 0.01:
                sig = "**"
            elif p_value < 0.05:
                sig = "*"
            elif p_value < 0.1:
                sig = "."
            else:
                sig = " "
            append(
                f"{name:20s}{coef['estimate']:12.4f}{coef['std_error']:12.4f}"
                f"{coef['z_value']:12.4f}{p_value:12.6f}    {sig}"
            )

        lines.append("---")
        lines.append("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1")

        lines.append(f"\nLog-likelihood: {result['log_likelihood']:.4f}")
        lines.append(
            f"Penalized log-likelihood: {result['penalized_log_likelihood']:.4f}"
        )
        lines.append(f"AIC: {result['aic']:.4f}")
        print("\n".join(lines))

    def _format_result(
        self,
        *,
        method: str,
        coef_names: list[str],
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
        ordered_names: list[str],
        alpha: float,
        convergence: int,
    ) -> dict:
        coefficients = {
            name: {
                "estimate": float(beta),
                "std_error": float(std_error),
                "z_value": float(z_stat),
                "p_value": float(p_value),
            }
            for name, beta, std_error, z_stat, p_value in zip(
                coef_names, beta_hat, se, z_stats, p_values
            )
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
