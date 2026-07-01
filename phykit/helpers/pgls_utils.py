"""
Shared PGLS (Phylogenetic Generalized Least Squares) utilities.

Provides reusable functions for:
- Pagel's lambda estimation via ML
- Concentrated PGLS log-likelihood
- GLS model fitting
- Lambda upper bound computation

Used by phylogenetic_regression, phylo_path, phylogenetic_signal,
phylogenetic_ordination, fit_continuous, and other comparative methods.
"""
from __future__ import annotations

from ..errors import PhykitUserError


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()
_CHO_FACTOR = None
_CHO_SOLVE = None
_MINIMIZE_SCALAR = None


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


def minimize_scalar(*args, **kwargs):
    global _MINIMIZE_SCALAR

    if _MINIMIZE_SCALAR is None:
        from scipy.optimize import minimize_scalar as _minimize_scalar

        _MINIMIZE_SCALAR = _minimize_scalar

    return _MINIMIZE_SCALAR(*args, **kwargs)


def _max_lambda_fallback(tree) -> float:
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


def max_lambda(tree) -> float:
    """Compute the upper bound for Pagel's lambda.

    For ultrametric trees, returns max_tip_height / max_parent_height.
    For non-ultrametric trees, returns 1.0.
    """
    direct_result = _max_lambda_direct(tree)
    if direct_result is not None:
        return direct_result

    try:
        tips = tree.get_terminals()
        root = tree.root
        depths = tree.depths()
        root_depth = depths[root]
        tip_heights = np.array(
            [depths[tip] - root_depth for tip in tips],
            dtype=np.float64,
        )
        max_tip_height = float(np.max(tip_heights))
        min_tip_height = float(np.min(tip_heights))

        is_ultrametric = (max_tip_height - min_tip_height) / max_tip_height < 1e-6

        if not is_ultrametric:
            return 1.0

        max_parent_height = 0.0
        for clade in tree.find_clades(order="level"):
            if clade == root:
                continue
            node_height = depths[clade] - root_depth
            parent_height = node_height - (clade.branch_length or 0.0)
            if parent_height > max_parent_height:
                max_parent_height = parent_height

        if max_parent_height == 0.0:
            return 1.0

        return max_tip_height / max_parent_height
    except (AttributeError, KeyError, TypeError, ValueError):
        return _max_lambda_fallback(tree)


def _max_lambda_direct(tree):
    try:
        root = tree.root
        root.clades
    except AttributeError:
        return None

    max_tip_height = -1.0
    min_tip_height = float("inf")
    max_parent_height = 0.0
    stack = [(root, 0.0)]
    try:
        pop = stack.pop
        append = stack.append
        while stack:
            clade, height = pop()
            children = clade.clades
            if children:
                if height > max_parent_height:
                    max_parent_height = height
                for child in children:
                    child_height = height + (child.branch_length or 0.0)
                    append((child, child_height))
            else:
                if height > max_tip_height:
                    max_tip_height = height
                if height < min_tip_height:
                    min_tip_height = height
    except AttributeError:
        return None

    is_ultrametric = (max_tip_height - min_tip_height) / max_tip_height < 1e-6
    if not is_ultrametric:
        return 1.0

    if max_parent_height == 0.0:
        return 1.0

    return max_tip_height / max_parent_height


def _pgls_log_likelihood_inverse(
    y: np.ndarray, X: np.ndarray, C: np.ndarray
) -> float:
    n = len(y)
    try:
        C_inv = np.linalg.inv(C)
        XtCiX = X.T @ C_inv @ X
        XtCiX_inv = np.linalg.inv(XtCiX)
    except np.linalg.LinAlgError:
        return -1e20

    beta_hat = XtCiX_inv @ (X.T @ (C_inv @ y))
    e = y - X @ beta_hat
    sigma2_ml = float(e @ C_inv @ e) / n

    sign, logdet_C = np.linalg.slogdet(C)
    if sign <= 0 or sigma2_ml <= 0:
        return -1e20

    ll = -0.5 * (
        n * np.log(2 * np.pi) + n * np.log(sigma2_ml) + logdet_C + n
    )
    return float(ll)


def pgls_log_likelihood(
    y: np.ndarray, X: np.ndarray, C: np.ndarray
) -> float:
    """Concentrated log-likelihood with beta and sigma^2 profiled out.

    Parameters
    ----------
    y : response vector (n,)
    X : design matrix (n, p)
    C : phylogenetic VCV matrix (n, n), possibly lambda-transformed
    """
    n = len(y)
    try:
        factor = cho_factor(C, lower=True, check_finite=False)
        rhs = np.empty((n, X.shape[1] + 1), dtype=np.result_type(y, X, C))
        rhs[:, :X.shape[1]] = X
        rhs[:, X.shape[1]] = y
        solved = cho_solve(factor, rhs, check_finite=False)
        C_inv_X = solved[:, :X.shape[1]]
        C_inv_y = solved[:, X.shape[1]]
        XtCiX = X.T @ C_inv_X
        beta_hat = np.linalg.solve(XtCiX, X.T @ C_inv_y)
        e = y - X @ beta_hat
        C_inv_e = C_inv_y - C_inv_X @ beta_hat
        sigma2_ml = float(e @ C_inv_e) / n

        diag = np.diag(factor[0])
        logdet_C = 2.0 * float(np.log(diag).sum())
        if sigma2_ml <= 0:
            return -1e20

        ll = -0.5 * (
            n * np.log(2 * np.pi) + n * np.log(sigma2_ml) + logdet_C + n
        )
        return float(ll)
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return _pgls_log_likelihood_inverse(y, X, C)


def estimate_lambda(
    y: np.ndarray,
    X: np.ndarray,
    vcv: np.ndarray,
    max_lam: float = 1.0,
) -> tuple[float, float]:
    """Optimize Pagel's lambda via ML using multi-interval bounded search.

    Parameters
    ----------
    y : response vector (n,)
    X : design matrix (n, p)
    vcv : phylogenetic VCV matrix (n, n)
    max_lam : upper bound for lambda (default 1.0)

    Returns
    -------
    (lambda_hat, log_likelihood_at_lambda)
    """
    diag_vals = vcv.diagonal().copy()
    diag_step = vcv.shape[0] + 1
    niter = 10

    def neg_ll(lam):
        C_lam = vcv * lam
        C_lam.ravel()[::diag_step] = diag_vals
        try:
            ll = pgls_log_likelihood(y, X, C_lam)
            return -ll
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return 1e10

    bounds_lo = np.linspace(0, max_lam - max_lam / niter, niter)
    bounds_hi = np.linspace(max_lam / niter, max_lam, niter)

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
    C_fitted.ravel()[::diag_step] = diag_vals
    ll_fitted = pgls_log_likelihood(y, X, C_fitted)

    return float(lambda_hat), float(ll_fitted)


def fit_gls(
    y: np.ndarray, X: np.ndarray, C_inv: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """GLS estimation: beta_hat = (X' C_inv X)^{-1} X' C_inv y.

    Parameters
    ----------
    y : response vector (n,)
    X : design matrix (n, p)
    C_inv : inverse of phylogenetic VCV matrix (n, n)

    Returns
    -------
    (beta_hat, residuals, sigma2_reml, var_beta)
    """
    n, k_plus_1 = X.shape
    rhs = np.empty((n, k_plus_1 + 1), dtype=np.result_type(X, y))
    rhs[:, :k_plus_1] = X
    rhs[:, k_plus_1] = y
    C_inv_rhs = C_inv @ rhs
    C_inv_X = C_inv_rhs[:, :k_plus_1]
    C_inv_y = C_inv_rhs[:, k_plus_1]
    XtCiX = X.T @ C_inv_X
    try:
        XtCiX_inv = np.linalg.solve(
            XtCiX,
            np.eye(k_plus_1, dtype=XtCiX.dtype),
        )
    except np.linalg.LinAlgError:
        raise PhykitUserError(
            [
                "Singular design matrix: cannot estimate coefficients.",
                "Check that predictors are not collinear.",
            ],
            code=2,
        )

    beta_hat = XtCiX_inv @ (X.T @ C_inv_y)
    residuals = y - X @ beta_hat
    C_inv_residuals = C_inv_y - C_inv_X @ beta_hat

    df_resid = n - k_plus_1
    sigma2 = float(residuals @ C_inv_residuals) / max(df_resid, 1)

    var_beta = sigma2 * XtCiX_inv

    return beta_hat, residuals, sigma2, var_beta


def apply_lambda(vcv: np.ndarray, lambda_val: float) -> np.ndarray:
    """Apply Pagel's lambda to a VCV matrix (scale off-diagonals, keep diagonal)."""
    diag_vals = np.diag(vcv).copy()
    vcv_lam = vcv * lambda_val
    np.fill_diagonal(vcv_lam, diag_vals)
    return vcv_lam
