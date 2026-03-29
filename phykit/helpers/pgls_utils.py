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
from typing import Tuple

import numpy as np
from scipy.optimize import minimize_scalar

from ..errors import PhykitUserError


def max_lambda(tree) -> float:
    """Compute the upper bound for Pagel's lambda.

    For ultrametric trees, returns max_tip_height / max_parent_height.
    For non-ultrametric trees, returns 1.0.
    """
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
        C_inv = np.linalg.inv(C)
        XtCiX = X.T @ C_inv @ X
        XtCiX_inv = np.linalg.inv(XtCiX)
    except np.linalg.LinAlgError:
        return -1e20

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


def estimate_lambda(
    y: np.ndarray,
    X: np.ndarray,
    vcv: np.ndarray,
    max_lam: float = 1.0,
) -> Tuple[float, float]:
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
    diag_vals = np.diag(vcv).copy()
    niter = 10

    def neg_ll(lam):
        C_lam = vcv * lam
        np.fill_diagonal(C_lam, diag_vals)
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
    np.fill_diagonal(C_fitted, diag_vals)
    ll_fitted = pgls_log_likelihood(y, X, C_fitted)

    return float(lambda_hat), float(ll_fitted)


def fit_gls(
    y: np.ndarray, X: np.ndarray, C_inv: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, float, np.ndarray]:
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

    df_resid = n - k_plus_1
    sigma2 = float(residuals @ C_inv @ residuals) / max(df_resid, 1)

    var_beta = sigma2 * XtCiX_inv

    return beta_hat, residuals, sigma2, var_beta


def apply_lambda(vcv: np.ndarray, lambda_val: float) -> np.ndarray:
    """Apply Pagel's lambda to a VCV matrix (scale off-diagonals, keep diagonal)."""
    diag_vals = np.diag(vcv).copy()
    vcv_lam = vcv * lambda_val
    np.fill_diagonal(vcv_lam, diag_vals)
    return vcv_lam
