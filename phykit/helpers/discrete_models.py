"""
Shared utilities for discrete trait evolution models.

Provides Q-matrix construction, Felsenstein pruning (likelihood computation),
maximum-likelihood Q-matrix fitting, and discrete trait data parsing. Used by
stochastic_character_map, ancestral_reconstruction, and fit_discrete.
"""
import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

from ..errors import PhykitUserError


VALID_DISCRETE_MODELS = frozenset(["ER", "SYM", "ARD"])


def count_params(k: int, model: str) -> int:
    """Return the number of free rate parameters for a discrete model."""
    if model == "ER":
        return 1
    elif model == "SYM":
        return k * (k - 1) // 2
    elif model == "ARD":
        return k * (k - 1)
    else:
        raise PhykitUserError(
            [f"Unknown model '{model}'. Use ER, SYM, or ARD."], code=2,
        )


def build_q_matrix(params: np.ndarray, k: int, model: str) -> np.ndarray:
    """Build a Q-matrix from parameters for ER, SYM, or ARD models.

    Rows sum to zero (standard continuous-time Markov chain convention).
    """
    Q = np.zeros((k, k))
    if model == "ER":
        rate = params[0]
        Q[:] = rate
        np.fill_diagonal(Q, 0.0)
    elif model == "SYM":
        idx = 0
        for i in range(k):
            for j in range(i + 1, k):
                Q[i, j] = params[idx]
                Q[j, i] = params[idx]
                idx += 1
    elif model == "ARD":
        idx = 0
        for i in range(k):
            for j in range(k):
                if i != j:
                    Q[i, j] = params[idx]
                    idx += 1
    for i in range(k):
        Q[i, i] = -np.sum(Q[i, :])
    return Q


def matrix_exp(Q: np.ndarray, t: float) -> np.ndarray:
    """Compute the matrix exponential P = exp(Q * t)."""
    return expm(Q * t)


def felsenstein_pruning(
    tree, tip_states: Dict[str, str], Q: np.ndarray,
    pi: np.ndarray, states: List[str]
) -> Tuple[Dict, float]:
    """Postorder traversal computing conditional likelihoods and log-likelihood.

    Returns (cond_liks, loglik) where cond_liks maps clade id to
    a k-length likelihood vector.
    """
    k = len(states)
    state_idx = {s: i for i, s in enumerate(states)}
    cond_liks = {}

    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            lik = np.zeros(k)
            if clade.name in tip_states:
                lik[state_idx[tip_states[clade.name]]] = 1.0
            cond_liks[id(clade)] = lik
        else:
            lik = np.ones(k)
            for child in clade.clades:
                t = child.branch_length if child.branch_length else 1e-8
                P = matrix_exp(Q, t)
                child_lik = cond_liks[id(child)]
                lik *= P @ child_lik
            cond_liks[id(clade)] = lik

    root_lik = cond_liks[id(tree.root)]
    total_lik = np.sum(pi * root_lik)
    if total_lik <= 0:
        loglik = -1e20
    else:
        loglik = np.log(total_lik)

    return cond_liks, loglik


def fit_q_matrix(
    tree, tip_states: Dict[str, str],
    states: List[str], model: str
) -> Tuple[np.ndarray, float]:
    """Fit Q-matrix parameters via maximum likelihood.

    Uses multi-start optimization with L-BFGS-B and Nelder-Mead,
    followed by a refinement step.

    Returns (Q_matrix, log_likelihood).
    """
    k = len(states)
    n_params = count_params(k, model)
    pi = np.ones(k) / k

    def neg_loglik(params):
        Q = build_q_matrix(np.abs(params), k, model)
        _, ll = felsenstein_pruning(tree, tip_states, Q, pi, states)
        return -ll

    bounds = [(1e-8, 100.0)] * n_params

    starting_values = [0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0]
    best_negll = np.inf
    best_params = np.ones(n_params) * 0.1

    for sv in starting_values:
        x0 = np.ones(n_params) * sv
        for method in ["L-BFGS-B", "Nelder-Mead"]:
            try:
                kwargs = {"method": method}
                if method == "L-BFGS-B":
                    kwargs["bounds"] = bounds
                result = minimize(neg_loglik, x0, **kwargs)
                if result.fun < best_negll:
                    best_negll = result.fun
                    best_params = np.abs(result.x)
            except (ValueError, np.linalg.LinAlgError):
                continue

    # Refine best result
    try:
        result = minimize(
            neg_loglik, best_params, method="Nelder-Mead",
            options={"maxiter": 10000, "xatol": 1e-10, "fatol": 1e-10},
        )
        if result.fun < best_negll:
            best_params = np.abs(result.x)
    except (ValueError, np.linalg.LinAlgError):
        pass

    Q = build_q_matrix(best_params, k, model)
    _, loglik = felsenstein_pruning(tree, tip_states, Q, pi, states)

    return Q, loglik


def parse_discrete_traits(
    path: str, tree_tips: List[str], trait_column: str = None
) -> Dict[str, str]:
    """Parse discrete trait data from a TSV file.

    If trait_column is None: expects 2-column format (taxon<tab>state),
    no header row.
    If trait_column is given: expects multi-column with header row,
    extracts the named column.

    Returns {taxon: state} dict for the intersection of file taxa and
    tree_tips. Requires at least 3 shared taxa.
    """
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

    if trait_column is not None:
        return _parse_multi_column(data_lines, path, tree_tips, trait_column)
    else:
        return _parse_two_column(data_lines, path, tree_tips)


def _parse_multi_column(
    data_lines: List[str], path: str, tree_tips: List[str], trait_column: str
) -> Dict[str, str]:
    if len(data_lines) < 2:
        raise PhykitUserError(
            ["Trait file must have a header row and at least one data row."],
            code=2,
        )

    header_parts = data_lines[0].split("\t")
    if len(header_parts) < 2:
        raise PhykitUserError(
            ["Header must have at least 2 columns (taxon + at least 1 trait)."],
            code=2,
        )

    trait_names = header_parts[1:]
    if trait_column not in trait_names:
        raise PhykitUserError(
            [
                f"Column '{trait_column}' not found in trait file.",
                f"Available columns: {', '.join(trait_names)}",
            ],
            code=2,
        )

    col_idx = trait_names.index(trait_column)
    n_cols = len(header_parts)

    traits = {}
    for line_idx, line in enumerate(data_lines[1:], 2):
        parts = line.split("\t")
        if len(parts) != n_cols:
            raise PhykitUserError(
                [f"Line {line_idx} has {len(parts)} columns; expected {n_cols}."],
                code=2,
            )
        taxon = parts[0]
        value = parts[1 + col_idx].strip()
        if not value:
            raise PhykitUserError(
                [f"Missing trait value for taxon '{taxon}' on line {line_idx}."],
                code=2,
            )
        traits[taxon] = value

    return _validate_shared_taxa(traits, tree_tips)


def _parse_two_column(
    data_lines: List[str], path: str, tree_tips: List[str]
) -> Dict[str, str]:
    traits = {}
    for line_idx, line in enumerate(data_lines, 1):
        parts = line.split("\t")
        if len(parts) != 2:
            raise PhykitUserError(
                [f"Line {line_idx} has {len(parts)} columns; expected 2 (taxon, state)."],
                code=2,
            )
        traits[parts[0]] = parts[1].strip()

    return _validate_shared_taxa(traits, tree_tips)


def _validate_shared_taxa(
    traits: Dict[str, str], tree_tips: List[str]
) -> Dict[str, str]:
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
