"""
Shared utilities for discrete trait evolution models.

Provides Q-matrix construction, Felsenstein pruning (likelihood computation),
maximum-likelihood Q-matrix fitting, and discrete trait data parsing. Used by
stochastic_character_map, ancestral_reconstruction, and fit_discrete.
"""
from __future__ import annotations

import math
import sys

from ..errors import PhykitUserError


VALID_DISCRETE_MODELS = frozenset(["ER", "SYM", "ARD"])
_EXPM = None
_MINIMIZE = None
_MINIMIZE_SCALAR = None


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


def expm(*args, **kwargs):
    global _EXPM

    if _EXPM is None:
        from scipy.linalg import expm as _expm

        _EXPM = _expm

    return _EXPM(*args, **kwargs)


def minimize(*args, **kwargs):
    global _MINIMIZE

    if _MINIMIZE is None:
        from scipy.optimize import minimize as _minimize

        _MINIMIZE = _minimize

    return _MINIMIZE(*args, **kwargs)


def minimize_scalar(*args, **kwargs):
    global _MINIMIZE_SCALAR

    if _MINIMIZE_SCALAR is None:
        from scipy.optimize import minimize_scalar as _minimize_scalar

        _MINIMIZE_SCALAR = _minimize_scalar

    return _MINIMIZE_SCALAR(*args, **kwargs)


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


_SYM_INDEX_CACHE = {}
_ARD_MASK_CACHE = {}


def _sym_offdiag_indices(k: int):
    indices = _SYM_INDEX_CACHE.get(k)
    if indices is None:
        indices = np.triu_indices(k, 1)
        _SYM_INDEX_CACHE[k] = indices
    return indices


def _ard_offdiag_mask(k: int):
    mask = _ARD_MASK_CACHE.get(k)
    if mask is None:
        mask = ~np.eye(k, dtype=bool)
        _ARD_MASK_CACHE[k] = mask
    return mask


def build_q_matrix(params: np.ndarray, k: int, model: str) -> np.ndarray:
    """Build a Q-matrix from parameters for ER, SYM, or ARD models.

    Rows sum to zero (standard continuous-time Markov chain convention).
    """
    Q = np.zeros((k, k))
    if model == "ER":
        rate = params[0]
        Q[:] = rate
        np.fill_diagonal(Q, -rate * (k - 1))
    elif model == "SYM":
        if k <= 3:
            idx = 0
            for i in range(k):
                for j in range(i + 1, k):
                    Q[i, j] = params[idx]
                    Q[j, i] = params[idx]
                    idx += 1
        else:
            rows, cols = _sym_offdiag_indices(k)
            Q[rows, cols] = params
            Q[cols, rows] = params
    elif model == "ARD":
        Q[_ard_offdiag_mask(k)] = params
    if model != "ER":
        np.fill_diagonal(Q, -Q.sum(axis=1))
    return Q


def matrix_exp(Q: np.ndarray, t: float) -> np.ndarray:
    """Compute the matrix exponential P = exp(Q * t)."""
    two_state_transition = _matrix_exp_two_state(Q, t)
    if two_state_transition is not None:
        return two_state_transition
    return expm(Q * t)


def _matrix_exp_two_state(Q: np.ndarray, t: float):
    rates = _two_state_ctmc_rates(Q)
    if rates is None:
        return None

    rate01, rate10 = rates
    total_rate = rate01 + rate10
    if total_rate == 0.0:
        return np.eye(2, dtype=float)

    decay = math.exp(-total_rate * t)
    stay0 = (rate10 / total_rate) + (rate01 / total_rate) * decay
    to1 = (rate01 / total_rate) * (1.0 - decay)
    to0 = (rate10 / total_rate) * (1.0 - decay)
    stay1 = (rate01 / total_rate) + (rate10 / total_rate) * decay
    return np.array([[stay0, to1], [to0, stay1]], dtype=float)


def felsenstein_pruning(
    tree, tip_states: dict[str, str], Q: np.ndarray,
    pi: np.ndarray, states: list[str]
) -> tuple[dict, float]:
    """Postorder traversal computing conditional likelihoods and log-likelihood.

    Returns (cond_liks, loglik) where cond_liks maps clade id to
    a k-length likelihood vector.
    """
    k = len(states)
    state_idx = {s: i for i, s in enumerate(states)}
    cond_liks = {}
    transition_cache = {}

    postorder = _postorder_clades_direct(tree)
    if postorder is None:
        postorder = tree.find_clades(order="postorder")

    for clade in postorder:
        if clade.is_terminal():
            lik = np.zeros(k)
            if clade.name in tip_states:
                lik[state_idx[tip_states[clade.name]]] = 1.0
            cond_liks[id(clade)] = lik
        else:
            lik = np.ones(k)
            for child in clade.clades:
                t = child.branch_length if child.branch_length else 1e-8
                P = transition_cache.get(t)
                if P is None:
                    P = matrix_exp(Q, t)
                    transition_cache[t] = P
                child_lik = cond_liks[id(child)]
                lik *= P @ child_lik
            cond_liks[id(clade)] = lik

    root_lik = cond_liks[id(tree.root)]
    total_lik = float(np.dot(pi, root_lik))
    if total_lik <= 0:
        loglik = -1e20
    else:
        loglik = math.log(total_lik)

    return cond_liks, loglik


def _prepare_felsenstein_context(tree, tip_states: dict[str, str], states: list[str]):
    """Precompute tree/state metadata reused by repeated pruning evaluations."""
    state_idx = {s: i for i, s in enumerate(states)}
    postorder = _postorder_clades_direct(tree)
    if postorder is None:
        postorder = list(tree.find_clades(order="postorder"))
    node_index = {id(clade): i for i, clade in enumerate(postorder)}

    child_indices = []
    branch_lengths = []
    tip_state_indices = np.full(len(postorder), -1, dtype=np.intp)
    internal_entries = []

    for idx, clade in enumerate(postorder):
        if clade.is_terminal():
            if clade.name in tip_states:
                tip_state_indices[idx] = state_idx[tip_states[clade.name]]
            child_indices.append(())
            branch_lengths.append(())
        else:
            child_indices.append(
                tuple(node_index[id(child)] for child in clade.clades)
            )
            branch_lengths.append(
                tuple(
                    child.branch_length if child.branch_length else 1e-8
                    for child in clade.clades
                )
            )
            internal_entries.append((idx, child_indices[-1], branch_lengths[-1]))

    tip_liks = np.zeros((len(postorder), len(states)), dtype=float)
    observed_tip_rows = np.nonzero(tip_state_indices >= 0)[0]
    tip_liks[
        observed_tip_rows,
        tip_state_indices[observed_tip_rows],
    ] = 1.0

    return {
        "child_indices": child_indices,
        "branch_lengths": branch_lengths,
        "tip_state_indices": tip_state_indices,
        "tip_liks": tip_liks,
        "internal_entries": internal_entries,
        "root_index": node_index[id(tree.root)],
        "n_states": len(states),
    }


def _postorder_clades_direct(tree):
    try:
        root = tree.root
        root.clades
    except AttributeError:
        return None

    clades = []
    stack = [root]
    try:
        while stack:
            clade = stack.pop()
            clades.append(clade)
            stack.extend(clade.clades)
    except AttributeError:
        return None
    clades.reverse()
    return clades


def _felsenstein_loglik_prepared(context, Q: np.ndarray, pi: np.ndarray) -> float:
    k = context["n_states"]
    tip_liks = context.get("tip_liks")
    internal_entries = context.get("internal_entries")
    if tip_liks is None or internal_entries is None:
        tip_state_indices = context["tip_state_indices"]
        child_indices = context["child_indices"]
        branch_lengths = context["branch_lengths"]
        cond_liks = np.zeros((len(tip_state_indices), k), dtype=float)
        internal_entries = []

        for idx in range(len(tip_state_indices)):
            state_idx = tip_state_indices[idx]
            if state_idx >= 0:
                cond_liks[idx, state_idx] = 1.0
                continue

            children = child_indices[idx]
            if children:
                internal_entries.append((idx, children, branch_lengths[idx]))
    else:
        cond_liks = tip_liks.copy()
    transition_cache = {}

    two_state_rates = _two_state_ctmc_rates(Q) if k == 2 else None
    if two_state_rates is not None:
        rate01, rate10 = two_state_rates
        total_rate = rate01 + rate10
        for idx, children, lengths in internal_entries:
            lik0 = 1.0
            lik1 = 1.0
            for child_idx, t in zip(children, lengths):
                transition = transition_cache.get(t)
                if transition is None:
                    if total_rate == 0.0:
                        transition = (1.0, 0.0, 0.0, 1.0)
                    else:
                        decay = math.exp(-total_rate * t)
                        transition = (
                            (rate10 / total_rate)
                            + (rate01 / total_rate) * decay,
                            (rate01 / total_rate) * (1.0 - decay),
                            (rate10 / total_rate) * (1.0 - decay),
                            (rate01 / total_rate)
                            + (rate10 / total_rate) * decay,
                        )
                    transition_cache[t] = transition
                p00, p01, p10, p11 = transition
                child_lik = cond_liks[child_idx]
                child0 = child_lik[0]
                child1 = child_lik[1]
                lik0 *= (p00 * child0) + (p01 * child1)
                lik1 *= (p10 * child0) + (p11 * child1)
            cond_liks[idx, 0] = lik0
            cond_liks[idx, 1] = lik1

        total_lik = (pi[0] * cond_liks[context["root_index"], 0]) + (
            pi[1] * cond_liks[context["root_index"], 1]
        )
        if total_lik <= 0:
            return -1e20
        return math.log(total_lik)

    for idx, children, lengths in internal_entries:
        lik = np.ones(k)
        for child_idx, t in zip(children, lengths):
            P = transition_cache.get(t)
            if P is None:
                P = matrix_exp(Q, t)
                transition_cache[t] = P
            lik *= P @ cond_liks[child_idx]
        cond_liks[idx] = lik

    total_lik = float(np.dot(pi, cond_liks[context["root_index"]]))
    if total_lik <= 0:
        return -1e20
    return math.log(total_lik)


def _felsenstein_loglik_two_state_er_rate(
    context, rate: float, pi: np.ndarray
) -> float:
    tip_liks = context.get("tip_liks")
    internal_entries = context.get("internal_entries")
    if tip_liks is None or internal_entries is None:
        Q = build_q_matrix(np.array([rate]), 2, "ER")
        return _felsenstein_loglik_prepared(context, Q, pi)

    cond_liks = tip_liks.copy()
    transition_cache = {}
    for idx, children, lengths in internal_entries:
        lik0 = 1.0
        lik1 = 1.0
        for child_idx, t in zip(children, lengths):
            transition = transition_cache.get(t)
            if transition is None:
                decay = math.exp(-2.0 * rate * t)
                same = 0.5 + (0.5 * decay)
                different = 0.5 - (0.5 * decay)
                transition = (same, different)
                transition_cache[t] = transition
            else:
                same, different = transition
            child_lik = cond_liks[child_idx]
            child0 = child_lik[0]
            child1 = child_lik[1]
            lik0 *= (same * child0) + (different * child1)
            lik1 *= (different * child0) + (same * child1)
        cond_liks[idx, 0] = lik0
        cond_liks[idx, 1] = lik1

    root_lik = cond_liks[context["root_index"]]
    total_lik = (pi[0] * root_lik[0]) + (pi[1] * root_lik[1])
    if total_lik <= 0:
        return -1e20
    return math.log(total_lik)


def _two_state_ctmc_rates(Q):
    try:
        if Q.shape[0] != 2 or Q.shape[1] != 2:
            return None
    except (AttributeError, IndexError):
        return None

    rate01 = float(Q[0, 1])
    rate10 = float(Q[1, 0])
    if rate01 < 0.0 or rate10 < 0.0:
        return None
    if Q[0, 0] != -rate01 or Q[1, 1] != -rate10:
        return None
    return rate01, rate10


def _fit_two_state_er_rate(pruning_context, pi: np.ndarray):
    def neg_loglik_rate(rate):
        rate = float(abs(rate))
        ll = _felsenstein_loglik_two_state_er_rate(
            pruning_context,
            rate,
            pi,
        )
        return -ll

    result = minimize_scalar(
        neg_loglik_rate,
        bounds=(1e-8, 100.0),
        method="bounded",
        options={"xatol": 1e-10, "maxiter": 500},
    )
    rate = abs(float(result.x))
    Q = build_q_matrix(np.array([rate]), 2, "ER")
    loglik = _felsenstein_loglik_two_state_er_rate(pruning_context, rate, pi)
    return Q, loglik


def fit_q_matrix(
    tree, tip_states: dict[str, str],
    states: list[str], model: str, pruning_context=None
) -> tuple[np.ndarray, float]:
    """Fit Q-matrix parameters via maximum likelihood.

    Uses multi-start optimization with L-BFGS-B and Nelder-Mead,
    followed by a refinement step.

    Returns (Q_matrix, log_likelihood).
    """
    k = len(states)
    n_params = count_params(k, model)
    pi = np.ones(k) / k
    if pruning_context is None:
        pruning_context = _prepare_felsenstein_context(tree, tip_states, states)

    if model == "ER" and k == 2:
        return _fit_two_state_er_rate(pruning_context, pi)

    def neg_loglik(params):
        Q = build_q_matrix(np.abs(params), k, model)
        ll = _felsenstein_loglik_prepared(pruning_context, Q, pi)
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
    loglik = _felsenstein_loglik_prepared(pruning_context, Q, pi)

    return Q, loglik


def parse_discrete_traits(
    path: str, tree_tips: list[str], trait_column: str = None
) -> dict[str, str]:
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
            if trait_column is None:
                return _parse_two_column_stream(f, tree_tips)
            return _parse_multi_column_stream(f, path, tree_tips, trait_column)
    except FileNotFoundError:
        raise PhykitUserError(
            [
                f"{path} corresponds to no such file or directory.",
                "Please check filename and pathing",
            ],
            code=2,
        )


def _parse_multi_column(
    data_lines: list[str], path: str, tree_tips: list[str], trait_column: str
) -> dict[str, str]:
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
    target_idx = 1 + col_idx
    max_splits = target_idx + 1

    traits = {}
    for line_idx, line in enumerate(data_lines[1:], 2):
        column_count = line.count("\t") + 1
        if column_count != n_cols:
            raise PhykitUserError(
                [f"Line {line_idx} has {column_count} columns; expected {n_cols}."],
                code=2,
            )
        parts = line.split("\t", max_splits)
        taxon = parts[0]
        value = parts[target_idx].strip()
        if not value:
            raise PhykitUserError(
                [f"Missing trait value for taxon '{taxon}' on line {line_idx}."],
                code=2,
            )
        traits[taxon] = value

    return _validate_shared_taxa(traits, tree_tips)


def _parse_two_column(
    data_lines: list[str], path: str, tree_tips: list[str]
) -> dict[str, str]:
    traits = {}
    for line_idx, line in enumerate(data_lines, 1):
        taxon, sep, state = line.partition("\t")
        if not sep or "\t" in state:
            column_count = line.count("\t") + 1
            raise PhykitUserError(
                [f"Line {line_idx} has {column_count} columns; expected 2 (taxon, state)."],
                code=2,
            )
        traits[taxon] = state.strip()

    return _validate_shared_taxa(traits, tree_tips)


def _parse_two_column_stream(lines, tree_tips: list[str]) -> dict[str, str]:
    traits = {}
    line_idx = 1
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        taxon, sep, state = line.partition("\t")
        if not sep or "\t" in state:
            column_count = line.count("\t") + 1
            raise PhykitUserError(
                [f"Line {line_idx} has {column_count} columns; expected 2 (taxon, state)."],
                code=2,
            )
        traits[taxon] = state.strip()
        line_idx += 1

    return _validate_shared_taxa(traits, tree_tips)


def _parse_multi_column_stream(
    lines, path: str, tree_tips: list[str], trait_column: str
) -> dict[str, str]:
    header_parts = None
    col_idx = None
    n_cols = 0
    traits = {}
    line_idx = 0

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        line_idx += 1

        if header_parts is None:
            header_parts = line.split("\t")
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
            target_idx = 1 + col_idx
            max_splits = target_idx + 1
            continue

        column_count = line.count("\t") + 1
        if column_count != n_cols:
            raise PhykitUserError(
                [f"Line {line_idx} has {column_count} columns; expected {n_cols}."],
                code=2,
            )
        parts = line.split("\t", max_splits)
        taxon = parts[0]
        value = parts[target_idx].strip()
        if not value:
            raise PhykitUserError(
                [f"Missing trait value for taxon '{taxon}' on line {line_idx}."],
                code=2,
            )
        traits[taxon] = value

    if header_parts is None or line_idx < 2:
        raise PhykitUserError(
            ["Trait file must have a header row and at least one data row."],
            code=2,
        )

    return _validate_shared_taxa(traits, tree_tips)


def _validate_shared_taxa(
    traits: dict[str, str], tree_tips: list[str]
) -> dict[str, str]:
    tree_tip_set = set(tree_tips)
    if (
        len(tree_tip_set) >= 3
        and len(tree_tip_set) == len(traits)
        and tree_tip_set == traits.keys()
    ):
        return traits

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
