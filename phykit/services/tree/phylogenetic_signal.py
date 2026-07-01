from __future__ import annotations

import math
import sys

from .base import Tree
from ...errors import PhykitUserError


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()
_CHO_FACTOR = None
_CHO_SOLVE = None
_MINIMIZE_SCALAR = None


def _permutation_p_value_ge(permutations: np.ndarray, observed: float) -> float:
    if permutations.size == 0:
        return float("nan")
    return float(np.count_nonzero(permutations >= observed) / permutations.size)


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def parse_multi_trait_file(*args, **kwargs):
    from ...helpers.trait_parsing import (
        parse_multi_trait_file as _parse_multi_trait_file,
    )

    return _parse_multi_trait_file(*args, **kwargs)


def trait_matrix_from_rows(*args, **kwargs):
    from ...helpers.trait_parsing import (
        trait_matrix_from_rows as _trait_matrix_from_rows,
    )

    return _trait_matrix_from_rows(*args, **kwargs)


def subset_traits_to_ordered_shared_taxa(*args, **kwargs):
    from ...helpers.trait_parsing import (
        subset_traits_to_ordered_shared_taxa as _subset_traits_to_ordered_shared_taxa,
    )

    return _subset_traits_to_ordered_shared_taxa(*args, **kwargs)


def _chi2_sf_df1(chi2_stat: float) -> float:
    return math.erfc(math.sqrt(chi2_stat / 2.0))


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


def _inverse_from_spd_matrix(matrix: np.ndarray) -> np.ndarray:
    try:
        factor = cho_factor(matrix, lower=True, check_finite=False)
        return cho_solve(
            factor,
            np.eye(matrix.shape[0]),
            check_finite=False,
        )
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return np.linalg.inv(matrix)


class PhylogeneticSignal(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.method = parsed["method"]
        self.permutations = parsed["permutations"]
        self.json_output = parsed["json_output"]
        self.gene_trees_path = parsed["gene_trees_path"]
        self.multivariate = parsed["multivariate"]

    def run(self) -> None:
        from .vcv_utils import build_vcv_matrix, build_discordance_vcv, parse_gene_trees

        tree = self.read_tree_file_unmodified()
        self.validate_tree(tree, min_tips=3, require_branch_lengths=True, context="phylogenetic signal analysis")

        tree_tips = self.get_tip_names_from_tree(tree)

        if self.multivariate:
            if self.method == "lambda":
                raise PhykitUserError(
                    [
                        "K_mult (--multivariate) only works with Blomberg's K framework.",
                        "Pagel's lambda is not supported for multivariate data.",
                        "Use --method blombergs_k (the default) instead.",
                    ],
                    code=2,
                )

            trait_names, traits_multi = parse_multi_trait_file(
                self.trait_data_path, tree_tips
            )
            ordered_names = sorted(traits_multi.keys())

            if self.gene_trees_path:
                gene_trees = parse_gene_trees(self.gene_trees_path)
                vcv, vcv_meta = build_discordance_vcv(tree, gene_trees, ordered_names)
                shared = vcv_meta["shared_taxa"]
                traits_multi, ordered_names = subset_traits_to_ordered_shared_taxa(
                    traits_multi, ordered_names, shared
                )
            else:
                vcv = build_vcv_matrix(tree, ordered_names)
                vcv_meta = None

            Y = trait_matrix_from_rows(traits_multi, ordered_names)

            result = self._kmult(Y, vcv, self.permutations)
            if vcv_meta is not None:
                result["vcv_metadata"] = vcv_meta
            if self.json_output:
                print_json(result)
                return
            print(
                f"{round(result['K_mult'], 4)}\t"
                f"{round(result['p_value'], 4)}\t"
                f"{result['n_traits']}\t"
                f"{result['permutations']}"
            )
            return

        traits = self._parse_trait_file(self.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())

        if self.gene_trees_path:
            gene_trees = parse_gene_trees(self.gene_trees_path)
            vcv, vcv_meta = build_discordance_vcv(tree, gene_trees, ordered_names)
            # Subset traits to shared taxa if needed
            shared = vcv_meta["shared_taxa"]
            traits, ordered_names = subset_traits_to_ordered_shared_taxa(
                traits, ordered_names, shared
            )
        else:
            vcv = build_vcv_matrix(tree, ordered_names)
            vcv_meta = None

        x = np.array([traits[name] for name in ordered_names])

        r2_phylo = self._compute_r2_phylo(x, vcv)

        if self.method == "blombergs_k":
            result = self._blombergs_k(x, vcv, self.permutations)
            result["r_squared_phylo"] = r2_phylo
            if vcv_meta is not None:
                result["vcv_metadata"] = vcv_meta
            if self.json_output:
                print_json(result)
                return
            print(
                f"{round(result['K'], 4)}\t"
                f"{round(result['p_value'], 4)}\t"
                f"{round(result['r_squared_phylo'], 4)}"
            )
        elif self.method == "lambda":
            from ...helpers.pgls_utils import max_lambda as compute_max_lambda

            max_lam = compute_max_lambda(tree)
            result = self._pagels_lambda(x, vcv, max_lam)
            result["r_squared_phylo"] = r2_phylo
            if vcv_meta is not None:
                result["vcv_metadata"] = vcv_meta
            if self.json_output:
                print_json(result)
                return
            print(
                f"{round(result['lambda'], 4)}\t"
                f"{round(result['log_likelihood'], 4)}\t"
                f"{round(result['p_value'], 4)}\t"
                f"{round(result['r_squared_phylo'], 4)}"
            )

    def process_args(self, args) -> dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            method=getattr(args, "method", "blombergs_k"),
            permutations=getattr(args, "permutations", 1000),
            json_output=getattr(args, "json", False),
            gene_trees_path=getattr(args, "gene_trees", None),
            multivariate=getattr(args, "multivariate", False),
        )

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

    def _kmult(
        self, Y: np.ndarray, vcv: np.ndarray, n_perm: int
    ) -> dict:
        """Multivariate K (Adams 2014).

        Y: n x p matrix of trait values (n taxa, p traits)
        vcv: n x n phylogenetic VCV matrix

        K_mult = observed_ratio / expected_ratio

        Where:
        - observed_ratio = MSE_obs / MSE_phylo
        - MSE_obs = sum of squared distances from each species to the
                    phylogenetic mean (using Euclidean distance in trait space)
        - MSE_phylo = sum of squared phylogenetically weighted distances
                      (using C_inv weighting)
        - expected_ratio = same formula as univariate K but uses trace
                          operations on the multivariate analogue
        """
        n, p = Y.shape
        ones = np.ones(n)

        C = vcv.copy()
        C_inv = _inverse_from_spd_matrix(C)
        sum_C_inv = float(ones @ C_inv @ ones)

        # GLS estimate of phylogenetic mean (p-dimensional vector)
        # a_hat = (1^T C^{-1} 1)^{-1} * (1^T C^{-1} Y)
        a_hat = (ones @ C_inv @ Y) / sum_C_inv  # shape: (p,)

        # Residual matrix: each row is (Y_i - a_hat)
        E = Y - a_hat  # n x p

        # MSE_obs: sum of squared Euclidean distances from each species to phylo mean
        # = trace(E^T @ E)
        MSE_obs = np.trace(E.T @ E)

        # MSE_phylo: phylogenetically-weighted version
        # = trace(E^T @ C_inv @ E)
        MSE_phylo = np.trace(E.T @ C_inv @ E)

        # Observed ratio
        if MSE_phylo == 0:
            observed_ratio = 0.0
        else:
            observed_ratio = MSE_obs / MSE_phylo

        # Expected ratio under BM (same as univariate, from Blomberg et al. 2003)
        expected_ratio = (np.trace(C) - n / sum_C_inv) / (n - 1)

        # K_mult
        if expected_ratio == 0:
            K_mult = 0.0
        else:
            K_mult = observed_ratio / expected_ratio

        # Permutation test: shuffle rows of Y, recompute K_mult
        k_perm = self._kmult_permutations(
            Y,
            C_inv,
            sum_C_inv,
            expected_ratio,
            n_perm,
        )

        p_value = _permutation_p_value_ge(k_perm, K_mult)

        return dict(
            K_mult=float(K_mult),
            p_value=p_value,
            permutations=n_perm,
            n_traits=p,
        )

    @staticmethod
    def _kmult_permutations(
        Y: np.ndarray,
        C_inv: np.ndarray,
        sum_C_inv: float,
        expected_ratio: float,
        n_perm: int,
        batch_size: int = 64,
    ) -> np.ndarray:
        rng = np.random.default_rng(seed=42)
        n, p = Y.shape
        k_perm = np.empty(n_perm, dtype=np.float64)
        if expected_ratio == 0:
            k_perm.fill(0.0)
            return k_perm

        weights = np.sum(C_inv, axis=0)

        for start in range(0, n_perm, batch_size):
            stop = min(start + batch_size, n_perm)
            batch_len = stop - start
            perms = np.empty((batch_len, n, p), dtype=np.float64)
            for row in range(batch_len):
                # Permute rows (species) of Y, keeping trait correlations intact.
                perms[row] = Y[rng.permutation(n)]

            a_perm = (
                np.einsum("n,bnp->bp", weights, perms, optimize=True)
                / sum_C_inv
            )
            E_perm = perms - a_perm[:, None, :]
            mse_obs = np.einsum("bnp,bnp->b", E_perm, E_perm, optimize=True)
            C_inv_E = np.einsum("ij,bjp->bip", C_inv, E_perm, optimize=True)
            mse_phylo = np.einsum(
                "bip,bip->b", E_perm, C_inv_E, optimize=True
            )
            obs_ratio = np.divide(
                mse_obs,
                mse_phylo,
                out=np.zeros_like(mse_obs),
                where=mse_phylo != 0,
            )
            k_perm[start:stop] = obs_ratio / expected_ratio

        return k_perm

    def _build_vcv_matrix(
        self, tree, ordered_names: list[str]
    ) -> np.ndarray:
        from .vcv_utils import build_vcv_matrix
        return build_vcv_matrix(tree, ordered_names)

    def _log_likelihood(
        self, x: np.ndarray, C: np.ndarray
    ) -> tuple[float, float]:
        """Concentrated log-likelihood with sigma^2 profiled out.

        Matches phytools::phylosig(method="lambda") internals.
        """
        try:
            n = len(x)
            ones = np.ones(n)
            factor = cho_factor(C, lower=True, check_finite=False)

            solve_rhs = np.empty((n, 2), dtype=np.result_type(x, C))
            solve_rhs[:, 0] = 1.0
            solve_rhs[:, 1] = x
            solved = cho_solve(factor, solve_rhs, check_finite=False)
            C_inv_ones = solved[:, 0]
            C_inv_x = solved[:, 1]

            # GLS estimate of phylogenetic mean
            a_hat = float((ones @ C_inv_x) / (ones @ C_inv_ones))

            # Residuals
            e = x - a_hat
            C_inv_e = C_inv_x - a_hat * C_inv_ones

            # MLE of sigma^2 (profiled out)
            sig2 = float(e @ C_inv_e / n)
            if sig2 <= 0:
                return self._log_likelihood_inverse(x, C)

            # Concentrated log-likelihood
            logdet_C = 2.0 * float(np.log(np.diag(factor[0])).sum())
            logdet_sig2C = n * np.log(sig2) + logdet_C
            ll = -0.5 * (n * np.log(2 * np.pi) + logdet_sig2C + n)
            return ll, a_hat
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._log_likelihood_inverse(x, C)

    def _log_likelihood_inverse(
        self, x: np.ndarray, C: np.ndarray
    ) -> tuple[float, float]:
        """Inverse-based likelihood implementation retained as a fallback."""
        n = len(x)
        ones = np.ones(n)

        C_inv = np.linalg.inv(C)

        # GLS estimate of phylogenetic mean
        a_hat = float((ones @ C_inv @ x) / (ones @ C_inv @ ones))

        # Residuals
        e = x - a_hat

        # MLE of sigma^2 (profiled out)
        sig2 = float(e @ C_inv @ e / n)

        # Concentrated log-likelihood
        sign, logdet_C = np.linalg.slogdet(C)
        # logdet(sig2 * C) = n * log(sig2) + logdet(C)
        logdet_sig2C = n * np.log(sig2) + logdet_C
        ll = -0.5 * (n * np.log(2 * np.pi) + logdet_sig2C + n)

        return ll, a_hat

    def _compute_r2_phylo(self, x: np.ndarray, vcv: np.ndarray) -> float:
        """Compute R²_phylo = 1 - (σ²_BM / σ²_WN).

        σ²_BM is the MLE of trait variance under Brownian Motion.
        σ²_WN is the MLE of trait variance under White Noise (= sample variance).
        """
        try:
            n = len(x)
            factor = cho_factor(vcv, lower=True, check_finite=False)
            ones = np.ones(n)
            solve_rhs = np.empty((n, 2), dtype=np.result_type(x, vcv))
            solve_rhs[:, 0] = 1.0
            solve_rhs[:, 1] = x
            solved = cho_solve(factor, solve_rhs, check_finite=False)
            C_inv_ones = solved[:, 0]
            C_inv_x = solved[:, 1]
            a_hat = float((ones @ C_inv_x) / (ones @ C_inv_ones))
            e = x - a_hat
            C_inv_e = C_inv_x - C_inv_ones * a_hat
            sig2_bm = float(e @ C_inv_e) / n
            sig2_wn = float(x.var())  # numpy uses /n by default
            if sig2_wn == 0:
                return float("nan")
            return 1.0 - sig2_bm / sig2_wn
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._compute_r2_phylo_inverse(x, vcv)

    def _compute_r2_phylo_inverse(self, x: np.ndarray, vcv: np.ndarray) -> float:
        """Inverse-based R²_phylo implementation retained as a fallback."""
        n = len(x)
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        a_hat = float((ones @ C_inv @ x) / (ones @ C_inv @ ones))
        e = x - a_hat
        sig2_bm = float(e @ C_inv @ e) / n
        sig2_wn = float(x.var())  # numpy uses /n by default
        if sig2_wn == 0:
            return float("nan")
        return 1.0 - sig2_bm / sig2_wn

    def _blombergs_k(
        self, x: np.ndarray, vcv: np.ndarray, n_perm: int
    ) -> dict:
        """Blomberg's K following phytools::phylosig(method="K").

        K = (e'e / e'C_inv'e) / ((trace(C) - n/sum(C_inv)) / (n-1))
        Uses raw sums (not MSE/MSE0 with sample-size normalization).
        """
        n = len(x)
        ones = np.ones(n)
        C = vcv.copy()

        C_inv = _inverse_from_spd_matrix(C)
        sum_C_inv = float(ones @ C_inv @ ones)

        # GLS estimate of phylogenetic mean
        a_hat = float((ones @ C_inv @ x) / sum_C_inv)

        # Residuals
        e = x - a_hat

        # Observed ratio (raw sums, matching phytools)
        observed_ratio = float(e @ e) / float(e @ C_inv @ e)

        # Expected ratio under BM (Blomberg et al. 2003)
        expected_ratio = (np.trace(C) - n / sum_C_inv) / (n - 1)

        # K statistic
        K = observed_ratio / expected_ratio

        k_perm = self._blombergs_k_permutations(
            x,
            C_inv,
            sum_C_inv,
            expected_ratio,
            n_perm,
        )

        p_value = _permutation_p_value_ge(k_perm, K)

        return dict(
            K=float(K),
            p_value=p_value,
            permutations=n_perm,
        )

    def _blombergs_k_permutations(
        self,
        x: np.ndarray,
        C_inv: np.ndarray,
        sum_C_inv: float,
        expected_ratio: float,
        n_perm: int,
        batch_size: int = 128,
    ) -> np.ndarray:
        rng = np.random.default_rng(seed=42)
        n = len(x)
        weights = np.sum(C_inv, axis=0)
        k_perm = np.empty(n_perm, dtype=np.float64)

        for start in range(0, n_perm, batch_size):
            stop = min(start + batch_size, n_perm)
            batch_len = stop - start
            perms = np.empty((batch_len, n), dtype=np.float64)
            for row in range(batch_len):
                perms[row] = rng.permutation(x)

            a_perm = (perms @ weights) / sum_C_inv
            E_perm = perms - a_perm[:, None]
            numerator = np.einsum("ij,ij->i", E_perm, E_perm)
            C_inv_E = E_perm @ C_inv
            denominator = np.einsum("ij,ij->i", E_perm, C_inv_E)
            k_perm[start:stop] = (numerator / denominator) / expected_ratio

        return k_perm

    def _pagels_lambda(
        self, x: np.ndarray, vcv: np.ndarray, max_lambda: float = 1.0
    ) -> dict:
        n = len(x)
        diag_vals = vcv.diagonal().copy()
        diag_step = vcv.shape[0] + 1
        niter = 10

        def neg_ll(lam):
            C_lam = vcv * lam
            C_lam.ravel()[::diag_step] = diag_vals
            try:
                ll, _ = self._log_likelihood(x, C_lam)
                return -ll
            except (np.linalg.LinAlgError, FloatingPointError, ValueError):
                return 1e10

        # Multi-interval optimization matching phytools (niter=10 sub-intervals)
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

        # Log-likelihood at fitted lambda
        C_fitted = vcv * lambda_hat
        C_fitted.ravel()[::diag_step] = diag_vals
        ll_fitted, _ = self._log_likelihood(x, C_fitted)

        # Log-likelihood at lambda = 0 (no phylogenetic signal)
        C_zero = vcv * 0.0
        C_zero.ravel()[::diag_step] = diag_vals
        ll_zero, _ = self._log_likelihood(x, C_zero)

        # Likelihood ratio test (1 df)
        lr_stat = 2.0 * (ll_fitted - ll_zero)
        lr_stat = max(lr_stat, 0.0)
        p_value = _chi2_sf_df1(lr_stat)

        return dict(
            **{"lambda": float(lambda_hat)},
            log_likelihood=float(ll_fitted),
            p_value=p_value,
        )
