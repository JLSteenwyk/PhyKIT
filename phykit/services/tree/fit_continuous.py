from __future__ import annotations

import math
import sys

from .base import Tree
from ...errors import PhykitUserError

ALL_MODELS = ["BM", "OU", "EB", "Lambda", "Delta", "Kappa", "White"]


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def subset_traits_to_ordered_shared_taxa(*args, **kwargs):
    from ...helpers.trait_parsing import (
        subset_traits_to_ordered_shared_taxa as _subset_traits_to_ordered_shared_taxa,
    )

    return _subset_traits_to_ordered_shared_taxa(*args, **kwargs)


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


class FitContinuous(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.selected_models = parsed["models"]
        self.json_output = parsed["json_output"]
        self.gene_trees_path = parsed["gene_trees_path"]

    def run(self) -> None:
        from .vcv_utils import build_vcv_matrix, build_discordance_vcv, parse_gene_trees

        tree = self.read_tree_file_unmodified()
        self.validate_tree(tree, min_tips=3, require_branch_lengths=True, context="model fitting")

        tree_tips = self.get_tip_names_from_tree(tree)
        traits = self._parse_trait_file(self.trait_data_path, tree_tips)

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

        x = np.array([traits[name] for name in ordered_names])
        n = len(x)

        # Precompute helpers needed by tree-transformation models
        parent_map = self._build_parent_map(tree)
        paths = self._build_root_to_tip_paths(tree, ordered_names, parent_map)
        if "Lambda" in self.selected_models:
            from ...helpers.pgls_utils import max_lambda as compute_max_lambda

            max_lam = compute_max_lambda(tree)
        else:
            max_lam = 1.0
        tree_height = float(np.max(np.diag(vcv)))

        results = []
        for model_name in self.selected_models:
            res = self._fit_model(
                model_name, x, vcv, tree, ordered_names,
                paths, max_lam, tree_height,
            )
            results.append(res)

        results = self._compute_model_comparison(results, n)

        # Always fit White for R² baseline
        white_sig2 = None
        for r in results:
            if r["model"] == "White":
                white_sig2 = r["sigma2"]
                break
        if white_sig2 is None:
            # White not in selected models — fit silently
            white_result = self._fit_white(x)
            white_sig2 = white_result["sigma2"]

        # Add R² to each model
        for r in results:
            if white_sig2 > 0:
                r["r_squared"] = 1.0 - r["sigma2"] / white_sig2
            else:
                r["r_squared"] = float("nan")

        if self.json_output:
            self._print_json_output(results, n, vcv_meta)
        else:
            self._print_text_output(results, n)

    def process_args(self, args) -> dict:
        models = ALL_MODELS[:]
        if hasattr(args, "models") and args.models:
            requested = [m.strip() for m in args.models.split(",")]
            valid = {m.lower(): m for m in ALL_MODELS}
            models = []
            for m in requested:
                canonical = valid.get(m.lower())
                if canonical is None:
                    raise PhykitUserError(
                        [
                            f"Unknown model '{m}'.",
                            f"Available models: {', '.join(ALL_MODELS)}",
                        ],
                        code=2,
                    )
                models.append(canonical)
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            models=models,
            json_output=getattr(args, "json", False),
            gene_trees_path=getattr(args, "gene_trees", None),
        )

    # ── Tree & trait parsing ─────────────────────────────────────────

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

    # ── VCV construction ─────────────────────────────────────────────

    def _build_vcv_matrix(
        self, tree, ordered_names: list[str]
    ) -> np.ndarray:
        from .vcv_utils import build_vcv_matrix
        return build_vcv_matrix(tree, ordered_names)

    def _build_parent_map(self, tree) -> dict:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            parent_map = {}
            for clade in tree.find_clades(order="preorder"):
                for child in clade.clades:
                    parent_map[id(child)] = clade
            return parent_map

        parent_map = {}
        stack = [root]
        try:
            pop = stack.pop
            extend = stack.extend
            while stack:
                clade = pop()
                children = clade.clades
                for child in children:
                    parent_map[id(child)] = clade
                if children:
                    extend(children)
        except AttributeError:
            parent_map = {}
            for clade in tree.find_clades(order="preorder"):
                for child in clade.clades:
                    parent_map[id(child)] = clade
        return parent_map

    def _build_root_to_tip_paths(
        self, tree, ordered_names: list[str], parent_map: dict
    ) -> dict[str, list]:
        ordered_name_set = set(ordered_names)
        try:
            root = tree.root
            root.clades
        except AttributeError:
            tip_map = {}
            for tip in tree.get_terminals():
                if tip.name in ordered_name_set:
                    tip_map[tip.name] = tip
        else:
            tip_map = {}
            stack = [root]
            while stack:
                clade = stack.pop()
                children = clade.clades
                if children:
                    stack.extend(reversed(children))
                elif clade.name in ordered_name_set:
                    tip_map[clade.name] = clade

        paths = {}
        for name in ordered_names:
            tip = tip_map[name]
            path = []
            current = tip
            while id(current) in parent_map:
                bl = current.branch_length if current.branch_length else 0.0
                path.append((id(current), bl))
                current = parent_map[id(current)]
            path.reverse()
            paths[name] = path

        return paths

    # ── Concentrated log-likelihood ──────────────────────────────────

    def _concentrated_ll(
        self, x: np.ndarray, C: np.ndarray
    ) -> tuple[float, float, float]:
        """Concentrated log-likelihood with sigma^2 profiled out.

        Returns (log_likelihood, sigma2, z0).
        """
        try:
            return self._concentrated_ll_cholesky(x, C)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._concentrated_ll_inverse(x, C)

    def _concentrated_ll_cholesky(
        self, x: np.ndarray, C: np.ndarray
    ) -> tuple[float, float, float]:
        n = len(x)
        ones = np.ones(n)
        factor = cho_factor(C, lower=True, check_finite=False)
        solve_rhs = np.empty((n, 2), dtype=np.result_type(x, C))
        solve_rhs[:, 0] = 1.0
        solve_rhs[:, 1] = x
        solved = cho_solve(factor, solve_rhs, check_finite=False)
        C_inv_ones = solved[:, 0]
        C_inv_x = solved[:, 1]

        denom = float(ones @ C_inv_ones)
        if abs(denom) < 1e-300:
            return float("-inf"), 0.0, 0.0
        z0 = float(ones @ C_inv_x) / denom

        e = x - z0
        C_inv_e = C_inv_x - z0 * C_inv_ones
        sig2 = float(e @ C_inv_e) / n

        if sig2 <= 0:
            return float("-inf"), 0.0, z0

        logdet_C = 2.0 * float(np.sum(np.log(np.diag(factor[0]))))
        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet_C + n)

        return float(ll), float(sig2), float(z0)

    def _concentrated_ll_inverse(
        self, x: np.ndarray, C: np.ndarray
    ) -> tuple[float, float, float]:
        n = len(x)
        ones = np.ones(n)

        try:
            C_inv = np.linalg.inv(C)
        except np.linalg.LinAlgError:
            return float("-inf"), 0.0, 0.0

        # GLS estimate of ancestral state
        denom = float(ones @ C_inv @ ones)
        if abs(denom) < 1e-300:
            return float("-inf"), 0.0, 0.0
        z0 = float(ones @ C_inv @ x) / denom

        # Residuals
        e = x - z0

        # MLE sigma^2
        sig2 = float(e @ C_inv @ e) / n

        if sig2 <= 0:
            return float("-inf"), 0.0, z0

        # Concentrated log-likelihood
        sign, logdet_C = np.linalg.slogdet(C)
        if sign <= 0:
            return float("-inf"), sig2, z0

        ll = -0.5 * (n * np.log(2 * np.pi) + n * np.log(sig2) + logdet_C + n)

        return float(ll), float(sig2), float(z0)

    # ── VCV transformations ──────────────────────────────────────────

    def _vcv_ou(self, C: np.ndarray, alpha: float) -> np.ndarray:
        """OU VCV using Martins & Hansen (1997) formula."""
        if alpha < 1e-10:
            return C.copy()

        diag_vals = np.diag(C)
        unique_path = diag_vals[:, None] + diag_vals[None, :] - 2.0 * C
        return (
            np.exp(-alpha * unique_path)
            * (1.0 - np.exp(-2.0 * alpha * C))
            / (2.0 * alpha)
        )

    def _vcv_lambda(self, C: np.ndarray, lam: float) -> np.ndarray:
        """Pagel's lambda VCV transformation."""
        diag_vals = np.diag(C).copy()
        C_lam = C * lam
        np.fill_diagonal(C_lam, diag_vals)
        return C_lam

    def _build_transformed_vcv(
        self, ordered_names: list[str], paths: dict[str, list],
        transform_fn,
    ) -> np.ndarray:
        """Build VCV from transformed branch lengths.

        transform_fn(bl, d_start, d_end) -> new branch length
        where d_start is distance from root to start of branch,
        d_end is distance from root to end of branch.
        """
        n = len(ordered_names)
        V = np.zeros((n, n))

        clade_indices = {}
        clade_lengths = {}
        for idx, name in enumerate(ordered_names):
            cumulative = 0.0
            for clade_id, bl in paths[name]:
                d_start = cumulative
                d_end = cumulative + bl
                clade_indices.setdefault(clade_id, []).append(idx)
                if clade_id not in clade_lengths:
                    clade_lengths[clade_id] = transform_fn(bl, d_start, d_end)
                cumulative = d_end

        for clade_id, indices in clade_indices.items():
            branch_length = clade_lengths[clade_id]
            if len(indices) == 1:
                V[indices[0], indices[0]] += branch_length
                continue
            idx = np.asarray(indices, dtype=np.intp)
            V[np.ix_(idx, idx)] += branch_length

        return V

    # ── Optimization ─────────────────────────────────────────────────

    def _optimize_parameter(self, neg_ll_func, bounds, niter=10):
        """Multi-interval bounded optimization."""
        lo_bound, hi_bound = bounds
        bounds_lo = np.linspace(lo_bound, hi_bound - (hi_bound - lo_bound) / niter, niter)
        bounds_hi = np.linspace(lo_bound + (hi_bound - lo_bound) / niter, hi_bound, niter)

        best_ll = -np.inf
        best_param = (lo_bound + hi_bound) / 2.0

        for lo, hi in zip(bounds_lo, bounds_hi):
            if lo >= hi:
                continue
            try:
                res = minimize_scalar(neg_ll_func, bounds=(lo, hi), method="bounded")
                ll_val = -res.fun
                if ll_val > best_ll:
                    best_ll = ll_val
                    best_param = res.x
            except (ValueError, RuntimeError):
                continue

        return best_param, best_ll

    # ── Model fitting ────────────────────────────────────────────────

    def _fit_model(
        self, name: str, x: np.ndarray, C: np.ndarray,
        tree, ordered_names: list[str],
        paths: dict[str, list], max_lam: float, tree_height: float,
    ) -> dict:
        if name == "BM":
            return self._fit_bm(x, C)
        elif name == "OU":
            return self._fit_ou(x, C, tree_height)
        elif name == "EB":
            return self._fit_eb(x, ordered_names, paths, tree_height)
        elif name == "Lambda":
            return self._fit_lambda(x, C, max_lam)
        elif name == "Delta":
            return self._fit_delta(x, ordered_names, paths)
        elif name == "Kappa":
            return self._fit_kappa(x, ordered_names, paths)
        elif name == "White":
            return self._fit_white(x)
        else:
            raise PhykitUserError([f"Unknown model: {name}"], code=2)

    def _fit_bm(self, x: np.ndarray, C: np.ndarray) -> dict:
        ll, sig2, z0 = self._concentrated_ll(x, C)
        return dict(
            model="BM", param_name=None, param_value=None,
            sigma2=sig2, z0=z0, log_likelihood=ll, k_params=2,
        )

    def _fit_ou(self, x: np.ndarray, C: np.ndarray, tree_height: float) -> dict:
        upper = 100.0 / tree_height if tree_height > 0 else 100.0

        def neg_ll(alpha):
            V = self._vcv_ou(C, alpha)
            ll, _, _ = self._concentrated_ll(x, V)
            return -ll

        alpha_hat, _ = self._optimize_parameter(neg_ll, (1e-8, upper))
        V = self._vcv_ou(C, alpha_hat)
        ll, sig2, z0 = self._concentrated_ll(x, V)
        return dict(
            model="OU", param_name="alpha", param_value=float(alpha_hat),
            sigma2=sig2, z0=z0, log_likelihood=ll, k_params=3,
        )

    def _fit_eb(
        self, x: np.ndarray, ordered_names: list[str],
        paths: dict[str, list], tree_height: float,
    ) -> dict:
        upper = 10.0 / tree_height if tree_height > 0 else 10.0
        lower = -10.0 / tree_height if tree_height > 0 else -10.0

        def neg_ll(a):
            def transform(bl, d_start, d_end):
                if abs(a) < 1e-10:
                    return bl
                return (np.exp(a * d_end) - np.exp(a * d_start)) / a

            V = self._build_transformed_vcv(ordered_names, paths, transform)
            ll, _, _ = self._concentrated_ll(x, V)
            return -ll

        a_hat, _ = self._optimize_parameter(neg_ll, (lower, upper))

        def transform_final(bl, d_start, d_end):
            if abs(a_hat) < 1e-10:
                return bl
            return (np.exp(a_hat * d_end) - np.exp(a_hat * d_start)) / a_hat

        V = self._build_transformed_vcv(ordered_names, paths, transform_final)
        ll, sig2, z0 = self._concentrated_ll(x, V)
        return dict(
            model="EB", param_name="a", param_value=float(a_hat),
            sigma2=sig2, z0=z0, log_likelihood=ll, k_params=3,
        )

    def _fit_lambda(
        self, x: np.ndarray, C: np.ndarray, max_lam: float,
    ) -> dict:
        diag_vals = np.diag(C).copy()

        def neg_ll(lam):
            C_lam = C * lam
            np.fill_diagonal(C_lam, diag_vals)
            ll, _, _ = self._concentrated_ll(x, C_lam)
            return -ll

        lam_hat, _ = self._optimize_parameter(neg_ll, (0.0, max_lam))
        C_fitted = C * lam_hat
        np.fill_diagonal(C_fitted, diag_vals)
        ll, sig2, z0 = self._concentrated_ll(x, C_fitted)
        return dict(
            model="Lambda", param_name="lambda", param_value=float(lam_hat),
            sigma2=sig2, z0=z0, log_likelihood=ll, k_params=3,
        )

    def _fit_delta(
        self, x: np.ndarray, ordered_names: list[str],
        paths: dict[str, list],
    ) -> dict:
        def neg_ll(delta):
            def transform(bl, d_start, d_end):
                return d_end ** delta - d_start ** delta

            V = self._build_transformed_vcv(ordered_names, paths, transform)
            ll, _, _ = self._concentrated_ll(x, V)
            return -ll

        delta_hat, _ = self._optimize_parameter(neg_ll, (0.01, 3.0))

        def transform_final(bl, d_start, d_end):
            return d_end ** delta_hat - d_start ** delta_hat

        V = self._build_transformed_vcv(ordered_names, paths, transform_final)
        ll, sig2, z0 = self._concentrated_ll(x, V)
        return dict(
            model="Delta", param_name="delta", param_value=float(delta_hat),
            sigma2=sig2, z0=z0, log_likelihood=ll, k_params=3,
        )

    def _fit_kappa(
        self, x: np.ndarray, ordered_names: list[str],
        paths: dict[str, list],
    ) -> dict:
        def neg_ll(kappa):
            def transform(bl, d_start, d_end):
                return bl ** kappa

            V = self._build_transformed_vcv(ordered_names, paths, transform)
            ll, _, _ = self._concentrated_ll(x, V)
            return -ll

        kappa_hat, _ = self._optimize_parameter(neg_ll, (0.01, 3.0))

        def transform_final(bl, d_start, d_end):
            return bl ** kappa_hat

        V = self._build_transformed_vcv(ordered_names, paths, transform_final)
        ll, sig2, z0 = self._concentrated_ll(x, V)
        return dict(
            model="Kappa", param_name="kappa", param_value=float(kappa_hat),
            sigma2=sig2, z0=z0, log_likelihood=ll, k_params=3,
        )

    def _fit_white(self, x: np.ndarray) -> dict:
        n = len(x)
        C = np.eye(n)
        ll, sig2, z0 = self._concentrated_ll(x, C)
        return dict(
            model="White", param_name=None, param_value=None,
            sigma2=sig2, z0=z0, log_likelihood=ll, k_params=2,
        )

    # ── Model comparison ─────────────────────────────────────────────

    def _compute_model_comparison(self, results: list[dict], n: int) -> list[dict]:
        for r in results:
            k = r["k_params"]
            ll = r["log_likelihood"]
            r["aic"] = -2.0 * ll + 2.0 * k
            r["bic"] = -2.0 * ll + k * np.log(n)

        # Sort by AIC
        results.sort(key=lambda r: r["aic"])

        min_aic = results[0]["aic"]
        min_bic = min(r["bic"] for r in results)

        for r in results:
            r["delta_aic"] = r["aic"] - min_aic
            r["delta_bic"] = r["bic"] - min_bic

        # AIC weights
        raw_weights = [math.exp(-0.5 * r["delta_aic"]) for r in results]
        total = sum(raw_weights)
        if total > 0.0:
            for r, w in zip(results, raw_weights):
                r["aic_weight"] = w / total
        else:
            for r in results:
                r["aic_weight"] = 0.0

        return results

    # ── Output ───────────────────────────────────────────────────────

    def _print_text_output(self, results: list[dict], n: int) -> None:
        header = (
            f"{'Model':<12}{'Param':<10}{'Value':<11}"
            f"{'Sigma2':<10}{'z0':<10}{'LL':<11}"
            f"{'AIC':<9}{'dAIC':<9}{'AICw':<9}"
            f"{'BIC':<9}{'dBIC':<9}"
            f"{'R2':<7}"
        )

        lines = [
            "Model Comparison (fitContinuous)",
            f"\nNumber of tips: {n}\n",
            header,
        ]
        best_bic = results[0]
        best_bic_value = best_bic["bic"]
        row_format = (
            "%-12s%-10s%-11s"
            "%-10.4f%-10.4f%-11.3f"
            "%-9.2f%-9.2f%-9.3f"
            "%-9.2f%-9.2f%-7.3f"
        )
        for r in results:
            param_name = r["param_name"] if r["param_name"] else "-"
            if r["param_value"] is not None:
                param_val = f"{r['param_value']:.4f}"
            else:
                param_val = "-"
            bic = r["bic"]
            if bic < best_bic_value:
                best_bic = r
                best_bic_value = bic
            lines.append(
                row_format % (
                    r["model"],
                    param_name,
                    param_val,
                    r["sigma2"],
                    r["z0"],
                    r["log_likelihood"],
                    r["aic"],
                    r["delta_aic"],
                    r["aic_weight"],
                    bic,
                    r["delta_bic"],
                    r["r_squared"],
                )
            )

        best_aic = results[0]["model"]
        lines.append(f"\nBest model (AIC): {best_aic}")
        lines.append(f"Best model (BIC): {best_bic['model']}")
        print("\n".join(lines))

    def _print_json_output(self, results: list[dict], n: int, vcv_meta=None) -> None:
        models = {}
        for r in results:
            models[r["model"]] = dict(
                param_name=r["param_name"],
                param_value=r["param_value"],
                sigma2=r["sigma2"],
                z0=r["z0"],
                log_likelihood=r["log_likelihood"],
                aic=r["aic"],
                delta_aic=r["delta_aic"],
                aic_weight=r["aic_weight"],
                bic=r["bic"],
                delta_bic=r["delta_bic"],
                k_params=r["k_params"],
                r_squared=r["r_squared"],
            )

        best_aic = results[0]["model"]
        best_bic = min(results, key=lambda r: r["bic"])["model"]

        payload = dict(
            n_tips=n,
            models=models,
            best_model_aic=best_aic,
            best_model_bic=best_bic,
        )
        if vcv_meta is not None:
            payload["vcv_metadata"] = vcv_meta
        print_json(payload, sort_keys=False)
