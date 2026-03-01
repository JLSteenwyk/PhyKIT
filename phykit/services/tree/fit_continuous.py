import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize_scalar

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError

ALL_MODELS = ["BM", "OU", "EB", "Lambda", "Delta", "Kappa", "White"]


class FitContinuous(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.selected_models = parsed["models"]
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file()
        self._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)
        traits = self._parse_trait_file(self.trait_data_path, tree_tips)

        ordered_names = sorted(traits.keys())
        x = np.array([traits[name] for name in ordered_names])
        n = len(x)

        vcv = self._build_vcv_matrix(tree, ordered_names)

        # Precompute helpers needed by tree-transformation models
        parent_map = self._build_parent_map(tree)
        paths = self._build_root_to_tip_paths(tree, ordered_names, parent_map)
        max_lam = self._max_lambda(tree) if "Lambda" in self.selected_models else 1.0
        tree_height = float(np.max(np.diag(vcv)))

        results = []
        for model_name in self.selected_models:
            res = self._fit_model(
                model_name, x, vcv, tree, ordered_names,
                paths, max_lam, tree_height,
            )
            results.append(res)

        results = self._compute_model_comparison(results, n)

        if self.json_output:
            self._print_json_output(results, n)
        else:
            self._print_text_output(results, n)

    def process_args(self, args) -> Dict:
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
        )

    # ── Tree & trait parsing ─────────────────────────────────────────

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for model fitting."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                raise PhykitUserError(
                    ["All branches in the tree must have lengths."],
                    code=2,
                )

    def _parse_trait_file(
        self, path: str, tree_tips: List[str]
    ) -> Dict[str, float]:
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

        traits = {}
        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise PhykitUserError(
                    [
                        f"Line {line_num} in trait file has {len(parts)} columns; expected 2.",
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

    # ── VCV construction ─────────────────────────────────────────────

    def _build_vcv_matrix(
        self, tree, ordered_names: List[str]
    ) -> np.ndarray:
        n = len(ordered_names)
        vcv = np.zeros((n, n))

        root_to_tip = {}
        for name in ordered_names:
            root_to_tip[name] = tree.distance(tree.root, name)

        for i in range(n):
            for j in range(i, n):
                if i == j:
                    vcv[i, j] = root_to_tip[ordered_names[i]]
                else:
                    d_ij = tree.distance(ordered_names[i], ordered_names[j])
                    shared_path = (
                        root_to_tip[ordered_names[i]]
                        + root_to_tip[ordered_names[j]]
                        - d_ij
                    ) / 2.0
                    vcv[i, j] = shared_path
                    vcv[j, i] = shared_path

        return vcv

    def _build_parent_map(self, tree) -> Dict:
        parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                parent_map[id(child)] = clade
        return parent_map

    def _build_root_to_tip_paths(
        self, tree, ordered_names: List[str], parent_map: Dict
    ) -> Dict[str, List]:
        tip_map = {}
        for tip in tree.get_terminals():
            if tip.name in ordered_names:
                tip_map[tip.name] = tip

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
    ) -> Tuple[float, float, float]:
        """Concentrated log-likelihood with sigma^2 profiled out.

        Returns (log_likelihood, sigma2, z0).
        """
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

        n = C.shape[0]
        V = np.zeros((n, n))
        for i in range(n):
            for j in range(i, n):
                s_ij = C[i, j]  # shared path (MRCA depth)
                d_i = C[i, i] - s_ij  # unique path for tip i
                d_j = C[j, j] - s_ij  # unique path for tip j
                val = np.exp(-alpha * (d_i + d_j)) * (
                    1.0 - np.exp(-2.0 * alpha * s_ij)
                ) / (2.0 * alpha)
                V[i, j] = val
                V[j, i] = val
        return V

    def _vcv_lambda(self, C: np.ndarray, lam: float) -> np.ndarray:
        """Pagel's lambda VCV transformation."""
        diag_vals = np.diag(C).copy()
        C_lam = C * lam
        np.fill_diagonal(C_lam, diag_vals)
        return C_lam

    def _build_transformed_vcv(
        self, ordered_names: List[str], paths: Dict[str, List],
        transform_fn,
    ) -> np.ndarray:
        """Build VCV from transformed branch lengths.

        transform_fn(bl, d_start, d_end) -> new branch length
        where d_start is distance from root to start of branch,
        d_end is distance from root to end of branch.
        """
        n = len(ordered_names)
        V = np.zeros((n, n))

        # Precompute transformed paths with cumulative distances
        transformed_paths = {}
        for name in ordered_names:
            raw_path = paths[name]
            t_path = []
            cum = 0.0
            for clade_id, bl in raw_path:
                d_start = cum
                d_end = cum + bl
                new_bl = transform_fn(bl, d_start, d_end)
                t_path.append((clade_id, new_bl))
                cum += bl
            transformed_paths[name] = t_path

        for i in range(n):
            path_i = transformed_paths[ordered_names[i]]
            # Diagonal: sum of all transformed branch lengths
            V[i, i] = sum(bl for _, bl in path_i)

            for j in range(i + 1, n):
                path_j = transformed_paths[ordered_names[j]]
                # Shared prefix
                shared = 0.0
                min_len = min(len(path_i), len(path_j))
                for s in range(min_len):
                    if path_i[s][0] == path_j[s][0]:
                        shared += path_i[s][1]
                    else:
                        break
                V[i, j] = shared
                V[j, i] = shared

        return V

    # ── Max lambda helper ────────────────────────────────────────────

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
        tree, ordered_names: List[str],
        paths: Dict[str, List], max_lam: float, tree_height: float,
    ) -> Dict:
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

    def _fit_bm(self, x: np.ndarray, C: np.ndarray) -> Dict:
        ll, sig2, z0 = self._concentrated_ll(x, C)
        return dict(
            model="BM", param_name=None, param_value=None,
            sigma2=sig2, z0=z0, log_likelihood=ll, k_params=2,
        )

    def _fit_ou(self, x: np.ndarray, C: np.ndarray, tree_height: float) -> Dict:
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
        self, x: np.ndarray, ordered_names: List[str],
        paths: Dict[str, List], tree_height: float,
    ) -> Dict:
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
    ) -> Dict:
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
        self, x: np.ndarray, ordered_names: List[str],
        paths: Dict[str, List],
    ) -> Dict:
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
        self, x: np.ndarray, ordered_names: List[str],
        paths: Dict[str, List],
    ) -> Dict:
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

    def _fit_white(self, x: np.ndarray) -> Dict:
        n = len(x)
        C = np.eye(n)
        ll, sig2, z0 = self._concentrated_ll(x, C)
        return dict(
            model="White", param_name=None, param_value=None,
            sigma2=sig2, z0=z0, log_likelihood=ll, k_params=2,
        )

    # ── Model comparison ─────────────────────────────────────────────

    def _compute_model_comparison(self, results: List[Dict], n: int) -> List[Dict]:
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
        delta_aics = np.array([r["delta_aic"] for r in results])
        raw_weights = np.exp(-0.5 * delta_aics)
        total = raw_weights.sum()
        aic_weights = raw_weights / total if total > 0 else raw_weights

        for r, w in zip(results, aic_weights):
            r["aic_weight"] = float(w)

        return results

    # ── Output ───────────────────────────────────────────────────────

    def _print_text_output(self, results: List[Dict], n: int) -> None:
        print("Model Comparison (fitContinuous)")
        print(f"\nNumber of tips: {n}\n")

        header = (
            f"{'Model':<12}{'Param':<10}{'Value':<11}"
            f"{'Sigma2':<10}{'z0':<10}{'LL':<11}"
            f"{'AIC':<9}{'dAIC':<9}{'AICw':<9}"
            f"{'BIC':<9}{'dBIC':<9}"
        )
        print(header)

        for r in results:
            param_name = r["param_name"] if r["param_name"] else "-"
            if r["param_value"] is not None:
                param_val = f"{r['param_value']:.4f}"
            else:
                param_val = "-"
            print(
                f"{r['model']:<12}{param_name:<10}{param_val:<11}"
                f"{r['sigma2']:<10.4f}{r['z0']:<10.4f}{r['log_likelihood']:<11.3f}"
                f"{r['aic']:<9.2f}{r['delta_aic']:<9.2f}{r['aic_weight']:<9.3f}"
                f"{r['bic']:<9.2f}{r['delta_bic']:<9.2f}"
            )

        best_aic = results[0]["model"]
        best_bic = min(results, key=lambda r: r["bic"])["model"]
        print(f"\nBest model (AIC): {best_aic}")
        print(f"Best model (BIC): {best_bic}")

    def _print_json_output(self, results: List[Dict], n: int) -> None:
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
            )

        best_aic = results[0]["model"]
        best_bic = min(results, key=lambda r: r["bic"])["model"]

        payload = dict(
            n_tips=n,
            models=models,
            best_model_aic=best_aic,
            best_model_bic=best_bic,
        )
        print_json(payload, sort_keys=False)
