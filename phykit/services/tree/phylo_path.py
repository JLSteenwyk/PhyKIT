"""
Phylogenetic path analysis (von Hardenberg & Gonzalez-Voyer 2013).

Compares competing causal DAGs using d-separation tests via PGLS,
ranks models by CICc, and estimates (model-averaged) path coefficients.
"""
from __future__ import annotations

import math
import sys
from collections import defaultdict

from .base import Tree
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def trait_matrix_from_rows(*args, **kwargs):
    from ...helpers.trait_parsing import (
        trait_matrix_from_rows as _trait_matrix_from_rows,
    )

    return _trait_matrix_from_rows(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


_CHDTRC = None
_STDTR = None


def cho_factor(*args, **kwargs):
    from scipy.linalg import cho_factor as _cho_factor

    return _cho_factor(*args, **kwargs)


def cho_solve(*args, **kwargs):
    from scipy.linalg import cho_solve as _cho_solve

    return _cho_solve(*args, **kwargs)


def _chi2_sf(value: float, df: int) -> float:
    if df > 0 and df % 2 == 0:
        if value <= 0:
            return 1.0
        half_value = value / 2.0
        exp_term = math.exp(-half_value)
        if df == 2:
            return exp_term
        if df == 4:
            return exp_term * (1.0 + half_value)
        if df == 6:
            return exp_term * (1.0 + half_value + half_value * half_value / 2.0)
        term = 1.0
        total = 1.0
        try:
            for j in range(1, df // 2):
                term *= half_value / j
                total += term
            return exp_term * total
        except OverflowError:
            pass

    global _CHDTRC

    if _CHDTRC is None:
        from scipy.special import chdtrc as _chdtrc

        _CHDTRC = _chdtrc

    return float(_CHDTRC(df, value))


def _fishers_c_from_p_values(p_values: list[float]) -> float:
    return -2.0 * sum(math.log(p_value) for p_value in p_values)


def _t_two_tailed_p_value(t_stat: float, df: int) -> float:
    global _STDTR

    if _STDTR is None:
        from scipy.special import stdtr as _stdtr

        _STDTR = _stdtr

    return float(2.0 * _STDTR(df, -abs(t_stat)))


class PhyloPath(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.models_path = parsed["models_path"]
        self.best_only = parsed["best_only"]
        self.plot_output = parsed["plot_output"]
        self.csv_output = parsed["csv_output"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.traits,
            models_path=args.models,
            best_only=getattr(args, "best_only", False),
            plot_output=getattr(args, "plot_output", None),
            csv_output=getattr(args, "csv", None),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        from .vcv_utils import build_vcv_matrix

        tree = self.read_tree_file_unmodified()
        if self._needs_default_branch_lengths(tree):
            tree = self._fast_copy(tree)
        self.validate_tree(tree, min_tips=4, assign_default_branch_length=1e-8, context="path analysis")

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_names, traits = self._parse_trait_file(
            self.trait_data_path, tree_tips
        )
        models = self._parse_models_file(self.models_path)

        trait_name_to_idx = {}
        for idx, name in enumerate(trait_names):
            if name not in trait_name_to_idx:
                trait_name_to_idx[name] = idx

        # Validate all variables in models exist in trait file
        all_vars = set()
        for name, edges in models.items():
            for src, dst in edges:
                all_vars.add(src)
                all_vars.add(dst)
        for var in all_vars:
            if var not in trait_name_to_idx:
                raise PhykitUserError(
                    [
                        f"Variable '{var}' in model file not found in trait file.",
                        f"Available columns: {', '.join(trait_names)}",
                    ],
                    code=2,
                )

        variables = sorted(all_vars)

        # Build ordered data
        ordered_names = sorted(traits.keys())
        n = len(ordered_names)

        # Z-score standardize traits
        trait_matrix = trait_matrix_from_rows(traits, ordered_names)
        z_data = {}
        for var in variables:
            vals = trait_matrix[:, trait_name_to_idx[var]].copy()
            mean = float(np.mean(vals))
            sd = float(np.std(vals, ddof=0))
            vals -= mean
            if sd > 0:
                vals /= sd
            z_data[var] = vals

        # Build VCV + lambda estimation
        vcv = build_vcv_matrix(tree, ordered_names)

        # Evaluate each model
        model_results = {}
        for model_name, edges in models.items():
            adj = self._build_adjacency(variables, edges)
            if not self._is_dag(adj, variables):
                raise PhykitUserError(
                    [f"Model '{model_name}' contains a cycle and is not a DAG."],
                    code=2,
                )
            basis = self._basis_set(adj, variables)
            k = len(basis)  # number of independence claims

            # Test each d-separation claim
            p_values = []
            for (i_var, j_var, cond_set) in basis:
                p_val = self._dsep_test(
                    j_var, i_var, list(cond_set), z_data, vcv, n
                )
                p_values.append(max(p_val, sys.float_info.epsilon))

            # Fisher's C
            if k > 0:
                C = _fishers_c_from_p_values(p_values)
                model_p = _chi2_sf(C, 2 * k)
            else:
                C = 0.0
                model_p = 1.0

            # CICc
            n_edges = len(edges)
            n_vars = len(variables)
            q = n_vars + n_edges
            if n - 1 - q > 0:
                CICc = C + 2.0 * q * (n / (n - 1.0 - q))
            else:
                CICc = C + 2.0 * q

            model_results[model_name] = {
                "edges": edges,
                "k": k,
                "q": q,
                "C": C,
                "CICc": CICc,
                "p": model_p,
                "p_values": p_values,
            }

        # Rank and compute weights
        min_CICc = min(r["CICc"] for r in model_results.values())
        for r in model_results.values():
            r["delta"] = r["CICc"] - min_CICc
            r["rel_lik"] = np.exp(-0.5 * r["delta"])

        sum_lik = sum(r["rel_lik"] for r in model_results.values())
        for r in model_results.values():
            r["weight"] = r["rel_lik"] / sum_lik if sum_lik > 0 else 0.0

        ranked = sorted(model_results.items(), key=lambda x: x[1]["CICc"])

        # Estimate path coefficients for each model
        for model_name, r in model_results.items():
            r["coefficients"] = self._estimate_path_coefficients(
                r["edges"], variables, z_data, vcv, n
            )

        # Model averaging or best only
        if self.best_only:
            best_name = ranked[0][0]
            avg_coefs = model_results[best_name]["coefficients"]
            coef_label = f"best model: {best_name}"
        else:
            avg_coefs = self._model_average(model_results, variables)
            coef_label = "model-averaged, conditional"

        # Output
        self._print_text(ranked, avg_coefs, coef_label, n)

        if self.json_output:
            self._print_json(ranked, avg_coefs, coef_label, variables, n)

        if self.csv_output:
            self._write_csv(ranked, avg_coefs, self.csv_output)

        if self.plot_output:
            best_edges = ranked[0][1]["edges"] if self.best_only else None
            self._plot_dag(
                avg_coefs, variables, self.plot_output, best_edges
            )

    @staticmethod
    def _needs_default_branch_lengths(tree) -> bool:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return True

        stack = [root]
        pop = stack.pop
        extend = stack.extend
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if clade is not root and clade.branch_length is None:
                    return True
                if children:
                    extend(children)
        except AttributeError:
            return True

        return False

    # ---- File parsing ----

    def _parse_trait_file(
        self, path: str, tree_tips: list[str]
    ) -> tuple[list[str], dict[str, list[float]]]:
        try:
            with open(path) as f:
                header_parts = None
                for line in f:
                    stripped = line.strip()
                    if stripped and not stripped.startswith("#"):
                        header_parts = stripped.split("\t")
                        break

                if header_parts is None:
                    raise PhykitUserError(
                        ["Trait file must have a header row and at least one data row."],
                        code=2,
                    )

                trait_names = header_parts[1:]

                traits = {}
                line_idx = 2
                saw_data = False
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue
                    saw_data = True
                    parts = stripped.split("\t")
                    taxon = parts[0]
                    value_parts = parts[1:]
                    try:
                        values = [float(val_str) for val_str in value_parts]
                    except ValueError:
                        values = []
                        for val_str in value_parts:
                            try:
                                values.append(float(val_str))
                            except ValueError:
                                raise PhykitUserError(
                                    [
                                        f"Non-numeric value '{val_str}' for taxon "
                                        f"'{taxon}' on line {line_idx}.",
                                    ],
                                    code=2,
                                )
                    traits[taxon] = values
                    line_idx += 1

                if not saw_data:
                    raise PhykitUserError(
                        ["Trait file must have a header row and at least one data row."],
                        code=2,
                    )
        except FileNotFoundError:
            raise PhykitUserError(
                [f"{path} corresponds to no such file or directory."],
                code=2,
            )

        tree_tip_set = set(tree_tips)
        if (
            len(tree_tip_set) >= 4
            and len(tree_tip_set) == len(traits)
            and tree_tip_set == traits.keys()
        ):
            return trait_names, traits

        shared = tree_tip_set & set(traits)
        if len(shared) < 4:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa. "
                    "At least 4 are required for path analysis.",
                ],
                code=2,
            )

        filtered = {t: traits[t] for t in shared}
        return trait_names, filtered

    def _parse_models_file(self, path: str) -> dict[str, list[tuple[str, str]]]:
        try:
            models = {}
            with open(path) as f:
                for line_num, line in enumerate(f, 1):
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue
                    if ":" not in stripped:
                        raise PhykitUserError(
                            [
                                f"Line {line_num} in model file is malformed.",
                                "Expected format: model_name: A->B, B->C, ...",
                            ],
                            code=2,
                        )
                    name, edge_str = stripped.split(":", 1)
                    name = name.strip()
                    edges = []
                    for part in edge_str.split(","):
                        part = part.strip()
                        if not part:
                            continue
                        if "->" not in part:
                            raise PhykitUserError(
                                [
                                    f"Invalid edge '{part}' in model '{name}'.",
                                    "Edges must use '->' syntax (e.g., A->B).",
                                ],
                                code=2,
                            )
                        src, dst = part.split("->", 1)
                        edges.append((src.strip(), dst.strip()))
                    if not edges:
                        raise PhykitUserError(
                            [f"Model '{name}' has no edges."], code=2
                        )
                    models[name] = edges
        except FileNotFoundError:
            raise PhykitUserError(
                [f"{path} corresponds to no such file or directory."],
                code=2,
            )

        if len(models) < 2:
            raise PhykitUserError(
                ["At least 2 candidate models are required for comparison."],
                code=2,
            )

        return models

    # ---- DAG operations ----

    @staticmethod
    def _build_adjacency(
        variables: list[str], edges: list[tuple[str, str]]
    ) -> dict[str, set]:
        adj = {v: set() for v in variables}
        for src, dst in edges:
            adj[src].add(dst)
        return adj

    @staticmethod
    def _is_dag(adj: dict[str, set], variables: list[str]) -> bool:
        """Check for cycles using topological sort (Kahn's algorithm)."""
        in_degree = {v: 0 for v in variables}
        for v in variables:
            for child in adj[v]:
                in_degree[child] += 1

        queue = [v for v in variables if in_degree[v] == 0]
        visited = 0
        queue_idx = 0
        while queue_idx < len(queue):
            node = queue[queue_idx]
            queue_idx += 1
            visited += 1
            for child in adj[node]:
                in_degree[child] -= 1
                if in_degree[child] == 0:
                    queue.append(child)

        return visited == len(variables)

    @staticmethod
    def _topological_order(
        adj: dict[str, set], variables: list[str]
    ) -> list[str]:
        in_degree = {v: 0 for v in variables}
        for v in variables:
            for child in adj[v]:
                in_degree[child] += 1
        queue = sorted([v for v in variables if in_degree[v] == 0])
        order = []
        queue_idx = 0
        while queue_idx < len(queue):
            node = queue[queue_idx]
            queue_idx += 1
            order.append(node)
            for child in sorted(adj[node]):
                in_degree[child] -= 1
                if in_degree[child] == 0:
                    queue.append(child)
        return order

    @staticmethod
    def _get_parents(
        adj: dict[str, set], variables: list[str], node: str
    ) -> set:
        return {v for v in variables if node in adj[v]}

    def _basis_set(
        self, adj: dict[str, set], variables: list[str]
    ) -> list[tuple[str, str, set]]:
        """Derive the basis set of d-separation statements (Shipley 2000)."""
        order = self._topological_order(adj, variables)

        # Build existing edges and parent sets in one pass over the model graph.
        adjacent = set()
        parents = {v: set() for v in variables}
        for src in variables:
            for dst in adj[src]:
                adjacent.add((src, dst))
                adjacent.add((dst, src))
                parents[dst].add(src)

        basis = []
        for i, vi in enumerate(order):
            parents_i = parents[vi]
            for j in range(i + 1, len(order)):
                vj = order[j]
                # Check non-adjacent
                if (vi, vj) in adjacent:
                    continue
                # Conditioning set: parents of both, excluding i and j
                cond = (parents_i | parents[vj]) - {vi, vj}
                basis.append((vi, vj, cond))

        return basis

    # ---- PGLS with lambda ----

    def _dsep_test(
        self, response: str, predictor: str,
        conditioning: list[str], z_data: dict[str, np.ndarray],
        vcv: np.ndarray, n: int,
    ) -> float:
        """Test independence of predictor and response via PGLS. Returns p-value."""
        from ...helpers.pgls_utils import apply_lambda, estimate_lambda

        y = z_data[response]
        predictors = conditioning + [predictor]
        X = self._design_matrix_from_z(z_data, predictors, n)

        lambda_val, _ = estimate_lambda(y, X, vcv)
        vcv_lam = apply_lambda(vcv, lambda_val)

        try:
            beta_hat, residuals, sigma2, var_beta = self._fit_gls_from_vcv(
                y, X, vcv_lam
            )
        except (SystemExit, np.linalg.LinAlgError):
            return 1.0

        df_resid = n - X.shape[1]
        if df_resid <= 0:
            return 1.0

        pred_idx = X.shape[1] - 1
        se = np.sqrt(max(var_beta[pred_idx, pred_idx], 0.0))
        if se <= 0:
            return 1.0
        t_stat = beta_hat[pred_idx] / se
        p_val = _t_two_tailed_p_value(t_stat, df_resid)

        return float(p_val)

    @staticmethod
    def _design_matrix_from_z(
        z_data: dict[str, np.ndarray],
        predictors: list[str],
        n: int,
    ) -> np.ndarray:
        X = np.empty((n, len(predictors) + 1), dtype=np.float64)
        X[:, 0] = 1.0
        for idx, predictor in enumerate(predictors, 1):
            X[:, idx] = z_data[predictor]
        return X

    def _pgls_fit(
        self, y: np.ndarray, X: np.ndarray,
        vcv: np.ndarray, n: int,
    ) -> tuple[np.ndarray, np.ndarray, float]:
        """Fit PGLS with lambda. Returns (beta, se, sigma2)."""
        from ...helpers.pgls_utils import apply_lambda, estimate_lambda

        lambda_val, _ = estimate_lambda(y, X, vcv)
        vcv_lam = apply_lambda(vcv, lambda_val)

        beta_hat, residuals, sigma2, var_beta = self._fit_gls_from_vcv(
            y, X, vcv_lam
        )
        se = np.sqrt(np.maximum(np.diag(var_beta), 0.0))

        return beta_hat, se, sigma2

    def _fit_gls_from_vcv(
        self, y: np.ndarray, X: np.ndarray, vcv: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
        try:
            factor = cho_factor(vcv, lower=True, check_finite=False)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._fit_gls_from_vcv_inverse(y, X, vcv)

        n_predictors = X.shape[1]
        solve_rhs = np.empty(
            (len(y), n_predictors + 1),
            dtype=np.result_type(X, y, vcv),
        )
        solve_rhs[:, :n_predictors] = X
        solve_rhs[:, n_predictors] = y
        solved = cho_solve(factor, solve_rhs, check_finite=False)
        C_inv_X = solved[:, :n_predictors]
        C_inv_y = solved[:, n_predictors]

        XtCiX = X.T @ C_inv_X
        try:
            normal_rhs = np.empty(
                (n_predictors, n_predictors + 1),
                dtype=XtCiX.dtype,
            )
            normal_rhs[:, 0] = X.T @ C_inv_y
            normal_rhs[:, 1:] = np.eye(n_predictors, dtype=XtCiX.dtype)
            normal_solved = np.linalg.solve(XtCiX, normal_rhs)
        except np.linalg.LinAlgError:
            raise PhykitUserError(
                [
                    "Singular design matrix: cannot estimate coefficients.",
                    "Check that predictors are not collinear.",
                ],
                code=2,
            )

        beta_hat = normal_solved[:, 0]
        XtCiX_inv = normal_solved[:, 1:]
        residuals = y - X @ beta_hat
        C_inv_residuals = C_inv_y - C_inv_X @ beta_hat

        df_resid = len(y) - X.shape[1]
        sigma2 = float(residuals @ C_inv_residuals) / max(df_resid, 1)
        var_beta = sigma2 * XtCiX_inv

        return beta_hat, residuals, sigma2, var_beta

    def _fit_gls_from_vcv_inverse(
        self, y: np.ndarray, X: np.ndarray, vcv: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
        from ...helpers.pgls_utils import fit_gls

        C_inv = np.linalg.inv(vcv)
        return fit_gls(y, X, C_inv)

    # ---- Path coefficient estimation ----

    def _estimate_path_coefficients(
        self, edges: list[tuple[str, str]],
        variables: list[str],
        z_data: dict[str, np.ndarray],
        vcv: np.ndarray, n: int,
    ) -> dict[str, dict]:
        """For each edge i->j, fit PGLS j ~ parents(j) and extract coef for i."""
        adj = self._build_adjacency(variables, edges)
        coefficients = {}

        # Group edges by response variable
        response_parents = defaultdict(list)
        for src, dst in edges:
            response_parents[dst].append(src)

        for response, parents in response_parents.items():
            y = z_data[response]
            X = self._design_matrix_from_z(z_data, parents, n)
            beta, se, sigma2 = self._pgls_fit(y, X, vcv, n)

            df_resid = n - X.shape[1]
            for idx, parent in enumerate(parents):
                coef_idx = idx + 1  # skip intercept
                t_stat = (
                    beta[coef_idx] / se[coef_idx]
                    if se[coef_idx] > 0 else 0.0
                )
                p_val = (
                    _t_two_tailed_p_value(t_stat, max(df_resid, 1))
                    if se[coef_idx] > 0 else 1.0
                )
                edge_key = f"{parent} -> {response}"
                coefficients[edge_key] = {
                    "coef": float(beta[coef_idx]),
                    "se": float(se[coef_idx]),
                    "p": float(p_val),
                }

        return coefficients

    # ---- Model averaging ----

    def _model_average(
        self, model_results: dict, variables: list[str]
    ) -> dict[str, dict]:
        """Conditional model averaging of path coefficients."""
        edge_entries = defaultdict(list)
        for r in model_results.values():
            weight = r["weight"]
            for edge, vals in r["coefficients"].items():
                edge_entries[edge].append((vals["coef"], vals["se"], weight))

        avg_coefs = {}
        for edge in sorted(edge_entries):
            entries = edge_entries[edge]
            w_sum = sum(weight for _, _, weight in entries)
            if w_sum <= 0:
                continue

            avg_coef = (
                sum(weight * coef for coef, _, weight in entries) / w_sum
            )
            # Combined SE: within-model + between-model variance
            avg_se = math.sqrt(
                sum(
                    weight * (se * se + (coef - avg_coef) ** 2)
                    for coef, se, weight in entries
                ) / w_sum
            )

            avg_coefs[edge] = {
                "coef": float(avg_coef),
                "se": float(avg_se),
                "lower": float(avg_coef - 1.96 * avg_se),
                "upper": float(avg_coef + 1.96 * avg_se),
            }

        return avg_coefs

    # ---- Output ----

    def _print_text(self, ranked, avg_coefs, coef_label, n) -> None:
        try:
            lines = [
                f"Phylogenetic Path Analysis (n = {n} taxa)",
                "",
                "Model comparison (ranked by CICc):",
                f"  {'Model':<20s} {'k':>4s} {'q':>4s} "
                f"{'C':>10s} {'CICc':>10s} {'delta':>8s} "
                f"{'weight':>8s} {'p':>8s}",
                "  " + "-" * 74,
            ]
            append = lines.append
            model_row = (
                "  %-20s %4d %4d %10.4f %10.4f "
                "%8.4f %8.4f %8.4f"
            )
            for name, r in ranked:
                append(
                    model_row % (
                        name,
                        r["k"],
                        r["q"],
                        r["C"],
                        r["CICc"],
                        r["delta"],
                        r["weight"],
                        r["p"],
                    )
                )

            best_name = ranked[0][0]
            append("")
            append(
                f"Best model: {best_name} "
                f"(CICc weight = {ranked[0][1]['weight']:.4f})"
            )

            append("")
            append(f"Path coefficients ({coef_label}):")
            append(
                f"  {'Path':<25s} {'coef':>10s} {'SE':>10s} "
                f"{'lower':>10s} {'upper':>10s}"
            )
            append("  " + "-" * 67)
            coef_row = "  %-25s %10.4f %10.4f %10.4f %10.4f %s"
            for edge, vals in sorted(avg_coefs.items()):
                coef = vals["coef"]
                se = vals["se"]
                lower = vals.get("lower", coef - 1.96 * se)
                upper = vals.get("upper", coef + 1.96 * se)
                sig = "*" if vals.get("p", 0) < 0.05 else " "
                if "p" in vals:
                    sig = (
                        "***" if vals["p"] < 0.001
                        else "**" if vals["p"] < 0.01
                        else "*" if vals["p"] < 0.05
                        else ""
                    )
                append(
                    coef_row % (edge, coef, se, lower, upper, sig)
                )
            print("\n".join(lines))

        except BrokenPipeError:
            return

    def _print_json(self, ranked, avg_coefs, coef_label, variables, n):
        models = []
        for name, r in ranked:
            models.append({
                "model": name,
                "k": r["k"],
                "q": r["q"],
                "C": round(r["C"], 4),
                "CICc": round(r["CICc"], 4),
                "delta_CICc": round(r["delta"], 4),
                "weight": round(r["weight"], 4),
                "p": round(r["p"], 4),
                "edges": [f"{s}->{d}" for s, d in r["edges"]],
                "coefficients": {
                    k: {kk: round(vv, 4) for kk, vv in v.items()}
                    for k, v in r["coefficients"].items()
                },
            })

        coefs = {}
        for edge, vals in sorted(avg_coefs.items()):
            coefs[edge] = {k: round(v, 4) for k, v in vals.items()}

        payload = {
            "n_taxa": n,
            "variables": variables,
            "n_models": len(ranked),
            "best_model": ranked[0][0],
            "model_comparison": models,
            "path_coefficients": coefs,
            "coefficient_type": coef_label,
        }
        print_json(payload, sort_keys=False)

    def _write_csv(self, ranked, avg_coefs, csv_path):
        with open(csv_path, "w") as f:
            f.write("model,k,q,C,CICc,delta_CICc,weight,p\n")
            for name, r in ranked:
                f.write(
                    f"{name},{r['k']},{r['q']},"
                    f"{r['C']:.4f},{r['CICc']:.4f},"
                    f"{r['delta']:.4f},{r['weight']:.4f},{r['p']:.4f}\n"
                )
            f.write("\npath,coef,se,lower,upper\n")
            for edge, vals in sorted(avg_coefs.items()):
                lower = vals.get("lower", vals["coef"] - 1.96 * vals["se"])
                upper = vals.get("upper", vals["coef"] + 1.96 * vals["se"])
                f.write(
                    f"{edge},{vals['coef']:.4f},{vals['se']:.4f},"
                    f"{lower:.4f},{upper:.4f}\n"
                )

    def _plot_dag(
        self, avg_coefs, variables, output_path, best_edges=None,
    ) -> None:
        """Plot DAG with path coefficients on edges."""
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.collections import PatchCollection
            from matplotlib.patches import FancyArrowPatch, Circle
        except ImportError:
            print("matplotlib is required for plotting.")
            return

        config = self.plot_config
        config.resolve(n_rows=len(variables), n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        # Determine edges from avg_coefs
        edges = []
        for edge_key, vals in avg_coefs.items():
            parts = edge_key.split(" -> ")
            if len(parts) == 2:
                edges.append((parts[0], parts[1], vals["coef"], vals["se"]))

        # Simple circular layout for nodes
        n_vars = len(variables)
        angles = np.linspace(0, 2 * np.pi, n_vars, endpoint=False)
        # Rotate so first node is at top
        angles = angles + np.pi / 2
        radius = 1.5
        pos = {}
        for i, var in enumerate(variables):
            pos[var] = (radius * np.cos(angles[i]), radius * np.sin(angles[i]))

        # Draw nodes
        node_radius = 0.25
        node_patches = [
            Circle((x, y), node_radius)
            for x, y in pos.values()
        ]
        if node_patches:
            ax.add_collection(
                PatchCollection(
                    node_patches,
                    facecolor="white",
                    edgecolor="black",
                    linewidths=1.5,
                    zorder=3,
                )
            )
        for var, (x, y) in pos.items():
            ax.text(
                x, y, var, ha="center", va="center",
                fontsize=9, fontweight="bold", zorder=4,
            )

        # Draw edges
        for src, dst, coef, se in edges:
            x0, y0 = pos[src]
            x1, y1 = pos[dst]

            # Shorten to node boundary
            dx = x1 - x0
            dy = y1 - y0
            dist = np.sqrt(dx**2 + dy**2)
            if dist == 0:
                continue
            ux, uy = dx / dist, dy / dist
            x0s = x0 + ux * (node_radius + 0.05)
            y0s = y0 + uy * (node_radius + 0.05)
            x1s = x1 - ux * (node_radius + 0.05)
            y1s = y1 - uy * (node_radius + 0.05)

            color = "#2b8cbe" if coef >= 0 else "#d62728"
            lw = max(1.0, min(4.0, abs(coef) * 4))

            arrow = FancyArrowPatch(
                (x0s, y0s), (x1s, y1s),
                arrowstyle="->,head_width=6,head_length=6",
                color=color, linewidth=lw,
                connectionstyle="arc3,rad=0.1",
                zorder=2,
            )
            ax.add_patch(arrow)

            # Label
            mx = (x0s + x1s) / 2 + uy * 0.15
            my = (y0s + y1s) / 2 - ux * 0.15
            sig = ""
            ax.text(
                mx, my, f"{coef:.2f}{sig}",
                ha="center", va="center", fontsize=8,
                color=color, fontweight="bold",
                bbox=dict(
                    facecolor="white", edgecolor="none",
                    alpha=0.8, pad=1,
                ),
                zorder=5,
            )

        ax.set_aspect("equal")
        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(-2.5, 2.5)
        ax.axis("off")

        if config.show_title:
            ax.set_title(
                config.title or "Phylogenetic Path Analysis",
                fontsize=config.title_fontsize,
            )

        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Plot saved: {output_path}")
