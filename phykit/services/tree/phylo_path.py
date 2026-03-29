"""
Phylogenetic path analysis (von Hardenberg & Gonzalez-Voyer 2013).

Compares competing causal DAGs using d-separation tests via PGLS,
ranks models by CICc, and estimates (model-averaged) path coefficients.
"""
import sys
from collections import defaultdict
from typing import Dict, List, Tuple

import numpy as np
from scipy.stats import chi2, t as t_dist

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig
from ...helpers.pgls_utils import (
    estimate_lambda,
    pgls_log_likelihood,
    fit_gls,
    apply_lambda,
)
from ...errors import PhykitUserError


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

    def process_args(self, args) -> Dict:
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

        tree = self.read_tree_file()
        self._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_names, traits = self._parse_trait_file(
            self.trait_data_path, tree_tips
        )
        models = self._parse_models_file(self.models_path)

        # Validate all variables in models exist in trait file
        all_vars = set()
        for name, edges in models.items():
            for src, dst in edges:
                all_vars.add(src)
                all_vars.add(dst)
        for var in all_vars:
            if var not in trait_names:
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
        raw_data = {}
        for var in variables:
            idx = trait_names.index(var)
            vals = np.array([traits[name][idx] for name in ordered_names])
            raw_data[var] = vals

        z_data = {}
        for var in variables:
            vals = raw_data[var]
            sd = np.std(vals, ddof=0)
            if sd > 0:
                z_data[var] = (vals - np.mean(vals)) / sd
            else:
                z_data[var] = vals - np.mean(vals)

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
                C = -2.0 * sum(np.log(p) for p in p_values)
                model_p = 1.0 - chi2.cdf(C, 2 * k)
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

    # ---- Tree validation ----

    def _validate_tree(self, tree) -> None:
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                clade.branch_length = 1e-8

    # ---- File parsing ----

    def _parse_trait_file(
        self, path: str, tree_tips: List[str]
    ) -> Tuple[List[str], Dict[str, List[float]]]:
        try:
            with open(path) as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise PhykitUserError(
                [f"{path} corresponds to no such file or directory."],
                code=2,
            )

        data_lines = [
            l.strip() for l in lines
            if l.strip() and not l.strip().startswith("#")
        ]
        if len(data_lines) < 2:
            raise PhykitUserError(
                ["Trait file must have a header row and at least one data row."],
                code=2,
            )

        header_parts = data_lines[0].split("\t")
        trait_names = header_parts[1:]

        traits = {}
        for line_idx, line in enumerate(data_lines[1:], 2):
            parts = line.split("\t")
            taxon = parts[0]
            values = []
            for i, val_str in enumerate(parts[1:]):
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

        tree_tip_set = set(tree_tips)
        shared = tree_tip_set & set(traits.keys())
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

    def _parse_models_file(self, path: str) -> Dict[str, List[Tuple[str, str]]]:
        try:
            with open(path) as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise PhykitUserError(
                [f"{path} corresponds to no such file or directory."],
                code=2,
            )

        models = {}
        for line_num, line in enumerate(lines, 1):
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

        if len(models) < 2:
            raise PhykitUserError(
                ["At least 2 candidate models are required for comparison."],
                code=2,
            )

        return models

    # ---- DAG operations ----

    @staticmethod
    def _build_adjacency(
        variables: List[str], edges: List[Tuple[str, str]]
    ) -> Dict[str, set]:
        adj = {v: set() for v in variables}
        for src, dst in edges:
            adj[src].add(dst)
        return adj

    @staticmethod
    def _is_dag(adj: Dict[str, set], variables: List[str]) -> bool:
        """Check for cycles using topological sort (Kahn's algorithm)."""
        in_degree = {v: 0 for v in variables}
        for v in variables:
            for child in adj[v]:
                in_degree[child] += 1

        queue = [v for v in variables if in_degree[v] == 0]
        visited = 0
        while queue:
            node = queue.pop(0)
            visited += 1
            for child in adj[node]:
                in_degree[child] -= 1
                if in_degree[child] == 0:
                    queue.append(child)

        return visited == len(variables)

    @staticmethod
    def _topological_order(
        adj: Dict[str, set], variables: List[str]
    ) -> List[str]:
        in_degree = {v: 0 for v in variables}
        for v in variables:
            for child in adj[v]:
                in_degree[child] += 1
        queue = sorted([v for v in variables if in_degree[v] == 0])
        order = []
        while queue:
            node = queue.pop(0)
            order.append(node)
            for child in sorted(adj[node]):
                in_degree[child] -= 1
                if in_degree[child] == 0:
                    queue.append(child)
        return order

    @staticmethod
    def _get_parents(
        adj: Dict[str, set], variables: List[str], node: str
    ) -> set:
        return {v for v in variables if node in adj[v]}

    def _basis_set(
        self, adj: Dict[str, set], variables: List[str]
    ) -> List[Tuple[str, str, set]]:
        """Derive the basis set of d-separation statements (Shipley 2000)."""
        order = self._topological_order(adj, variables)
        order_idx = {v: i for i, v in enumerate(order)}

        # Build set of existing edges (both directions for adjacency check)
        adjacent = set()
        for src in variables:
            for dst in adj[src]:
                adjacent.add((src, dst))
                adjacent.add((dst, src))

        basis = []
        for i, vi in enumerate(order):
            for j, vj in enumerate(order):
                if j <= i:
                    continue
                # Check non-adjacent
                if (vi, vj) in adjacent:
                    continue
                # Conditioning set: parents of both, excluding i and j
                parents_i = self._get_parents(adj, variables, vi)
                parents_j = self._get_parents(adj, variables, vj)
                cond = (parents_i | parents_j) - {vi, vj}
                basis.append((vi, vj, cond))

        return basis

    # ---- PGLS with lambda ----

    def _dsep_test(
        self, response: str, predictor: str,
        conditioning: List[str], z_data: Dict[str, np.ndarray],
        vcv: np.ndarray, n: int,
    ) -> float:
        """Test independence of predictor and response via PGLS. Returns p-value."""
        y = z_data[response]
        predictors = conditioning + [predictor]
        X = np.column_stack(
            [np.ones(n)] + [z_data[p] for p in predictors]
        )

        lambda_val, _ = estimate_lambda(y, X, vcv)
        vcv_lam = apply_lambda(vcv, lambda_val)

        try:
            C_inv = np.linalg.inv(vcv_lam)
        except np.linalg.LinAlgError:
            return 1.0

        try:
            beta_hat, residuals, sigma2, var_beta = fit_gls(y, X, C_inv)
        except SystemExit:
            return 1.0

        df_resid = n - X.shape[1]
        if df_resid <= 0:
            return 1.0

        pred_idx = X.shape[1] - 1
        se = np.sqrt(max(var_beta[pred_idx, pred_idx], 0.0))
        if se <= 0:
            return 1.0
        t_stat = beta_hat[pred_idx] / se
        p_val = 2.0 * t_dist.sf(abs(t_stat), df=df_resid)

        return float(p_val)

    def _pgls_fit(
        self, y: np.ndarray, X: np.ndarray,
        vcv: np.ndarray, n: int,
    ) -> Tuple[np.ndarray, np.ndarray, float]:
        """Fit PGLS with lambda. Returns (beta, se, sigma2)."""
        lambda_val, _ = estimate_lambda(y, X, vcv)
        vcv_lam = apply_lambda(vcv, lambda_val)

        C_inv = np.linalg.inv(vcv_lam)
        beta_hat, residuals, sigma2, var_beta = fit_gls(y, X, C_inv)
        se = np.sqrt(np.maximum(np.diag(var_beta), 0.0))

        return beta_hat, se, sigma2

    # ---- Path coefficient estimation ----

    def _estimate_path_coefficients(
        self, edges: List[Tuple[str, str]],
        variables: List[str],
        z_data: Dict[str, np.ndarray],
        vcv: np.ndarray, n: int,
    ) -> Dict[str, Dict]:
        """For each edge i->j, fit PGLS j ~ parents(j) and extract coef for i."""
        adj = self._build_adjacency(variables, edges)
        coefficients = {}

        # Group edges by response variable
        response_parents = defaultdict(list)
        for src, dst in edges:
            response_parents[dst].append(src)

        for response, parents in response_parents.items():
            y = z_data[response]
            X = np.column_stack(
                [np.ones(n)] + [z_data[p] for p in parents]
            )
            beta, se, sigma2 = self._pgls_fit(y, X, vcv, n)

            df_resid = n - X.shape[1]
            for idx, parent in enumerate(parents):
                coef_idx = idx + 1  # skip intercept
                t_stat = (
                    beta[coef_idx] / se[coef_idx]
                    if se[coef_idx] > 0 else 0.0
                )
                p_val = (
                    2.0 * t_dist.sf(abs(t_stat), df=max(df_resid, 1))
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
        self, model_results: Dict, variables: List[str]
    ) -> Dict[str, Dict]:
        """Conditional model averaging of path coefficients."""
        all_edges = set()
        for r in model_results.values():
            all_edges.update(r["coefficients"].keys())

        avg_coefs = {}
        for edge in sorted(all_edges):
            coefs = []
            ses = []
            weights = []
            for r in model_results.values():
                if edge in r["coefficients"]:
                    coefs.append(r["coefficients"][edge]["coef"])
                    ses.append(r["coefficients"][edge]["se"])
                    weights.append(r["weight"])

            if not coefs:
                continue

            w = np.array(weights)
            w_sum = w.sum()
            if w_sum <= 0:
                continue
            w_norm = w / w_sum

            c = np.array(coefs)
            s = np.array(ses)

            avg_coef = float(np.sum(w_norm * c))
            # Combined SE: within-model + between-model variance
            avg_se = float(np.sqrt(
                np.sum(w_norm * (s**2 + (c - avg_coef)**2))
            ))

            avg_coefs[edge] = {
                "coef": avg_coef,
                "se": avg_se,
                "lower": avg_coef - 1.96 * avg_se,
                "upper": avg_coef + 1.96 * avg_se,
            }

        return avg_coefs

    # ---- Output ----

    def _print_text(self, ranked, avg_coefs, coef_label, n) -> None:
        try:
            print(f"Phylogenetic Path Analysis (n = {n} taxa)")
            print(f"\nModel comparison (ranked by CICc):")
            print(
                f"  {'Model':<20s} {'k':>4s} {'q':>4s} "
                f"{'C':>10s} {'CICc':>10s} {'delta':>8s} "
                f"{'weight':>8s} {'p':>8s}"
            )
            print("  " + "-" * 74)
            for name, r in ranked:
                print(
                    f"  {name:<20s} {r['k']:>4d} {r['q']:>4d} "
                    f"{r['C']:>10.4f} {r['CICc']:>10.4f} "
                    f"{r['delta']:>8.4f} {r['weight']:>8.4f} "
                    f"{r['p']:>8.4f}"
                )

            best_name = ranked[0][0]
            print(f"\nBest model: {best_name} "
                  f"(CICc weight = {ranked[0][1]['weight']:.4f})")

            print(f"\nPath coefficients ({coef_label}):")
            print(
                f"  {'Path':<25s} {'coef':>10s} {'SE':>10s} "
                f"{'lower':>10s} {'upper':>10s}"
            )
            print("  " + "-" * 67)
            for edge, vals in sorted(avg_coefs.items()):
                lower = vals.get("lower", vals["coef"] - 1.96 * vals["se"])
                upper = vals.get("upper", vals["coef"] + 1.96 * vals["se"])
                sig = "*" if vals.get("p", 0) < 0.05 else " "
                if "p" in vals:
                    sig = (
                        "***" if vals["p"] < 0.001
                        else "**" if vals["p"] < 0.01
                        else "*" if vals["p"] < 0.05
                        else ""
                    )
                print(
                    f"  {edge:<25s} {vals['coef']:>10.4f} "
                    f"{vals['se']:>10.4f} {lower:>10.4f} "
                    f"{upper:>10.4f} {sig}"
                )

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
        for var, (x, y) in pos.items():
            circle = Circle(
                (x, y), node_radius, facecolor="white",
                edgecolor="black", linewidth=1.5, zorder=3,
            )
            ax.add_patch(circle)
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
