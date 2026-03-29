import copy
import os
import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize_scalar

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.pgls_utils import max_lambda as compute_max_lambda
from ...helpers.plot_config import PlotConfig
from ...errors import PhykitUserError


class PhylogeneticOrdination(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.method = parsed["method"]
        self.correction = parsed["correction"]
        self.mode = parsed["mode"]
        self.n_components = parsed["n_components"]
        self.perplexity = parsed["perplexity"]
        self.n_neighbors = parsed["n_neighbors"]
        self.min_dist = parsed["min_dist"]
        self.seed = parsed["seed"]
        self.json_output = parsed["json_output"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]
        self.plot_tree = parsed["plot_tree"]
        self.color_by = parsed["color_by"]
        self.tree_color_by = parsed["tree_color_by"]
        self.gene_trees_path = parsed["gene_trees_path"]
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        from .vcv_utils import build_vcv_matrix, build_discordance_vcv, parse_gene_trees

        tree = self.read_tree_file()
        self._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_names, traits = self._parse_multi_trait_file(
            self.trait_data_path, tree_tips
        )

        ordered_names = sorted(traits.keys())

        if self.gene_trees_path:
            gene_trees = parse_gene_trees(self.gene_trees_path)
            vcv, vcv_meta = build_discordance_vcv(tree, gene_trees, ordered_names)
            shared = vcv_meta["shared_taxa"]
            if set(shared) != set(ordered_names):
                traits = {k: traits[k] for k in shared}
                ordered_names = shared
        else:
            vcv = build_vcv_matrix(tree, ordered_names)
            vcv_meta = None

        n = len(ordered_names)
        p = len(trait_names)

        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])

        lambda_val = None
        log_likelihood = None

        if self.correction == "lambda":
            max_lam = compute_max_lambda(tree)
            lambda_val, log_likelihood = self._multi_trait_lambda(Y, vcv, max_lam)
            diag_vals = np.diag(vcv).copy()
            vcv = vcv * lambda_val
            np.fill_diagonal(vcv, diag_vals)

        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)

        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])

        Z = Y - np.outer(ones, a_hat)

        if self.method == "pca":
            self._run_pca(
                Z, C_inv, n, p, tree, ordered_names, trait_names, Y,
                lambda_val, log_likelihood, vcv_meta,
            )
        else:  # tsne or umap
            self._run_dimreduce(
                Z, n, tree, ordered_names, trait_names, Y,
                lambda_val, log_likelihood, vcv_meta,
            )

    def _run_pca(
        self, Z, C_inv, n, p, tree, ordered_names, trait_names, Y,
        lambda_val, log_likelihood, vcv_meta=None,
    ) -> None:
        R = (Z.T @ C_inv @ Z) / (n - 1)

        Z_std = None
        if self.mode == "corr":
            D = np.sqrt(np.diag(R))
            D_inv = np.diag(1.0 / D)
            R_corr = D_inv @ R @ D_inv
            eigenvalues, eigenvectors = np.linalg.eigh(R_corr)
        else:
            eigenvalues, eigenvectors = np.linalg.eigh(R)

        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]

        pc_labels = [f"PC{i+1}" for i in range(p)]

        total_var = np.sum(eigenvalues)
        proportions = eigenvalues / total_var

        if self.mode == "corr":
            Z_std = Z @ D_inv
            scores = Z_std @ eigenvectors
        else:
            scores = Z @ eigenvectors

        if self.plot:
            self._plot_pca(
                scores, ordered_names, pc_labels, proportions,
                tree=tree, eigenvectors=eigenvectors, Z=Z, Z_std=Z_std,
                trait_names=trait_names, Y=Y,
            )

        result = self._format_pca_result(
            eigenvalues, proportions, eigenvectors, scores,
            trait_names, ordered_names, pc_labels,
            lambda_val, log_likelihood,
        )

        if self.json_output:
            if self.plot:
                result["plot_output"] = self.plot_output
            if vcv_meta is not None:
                result["vcv_metadata"] = vcv_meta
            print_json(result)
        else:
            self._print_pca_text_output(
                eigenvalues, proportions, eigenvectors, scores,
                trait_names, ordered_names, pc_labels,
                lambda_val, log_likelihood,
            )
            if self.plot:
                print(f"\nSaved PCA plot: {self.plot_output}")

    def _run_dimreduce(
        self, Z, n, tree, ordered_names, trait_names, Y,
        lambda_val, log_likelihood, vcv_meta=None,
    ) -> None:
        if self.method == "tsne":
            embedding, params = self._embed_tsne(Z, n)
        else:
            embedding, params = self._embed_umap(Z, n)

        if self.plot:
            self._plot_dimreduce(
                embedding, ordered_names,
                tree=tree, Z=Z,
                trait_names=trait_names, Y=Y,
            )

        result = self._format_dimreduce_result(
            embedding, ordered_names, params,
            lambda_val, log_likelihood,
        )

        if self.json_output:
            if self.plot:
                result["plot_output"] = self.plot_output
            if vcv_meta is not None:
                result["vcv_metadata"] = vcv_meta
            print_json(result)
        else:
            self._print_dimreduce_text_output(
                embedding, ordered_names, params,
                lambda_val, log_likelihood,
            )
            if self.plot:
                print(f"\nSaved plot: {self.plot_output}")

    def process_args(self, args) -> Dict[str, str]:
        method = getattr(args, "method", "pca")
        plot_tree = getattr(args, "plot_tree", False)
        no_plot_tree = getattr(args, "no_plot_tree", False)
        # For tsne/umap, default to showing the tree unless --no-plot-tree
        if method in ("tsne", "umap") and not no_plot_tree:
            plot_tree = True
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            method=method,
            correction=getattr(args, "correction", "BM"),
            mode=getattr(args, "mode", "cov"),
            n_components=getattr(args, "n_components", 2),
            perplexity=getattr(args, "perplexity", None),
            n_neighbors=getattr(args, "n_neighbors", None),
            min_dist=getattr(args, "min_dist", 0.1),
            seed=getattr(args, "seed", None),
            json_output=getattr(args, "json", False),
            plot=getattr(args, "plot", False),
            plot_output=getattr(args, "plot_output", "phylo_ordination_plot.png"),
            plot_tree=plot_tree,
            color_by=getattr(args, "color_by", None),
            tree_color_by=getattr(args, "tree_color_by", None),
            gene_trees_path=getattr(args, "gene_trees", None),
            plot_config=PlotConfig.from_args(args),
        )

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for phylogenetic ordination."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                raise PhykitUserError(
                    ["All branches in the tree must have lengths."],
                    code=2,
                )

    def _parse_multi_trait_file(
        self, path: str, tree_tips: List[str]
    ) -> Tuple[List[str], Dict[str, List[float]]]:
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

        if len(data_lines) < 2:
            raise PhykitUserError(
                [
                    "Multi-trait file must have a header row and at least one data row.",
                ],
                code=2,
            )

        header_parts = data_lines[0].split("\t")
        n_cols = len(header_parts)
        if n_cols < 2:
            raise PhykitUserError(
                [
                    "Header must have at least 2 columns (taxon + at least 1 trait).",
                ],
                code=2,
            )
        trait_names = header_parts[1:]

        traits = {}
        for line_idx, line in enumerate(data_lines[1:], 2):
            parts = line.split("\t")
            if len(parts) != n_cols:
                raise PhykitUserError(
                    [
                        f"Line {line_idx} has {len(parts)} columns; expected {n_cols}.",
                        f"Each line should have: taxon_name<tab>{'<tab>'.join(['trait'] * len(trait_names))}",
                    ],
                    code=2,
                )
            taxon = parts[0]
            values = []
            for i, val_str in enumerate(parts[1:]):
                try:
                    values.append(float(val_str))
                except ValueError:
                    raise PhykitUserError(
                        [
                            f"Non-numeric trait value '{val_str}' for taxon '{taxon}' "
                            f"(trait '{trait_names[i]}') on line {line_idx}.",
                        ],
                        code=2,
                    )
            traits[taxon] = values

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

        filtered = {taxon: traits[taxon] for taxon in shared}
        return trait_names, filtered

    def _build_vcv_matrix(
        self, tree, ordered_names: List[str]
    ) -> np.ndarray:
        from .vcv_utils import build_vcv_matrix
        return build_vcv_matrix(tree, ordered_names)

    def _multi_trait_log_likelihood(
        self, Y: np.ndarray, C: np.ndarray
    ) -> float:
        n, p = Y.shape
        ones = np.ones(n)

        C_inv = np.linalg.inv(C)

        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])

        Z = Y - np.outer(ones, a_hat)

        R_mle = (Z.T @ C_inv @ Z) / n

        sign_C, logdet_C = np.linalg.slogdet(C)
        sign_R, logdet_R = np.linalg.slogdet(R_mle)

        if sign_C <= 0 or sign_R <= 0:
            return -1e20

        ll = -0.5 * (n * p * np.log(2 * np.pi) + p * logdet_C + n * logdet_R + n * p)
        return float(ll)

    def _multi_trait_lambda(
        self, Y: np.ndarray, vcv: np.ndarray, max_lambda: float
    ) -> Tuple[float, float]:
        diag_vals = np.diag(vcv).copy()
        niter = 10

        def neg_ll(lam):
            C_lam = vcv * lam
            np.fill_diagonal(C_lam, diag_vals)
            try:
                ll = self._multi_trait_log_likelihood(Y, C_lam)
                return -ll
            except (np.linalg.LinAlgError, FloatingPointError, ValueError):
                return 1e10

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

        C_fitted = vcv * lambda_hat
        np.fill_diagonal(C_fitted, diag_vals)
        ll_fitted = self._multi_trait_log_likelihood(Y, C_fitted)

        return float(lambda_hat), float(ll_fitted)

    def _embed_tsne(self, Z: np.ndarray, n: int) -> Tuple[np.ndarray, Dict]:
        from sklearn.manifold import TSNE

        if self.perplexity is not None:
            perplexity = self.perplexity
        else:
            perplexity = min(30.0, (n - 1) / 3.0)

        if perplexity < 1:
            raise PhykitUserError(
                [
                    f"t-SNE requires perplexity >= 1, but auto-computed {perplexity:.2f}.",
                    "Need at least 4 taxa for t-SNE.",
                ],
                code=2,
            )

        tsne = TSNE(
            n_components=self.n_components,
            perplexity=perplexity,
            random_state=self.seed,
        )
        embedding = tsne.fit_transform(Z)

        params = {
            "perplexity": round(perplexity, 2),
        }
        if self.seed is not None:
            params["seed"] = self.seed

        return embedding, params

    def _embed_umap(self, Z: np.ndarray, n: int) -> Tuple[np.ndarray, Dict]:
        import umap

        if self.n_neighbors is not None:
            n_neighbors = self.n_neighbors
        else:
            n_neighbors = min(15, n - 1)

        if n_neighbors < 2:
            raise PhykitUserError(
                [
                    f"UMAP requires n_neighbors >= 2, but auto-computed {n_neighbors}.",
                    "Need at least 3 taxa for UMAP.",
                ],
                code=2,
            )

        reducer = umap.UMAP(
            n_components=self.n_components,
            n_neighbors=n_neighbors,
            min_dist=self.min_dist,
            random_state=self.seed,
        )
        embedding = reducer.fit_transform(Z)

        params = {
            "n_neighbors": n_neighbors,
            "min_dist": self.min_dist,
        }
        if self.seed is not None:
            params["seed"] = self.seed

        return embedding, params

    def _reconstruct_ancestral_scores(
        self, tree, scores: np.ndarray, ordered_names: List[str]
    ) -> Tuple[Dict, Dict]:
        tree_copy = copy.deepcopy(tree)

        tip_names_in_tree = [t.name for t in tree_copy.get_terminals()]
        tips_to_prune = [t for t in tip_names_in_tree if t not in ordered_names]
        if tips_to_prune:
            tree_copy = self.prune_tree_using_taxa_list(tree_copy, tips_to_prune)

        name_to_idx = {name: i for i, name in enumerate(ordered_names)}

        node_estimates = {}
        node_variances = {}
        node_distances = {}

        root = tree_copy.root
        for clade in tree_copy.find_clades(order="postorder"):
            if clade == root:
                node_distances[id(clade)] = 0.0
            else:
                node_distances[id(clade)] = tree_copy.distance(root, clade)

        for clade in tree_copy.find_clades(order="postorder"):
            if clade.is_terminal():
                name = clade.name
                if name in name_to_idx:
                    node_estimates[id(clade)] = scores[name_to_idx[name]]
                    node_variances[id(clade)] = 0.0
            else:
                children = clade.clades
                precisions = []
                weighted_scores = []
                for child in children:
                    child_id = id(child)
                    if child_id not in node_estimates:
                        continue
                    v_i = child.branch_length if child.branch_length else 0.0
                    child_var = node_variances[child_id]
                    denom = v_i + child_var
                    if denom == 0:
                        denom = 1e-10
                    prec = 1.0 / denom
                    precisions.append(prec)
                    weighted_scores.append(prec * node_estimates[child_id])

                if precisions:
                    total_prec = sum(precisions)
                    est = sum(weighted_scores) / total_prec
                    node_estimates[id(clade)] = est
                    node_variances[id(clade)] = 1.0 / total_prec

        return node_estimates, node_distances, tree_copy

    def _parse_color_by(
        self,
        color_by: str,
        trait_names: List[str],
        Y: np.ndarray,
        ordered_names: List[str],
    ) -> Tuple[np.ndarray, List[str], str]:
        if color_by in trait_names:
            col_idx = trait_names.index(color_by)
            return Y[:, col_idx], [], "continuous"

        if not os.path.isfile(color_by):
            raise PhykitUserError(
                [
                    f"--color-by '{color_by}' is not a trait column name "
                    f"({', '.join(trait_names)}) and is not a valid file path.",
                ],
                code=2,
            )

        name_to_idx = {name: i for i, name in enumerate(ordered_names)}
        values = [None] * len(ordered_names)

        with open(color_by) as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                parts = stripped.split("\t")
                if len(parts) < 2:
                    continue
                taxon = parts[0]
                if taxon in name_to_idx:
                    values[name_to_idx[taxon]] = parts[1]

        missing = [ordered_names[i] for i, v in enumerate(values) if v is None]
        if missing:
            raise PhykitUserError(
                [
                    f"--color-by file missing values for taxa: {', '.join(missing)}",
                ],
                code=2,
            )

        is_numeric = True
        for v in values:
            try:
                float(v)
            except (ValueError, TypeError):
                is_numeric = False
                break

        if is_numeric:
            return np.array([float(v) for v in values]), [], "continuous"
        else:
            categories = sorted(set(values))
            return np.array(values), categories, "discrete"

    # --- PCA output ---

    def _format_pca_result(
        self,
        eigenvalues, proportions, eigenvectors, scores,
        trait_names, taxon_names, pc_labels,
        lambda_val, log_likelihood,
    ) -> Dict:
        result = {
            "eigenvalues": {pc_labels[i]: float(eigenvalues[i]) for i in range(len(pc_labels))},
            "proportion_of_variance": {pc_labels[i]: float(proportions[i]) for i in range(len(pc_labels))},
            "loadings": {
                trait_names[j]: {pc_labels[i]: float(eigenvectors[j, i]) for i in range(len(pc_labels))}
                for j in range(len(trait_names))
            },
            "scores": {
                taxon_names[k]: {pc_labels[i]: float(scores[k, i]) for i in range(len(pc_labels))}
                for k in range(len(taxon_names))
            },
        }
        if lambda_val is not None:
            result["lambda"] = float(lambda_val)
            result["log_likelihood"] = float(log_likelihood)
        return result

    def _print_pca_text_output(
        self,
        eigenvalues, proportions, eigenvectors, scores,
        trait_names, taxon_names, pc_labels,
        lambda_val, log_likelihood,
    ) -> None:
        def fmt(v):
            return f"{v:.6f}"

        print("Eigenvalues:")
        print("\t" + "\t".join(pc_labels))
        print("eigenvalue\t" + "\t".join(fmt(v) for v in eigenvalues))
        print("proportion\t" + "\t".join(fmt(v) for v in proportions))

        print("\nLoadings:")
        print("\t" + "\t".join(pc_labels))
        for j, trait in enumerate(trait_names):
            row = "\t".join(fmt(eigenvectors[j, i]) for i in range(len(pc_labels)))
            print(f"{trait}\t{row}")

        print("\nScores:")
        print("\t" + "\t".join(pc_labels))
        for k, taxon in enumerate(taxon_names):
            row = "\t".join(fmt(scores[k, i]) for i in range(len(pc_labels)))
            print(f"{taxon}\t{row}")

        if lambda_val is not None:
            print(f"\nLambda: {fmt(lambda_val)}")
            print(f"Log-likelihood: {fmt(log_likelihood)}")

    # --- Dimreduce output ---

    def _format_dimreduce_result(
        self,
        embedding: np.ndarray,
        taxon_names: List[str],
        params: Dict,
        lambda_val,
        log_likelihood,
    ) -> Dict:
        dim_labels = [f"Dim{i+1}" for i in range(self.n_components)]
        result = {
            "method": self.method,
            "correction": self.correction,
            "n_components": self.n_components,
            "parameters": params,
            "embedding": {
                taxon_names[k]: {
                    dim_labels[i]: round(float(embedding[k, i]), 6)
                    for i in range(self.n_components)
                }
                for k in range(len(taxon_names))
            },
        }
        if lambda_val is not None:
            result["lambda"] = float(lambda_val)
            result["log_likelihood"] = float(log_likelihood)
        return result

    def _print_dimreduce_text_output(
        self,
        embedding: np.ndarray,
        taxon_names: List[str],
        params: Dict,
        lambda_val,
        log_likelihood,
    ) -> None:
        def fmt(v):
            return f"{v:.6f}"

        dim_labels = [f"Dim{i+1}" for i in range(self.n_components)]

        print(f"Method: {self.method}")
        print(f"Correction: {self.correction}")

        if "perplexity" in params:
            print(f"Perplexity: {params['perplexity']}")
        if "n_neighbors" in params:
            print(f"n_neighbors: {params['n_neighbors']}")
        if "min_dist" in params:
            print(f"min_dist: {params['min_dist']}")

        print(f"\nEmbedding:")
        print("\t" + "\t".join(dim_labels))
        for k, taxon in enumerate(taxon_names):
            row = "\t".join(fmt(embedding[k, i]) for i in range(self.n_components))
            print(f"{taxon}\t{row}")

        if lambda_val is not None:
            print(f"\nLambda: {fmt(lambda_val)}")
            print(f"Log-likelihood: {fmt(log_likelihood)}")

    # --- Plotting ---

    def _plot_pca(
        self,
        scores: np.ndarray,
        taxon_names: List[str],
        pc_labels: List[str],
        proportions: np.ndarray,
        tree=None,
        eigenvectors=None,
        Z=None,
        Z_std=None,
        trait_names=None,
        Y=None,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot. Install matplotlib and retry.")
            raise SystemExit(2)

        config = self.plot_config
        config.resolve(n_rows=None, n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        if self.plot_tree and tree is not None and eigenvectors is not None:
            data_for_anc = Z_std if Z_std is not None else Z
            self._draw_tree_overlay(
                ax, fig, tree, data_for_anc @ eigenvectors, taxon_names,
                trait_names, Y,
            )

        self._draw_points(ax, fig, scores, taxon_names, trait_names, Y)

        ax.axhline(0, color="#cccccc", linewidth=0.8, linestyle="--", zorder=1)
        ax.axvline(0, color="#cccccc", linewidth=0.8, linestyle="--", zorder=1)

        ax.set_xlabel(
            f"{pc_labels[0]} ({proportions[0]*100:.1f}% variance explained)"
        )
        ax.set_ylabel(
            f"{pc_labels[1]} ({proportions[1]*100:.1f}% variance explained)"
        )
        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def _plot_dimreduce(
        self,
        embedding: np.ndarray,
        taxon_names: List[str],
        tree=None,
        Z=None,
        trait_names=None,
        Y=None,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot. Install matplotlib and retry.")
            raise SystemExit(2)

        config = self.plot_config
        config.resolve(n_rows=None, n_cols=None)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        if self.plot_tree and tree is not None and Z is not None:
            self._draw_tree_overlay(ax, fig, tree, embedding, taxon_names, trait_names, Y)

        self._draw_points(ax, fig, embedding, taxon_names, trait_names, Y)

        if self.method == "umap":
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel("")
            ax.set_ylabel("")
        else:
            dim_labels = [f"Dim{i+1}" for i in range(self.n_components)]
            ax.set_xlabel(dim_labels[0])
            ax.set_ylabel(dim_labels[1])
        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def _draw_tree_overlay(
        self, ax, fig, tree, scores, taxon_names, trait_names, Y,
    ):
        from matplotlib.collections import LineCollection
        from matplotlib.colors import Normalize
        import matplotlib.pyplot as plt

        node_estimates, node_distances, tree_pruned = \
            self._reconstruct_ancestral_scores(tree, scores, taxon_names)

        # Determine edge coloring strategy
        tree_color_by = self.tree_color_by
        edge_values = None
        edge_label = "Distance from root"
        edge_cmap = "coolwarm"

        if tree_color_by is not None:
            # Reconstruct ancestral values for the chosen trait
            anc_trait_vals, _, _ = self._resolve_tree_color_trait(
                tree_color_by, trait_names, Y, taxon_names, tree,
            )
            if anc_trait_vals is not None:
                edge_values = anc_trait_vals
                edge_label = tree_color_by
                edge_cmap = "viridis"

        all_dists = list(node_distances.values())
        max_dist = max(all_dists) if all_dists else 1.0

        segments = []
        colors = []

        for clade in tree_pruned.find_clades(order="preorder"):
            parent_id = id(clade)
            if parent_id not in node_estimates:
                continue
            parent_scores = node_estimates[parent_id]

            for child in clade.clades:
                child_id = id(child)
                if child_id not in node_estimates:
                    continue
                child_scores = node_estimates[child_id]

                x0, y0 = parent_scores[0], parent_scores[1]
                x1, y1 = child_scores[0], child_scores[1]
                segments.append([(x0, y0), (x1, y1)])

                if edge_values is not None:
                    p_val = edge_values.get(parent_id, 0)
                    c_val = edge_values.get(child_id, 0)
                    colors.append((p_val + c_val) / 2.0)
                else:
                    parent_dist = node_distances.get(parent_id, 0)
                    child_dist = node_distances.get(child_id, 0)
                    colors.append((parent_dist + child_dist) / 2.0)

        if segments:
            if edge_values is not None:
                vmin = min(colors)
                vmax = max(colors)
            else:
                vmin = 0
                vmax = max_dist
            norm = Normalize(vmin=vmin, vmax=vmax)
            cmap = plt.get_cmap(edge_cmap)

            lc = LineCollection(
                segments,
                array=np.array(colors),
                cmap=cmap,
                norm=norm,
                linewidths=1.0,
                alpha=0.7,
                zorder=2,
            )
            ax.add_collection(lc)
            cbar = fig.colorbar(lc, ax=ax, pad=0.02, fraction=0.046)
            cbar.set_label(edge_label)

    def _resolve_tree_color_trait(
        self, tree_color_by, trait_names, Y, taxon_names, tree,
    ):
        """Reconstruct ancestral values for a trait to color tree edges."""
        # Get per-tip values
        if tree_color_by in (trait_names or []):
            col_idx = trait_names.index(tree_color_by)
            tip_vals = Y[:, col_idx]
        elif os.path.isfile(tree_color_by):
            name_to_idx = {name: i for i, name in enumerate(taxon_names)}
            tip_vals = np.zeros(len(taxon_names))
            with open(tree_color_by) as f:
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue
                    parts = stripped.split("\t")
                    if len(parts) < 2:
                        continue
                    taxon = parts[0]
                    if taxon in name_to_idx:
                        try:
                            tip_vals[name_to_idx[taxon]] = float(parts[1])
                        except ValueError:
                            return None, None, None
            # Check all taxa covered
            found = set()
            with open(tree_color_by) as f:
                for line in f:
                    stripped = line.strip()
                    if stripped and not stripped.startswith("#"):
                        parts = stripped.split("\t")
                        if parts[0] in name_to_idx:
                            found.add(parts[0])
            if len(found) < len(taxon_names):
                return None, None, None
        else:
            return None, None, None

        # Reconstruct ancestral trait values via the same post-order method
        # but on the 1D trait column reshaped to (n, 1)
        trait_col = tip_vals.reshape(-1, 1)
        node_estimates, node_distances, tree_pruned = \
            self._reconstruct_ancestral_scores(tree, trait_col, taxon_names)

        # Flatten: node_id -> scalar
        anc_vals = {}
        for node_id, est in node_estimates.items():
            anc_vals[node_id] = float(est[0])

        return anc_vals, node_distances, tree_pruned

    def _draw_points(self, ax, fig, coords, taxon_names, trait_names, Y):
        import matplotlib.pyplot as plt

        color_values = None
        color_categories = None
        color_kind = None
        if self.color_by and trait_names is not None and Y is not None:
            color_values, color_categories, color_kind = self._parse_color_by(
                self.color_by, trait_names, Y, taxon_names
            )

        if color_values is not None and color_kind == "continuous":
            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                s=40,
                alpha=0.8,
                c=color_values,
                cmap="viridis",
                edgecolors="white",
                linewidth=0.5,
                zorder=3,
            )
            cbar = fig.colorbar(scatter, ax=ax, pad=0.02, fraction=0.046)
            cbar.set_label(self.color_by if self.color_by in (trait_names or []) else "color value")
        elif color_values is not None and color_kind == "discrete":
            unique_cats = color_categories
            cat_colors = plt.get_cmap("tab10")
            cat_map = {cat: cat_colors(i / max(len(unique_cats) - 1, 1))
                       for i, cat in enumerate(unique_cats)}
            point_colors = [cat_map[v] for v in color_values]
            ax.scatter(
                coords[:, 0],
                coords[:, 1],
                s=40,
                alpha=0.8,
                c=point_colors,
                edgecolors="white",
                linewidth=0.5,
                zorder=3,
            )
            from matplotlib.lines import Line2D
            handles = [
                Line2D([0], [0], marker='o', color='w', markerfacecolor=cat_map[cat],
                       markersize=7, label=cat)
                for cat in unique_cats
            ]
            ax.legend(handles=handles, title=self.color_by, fontsize=7, title_fontsize=8)
        else:
            ax.scatter(
                coords[:, 0],
                coords[:, 1],
                s=40,
                alpha=0.8,
                color="#2b8cbe",
                edgecolors="white",
                linewidth=0.5,
                zorder=3,
            )

        for k, name in enumerate(taxon_names):
            ax.annotate(
                name,
                (coords[k, 0], coords[k, 1]),
                textcoords="offset points",
                xytext=(5, 5),
                fontsize=8,
            )
