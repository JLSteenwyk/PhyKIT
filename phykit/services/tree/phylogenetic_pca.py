import copy
import os
import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize_scalar

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class PhylogeneticPCA(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.method = parsed["method"]
        self.mode = parsed["mode"]
        self.json_output = parsed["json_output"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]
        self.plot_tree = parsed["plot_tree"]
        self.color_by = parsed["color_by"]

    def run(self) -> None:
        tree = self.read_tree_file()
        self._validate_tree(tree)

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_names, traits = self._parse_multi_trait_file(
            self.trait_data_path, tree_tips
        )

        ordered_names = sorted(traits.keys())
        n = len(ordered_names)
        p = len(trait_names)

        # Build data matrix Y (n x p)
        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])

        vcv = self._build_vcv_matrix(tree, ordered_names)

        lambda_val = None
        log_likelihood = None

        if self.method == "lambda":
            max_lam = self._max_lambda(tree)
            lambda_val, log_likelihood = self._multi_trait_lambda(Y, vcv, max_lam)
            # Transform VCV
            diag_vals = np.diag(vcv).copy()
            vcv = vcv * lambda_val
            np.fill_diagonal(vcv, diag_vals)

        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)

        # GLS mean vector: a_j = (1' C_inv x_j) / (1' C_inv 1)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])

        # GLS-centered data
        Z = Y - np.outer(ones, a_hat)

        # Evolutionary rate matrix
        R = (Z.T @ C_inv @ Z) / (n - 1)

        if self.mode == "corr":
            # Standardize R to correlation matrix
            D = np.sqrt(np.diag(R))
            D_inv = np.diag(1.0 / D)
            R_corr = D_inv @ R @ D_inv
            eigenvalues, eigenvectors = np.linalg.eigh(R_corr)
        else:
            eigenvalues, eigenvectors = np.linalg.eigh(R)

        # Sort descending
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]

        # PC labels
        pc_labels = [f"PC{i+1}" for i in range(p)]

        # Proportion of variance
        total_var = np.sum(eigenvalues)
        proportions = eigenvalues / total_var

        # Compute scores
        Z_std = None
        if self.mode == "corr":
            # Standardize Z columns by sqrt of rate matrix diagonal
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

        # Build result
        result = self._format_result(
            eigenvalues, proportions, eigenvectors, scores,
            trait_names, ordered_names, pc_labels,
            lambda_val, log_likelihood,
        )

        if self.json_output:
            if self.plot:
                result["plot_output"] = self.plot_output
            print_json(result)
        else:
            self._print_text_output(
                eigenvalues, proportions, eigenvectors, scores,
                trait_names, ordered_names, pc_labels,
                lambda_val, log_likelihood,
            )
            if self.plot:
                print(f"\nSaved PCA plot: {self.plot_output}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            method=getattr(args, "method", "BM"),
            mode=getattr(args, "mode", "cov"),
            json_output=getattr(args, "json", False),
            plot=getattr(args, "plot", False),
            plot_output=getattr(args, "plot_output", "phylogenetic_pca_plot.png"),
            plot_tree=getattr(args, "plot_tree", False),
            color_by=getattr(args, "color_by", None),
        )

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for phylogenetic PCA."],
                code=2,
            )
        for clade in tree.find_clades():
            if clade.branch_length is None and clade != tree.root:
                raise PhykitUserError(
                    ["All branches in the tree must have lengths."],
                    code=2,
                )

    def _reconstruct_ancestral_scores(
        self, tree, scores: np.ndarray, ordered_names: List[str]
    ) -> Tuple[Dict, Dict]:
        """ML ancestral reconstruction of PC scores at internal nodes via pruning.

        Returns node_estimates (clade -> score array) and node_distances
        (clade -> distance from root).
        """
        tree_copy = copy.deepcopy(tree)

        # Prune tree to only shared taxa
        tip_names_in_tree = [t.name for t in tree_copy.get_terminals()]
        tips_to_prune = [t for t in tip_names_in_tree if t not in ordered_names]
        if tips_to_prune:
            tree_copy = self.prune_tree_using_taxa_list(tree_copy, tips_to_prune)

        # Map tip names to score arrays
        name_to_idx = {name: i for i, name in enumerate(ordered_names)}

        node_estimates = {}
        node_variances = {}
        node_distances = {}

        # Compute distances from root for all clades
        root = tree_copy.root
        for clade in tree_copy.find_clades(order="postorder"):
            if clade == root:
                node_distances[id(clade)] = 0.0
            else:
                node_distances[id(clade)] = tree_copy.distance(root, clade)

        # Post-order traversal for ancestral reconstruction
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

        # Filter out comments and blank lines
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

        # First data line is the header
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

    def _multi_trait_log_likelihood(
        self, Y: np.ndarray, C: np.ndarray
    ) -> float:
        """Multivariate concentrated log-likelihood for joint lambda estimation.

        logL = -0.5 * (n*p*log(2pi) + p*logdet(C) + n*logdet(R_mle) + n*p)
        where R_mle = Z' C_inv Z / n
        """
        n, p = Y.shape
        ones = np.ones(n)

        C_inv = np.linalg.inv(C)

        # GLS mean vector
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])

        # GLS-centered data
        Z = Y - np.outer(ones, a_hat)

        # MLE rate matrix
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
        """Optimize single lambda jointly across all traits."""
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

        # Compute log-likelihood at fitted lambda
        C_fitted = vcv * lambda_hat
        np.fill_diagonal(C_fitted, diag_vals)
        ll_fitted = self._multi_trait_log_likelihood(Y, C_fitted)

        return float(lambda_hat), float(ll_fitted)

    def _format_result(
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

    def _print_text_output(
        self,
        eigenvalues, proportions, eigenvectors, scores,
        trait_names, taxon_names, pc_labels,
        lambda_val, log_likelihood,
    ) -> None:
        def fmt(v):
            return f"{v:.6f}"

        # Eigenvalues section
        print("Eigenvalues:")
        print("\t" + "\t".join(pc_labels))
        print("eigenvalue\t" + "\t".join(fmt(v) for v in eigenvalues))
        print("proportion\t" + "\t".join(fmt(v) for v in proportions))

        # Loadings section
        print("\nLoadings:")
        print("\t" + "\t".join(pc_labels))
        for j, trait in enumerate(trait_names):
            row = "\t".join(fmt(eigenvectors[j, i]) for i in range(len(pc_labels)))
            print(f"{trait}\t{row}")

        # Scores section
        print("\nScores:")
        print("\t" + "\t".join(pc_labels))
        for k, taxon in enumerate(taxon_names):
            row = "\t".join(fmt(scores[k, i]) for i in range(len(pc_labels)))
            print(f"{taxon}\t{row}")

        # Lambda section (if applicable)
        if lambda_val is not None:
            print(f"\nLambda: {fmt(lambda_val)}")
            print(f"Log-likelihood: {fmt(log_likelihood)}")

    def _parse_color_by(
        self,
        color_by: str,
        trait_names: List[str],
        Y: np.ndarray,
        ordered_names: List[str],
    ) -> Tuple[np.ndarray, List[str], str]:
        """Parse --color-by as column name or file path.

        Returns (values, labels, kind) where kind is 'continuous' or 'discrete'.
        values is an array parallel to ordered_names.
        labels is the list of category names for discrete, or empty for continuous.
        """
        # Check if color_by matches a column name
        if color_by in trait_names:
            col_idx = trait_names.index(color_by)
            return Y[:, col_idx], [], "continuous"

        # Try to open as a file
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

        # Check for missing taxa
        missing = [ordered_names[i] for i, v in enumerate(values) if v is None]
        if missing:
            raise PhykitUserError(
                [
                    f"--color-by file missing values for taxa: {', '.join(missing)}",
                ],
                code=2,
            )

        # Determine if continuous or discrete
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
            from matplotlib.collections import LineCollection
            from matplotlib.colors import Normalize
            import matplotlib.cm as cm
        except ImportError:
            print("matplotlib is required for --plot. Install matplotlib and retry.")
            raise SystemExit(2)

        fig, ax = plt.subplots(figsize=(7, 5))

        # Draw phylomorphospace edges if requested
        if self.plot_tree and tree is not None and eigenvectors is not None:
            # Project centered data into PC space for ancestral reconstruction
            # Use first 2 PCs
            data_for_anc = Z_std if Z_std is not None else Z
            node_estimates, node_distances, tree_pruned = \
                self._reconstruct_ancestral_scores(
                    tree, data_for_anc @ eigenvectors, taxon_names
                )

            # Collect all distances for normalization
            all_dists = [d for d in node_distances.values()]
            max_dist = max(all_dists) if all_dists else 1.0

            # Draw edges colored by time from root
            segments = []
            colors = []
            norm = Normalize(vmin=0, vmax=max_dist)
            cmap = plt.get_cmap("coolwarm")

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

                    # Color by midpoint distance from root
                    parent_dist = node_distances.get(parent_id, 0)
                    child_dist = node_distances.get(child_id, 0)
                    mid_dist = (parent_dist + child_dist) / 2.0
                    colors.append(mid_dist)

            if segments:
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
                cbar.set_label("Distance from root")

        # Determine tip point colors
        color_values = None
        color_categories = None
        color_kind = None
        if self.color_by and trait_names is not None and Y is not None:
            color_values, color_categories, color_kind = self._parse_color_by(
                self.color_by, trait_names, Y, taxon_names
            )

        if color_values is not None and color_kind == "continuous":
            scatter = ax.scatter(
                scores[:, 0],
                scores[:, 1],
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
                scores[:, 0],
                scores[:, 1],
                s=40,
                alpha=0.8,
                c=point_colors,
                edgecolors="white",
                linewidth=0.5,
                zorder=3,
            )
            # Add legend
            from matplotlib.lines import Line2D
            handles = [
                Line2D([0], [0], marker='o', color='w', markerfacecolor=cat_map[cat],
                       markersize=7, label=cat)
                for cat in unique_cats
            ]
            ax.legend(handles=handles, title=self.color_by, fontsize=7, title_fontsize=8)
        else:
            ax.scatter(
                scores[:, 0],
                scores[:, 1],
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
                (scores[k, 0], scores[k, 1]),
                textcoords="offset points",
                xytext=(5, 5),
                fontsize=8,
            )

        ax.axhline(0, color="#cccccc", linewidth=0.8, linestyle="--", zorder=1)
        ax.axvline(0, color="#cccccc", linewidth=0.8, linestyle="--", zorder=1)

        ax.set_xlabel(
            f"{pc_labels[0]} ({proportions[0]*100:.1f}% variance explained)"
        )
        ax.set_ylabel(
            f"{pc_labels[1]} ({proportions[1]*100:.1f}% variance explained)"
        )
        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=300, bbox_inches="tight")
        plt.close(fig)
