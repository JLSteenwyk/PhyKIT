from __future__ import annotations

import os

from .base import Tree
from ...errors import PhykitUserError


class _LazyNumpy:
    _module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = _np
            self._module = module

        value = getattr(module, name)
        setattr(self, name, value)
        return value


np = _LazyNumpy()


def _eigenvalue_total_variance(eigenvalues) -> float:
    return float(eigenvalues.sum())


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


def cho_factor(*args, **kwargs):
    from scipy.linalg import cho_factor as _cho_factor

    return _cho_factor(*args, **kwargs)


def cho_solve(*args, **kwargs):
    from scipy.linalg import cho_solve as _cho_solve

    return _cho_solve(*args, **kwargs)


def minimize_scalar(*args, **kwargs):
    from scipy.optimize import minimize_scalar as _minimize_scalar

    return _minimize_scalar(*args, **kwargs)


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
        from .vcv_utils import (
            build_discordance_vcv,
            build_vcv_matrix,
            parse_gene_trees,
        )

        tree = self.read_tree_file_unmodified()
        self.validate_tree(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="phylogenetic ordination",
        )

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(
            self.trait_data_path, tree_tips
        )

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

        n = len(ordered_names)
        p = len(trait_names)

        Y = trait_matrix_from_rows(traits, ordered_names)

        lambda_val = None
        log_likelihood = None

        if self.correction == "lambda":
            from ...helpers.pgls_utils import max_lambda as compute_max_lambda

            max_lam = compute_max_lambda(tree)
            lambda_val, log_likelihood = self._multi_trait_lambda(Y, vcv, max_lam)
            diag_vals = vcv.diagonal().copy()
            vcv = vcv * lambda_val
            np.fill_diagonal(vcv, diag_vals)

        Z, C_inv_Z = self._center_traits_by_vcv(
            Y, vcv, include_weighted=self.method == "pca"
        )

        if self.method == "pca":
            self._run_pca(
                Z, C_inv_Z, n, p, tree, ordered_names, trait_names, Y,
                lambda_val, log_likelihood, vcv_meta,
            )
        else:  # tsne or umap
            self._run_dimreduce(
                Z, n, tree, ordered_names, trait_names, Y,
                lambda_val, log_likelihood, vcv_meta,
            )

    def _run_pca(
        self, Z, C_inv_Z, n, p, tree, ordered_names, trait_names, Y,
        lambda_val, log_likelihood, vcv_meta=None,
    ) -> None:
        R = (Z.T @ C_inv_Z) / (n - 1)

        Z_std = None
        if self.mode == "corr":
            D_inv = 1.0 / np.sqrt(R.diagonal())
            R_corr = (R * D_inv[:, None]) * D_inv[None, :]
            eigenvalues, eigenvectors = np.linalg.eigh(R_corr)
        else:
            eigenvalues, eigenvectors = np.linalg.eigh(R)

        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]

        pc_labels = [f"PC{i+1}" for i in range(p)]

        total_var = _eigenvalue_total_variance(eigenvalues)
        proportions = eigenvalues / total_var

        if self.mode == "corr":
            Z_std = Z * D_inv
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

    def process_args(self, args) -> dict[str, str]:
        from ...helpers.plot_config import PlotConfig

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

    def _build_vcv_matrix(
        self, tree, ordered_names: list[str]
    ) -> np.ndarray:
        from .vcv_utils import build_vcv_matrix
        return build_vcv_matrix(tree, ordered_names)

    def _center_traits_by_vcv(
        self, Y: np.ndarray, C: np.ndarray, include_weighted: bool = False
    ) -> tuple[np.ndarray, np.ndarray]:
        try:
            return self._center_traits_by_vcv_cholesky(Y, C, include_weighted)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._center_traits_by_vcv_inverse(Y, C, include_weighted)

    def _center_traits_by_vcv_cholesky(
        self, Y: np.ndarray, C: np.ndarray, include_weighted: bool = False
    ) -> tuple[np.ndarray, np.ndarray]:
        n, p = Y.shape
        ones = np.ones(n)
        factor = cho_factor(C, lower=True, check_finite=False)

        solve_rhs = np.empty((n, p + 1), dtype=np.result_type(Y, C))
        solve_rhs[:, 0] = ones
        solve_rhs[:, 1:] = Y
        solved = cho_solve(factor, solve_rhs, check_finite=False)
        C_inv_ones = solved[:, 0]
        C_inv_Y = solved[:, 1:]

        denom = ones @ C_inv_ones
        a_hat = (ones @ C_inv_Y) / denom
        Z = Y - a_hat
        if not include_weighted:
            return Z, None

        C_inv_Z = C_inv_Y - C_inv_ones[:, None] * a_hat
        return Z, C_inv_Z

    def _center_traits_by_vcv_inverse(
        self, Y: np.ndarray, C: np.ndarray, include_weighted: bool = False
    ) -> tuple[np.ndarray, np.ndarray]:
        n, p = Y.shape
        ones = np.ones(n)
        C_inv = np.linalg.inv(C)

        C_inv_rhs = C_inv @ np.column_stack((ones, Y))
        C_inv_ones = C_inv_rhs[:, 0]
        C_inv_Y = C_inv_rhs[:, 1:]

        denom = ones @ C_inv_ones
        a_hat = (ones @ C_inv_Y) / denom
        Z = Y - a_hat
        if not include_weighted:
            return Z, None

        C_inv_Z = C_inv_Y - C_inv_ones[:, None] * a_hat
        return Z, C_inv_Z

    def _multi_trait_log_likelihood(
        self, Y: np.ndarray, C: np.ndarray
    ) -> float:
        try:
            return self._multi_trait_log_likelihood_cholesky(Y, C)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._multi_trait_log_likelihood_inverse(Y, C)

    def _multi_trait_log_likelihood_cholesky(
        self, Y: np.ndarray, C: np.ndarray
    ) -> float:
        n, p = Y.shape
        ones = np.ones(n)

        factor = cho_factor(C, lower=True, check_finite=False)
        solve_rhs = np.empty((n, p + 1), dtype=np.result_type(Y, C))
        solve_rhs[:, 0] = 1.0
        solve_rhs[:, 1:] = Y
        solved = cho_solve(factor, solve_rhs, check_finite=False)
        C_inv_ones = solved[:, 0]
        C_inv_Y = solved[:, 1:]

        denom = ones @ C_inv_ones
        a_hat = (ones @ C_inv_Y) / denom

        Z = Y - a_hat
        C_inv_Z = C_inv_Y - C_inv_ones[:, None] * a_hat

        R_mle = (Z.T @ C_inv_Z) / n

        logdet_C = 2.0 * float(np.log(np.diagonal(factor[0])).sum())
        sign_R, logdet_R = np.linalg.slogdet(R_mle)

        if sign_R <= 0:
            return -1e20

        ll = -0.5 * (n * p * np.log(2 * np.pi) + p * logdet_C + n * logdet_R + n * p)
        return float(ll)

    def _multi_trait_log_likelihood_inverse(
        self, Y: np.ndarray, C: np.ndarray
    ) -> float:
        n, p = Y.shape
        ones = np.ones(n)

        C_inv = np.linalg.inv(C)

        C_inv_rhs = C_inv @ np.column_stack((ones, Y))
        C_inv_ones = C_inv_rhs[:, 0]
        C_inv_Y = C_inv_rhs[:, 1:]

        denom = ones @ C_inv_ones
        a_hat = (ones @ C_inv_Y) / denom

        Z = Y - a_hat
        C_inv_Z = C_inv_Y - C_inv_ones[:, None] * a_hat

        R_mle = (Z.T @ C_inv_Z) / n

        sign_C, logdet_C = np.linalg.slogdet(C)
        sign_R, logdet_R = np.linalg.slogdet(R_mle)

        if sign_C <= 0 or sign_R <= 0:
            return -1e20

        ll = -0.5 * (n * p * np.log(2 * np.pi) + p * logdet_C + n * logdet_R + n * p)
        return float(ll)

    def _multi_trait_lambda(
        self, Y: np.ndarray, vcv: np.ndarray, max_lambda: float
    ) -> tuple[float, float]:
        diag_vals = vcv.diagonal().copy()
        diag_step = vcv.shape[0] + 1
        niter = 10

        def neg_ll(lam):
            C_lam = vcv * lam
            C_lam.ravel()[::diag_step] = diag_vals
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
        C_fitted.ravel()[::diag_step] = diag_vals
        ll_fitted = self._multi_trait_log_likelihood(Y, C_fitted)

        return float(lambda_hat), float(ll_fitted)

    def _embed_tsne(self, Z: np.ndarray, n: int) -> tuple[np.ndarray, dict]:
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

    def _embed_umap(self, Z: np.ndarray, n: int) -> tuple[np.ndarray, dict]:
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
        self, tree, scores: np.ndarray, ordered_names: list[str]
    ) -> tuple[dict, dict]:
        tip_names_in_tree = self.get_tip_names_from_tree(tree)
        tips_to_prune = self._tips_to_prune_for_ordered_names(
            tip_names_in_tree,
            ordered_names,
        )
        tree_for_analysis = self._fast_copy(tree) if tips_to_prune else tree
        if tips_to_prune:
            tree_for_analysis = self.prune_tree_using_taxa_list(
                tree_for_analysis, tips_to_prune
            )

        name_to_idx = {name: i for i, name in enumerate(ordered_names)}

        root = tree_for_analysis.root
        try:
            preorder_clades = []
            node_distances = {}
            stack = [(root, 0.0)]
            while stack:
                clade, distance = stack.pop()
                preorder_clades.append(clade)
                node_distances[id(clade)] = distance
                children = clade.clades
                if children:
                    for child in reversed(children):
                        branch_length = (
                            child.branch_length if child.branch_length else 0.0
                        )
                        stack.append((child, distance + branch_length))
        except AttributeError:
            node_distances = {}
            for clade in tree_for_analysis.find_clades(order="postorder"):
                if clade == root:
                    node_distances[id(clade)] = 0.0
                else:
                    node_distances[id(clade)] = tree_for_analysis.distance(root, clade)
            preorder_clades = list(tree_for_analysis.find_clades(order="preorder"))

        node_estimates = {}
        node_variances = {}

        for clade in reversed(preorder_clades):
            children = clade.clades
            if not children:
                name = clade.name
                if name in name_to_idx:
                    node_estimates[id(clade)] = scores[name_to_idx[name]]
                    node_variances[id(clade)] = 0.0
            else:
                total_prec = 0.0
                weighted_score = None
                for child in children:
                    child_id = id(child)
                    child_estimate = node_estimates.get(child_id)
                    if child_estimate is None:
                        continue
                    v_i = child.branch_length if child.branch_length else 0.0
                    child_var = node_variances[child_id]
                    denom = v_i + child_var
                    if denom == 0:
                        denom = 1e-10
                    prec = 1.0 / denom
                    if weighted_score is None:
                        weighted_score = prec * child_estimate
                    else:
                        weighted_score += prec * child_estimate
                    total_prec += prec

                if total_prec:
                    node_estimates[id(clade)] = weighted_score / total_prec
                    node_variances[id(clade)] = 1.0 / total_prec

        return node_estimates, node_distances, tree_for_analysis

    @staticmethod
    def _preorder_clades_direct(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        clades = []
        stack = [root]
        append_clade = clades.append
        append_stack = stack.append
        pop = stack.pop
        try:
            while stack:
                clade = pop()
                append_clade(clade)
                children = clade.clades
                child_count = len(children)
                if child_count == 2:
                    append_stack(children[1])
                    append_stack(children[0])
                else:
                    for index in range(child_count - 1, -1, -1):
                        append_stack(children[index])
        except AttributeError:
            return None
        return clades

    def _parse_color_by(
        self,
        color_by: str,
        trait_names: list[str],
        Y: np.ndarray,
        ordered_names: list[str],
    ) -> tuple[np.ndarray, list[str], str]:
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
                if not stripped or stripped[0] == "#":
                    continue
                taxon, sep, rest = stripped.partition("\t")
                if not sep:
                    continue
                idx = name_to_idx.get(taxon)
                if idx is not None:
                    values[idx] = rest.partition("\t")[0]

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
            return (
                np.fromiter(
                    (float(v) for v in values),
                    dtype=np.float64,
                    count=len(values),
                ),
                [],
                "continuous",
            )
        else:
            categories = sorted(set(values))
            return np.array(values), categories, "discrete"

    # --- PCA output ---

    def _format_pca_result(
        self,
        eigenvalues, proportions, eigenvectors, scores,
        trait_names, taxon_names, pc_labels,
        lambda_val, log_likelihood,
    ) -> dict:
        def format_rows(names, values):
            if len(pc_labels) == 2:
                pc1, pc2 = pc_labels
                col1 = values[:, 0].tolist()
                col2 = values[:, 1].tolist()
                return {
                    name: {
                        pc1: float(val1),
                        pc2: float(val2),
                    }
                    for name, val1, val2 in zip(names, col1, col2)
                }
            if len(pc_labels) == 3:
                pc1, pc2, pc3 = pc_labels
                col1 = values[:, 0].tolist()
                col2 = values[:, 1].tolist()
                col3 = values[:, 2].tolist()
                return {
                    name: {
                        pc1: float(val1),
                        pc2: float(val2),
                        pc3: float(val3),
                    }
                    for name, val1, val2, val3 in zip(
                        names,
                        col1,
                        col2,
                        col3,
                    )
                }
            if len(pc_labels) == 4:
                pc1, pc2, pc3, pc4 = pc_labels
                col1 = values[:, 0].tolist()
                col2 = values[:, 1].tolist()
                col3 = values[:, 2].tolist()
                col4 = values[:, 3].tolist()
                return {
                    name: {
                        pc1: float(val1),
                        pc2: float(val2),
                        pc3: float(val3),
                        pc4: float(val4),
                    }
                    for name, val1, val2, val3, val4 in zip(
                        names,
                        col1,
                        col2,
                        col3,
                        col4,
                    )
                }
            return {
                name: {
                    label: float(value)
                    for label, value in zip(pc_labels, row_values)
                }
                for name, row_values in zip(names, values)
            }

        result = {
            "eigenvalues": {
                label: float(value)
                for label, value in zip(pc_labels, eigenvalues)
            },
            "proportion_of_variance": {
                label: float(value)
                for label, value in zip(pc_labels, proportions)
            },
            "loadings": format_rows(trait_names, eigenvectors),
            "scores": format_rows(taxon_names, scores),
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

        lines = [
            "Eigenvalues:",
            "\t" + "\t".join(pc_labels),
            "eigenvalue\t" + "\t".join(fmt(v) for v in eigenvalues),
            "proportion\t" + "\t".join(fmt(v) for v in proportions),
            "",
            "Loadings:",
            "\t" + "\t".join(pc_labels),
        ]
        append = lines.append

        def append_rows(names, values):
            if len(pc_labels) == 2:
                col1 = values[:, 0].tolist()
                col2 = values[:, 1].tolist()
                for name, val1, val2 in zip(names, col1, col2):
                    append(f"{name}\t{val1:.6f}\t{val2:.6f}")
            elif len(pc_labels) == 3:
                col1 = values[:, 0].tolist()
                col2 = values[:, 1].tolist()
                col3 = values[:, 2].tolist()
                for name, val1, val2, val3 in zip(
                    names,
                    col1,
                    col2,
                    col3,
                ):
                    append(
                        f"{name}\t{val1:.6f}\t{val2:.6f}\t{val3:.6f}"
                    )
            elif len(pc_labels) == 4:
                col1 = values[:, 0].tolist()
                col2 = values[:, 1].tolist()
                col3 = values[:, 2].tolist()
                col4 = values[:, 3].tolist()
                for name, val1, val2, val3, val4 in zip(
                    names,
                    col1,
                    col2,
                    col3,
                    col4,
                ):
                    append(
                        f"{name}\t{val1:.6f}\t{val2:.6f}"
                        f"\t{val3:.6f}\t{val4:.6f}"
                    )
            else:
                for name, row_values in zip(names, values.tolist()):
                    row = "\t".join(fmt(value) for value in row_values)
                    append(f"{name}\t{row}")

        append_rows(trait_names, eigenvectors)

        append("")
        append("Scores:")
        append("\t" + "\t".join(pc_labels))
        append_rows(taxon_names, scores)

        if lambda_val is not None:
            append("")
            append(f"Lambda: {fmt(lambda_val)}")
            append(f"Log-likelihood: {fmt(log_likelihood)}")

        print("\n".join(lines))

    # --- Dimreduce output ---

    def _format_dimreduce_result(
        self,
        embedding: np.ndarray,
        taxon_names: list[str],
        params: dict,
        lambda_val,
        log_likelihood,
    ) -> dict:
        dim_labels = [f"Dim{i+1}" for i in range(self.n_components)]
        if self.n_components == 2:
            embedding_rows = {
                taxon: {
                    "Dim1": round(float(row[0]), 6),
                    "Dim2": round(float(row[1]), 6),
                }
                for taxon, row in zip(taxon_names, embedding)
            }
        elif self.n_components == 3:
            dim1_values = embedding[:, 0].tolist()
            dim2_values = embedding[:, 1].tolist()
            dim3_values = embedding[:, 2].tolist()
            embedding_rows = {
                taxon: {
                    "Dim1": round(float(dim1), 6),
                    "Dim2": round(float(dim2), 6),
                    "Dim3": round(float(dim3), 6),
                }
                for taxon, dim1, dim2, dim3 in zip(
                    taxon_names,
                    dim1_values,
                    dim2_values,
                    dim3_values,
                )
            }
        else:
            embedding_rows = {
                taxon: {
                    label: round(float(value), 6)
                    for label, value in zip(dim_labels, row)
                }
                for taxon, row in zip(taxon_names, embedding.tolist())
            }
        result = {
            "method": self.method,
            "correction": self.correction,
            "n_components": self.n_components,
            "parameters": params,
            "embedding": embedding_rows,
        }
        if lambda_val is not None:
            result["lambda"] = float(lambda_val)
            result["log_likelihood"] = float(log_likelihood)
        return result

    def _print_dimreduce_text_output(
        self,
        embedding: np.ndarray,
        taxon_names: list[str],
        params: dict,
        lambda_val,
        log_likelihood,
    ) -> None:
        def fmt(v):
            return f"{v:.6f}"

        dim_labels = [f"Dim{i+1}" for i in range(self.n_components)]

        lines = [
            f"Method: {self.method}",
            f"Correction: {self.correction}",
        ]
        append = lines.append

        if "perplexity" in params:
            append(f"Perplexity: {params['perplexity']}")
        if "n_neighbors" in params:
            append(f"n_neighbors: {params['n_neighbors']}")
        if "min_dist" in params:
            append(f"min_dist: {params['min_dist']}")

        append("")
        append("Embedding:")
        append("\t" + "\t".join(dim_labels))
        if self.n_components == 2:
            try:
                dim1_values = embedding[:, 0].tolist()
                dim2_values = embedding[:, 1].tolist()
            except (AttributeError, IndexError, TypeError):
                for taxon, values in zip(taxon_names, embedding):
                    append(f"{taxon}\t{values[0]:.6f}\t{values[1]:.6f}")
            else:
                for taxon, dim1, dim2 in zip(taxon_names, dim1_values, dim2_values):
                    append(f"{taxon}\t{dim1:.6f}\t{dim2:.6f}")
        elif self.n_components == 3:
            try:
                dim1_values = embedding[:, 0].tolist()
                dim2_values = embedding[:, 1].tolist()
                dim3_values = embedding[:, 2].tolist()
            except (AttributeError, IndexError, TypeError):
                for taxon, values in zip(taxon_names, embedding):
                    append(
                        f"{taxon}\t{values[0]:.6f}\t"
                        f"{values[1]:.6f}\t{values[2]:.6f}"
                    )
            else:
                for taxon, dim1, dim2, dim3 in zip(
                    taxon_names,
                    dim1_values,
                    dim2_values,
                    dim3_values,
                ):
                    append(f"{taxon}\t{dim1:.6f}\t{dim2:.6f}\t{dim3:.6f}")
        else:
            for taxon, values in zip(taxon_names, embedding.tolist()):
                row = "\t".join(fmt(value) for value in values)
                append(f"{taxon}\t{row}")

        if lambda_val is not None:
            append("")
            append(f"Lambda: {fmt(lambda_val)}")
            append(f"Log-likelihood: {fmt(log_likelihood)}")

        print("\n".join(lines))

    # --- Plotting ---

    def _plot_pca(
        self,
        scores: np.ndarray,
        taxon_names: list[str],
        pc_labels: list[str],
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
        taxon_names: list[str],
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

        max_dist = max(node_distances.values(), default=1.0)

        segments = []
        colors = []

        preorder_clades = self._preorder_clades_direct(tree_pruned)
        if preorder_clades is None:
            preorder_clades = tree_pruned.find_clades(order="preorder")

        for clade in preorder_clades:
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
            found = set()
            with open(tree_color_by) as f:
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped[0] == "#":
                        continue
                    taxon, sep, rest = stripped.partition("\t")
                    if not sep:
                        continue
                    idx = name_to_idx.get(taxon)
                    if idx is not None:
                        try:
                            tip_vals[idx] = float(rest.partition("\t")[0])
                        except ValueError:
                            return None, None, None
                        found.add(taxon)
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
