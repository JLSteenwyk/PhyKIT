import pickle
import sys
from io import StringIO
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from Bio import Phylo

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig
from ...errors import PhykitUserError


class SpectralDiscordance(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        tree_file = parsed["tree_file_path"]
        if tree_file is not None:
            super().__init__(tree_file_path=tree_file)
        else:
            self.tree_file_path = None
        self.gene_trees_path = parsed["gene_trees_path"]
        self.metric = parsed["metric"]
        self.n_clusters = parsed["n_clusters"]
        self.n_pcs = parsed["n_pcs"]
        self.top_loadings = parsed["top_loadings"]
        self.plot_output = parsed["plot_output"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=getattr(args, "tree", None),
            gene_trees_path=args.gene_trees,
            metric=getattr(args, "metric", "nrf"),
            n_clusters=getattr(args, "clusters", None),
            n_pcs=getattr(args, "n_pcs", None),
            top_loadings=getattr(args, "top_loadings", 5),
            plot_output=getattr(args, "plot", None),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def run(self) -> None:
        gene_trees = self._parse_gene_trees(self.gene_trees_path)

        species_tree = None
        if self.tree_file_path is not None:
            species_tree = self.read_tree_file()

        shared = self._get_shared_taxa(gene_trees, species_tree)

        X, bip_index, sp_flags = self._build_bipartition_matrix(
            gene_trees, shared, metric=self.metric, species_tree=species_tree
        )

        G = X.shape[0]
        if G < 5:
            print(
                f"Warning: only {G} gene trees; clustering may be unreliable.",
                file=sys.stderr,
            )

        scores, var_explained, loadings = self._run_pca(X)
        n_pcs = self.n_pcs if self.n_pcs else min(10, scores.shape[1])
        n_pcs = min(n_pcs, scores.shape[1])

        X_centered = X - X.mean(axis=0)
        labels, K, eigengaps = self._spectral_cluster(
            X_centered, n_clusters=self.n_clusters
        )

        top_loadings = self._get_top_loadings(
            loadings, bip_index, sp_flags, n_pcs, self.top_loadings
        )

        if self.plot_output:
            self._plot_scatter(scores, labels, K, self.plot_output + "_scatter.png")
            self._plot_eigengap(eigengaps, K, self.plot_output + "_eigengap.png")

        if self.json_output:
            result = self._format_json(
                scores, var_explained, top_loadings, labels, K,
                eigengaps, n_pcs, bip_index, sp_flags
            )
            print_json(result)
        else:
            self._print_text(
                scores, var_explained, top_loadings, labels, K,
                eigengaps, n_pcs, G, len(shared), len(bip_index)
            )

    def _bipartition_to_str(self, bip: frozenset) -> str:
        return "{" + ",".join(sorted(bip)) + "}"

    def _get_top_loadings(
        self, loadings, bip_index, sp_flags, n_pcs, top_n,
    ) -> Dict[str, list]:
        result = {}
        for pc in range(min(n_pcs, loadings.shape[0])):
            pc_loadings = loadings[pc]
            top_idx = np.argsort(np.abs(pc_loadings))[::-1][:top_n]
            entries = []
            for idx in top_idx:
                bip = bip_index[idx]
                entries.append({
                    "bipartition": self._bipartition_to_str(bip),
                    "loading": float(pc_loadings[idx]),
                    "in_species_tree": sp_flags.get(bip, False),
                })
            result[f"PC{pc + 1}"] = entries
        return result

    def _format_json(
        self, scores, var_explained, top_loadings, labels, K,
        eigengaps, n_pcs, bip_index, sp_flags,
    ) -> Dict:
        G = scores.shape[0]
        score_dict = {}
        for g in range(G):
            row = {}
            for pc in range(min(n_pcs, scores.shape[1])):
                row[f"PC{pc + 1}"] = float(scores[g, pc])
            row["cluster"] = int(labels[g])
            score_dict[f"gene_tree_{g + 1}"] = row

        return {
            "metric": self.metric,
            "n_gene_trees": G,
            "n_bipartitions": len(bip_index),
            "n_clusters": int(K),
            "scores": score_dict,
            "variance_explained": {
                f"PC{i + 1}": float(var_explained[i])
                for i in range(min(n_pcs, len(var_explained)))
            },
            "top_loadings": top_loadings,
            "cluster_assignments": [int(l) for l in labels],
            "eigengap_values": [float(g) for g in eigengaps],
        }

    def _print_text(
        self, scores, var_explained, top_loadings, labels, K,
        eigengaps, n_pcs, n_trees, n_taxa, n_bips,
    ) -> None:
        print("Spectral Discordance Decomposition")
        print(f"\nMetric: {self.metric}")
        print(f"Gene trees: {n_trees}")
        print(f"Shared taxa: {n_taxa}")
        print(f"Unique bipartitions: {n_bips}")
        print(f"Clusters detected: {K}")

        print("\nVariance explained:")
        for i in range(min(n_pcs, len(var_explained))):
            bar = "#" * int(var_explained[i] * 50)
            print(f"  PC{i + 1:<3d} {var_explained[i]:>8.4f}  {bar}")

        for pc_name, entries in top_loadings.items():
            print(f"\nTop loadings for {pc_name}:")
            for e in entries:
                flag = " *" if e["in_species_tree"] else ""
                print(f"  {e['loading']:>8.4f}  {e['bipartition']}{flag}")

        show_pcs = min(n_pcs, scores.shape[1], 5)
        header = f"  {'Gene tree':<14s}"
        for pc in range(show_pcs):
            header += f"{'PC' + str(pc + 1):>10s}"
        header += f"{'Cluster':>10s}"
        print(f"\nPC scores:")
        print(header)
        for g in range(scores.shape[0]):
            row = f"  {'gene_tree_' + str(g + 1):<14s}"
            for pc in range(show_pcs):
                row += f"{scores[g, pc]:>10.4f}"
            row += f"{labels[g]:>10d}"
            print(row)

    def _plot_scatter(self, scores, labels, K, output_path) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for plotting. Install matplotlib and retry.")
            raise SystemExit(2)

        config = self.plot_config
        config.resolve(n_rows=len(scores), n_cols=None)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        cmap = plt.get_cmap("tab10")
        for k in range(K):
            mask = labels == k
            ax.scatter(
                scores[mask, 0], scores[mask, 1],
                c=[cmap(k)], label=f"Cluster {k + 1}",
                s=60, edgecolors="black", linewidths=0.5, alpha=0.8,
            )
            idxs = np.where(mask)[0]
            for idx in idxs:
                ax.annotate(
                    str(idx + 1), (scores[idx, 0], scores[idx, 1]),
                    fontsize=7, ha="center", va="bottom",
                    xytext=(0, 4), textcoords="offset points",
                )
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        if config.show_title:
            ax.set_title(config.title or "Gene Tree Space — Spectral Discordance", fontsize=config.title_fontsize)
        ax.legend()
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)
        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved scatter plot: {output_path}")

    def _plot_eigengap(self, eigengaps, K, output_path) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for plotting. Install matplotlib and retry.")
            raise SystemExit(2)

        config = self.plot_config
        config.resolve(n_rows=len(eigengaps), n_cols=None)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        n = len(eigengaps)
        x = np.arange(1, n + 1)
        colors = ["tab:red" if i == K - 1 else "tab:blue" for i in range(n)]
        ax.bar(x, eigengaps, color=colors, edgecolor="black", linewidth=0.5)
        ax.set_xlabel("Eigenvalue index (i)")
        ax.set_ylabel(r"Eigengap ($\lambda_{i+1} - \lambda_i$)")
        if config.show_title:
            ax.set_title(config.title or f"Eigengap Heuristic — K = {K} clusters", fontsize=config.title_fontsize)
        ax.axvline(K, color="tab:red", linestyle="--", alpha=0.7, label=f"K = {K}")
        ax.legend()
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)
        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved eigengap plot: {output_path}")

    # ------------------------------------------------------------------
    # Task 3: Gene-tree loading and bipartition extraction
    # ------------------------------------------------------------------

    def _parse_gene_trees(self, path: str) -> list:
        """Read a file of gene trees (one Newick string per line, or
        file-of-filenames where each line is a path to a single-tree file).

        Returns a list of Bio.Phylo tree objects.
        Raises PhykitUserError if the file is not found or contains < 2 trees.
        """
        p = Path(path)
        if not p.exists():
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        trees: list = []
        with open(p) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                # Heuristic: if the line starts with '(' it is inline Newick
                if line.startswith("("):
                    tree = Phylo.read(StringIO(line), "newick")
                    trees.append(tree)
                else:
                    # Treat as a file path (file-of-filenames)
                    tree_path = Path(line)
                    if not tree_path.exists():
                        raise PhykitUserError(
                            [
                                f"{line} corresponds to no such file or directory.",
                                "Please check filename and pathing",
                            ],
                            code=2,
                        )
                    tree = Phylo.read(str(tree_path), "newick")
                    trees.append(tree)

        if len(trees) < 2:
            raise PhykitUserError(
                [
                    "At least 2 gene trees are required for spectral discordance analysis.",
                    f"Only {len(trees)} tree(s) found in {path}.",
                ],
                code=2,
            )

        return trees

    def _get_shared_taxa(self, gene_trees, species_tree=None) -> set:
        """Compute the intersection of tip names across all gene trees
        (and the species tree, if provided).

        Raises PhykitUserError if fewer than 4 taxa are shared.
        """
        if not gene_trees:
            raise PhykitUserError(
                ["No gene trees provided."],
                code=2,
            )

        shared = set(tip.name for tip in gene_trees[0].get_terminals())
        for gt in gene_trees[1:]:
            shared &= set(tip.name for tip in gt.get_terminals())

        if species_tree is not None:
            sp_taxa = set(tip.name for tip in species_tree.get_terminals())
            shared &= sp_taxa

        if len(shared) < 4:
            raise PhykitUserError(
                [
                    "Fewer than 4 shared taxa across gene trees.",
                    f"Only {len(shared)} shared taxa found.",
                    "At least 4 shared taxa are required for meaningful "
                    "bipartition analysis.",
                ],
                code=2,
            )

        return shared

    @staticmethod
    def _canonical_split(taxa_side, all_taxa) -> frozenset:
        """Return the canonical representation of a bipartition split.

        The canonical form is the smaller side of the bipartition.
        Ties (equal-sized sides) are broken lexicographically: the side
        whose sorted tuple is lexicographically smaller is chosen.
        """
        complement = all_taxa - taxa_side
        if len(taxa_side) < len(complement):
            return frozenset(taxa_side)
        elif len(taxa_side) > len(complement):
            return frozenset(complement)
        else:
            # Tie-breaking: pick the lexicographically smaller side
            if sorted(taxa_side) <= sorted(complement):
                return frozenset(taxa_side)
            else:
                return frozenset(complement)

    def _extract_splits(self, tree, all_taxa_fs) -> set:
        """Extract canonical bipartitions from a tree via postorder traversal.

        Skips trivial splits (single-tip on both sides) and the full-taxa split.
        Returns a set of frozensets.
        """
        splits = set()
        # Map terminal clades to their single-taxon sets
        clade_taxa: dict = {}

        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                if clade.name in all_taxa_fs:
                    clade_taxa[id(clade)] = frozenset({clade.name})
                else:
                    clade_taxa[id(clade)] = frozenset()
            else:
                taxa = frozenset()
                for child in clade.clades:
                    taxa = taxa | clade_taxa.get(id(child), frozenset())
                clade_taxa[id(clade)] = taxa

                # Skip polytomous nodes (>2 children = unresolved branching),
                # but allow trifurcating roots (standard unrooted Newick).
                n_children = len(clade.clades)
                if n_children > 2:
                    is_root = (clade == tree.root)
                    if not (is_root and n_children == 3):
                        continue

                # Skip trivial and full-taxa splits
                if len(taxa) <= 1:
                    continue
                complement_size = len(all_taxa_fs) - len(taxa)
                if complement_size <= 0:
                    continue
                # A split where one side has 1 taxon is trivial only if the
                # complement also has 1 taxon (i.e. 2-taxon tree). For larger
                # trees, a split of {A} vs rest is still informative only if
                # rest > 1, but we skip single-tip sides as they are trivial.
                # Actually, the spec says skip trivial (single-tip) splits:
                # that means skip if len(taxa)==1 (already handled above) OR
                # if complement is 1 taxon.
                if complement_size < 1:
                    continue

                canonical = self._canonical_split(taxa, all_taxa_fs)
                # Skip if canonical is single-tip (trivial)
                if len(canonical) <= 1:
                    continue
                # Skip full-taxa split
                if canonical == all_taxa_fs:
                    continue
                splits.add(canonical)

        return splits

    def _extract_splits_with_lengths(self, tree, all_taxa_fs) -> Dict[frozenset, float]:
        """Extract canonical bipartitions with their branch lengths.

        Same traversal as _extract_splits but returns a dict mapping
        canonical split -> branch length (for the wrf metric).
        """
        split_lengths: Dict[frozenset, float] = {}
        clade_taxa: dict = {}

        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                if clade.name in all_taxa_fs:
                    clade_taxa[id(clade)] = frozenset({clade.name})
                else:
                    clade_taxa[id(clade)] = frozenset()
            else:
                taxa = frozenset()
                for child in clade.clades:
                    taxa = taxa | clade_taxa.get(id(child), frozenset())
                clade_taxa[id(clade)] = taxa

                # Skip polytomous nodes (>2 children = unresolved branching),
                # but allow trifurcating roots (standard unrooted Newick).
                n_children = len(clade.clades)
                if n_children > 2:
                    is_root = (clade == tree.root)
                    if not (is_root and n_children == 3):
                        continue

                if len(taxa) <= 1:
                    continue
                complement_size = len(all_taxa_fs) - len(taxa)
                if complement_size <= 0:
                    continue

                canonical = self._canonical_split(taxa, all_taxa_fs)
                if len(canonical) <= 1:
                    continue
                if canonical == all_taxa_fs:
                    continue

                bl = clade.branch_length if clade.branch_length is not None else 0.0
                split_lengths[canonical] = bl

        return split_lengths

    def _build_bipartition_matrix(
        self,
        gene_trees,
        shared_taxa,
        metric="nrf",
        species_tree=None,
    ) -> Tuple[np.ndarray, List[frozenset], Dict[frozenset, bool]]:
        """Build the gene-tree x bipartition matrix.

        Parameters
        ----------
        gene_trees : list
            List of Bio.Phylo tree objects.
        shared_taxa : set
            Set of taxon names shared across all trees.
        metric : str
            "nrf" for binary (presence/absence) or "wrf" for branch-length
            weighted matrix.
        species_tree : Bio.Phylo tree or None
            If provided, bipartitions are flagged as concordant/discordant
            with the species tree.

        Returns
        -------
        X : np.ndarray of shape (G, B)
            The bipartition matrix.
        bip_index : list of frozenset
            Ordered list of bipartitions (columns of X).
        sp_flags : dict
            Maps each bipartition to True if it appears in the species tree,
            False otherwise. Empty dict if no species tree provided.
        """
        all_taxa_fs = frozenset(shared_taxa)

        # Prune gene trees to shared taxa and extract splits
        per_tree_data: list = []
        all_bipartitions: set = set()

        for gt in gene_trees:
            pruned = pickle.loads(pickle.dumps(gt, protocol=pickle.HIGHEST_PROTOCOL))
            tips_to_remove = [
                tip.name
                for tip in pruned.get_terminals()
                if tip.name not in shared_taxa
            ]
            for tip_name in tips_to_remove:
                pruned.prune(tip_name)

            if metric == "wrf":
                split_data = self._extract_splits_with_lengths(pruned, all_taxa_fs)
            else:
                split_data = self._extract_splits(pruned, all_taxa_fs)

            per_tree_data.append(split_data)

            if isinstance(split_data, dict):
                all_bipartitions.update(split_data.keys())
            else:
                all_bipartitions.update(split_data)

        # Stable sort: by (length of split, sorted taxon names)
        bip_index = sorted(
            all_bipartitions,
            key=lambda s: (len(s), sorted(s)),
        )

        # Build matrix
        n_genes = len(gene_trees)
        n_bips = len(bip_index)
        X = np.zeros((n_genes, n_bips), dtype=float)

        for i, tree_data in enumerate(per_tree_data):
            for j, bip in enumerate(bip_index):
                if metric == "wrf":
                    # tree_data is a dict: split -> branch_length
                    if bip in tree_data:
                        X[i, j] = tree_data[bip]
                else:
                    # tree_data is a set of splits
                    if bip in tree_data:
                        X[i, j] = 1.0

        # Species tree flags
        sp_flags: Dict[frozenset, bool] = {}
        if species_tree is not None:
            sp_pruned = pickle.loads(pickle.dumps(species_tree, protocol=pickle.HIGHEST_PROTOCOL))
            sp_tips_to_remove = [
                tip.name
                for tip in sp_pruned.get_terminals()
                if tip.name not in shared_taxa
            ]
            for tip_name in sp_tips_to_remove:
                sp_pruned.prune(tip_name)

            sp_splits = self._extract_splits(sp_pruned, all_taxa_fs)
            for bip in bip_index:
                sp_flags[bip] = bip in sp_splits

        return X, bip_index, sp_flags

    # ------------------------------------------------------------------
    # Task 4: PCA via SVD
    # ------------------------------------------------------------------

    def _run_pca(
        self, X: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """PCA via SVD on centered bipartition matrix.

        Returns (scores, variance_explained, loadings).
        scores: G x n_components
        variance_explained: n_components (sums to 1.0)
        loadings: n_components x B
        """
        X_centered = X - X.mean(axis=0)

        if np.allclose(X_centered, 0):
            n_comp = min(X.shape)
            return (
                np.zeros((X.shape[0], n_comp)),
                np.zeros(n_comp),
                np.zeros((n_comp, X.shape[1])),
            )

        U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
        scores = U * S
        total_var = np.sum(S ** 2)
        var_explained = S ** 2 / total_var if total_var > 0 else np.zeros_like(S)
        loadings = Vt

        return scores, var_explained, loadings

    # ------------------------------------------------------------------
    # Task 5: Spectral clustering with eigengap
    # ------------------------------------------------------------------

    def _spectral_cluster(
        self, X_centered: np.ndarray, n_clusters: int = None,
    ) -> Tuple[np.ndarray, int, np.ndarray]:
        """Spectral clustering with eigengap auto-K detection."""
        from scipy.spatial.distance import pdist, squareform

        G = X_centered.shape[0]

        dists = squareform(pdist(X_centered, metric="euclidean"))

        upper_tri = dists[np.triu_indices(G, k=1)]
        sigma = float(np.median(upper_tri))
        if sigma == 0:
            nonzero = upper_tri[upper_tri > 0]
            sigma = float(np.mean(nonzero)) if len(nonzero) > 0 else 1.0

        W = np.exp(-dists ** 2 / (2 * sigma ** 2))
        np.fill_diagonal(W, W.diagonal() + 1e-10)

        d = np.sum(W, axis=1)
        d_inv_sqrt = 1.0 / np.sqrt(d)
        D_inv_sqrt = np.diag(d_inv_sqrt)
        L_norm = np.eye(G) - D_inv_sqrt @ W @ D_inv_sqrt

        eigenvalues, eigenvectors = np.linalg.eigh(L_norm)

        max_k = min(20, G // 2)
        if max_k < 2:
            max_k = 2
        eigengaps = np.diff(eigenvalues[:max_k + 1])

        if n_clusters is not None:
            K = n_clusters
        else:
            K = int(np.argmax(eigengaps[1:max_k]) + 2)
            K = max(K, 2)

        vecs = eigenvectors[:, :K]
        row_norms = np.linalg.norm(vecs, axis=1, keepdims=True)
        row_norms[row_norms == 0] = 1.0
        vecs = vecs / row_norms

        labels = self._kmeans(vecs, K)

        return labels, K, eigengaps

    @staticmethod
    def _kmeans(X: np.ndarray, K: int, max_iter: int = 300) -> np.ndarray:
        """Simple K-means clustering (avoids sklearn dependency)."""
        G = X.shape[0]
        rng = np.random.RandomState(42)

        centers = [X[rng.randint(G)]]
        for _ in range(1, K):
            dists = np.array([
                min(np.sum((x - c) ** 2) for c in centers)
                for x in X
            ])
            probs = dists / dists.sum() if dists.sum() > 0 else np.ones(G) / G
            centers.append(X[rng.choice(G, p=probs)])
        centers = np.array(centers)

        for _ in range(max_iter):
            dists_to_centers = np.array([
                np.sum((X - c) ** 2, axis=1) for c in centers
            ]).T
            labels = np.argmin(dists_to_centers, axis=1)

            new_centers = np.array([
                X[labels == k].mean(axis=0) if np.any(labels == k) else centers[k]
                for k in range(K)
            ])
            if np.allclose(new_centers, centers):
                break
            centers = new_centers

        return labels
