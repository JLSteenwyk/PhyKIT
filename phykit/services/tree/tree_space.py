import math
import pickle
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from Bio import Phylo

from .base import Tree
from ...errors import PhykitUserError
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig


class TreeSpace(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(trees=parsed["trees"])
        self.output = parsed["output"]
        self.metric = parsed["metric"]
        self.method = parsed["method"]
        self.species_tree_path = parsed["species_tree"]
        self.k = parsed["k"]
        self.seed = parsed["seed"]
        self.json_output = parsed["json_output"]
        self.heatmap = parsed["heatmap"]
        self.distance_matrix_path = parsed["distance_matrix"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> Dict:
        return dict(
            trees=args.trees,
            output=args.output,
            metric=getattr(args, "metric", "rf"),
            method=getattr(args, "method", "mds"),
            species_tree=getattr(args, "species_tree", None),
            k=getattr(args, "k", None),
            seed=getattr(args, "seed", None),
            json_output=getattr(args, "json", False),
            heatmap=getattr(args, "heatmap", False),
            distance_matrix=getattr(args, "distance_matrix", None),
            plot_config=PlotConfig.from_args(args),
        )

    # ------------------------------------------------------------------
    # Tree parsing
    # ------------------------------------------------------------------

    def _parse_trees_from_source(self, trees_path: str) -> List:
        source = Path(trees_path)
        try:
            lines = source.read_text().splitlines()
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{trees_path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        cleaned = [
            line.strip()
            for line in lines
            if line.strip() and not line.strip().startswith("#")
        ]
        if not cleaned:
            raise PhykitUserError(["No trees found in input."], code=2)

        if all(line.startswith("(") for line in cleaned):
            trees = []
            try:
                for line in cleaned:
                    trees.append(Phylo.read(StringIO(line), "newick"))
            except Exception as exc:
                raise PhykitUserError(
                    [f"Failed to parse Newick trees: {exc}"], code=2
                )
            return trees

        trees = []
        for line in cleaned:
            tree_path = Path(line)
            if not tree_path.is_absolute():
                tree_path = source.parent / tree_path
            if not tree_path.exists():
                raise PhykitUserError(
                    [
                        f"{tree_path} corresponds to no such file or directory.",
                        "Please check filename and pathing",
                    ],
                    code=2,
                )
            trees.append(Phylo.read(str(tree_path), "newick"))

        return trees

    # ------------------------------------------------------------------
    # Taxa helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _tips(tree) -> Set[str]:
        return {tip.name for tip in tree.get_terminals()}

    @staticmethod
    def _prune_to_taxa(tree, taxa: Set[str]):
        remove = [
            tip.name for tip in tree.get_terminals() if tip.name not in taxa
        ]
        for tip in remove:
            tree.prune(tip)
        return tree

    def _get_shared_taxa(
        self, gene_trees: List, species_tree=None
    ) -> Set[str]:
        if not gene_trees:
            raise PhykitUserError(["No gene trees provided."], code=2)

        shared = self._tips(gene_trees[0])
        for gt in gene_trees[1:]:
            shared &= self._tips(gt)

        if species_tree is not None:
            shared &= self._tips(species_tree)

        if len(shared) < 4:
            raise PhykitUserError(
                [
                    "Fewer than 4 shared taxa across gene trees.",
                    f"Only {len(shared)} shared taxa found.",
                    "At least 4 shared taxa are required.",
                ],
                code=2,
            )

        return shared

    # ------------------------------------------------------------------
    # Bipartition extraction (used by RF)
    # ------------------------------------------------------------------

    @staticmethod
    def _canonical_split(taxa_side: frozenset, all_taxa: frozenset) -> frozenset:
        complement = all_taxa - taxa_side
        if len(taxa_side) < len(complement):
            return taxa_side
        if len(taxa_side) > len(complement):
            return complement
        if sorted(taxa_side) <= sorted(complement):
            return taxa_side
        return complement

    def _extract_splits(self, tree, all_taxa_fs: frozenset) -> Set[frozenset]:
        splits = set()
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

                n_children = len(clade.clades)
                if n_children > 2:
                    is_root = clade == tree.root
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
                splits.add(canonical)

        return splits

    # ------------------------------------------------------------------
    # Splits with branch lengths (used by KF)
    # ------------------------------------------------------------------

    def _extract_splits_with_lengths(
        self, tree, all_taxa_fs: frozenset
    ) -> Dict[frozenset, float]:
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

                n_children = len(clade.clades)
                if n_children > 2:
                    is_root = clade == tree.root
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

    # ------------------------------------------------------------------
    # Pairwise distance computation
    # ------------------------------------------------------------------

    def _compute_rf_distance(
        self,
        splits_i: Set[frozenset],
        splits_j: Set[frozenset],
        n_taxa: int,
    ) -> float:
        sym_diff = len(splits_i ^ splits_j)
        max_splits = 2 * (n_taxa - 3)
        if max_splits == 0:
            return 0.0
        return sym_diff / max_splits

    def _compute_kf_distance(
        self,
        splits_i: Dict[frozenset, float],
        splits_j: Dict[frozenset, float],
    ) -> float:
        all_splits = set(splits_i.keys()) | set(splits_j.keys())
        kf_squared = 0.0
        for split in all_splits:
            b0 = splits_i.get(split, 0.0)
            b1 = splits_j.get(split, 0.0)
            kf_squared += (b0 - b1) ** 2
        return math.sqrt(kf_squared)

    def _build_distance_matrix(
        self, trees: List, shared_taxa: Set[str], metric: str
    ) -> np.ndarray:
        n = len(trees)
        all_taxa_fs = frozenset(shared_taxa)

        # Prune trees and extract splits
        tree_data = []
        for gt in trees:
            pruned = pickle.loads(pickle.dumps(gt, protocol=pickle.HIGHEST_PROTOCOL))
            tips_to_remove = [
                tip.name
                for tip in pruned.get_terminals()
                if tip.name not in shared_taxa
            ]
            for tip_name in tips_to_remove:
                pruned.prune(tip_name)

            if metric == "kf":
                tree_data.append(
                    self._extract_splits_with_lengths(pruned, all_taxa_fs)
                )
            else:
                tree_data.append(self._extract_splits(pruned, all_taxa_fs))

        n_taxa = len(shared_taxa)
        dist = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                if metric == "kf":
                    d = self._compute_kf_distance(tree_data[i], tree_data[j])
                else:
                    d = self._compute_rf_distance(
                        tree_data[i], tree_data[j], n_taxa
                    )
                dist[i, j] = d
                dist[j, i] = d

        return dist

    # ------------------------------------------------------------------
    # Dimensionality reduction
    # ------------------------------------------------------------------

    def _reduce_dimensions(
        self, dist_matrix: np.ndarray, method: str, seed: Optional[int]
    ) -> np.ndarray:
        if method == "tsne":
            from sklearn.manifold import TSNE

            n_samples = dist_matrix.shape[0]
            perplexity = min(30.0, max(2.0, (n_samples - 1) / 3.0))
            reducer = TSNE(
                n_components=2,
                metric="precomputed",
                random_state=seed,
                perplexity=perplexity,
                init="random",
            )
            return reducer.fit_transform(dist_matrix)
        elif method == "umap":
            try:
                import umap
            except ImportError:
                raise PhykitUserError(
                    [
                        "umap-learn is required for UMAP dimensionality reduction.",
                        "Install it with: pip install umap-learn",
                    ],
                    code=2,
                )
            reducer = umap.UMAP(
                n_components=2, metric="precomputed", random_state=seed
            )
            return reducer.fit_transform(dist_matrix)
        else:
            # MDS (default)
            from sklearn.manifold import MDS

            reducer = MDS(
                n_components=2, dissimilarity="precomputed", random_state=seed
            )
            return reducer.fit_transform(dist_matrix)

    # ------------------------------------------------------------------
    # Spectral clustering with eigengap auto-K
    # ------------------------------------------------------------------

    def _auto_detect_k(self, affinity: np.ndarray) -> int:
        n = affinity.shape[0]
        d = np.sum(affinity, axis=1)
        d_inv_sqrt = 1.0 / np.sqrt(np.maximum(d, 1e-10))
        D_inv_sqrt = np.diag(d_inv_sqrt)
        L_norm = np.eye(n) - D_inv_sqrt @ affinity @ D_inv_sqrt

        eigenvalues = np.sort(np.linalg.eigvalsh(L_norm))

        max_k = min(10, n // 2)
        if max_k < 2:
            max_k = 2
        eigengaps = np.diff(eigenvalues[: max_k + 1])
        # Skip the first eigengap (index 0 -> trivial zero eigenvalue)
        if len(eigengaps) > 1:
            k = int(np.argmax(eigengaps[1:max_k]) + 2)
        else:
            k = 2
        k = max(2, min(k, n // 2))
        return k

    def _spectral_cluster(
        self, dist_matrix: np.ndarray, k: Optional[int], seed: Optional[int]
    ) -> Tuple[np.ndarray, int]:
        from sklearn.cluster import SpectralClustering

        # Convert distance to affinity
        upper_tri = dist_matrix[np.triu_indices(dist_matrix.shape[0], k=1)]
        med = float(np.median(upper_tri))
        if med == 0:
            nonzero = upper_tri[upper_tri > 0]
            med = float(np.mean(nonzero)) if len(nonzero) > 0 else 1.0

        affinity = np.exp(-(dist_matrix ** 2) / (2 * med ** 2))
        np.fill_diagonal(affinity, 1.0)

        if k is None:
            k = self._auto_detect_k(affinity)

        sc = SpectralClustering(
            n_clusters=k,
            affinity="precomputed",
            random_state=seed if seed is not None else 42,
        )
        labels = sc.fit_predict(affinity)

        return labels, k

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def _plot(
        self,
        coords: np.ndarray,
        labels: np.ndarray,
        k: int,
        species_tree_idx: Optional[int],
        output_path: str,
    ) -> None:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        config = self.plot_config
        config.resolve(n_rows=len(coords), n_cols=None)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        cmap = plt.get_cmap("tab10")
        for i in range(k):
            mask = labels == i
            # Exclude species tree point from gene tree clusters
            if species_tree_idx is not None:
                mask_indices = np.where(mask)[0]
                mask_indices = mask_indices[mask_indices != species_tree_idx]
                if len(mask_indices) > 0:
                    ax.scatter(
                        coords[mask_indices, 0],
                        coords[mask_indices, 1],
                        c=[cmap(i)],
                        label=f"Cluster {i + 1}",
                        s=30,
                        alpha=0.7,
                        edgecolors="black",
                        linewidths=0.3,
                    )
            else:
                if mask.any():
                    ax.scatter(
                        coords[mask, 0],
                        coords[mask, 1],
                        c=[cmap(i)],
                        label=f"Cluster {i + 1}",
                        s=30,
                        alpha=0.7,
                        edgecolors="black",
                        linewidths=0.3,
                    )

        if species_tree_idx is not None:
            ax.scatter(
                coords[species_tree_idx, 0],
                coords[species_tree_idx, 1],
                marker="*",
                s=200,
                c="black",
                zorder=10,
                label="Species tree",
                edgecolors="gold",
                linewidths=1,
            )

        ax.set_xlabel("Dim 1")
        ax.set_ylabel("Dim 2")
        ax.legend()

        config.apply_to_figure(
            fig,
            ax,
            f"Tree Space ({self.method.upper()}, {self.metric.upper()})",
            [],
        )

        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def _plot_heatmap(
        self,
        dist_matrix: np.ndarray,
        tree_labels: List[str],
        output_path: str,
    ) -> None:
        """Draw a clustered heatmap of pairwise tree distances."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from scipy.cluster.hierarchy import linkage, dendrogram

        config = self.plot_config
        n = len(tree_labels)
        config.resolve(n_rows=n, n_cols=n)

        # Hierarchical clustering for row/column ordering
        condensed = []
        for i in range(n):
            for j in range(i + 1, n):
                condensed.append(dist_matrix[i, j])
        condensed = np.array(condensed)

        Z = linkage(condensed, method="average")
        dendro = dendrogram(Z, no_plot=True)
        order = dendro["leaves"]

        # Reorder matrix and labels
        reordered = dist_matrix[np.ix_(order, order)]
        reordered_labels = [tree_labels[i] for i in order]

        # Create figure with dendrogram + heatmap
        fig = plt.figure(figsize=(config.fig_width, config.fig_height))

        # Dendrogram on top
        ax_dendro = fig.add_axes([0.10, 0.80, 0.65, 0.15])
        dendrogram(Z, ax=ax_dendro, leaf_rotation=90, leaf_font_size=0,
                   color_threshold=0, above_threshold_color="#333")
        ax_dendro.axis("off")

        # Heatmap
        ax_heat = fig.add_axes([0.10, 0.10, 0.65, 0.68])
        im = ax_heat.imshow(reordered, cmap="viridis_r", aspect="auto",
                            interpolation="nearest")

        # Labels
        label_fontsize = config.ylabel_fontsize if config.ylabel_fontsize and config.ylabel_fontsize > 0 else 6
        ax_heat.set_xticks(range(n))
        ax_heat.set_yticks(range(n))
        ax_heat.set_xticklabels(reordered_labels, rotation=90,
                                fontsize=label_fontsize)
        ax_heat.set_yticklabels(reordered_labels, fontsize=label_fontsize)

        # Colorbar
        ax_cbar = fig.add_axes([0.78, 0.10, 0.03, 0.68])
        cbar = fig.colorbar(im, cax=ax_cbar)
        cbar.set_label(f"{self.metric.upper()} distance", fontsize=9)

        if config.show_title:
            title = config.title or f"Tree Distance Heatmap ({self.metric.upper()})"
            fig.suptitle(title, fontsize=config.title_fontsize, y=0.98)

        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    @staticmethod
    def _write_distance_matrix(
        dist_matrix: np.ndarray,
        tree_labels: List[str],
        output_path: str,
    ) -> None:
        """Write pairwise distance matrix as CSV."""
        with open(output_path, "w") as f:
            f.write("," + ",".join(tree_labels) + "\n")
            for i, label in enumerate(tree_labels):
                row_vals = ",".join(f"{dist_matrix[i, j]:.6f}"
                                    for j in range(len(tree_labels)))
                f.write(f"{label},{row_vals}\n")

    # ------------------------------------------------------------------
    # Main run
    # ------------------------------------------------------------------

    def run(self) -> None:
        # 1. Parse gene trees
        gene_trees = self._parse_trees_from_source(self.trees)
        if len(gene_trees) < 3:
            raise PhykitUserError(
                [
                    "At least 3 gene trees are required for tree space analysis.",
                    f"Only {len(gene_trees)} tree(s) found.",
                ],
                code=2,
            )

        # 2. Optionally load species tree
        species_tree = None
        if self.species_tree_path is not None:
            species_tree = Phylo.read(self.species_tree_path, "newick")

        # 3. Get shared taxa
        shared_taxa = self._get_shared_taxa(gene_trees, species_tree)
        n_taxa = len(shared_taxa)

        # 4. Build list of all trees (gene trees + optional species tree)
        all_trees = list(gene_trees)
        species_tree_idx = None
        if species_tree is not None:
            species_tree_idx = len(all_trees)
            all_trees.append(species_tree)

        # 5. Compute pairwise distance matrix
        dist_matrix = self._build_distance_matrix(
            all_trees, shared_taxa, self.metric
        )

        # 6. Dimensionality reduction
        coords = self._reduce_dimensions(dist_matrix, self.method, self.seed)

        # 7. Spectral clustering (on gene trees only)
        if species_tree_idx is not None:
            gene_dist = dist_matrix[:species_tree_idx, :species_tree_idx]
        else:
            gene_dist = dist_matrix

        labels_gene, k = self._spectral_cluster(gene_dist, self.k, self.seed)

        # Assign species tree to nearest cluster
        species_tree_cluster = None
        if species_tree_idx is not None:
            # Full labels array including species tree
            all_labels = np.zeros(len(all_trees), dtype=int)
            all_labels[:species_tree_idx] = labels_gene
            # Assign species tree to the cluster of its nearest gene tree
            sp_dists = dist_matrix[species_tree_idx, :species_tree_idx]
            nearest_gene = int(np.argmin(sp_dists))
            all_labels[species_tree_idx] = labels_gene[nearest_gene]
            species_tree_cluster = int(labels_gene[nearest_gene]) + 1
        else:
            all_labels = labels_gene

        # 8. Write distance matrix if requested
        n_gene_trees = len(gene_trees)
        if self.distance_matrix_path:
            tree_labels = [f"tree_{i+1}" for i in range(n_gene_trees)]
            if species_tree_idx is not None:
                tree_labels.append("species_tree")
            self._write_distance_matrix(
                dist_matrix, tree_labels, self.distance_matrix_path
            )

        # 9. Plot
        if self.heatmap:
            tree_labels = [f"tree_{i+1}" for i in range(n_gene_trees)]
            if species_tree_idx is not None:
                tree_labels.append("species_tree")
            self._plot_heatmap(dist_matrix, tree_labels, self.output)
        else:
            self._plot(coords, all_labels, k, species_tree_idx, self.output)

        # 10. Build cluster info
        clusters_info = []
        for i in range(k):
            tree_indices = [
                int(idx)
                for idx in range(n_gene_trees)
                if labels_gene[idx] == i
            ]
            clusters_info.append(
                {
                    "id": i + 1,
                    "n_trees": len(tree_indices),
                    "tree_indices": tree_indices,
                }
            )

        # 10. Output
        k_method = "user-specified" if self.k is not None else "auto-detected"
        plot_mode = "heatmap" if self.heatmap else "scatter"

        if self.json_output:
            result = {
                "method": self.method,
                "metric": self.metric,
                "n_trees": n_gene_trees,
                "n_taxa": n_taxa,
                "n_clusters": k,
                "clusters": clusters_info,
                "coordinates": [
                    [round(float(coords[i, 0]), 6), round(float(coords[i, 1]), 6)]
                    for i in range(n_gene_trees)
                ],
                "output_file": self.output,
            }
            if species_tree_cluster is not None:
                result["species_tree_cluster"] = species_tree_cluster
            print_json(result)
        else:
            print("Tree Space Analysis")
            print(f"Method: {self.method.upper()}")
            print(f"Metric: {self.metric.upper()}")
            print(f"Gene trees: {n_gene_trees}")
            print(f"Shared taxa: {n_taxa}")
            print(f"Clusters: {k} ({k_method})")
            for ci in clusters_info:
                print(f"  Cluster {ci['id']}: {ci['n_trees']} trees")
            if species_tree_cluster is not None:
                print(f"Species tree: Cluster {species_tree_cluster}")
            print(f"Plot mode: {plot_mode}")
            print(f"Output: {self.output}")
            if self.distance_matrix_path:
                print(f"Distance matrix: {self.distance_matrix_path}")
