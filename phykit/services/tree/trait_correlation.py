"""
Phylogenetic trait correlation: compute phylogenetic correlations between
all pairs of traits and display them as a heatmap with significance
indicators.
"""
import sys
from typing import Dict, List, Tuple

import numpy as np
from scipy.stats import t as t_dist

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig
from ...errors import PhykitUserError


class TraitCorrelation(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.output_path = parsed["output_path"]
        self.alpha = parsed["alpha"]
        self.cluster = parsed["cluster"]
        self.gene_trees_path = parsed["gene_trees_path"]
        self.json_output = parsed["json_output"]
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
            vcv_type = "discordance"
        else:
            vcv = build_vcv_matrix(tree, ordered_names)
            vcv_meta = None
            vcv_type = "BM"

        n = len(ordered_names)
        p = len(trait_names)

        Y = np.array([[traits[name][j] for j in range(p)] for name in ordered_names])

        # GLS-centered data
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        a_hat = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        Y_centered = Y - a_hat  # broadcast subtraction

        # Phylogenetic covariance matrix
        phylo_cov = (Y_centered.T @ C_inv @ Y_centered) / (n - 1)

        # Convert to correlation matrix
        std_devs = np.sqrt(np.diag(phylo_cov))
        corr_matrix = phylo_cov / np.outer(std_devs, std_devs)
        np.fill_diagonal(corr_matrix, 1.0)

        # Compute p-values
        p_matrix = np.ones((p, p))
        for i in range(p):
            for j in range(i + 1, p):
                r = corr_matrix[i, j]
                if abs(r) >= 1.0 - 1e-10:
                    pval = 0.0
                else:
                    t_stat = r * np.sqrt((n - 2) / (1 - r ** 2))
                    pval = 2.0 * t_dist.sf(abs(t_stat), df=n - 2)
                p_matrix[i, j] = pval
                p_matrix[j, i] = pval

        # Plot heatmap
        self._plot_heatmap(
            corr_matrix, p_matrix, trait_names, self.output_path
        )

        # Output
        if self.json_output:
            self._print_json(
                n, p, trait_names, vcv_type, corr_matrix, p_matrix
            )
        else:
            self._print_text(
                n, p, trait_names, vcv_type, corr_matrix, p_matrix,
                self.trait_data_path,
            )

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            trait_data_path=args.trait_data,
            output_path=args.output,
            alpha=getattr(args, "alpha", 0.05),
            cluster=getattr(args, "cluster", False),
            gene_trees_path=getattr(args, "gene_trees", None),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def _validate_tree(self, tree) -> None:
        tips = list(tree.get_terminals())
        if len(tips) < 3:
            raise PhykitUserError(
                ["Tree must have at least 3 tips for trait correlation analysis."],
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
        if n_cols < 3:
            raise PhykitUserError(
                [
                    "Header must have at least 3 columns (taxon + at least 2 traits).",
                    "Trait correlation requires at least 2 traits.",
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

    def _significance_stars(self, pval: float) -> str:
        if pval < 0.001:
            return "***"
        elif pval < 0.01:
            return "**"
        elif pval < self.alpha:
            return "*"
        return ""

    def _plot_heatmap(
        self,
        corr_matrix: np.ndarray,
        p_matrix: np.ndarray,
        trait_names: List[str],
        output_path: str,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required. Install matplotlib and retry.")
            raise SystemExit(2)

        config = self.plot_config
        p = len(trait_names)
        config.resolve(n_rows=p, n_cols=p)

        if self.cluster and p >= 2:
            self._plot_clustered_heatmap(
                corr_matrix, p_matrix, trait_names, output_path, config, plt
            )
        else:
            self._plot_simple_heatmap(
                corr_matrix, p_matrix, trait_names, output_path, config, plt
            )

    def _plot_simple_heatmap(
        self, corr_matrix, p_matrix, trait_names, output_path, config, plt
    ) -> None:
        p = len(trait_names)
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))

        im = ax.imshow(
            corr_matrix, cmap="RdBu_r", vmin=-1, vmax=1,
            aspect="auto", interpolation="nearest",
        )

        # Add significance stars
        for i in range(p):
            for j in range(p):
                if i != j:
                    stars = self._significance_stars(p_matrix[i, j])
                    if stars:
                        ax.text(
                            j, i, stars,
                            ha="center", va="center",
                            fontsize=8, color="black",
                        )

        ax.set_xticks(np.arange(p))
        ax.set_yticks(np.arange(p))
        xlabel_fs = config.xlabel_fontsize if config.xlabel_fontsize and config.xlabel_fontsize > 0 else 9
        ylabel_fs = config.ylabel_fontsize if config.ylabel_fontsize and config.ylabel_fontsize > 0 else 9
        ax.set_xticklabels(trait_names, rotation=45, ha="right", fontsize=xlabel_fs)
        ax.set_yticklabels(trait_names, fontsize=ylabel_fs)

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Phylogenetic correlation (r)")

        if config.show_title:
            ax.set_title(
                config.title or "Phylogenetic Trait Correlation",
                fontsize=config.title_fontsize,
            )

        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def _plot_clustered_heatmap(
        self, corr_matrix, p_matrix, trait_names, output_path, config, plt
    ) -> None:
        from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
        from scipy.spatial.distance import squareform

        p = len(trait_names)

        # Distance: 1 - |r|, ensure symmetric for squareform
        dist_matrix = 1.0 - np.abs(corr_matrix)
        dist_matrix = (dist_matrix + dist_matrix.T) / 2.0
        np.fill_diagonal(dist_matrix, 0.0)
        condensed = squareform(dist_matrix, checks=False)
        Z = linkage(condensed, method="average")
        leaf_order = list(leaves_list(Z))

        # Reorder matrices
        corr_ordered = corr_matrix[np.ix_(leaf_order, leaf_order)]
        p_ordered = p_matrix[np.ix_(leaf_order, leaf_order)]
        names_ordered = [trait_names[i] for i in leaf_order]

        fig = plt.figure(figsize=(config.fig_width, config.fig_height))

        # Layout: dendrogram on top + left, heatmap center, colorbar right
        # Using add_axes for precise control
        dendro_top_ax = fig.add_axes([0.15, 0.82, 0.65, 0.12])
        dendro_left_ax = fig.add_axes([0.0, 0.10, 0.12, 0.70])
        heat_ax = fig.add_axes([0.15, 0.10, 0.65, 0.70])
        cbar_ax = fig.add_axes([0.83, 0.10, 0.03, 0.70])

        # Draw top dendrogram
        dendrogram(Z, ax=dendro_top_ax, no_labels=True, color_threshold=0)
        dendro_top_ax.axis("off")

        # Draw left dendrogram (rotated)
        dendrogram(Z, ax=dendro_left_ax, no_labels=True, color_threshold=0,
                   orientation="left")
        dendro_left_ax.axis("off")

        # Draw heatmap
        im = heat_ax.imshow(
            corr_ordered, cmap="RdBu_r", vmin=-1, vmax=1,
            aspect="auto", interpolation="nearest",
        )

        # Significance stars
        for i in range(p):
            for j in range(p):
                if i != j:
                    stars = self._significance_stars(p_ordered[i, j])
                    if stars:
                        heat_ax.text(
                            j, i, stars,
                            ha="center", va="center",
                            fontsize=8, color="black",
                        )

        heat_ax.set_xticks(np.arange(p))
        heat_ax.set_yticks(np.arange(p))
        xlabel_fs = config.xlabel_fontsize if config.xlabel_fontsize and config.xlabel_fontsize > 0 else 9
        ylabel_fs = config.ylabel_fontsize if config.ylabel_fontsize and config.ylabel_fontsize > 0 else 9
        heat_ax.set_xticklabels(names_ordered, rotation=45, ha="right", fontsize=xlabel_fs)
        heat_ax.set_yticklabels(names_ordered, fontsize=ylabel_fs)

        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label("Phylogenetic correlation (r)")

        if config.show_title:
            fig.suptitle(
                config.title or "Phylogenetic Trait Correlation",
                fontsize=config.title_fontsize,
                y=0.97,
            )

        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def _print_text(
        self, n, p, trait_names, vcv_type, corr_matrix, p_matrix,
        tree_path,
    ) -> None:
        print("Phylogenetic Trait Correlation")
        print(f"Tree: {tree_path}")
        print(f"Taxa: {n}")
        print(f"Traits: {p}")
        print(f"VCV: {vcv_type} ({'standard' if vcv_type == 'BM' else 'discordance-aware'})")
        print(f"Alpha: {self.alpha}")

        # Determine column width
        max_name_len = max(len(name) for name in trait_names)
        col_width = max(max_name_len, 12)

        # Header
        header = " " * (max_name_len + 2)
        for name in trait_names:
            header += name.rjust(col_width) + "  "
        print(f"\n{header}")

        # Rows
        for i, name in enumerate(trait_names):
            row = name.ljust(max_name_len + 2)
            for j in range(p):
                if i == j:
                    cell = "1.0000"
                else:
                    stars = self._significance_stars(p_matrix[i, j])
                    cell = f"{corr_matrix[i, j]:.4f}{stars}"
                row += cell.rjust(col_width) + "  "
            print(row)

        # Count significant pairs
        n_sig = 0
        n_total = 0
        for i in range(p):
            for j in range(i + 1, p):
                n_total += 1
                if p_matrix[i, j] < self.alpha:
                    n_sig += 1
        print(f"\nSignificant pairs (p < {self.alpha}): {n_sig} of {n_total}")
        print(f"Output: {self.output_path}")

    def _print_json(
        self, n, p, trait_names, vcv_type, corr_matrix, p_matrix
    ) -> None:
        sig_pairs = []
        for i in range(p):
            for j in range(i + 1, p):
                if p_matrix[i, j] < self.alpha:
                    sig_pairs.append({
                        "trait_i": trait_names[i],
                        "trait_j": trait_names[j],
                        "r": round(float(corr_matrix[i, j]), 6),
                        "p": round(float(p_matrix[i, j]), 6),
                    })

        payload = {
            "n_taxa": n,
            "n_traits": p,
            "trait_names": trait_names,
            "vcv_type": vcv_type,
            "alpha": self.alpha,
            "correlation_matrix": [
                [round(float(corr_matrix[i, j]), 6) for j in range(p)]
                for i in range(p)
            ],
            "p_value_matrix": [
                [round(float(p_matrix[i, j]), 6) for j in range(p)]
                for i in range(p)
            ],
            "significant_pairs": sig_pairs,
            "output_file": self.output_path,
        }
        print_json(payload)
