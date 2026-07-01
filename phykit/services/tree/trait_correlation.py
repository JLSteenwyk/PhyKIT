"""
Phylogenetic trait correlation: compute phylogenetic correlations between
all pairs of traits and display them as a heatmap with significance
indicators.
"""
from __future__ import annotations

from .base import Tree
from ...errors import PhykitUserError


_STDTR = None


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


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


def _t_two_tailed_p_values(t_stats: np.ndarray, df: int) -> np.ndarray:
    global _STDTR

    if _STDTR is None:
        from scipy.special import stdtr as _stdtr

        _STDTR = _stdtr

    return 2.0 * _STDTR(df, -np.abs(t_stats))


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
            context="trait correlation analysis",
        )

        tree_tips = self.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(
            self.trait_data_path, tree_tips, min_columns=3
        )

        ordered_names = sorted(traits.keys())

        if self.gene_trees_path:
            gene_trees = parse_gene_trees(self.gene_trees_path)
            vcv, vcv_meta = build_discordance_vcv(tree, gene_trees, ordered_names)
            shared = vcv_meta["shared_taxa"]
            traits, ordered_names = subset_traits_to_ordered_shared_taxa(
                traits, ordered_names, shared
            )
            vcv_type = "discordance"
        else:
            vcv = build_vcv_matrix(tree, ordered_names)
            vcv_meta = None
            vcv_type = "BM"

        n = len(ordered_names)
        p = len(trait_names)

        Y = trait_matrix_from_rows(traits, ordered_names)

        corr_matrix, p_matrix = self._compute_correlation_matrices(Y, vcv)

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

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

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

    def _compute_correlation_matrices(
        self, Y: np.ndarray, vcv: np.ndarray
    ):
        try:
            return self._compute_correlation_matrices_cholesky(Y, vcv)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            return self._compute_correlation_matrices_inverse(Y, vcv)

    def _compute_correlation_matrices_cholesky(
        self, Y: np.ndarray, vcv: np.ndarray
    ):
        n, p = Y.shape
        factor = cho_factor(vcv, lower=True, check_finite=False)
        ones = np.ones(n)
        solve_rhs = np.empty(
            (n, p + 1),
            dtype=np.result_type(Y, vcv),
        )
        solve_rhs[:, 0] = 1.0
        solve_rhs[:, 1:] = Y
        solved = cho_solve(factor, solve_rhs, check_finite=False)
        C_inv_ones = solved[:, 0]
        C_inv_Y = solved[:, 1:]
        denom = ones @ C_inv_ones
        a_hat = (ones @ C_inv_Y) / denom
        Y_centered = Y - a_hat
        C_inv_Y_centered = C_inv_Y - C_inv_ones[:, None] * a_hat
        phylo_cov = (Y_centered.T @ C_inv_Y_centered) / (n - 1)
        return self._correlation_and_p_values(phylo_cov, n, p)

    def _compute_correlation_matrices_inverse(
        self, Y: np.ndarray, vcv: np.ndarray
    ):
        n, p = Y.shape
        C_inv = np.linalg.inv(vcv)
        ones = np.ones(n)
        C_inv_ones = C_inv @ ones
        C_inv_Y = C_inv @ Y
        denom = ones @ C_inv_ones
        a_hat = (ones @ C_inv_Y) / denom
        Y_centered = Y - a_hat  # broadcast subtraction
        C_inv_Y_centered = C_inv_Y - C_inv_ones[:, None] * a_hat
        phylo_cov = (Y_centered.T @ C_inv_Y_centered) / (n - 1)
        return self._correlation_and_p_values(phylo_cov, n, p)

    def _correlation_and_p_values(
        self, phylo_cov: np.ndarray, n: int, p: int
    ):
        std_devs = np.sqrt(phylo_cov.diagonal())
        corr_matrix = phylo_cov / np.outer(std_devs, std_devs)
        np.fill_diagonal(corr_matrix, 1.0)

        p_matrix = np.ones((p, p))
        upper_rows, upper_cols = np.triu_indices(p, k=1)
        upper_r = corr_matrix[upper_rows, upper_cols]
        perfect = np.abs(upper_r) >= 1.0 - 1e-10
        if perfect.any():
            rows = upper_rows[perfect]
            cols = upper_cols[perfect]
            p_matrix[rows, cols] = 0.0
            p_matrix[cols, rows] = 0.0

        testable = ~perfect
        if testable.any():
            rows = upper_rows[testable]
            cols = upper_cols[testable]
            r = upper_r[testable]
            t_stat = r * np.sqrt((n - 2) / (1 - r ** 2))
            p_values = _t_two_tailed_p_values(t_stat, n - 2)
            p_matrix[rows, cols] = p_values
            p_matrix[cols, rows] = p_values

        return corr_matrix, p_matrix

    def _significance_stars(self, pval: float) -> str:
        if pval < 0.001:
            return "***"
        elif pval < 0.01:
            return "**"
        elif pval < self.alpha:
            return "*"
        return ""

    def _draw_significance_stars(self, ax, p_matrix: np.ndarray) -> None:
        if p_matrix.size == 0 or np.min(p_matrix) >= self.alpha:
            return

        flat_indices = np.flatnonzero(p_matrix < self.alpha)
        if flat_indices.size == 0:
            return

        rows, cols = np.unravel_index(flat_indices, p_matrix.shape)
        offdiag = rows != cols
        if not offdiag.any():
            return

        rows = rows[offdiag]
        cols = cols[offdiag]
        p_values = p_matrix[rows, cols]
        star_groups = [
            ("***", p_values < 0.001),
            ("**", (p_values >= 0.001) & (p_values < 0.01)),
            ("*", p_values >= 0.01),
        ]
        for stars, mask in star_groups:
            if mask.any():
                ax.scatter(
                    cols[mask],
                    rows[mask],
                    marker=f"${stars}$",
                    s=85,
                    c="black",
                    linewidths=0,
                    zorder=3,
                )

    def _plot_heatmap(
        self,
        corr_matrix: np.ndarray,
        p_matrix: np.ndarray,
        trait_names: list[str],
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

        self._draw_significance_stars(ax, p_matrix)

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
        # Leave more space on the left for the dendrogram + label gap
        dendro_top_ax = fig.add_axes([0.22, 0.82, 0.58, 0.12])
        dendro_left_ax = fig.add_axes([0.02, 0.10, 0.10, 0.70])
        heat_ax = fig.add_axes([0.22, 0.10, 0.58, 0.70])
        cbar_ax = fig.add_axes([0.83, 0.10, 0.03, 0.70])

        # Draw top dendrogram (black)
        dendrogram(Z, ax=dendro_top_ax, no_labels=True,
                   color_threshold=0, above_threshold_color="#333333")
        dendro_top_ax.axis("off")

        # Draw left dendrogram (black, rotated)
        dendrogram(Z, ax=dendro_left_ax, no_labels=True,
                   color_threshold=0, above_threshold_color="#333333",
                   orientation="left")
        dendro_left_ax.axis("off")

        # Draw heatmap
        im = heat_ax.imshow(
            corr_ordered, cmap="RdBu_r", vmin=-1, vmax=1,
            aspect="auto", interpolation="nearest",
        )

        self._draw_significance_stars(heat_ax, p_ordered)

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
        lines = [
            "Phylogenetic Trait Correlation",
            f"Tree: {tree_path}",
            f"Taxa: {n}",
            f"Traits: {p}",
            f"VCV: {vcv_type} ({'standard' if vcv_type == 'BM' else 'discordance-aware'})",
            f"Alpha: {self.alpha}",
        ]

        # Determine column width
        max_name_len = max(len(name) for name in trait_names)
        col_width = max(max_name_len, 12)
        suffix = "  "

        # Header
        header = (" " * (max_name_len + 2)) + "".join(
            name.rjust(col_width) + suffix for name in trait_names
        )
        lines.append(f"\n{header}")

        # Rows
        alpha = self.alpha
        row_label_width = max_name_len + 2
        for i, name in enumerate(trait_names):
            row_cells = []
            row_cells_append = row_cells.append
            corr_row = corr_matrix[i]
            p_row = p_matrix[i]
            for j, (corr_value, pval) in enumerate(zip(corr_row, p_row)):
                if i == j:
                    cell = "1.0000"
                else:
                    if pval < 0.001:
                        stars = "***"
                    elif pval < 0.01:
                        stars = "**"
                    elif pval < alpha:
                        stars = "*"
                    else:
                        stars = ""
                    cell = f"{corr_value:.4f}{stars}"
                row_cells_append(cell.rjust(col_width))
            lines.append(name.ljust(row_label_width) + suffix.join(row_cells) + suffix)

        n_total = p * (p - 1) // 2
        n_sig = self._count_significant_upper_triangle(p_matrix, self.alpha)
        lines.append(f"\nSignificant pairs (p < {self.alpha}): {n_sig} of {n_total}")
        lines.append(f"Output: {self.output_path}")
        print("\n".join(lines))

    @staticmethod
    def _count_significant_upper_triangle(p_matrix: np.ndarray, alpha: float) -> int:
        n_traits = p_matrix.shape[0]
        count = 0
        for row_idx in range(n_traits - 1):
            count += int(np.count_nonzero(p_matrix[row_idx, row_idx + 1:] < alpha))
        return count

    def _print_json(
        self, n, p, trait_names, vcv_type, corr_matrix, p_matrix
    ) -> None:
        corr_block = corr_matrix[:p, :p]
        p_block = p_matrix[:p, :p]
        sig_pairs = self._build_significant_pairs(
            trait_names, corr_block, p_block, self.alpha
        )

        payload = {
            "n_taxa": n,
            "n_traits": p,
            "trait_names": trait_names,
            "vcv_type": vcv_type,
            "alpha": self.alpha,
            "correlation_matrix": np.round(corr_block, 6).tolist(),
            "p_value_matrix": np.round(p_block, 6).tolist(),
            "significant_pairs": sig_pairs,
            "output_file": self.output_path,
        }
        print_json(payload)

    @staticmethod
    def _build_significant_pairs(
        trait_names, corr_block: np.ndarray, p_block: np.ndarray, alpha: float
    ) -> list[dict]:
        n_traits = len(trait_names)
        total_pairs = n_traits * (n_traits - 1) // 2
        significant_entries = int(np.count_nonzero(p_block < alpha))
        if significant_entries == 0:
            return []

        if significant_entries <= max(64, total_pairs // 20):
            pairs = []
            append = pairs.append
            for i in range(n_traits - 1):
                sig_offsets = np.flatnonzero(p_block[i, i + 1:] < alpha)
                if sig_offsets.size:
                    trait_i = trait_names[i]
                    for offset in sig_offsets:
                        j = i + 1 + int(offset)
                        append(
                            {
                                "trait_i": trait_i,
                                "trait_j": trait_names[j],
                                "r": round(float(corr_block[i, j]), 6),
                                "p": round(float(p_block[i, j]), 6),
                            }
                        )
            return pairs

        upper_rows, upper_cols = np.nonzero(np.triu(p_block < alpha, k=1))
        corr_values = corr_block[upper_rows, upper_cols]
        p_values = p_block[upper_rows, upper_cols]
        return [
            {
                "trait_i": trait_names[int(i)],
                "trait_j": trait_names[int(j)],
                "r": round(float(corr_value), 6),
                "p": round(float(p_value), 6),
            }
            for i, j, corr_value, p_value in zip(
                upper_rows, upper_cols, corr_values, p_values
            )
        ]
