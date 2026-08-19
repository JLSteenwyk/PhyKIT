"""Topology-robust evolutionary-rate covariation by distance projection."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .base import Tree
from .covarying_evolutionary_rates import _pearsonr, _zscore
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


@dataclass(frozen=True)
class ReferenceEdge:
    """One identifiable unrooted edge in the reference-tree basis."""

    split: tuple[str, ...]
    length: float

    @property
    def label(self) -> str:
        return ";".join(self.split)


@dataclass(frozen=True)
class ProjectionResult:
    """Projected edge lengths and goodness-of-fit diagnostics."""

    edge_lengths: np.ndarray
    nrmse: float
    weighted_sse: float
    max_absolute_residual: float


@dataclass(frozen=True)
class ProjectedRateAnalysis:
    """Reusable projected-rate vectors and their supporting diagnostics."""

    reference_tree: object
    shared_taxa: list[str]
    total_pair_count: int
    used_pair_count: int
    edges: list[ReferenceEdge]
    rows: list[dict]
    retained_edge_indices: np.ndarray
    standardized_zero: np.ndarray
    standardized_one: np.ndarray
    correlation: float
    p_value: float
    projection_zero: ProjectionResult
    projection_one: ProjectionResult


class ProjectedCovaryingRates(Tree):
    """Correlate gene rates after projecting distances onto a reference tree."""

    MAX_DESIGN_CELLS = 20_000_000

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        self._processed_args = parsed
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            tree1_file_path=parsed["tree1_file_path"],
            reference=parsed["reference"],
            verbose=parsed["verbose"],
        )
        self.weighting = parsed["weighting"]
        self.max_rate = parsed["max_rate"]
        self.max_pairs = parsed["max_pairs"]
        self.seed = parsed["seed"]
        self.json_output = parsed["json_output"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> dict:
        from ...helpers.plot_config import PlotConfig

        weighting = getattr(args, "weighting", "uniform")
        if weighting not in {"uniform", "inverse_reference"}:
            raise PhykitUserError(
                [
                    "--weighting must be either 'uniform' or "
                    "'inverse_reference'."
                ],
                code=2,
            )

        max_rate = float(getattr(args, "max_rate", 5.0))
        if not math.isfinite(max_rate) or max_rate <= 0.0:
            raise PhykitUserError(
                ["--max-rate must be a finite number greater than zero."],
                code=2,
            )

        max_pairs = int(getattr(args, "max_pairs", 50_000))
        if max_pairs < 1:
            raise PhykitUserError(
                ["--max-pairs must be at least 1."],
                code=2,
            )

        return {
            "tree_file_path": args.tree_zero,
            "tree1_file_path": args.tree_one,
            "reference": args.reference,
            "verbose": getattr(args, "verbose", False),
            "weighting": weighting,
            "max_rate": max_rate,
            "max_pairs": max_pairs,
            "seed": int(getattr(args, "seed", 0)),
            "json_output": getattr(args, "json", False),
            "plot": getattr(args, "plot", False),
            "plot_output": getattr(
                args,
                "plot_output",
                "projected_covarying_rates_plot.png",
            ),
            "plot_config": PlotConfig.from_args(args),
        }

    def run(self) -> None:
        analysis = self._analyze_projected_rates()

        if self.plot:
            self._plot_projected_rates(
                analysis.standardized_zero,
                analysis.standardized_one,
                analysis.correlation,
                analysis.p_value,
            )

        payload = self._result_payload(
            correlation=analysis.correlation,
            p_value=analysis.p_value,
            shared_taxa=analysis.shared_taxa,
            total_pair_count=analysis.total_pair_count,
            used_pair_count=analysis.used_pair_count,
            edges=analysis.edges,
            rows=analysis.rows,
            projection_zero=analysis.projection_zero,
            projection_one=analysis.projection_one,
        )
        if self.json_output:
            print_json(payload)
        else:
            self._print_text(payload)

    def _analyze_projected_rates(self) -> ProjectedRateAnalysis:
        tree_zero = self.read_tree_file_unmodified()
        tree_one = self.read_tree1_file_unmodified()
        tree_ref = self.read_reference_tree_file_unmodified()

        for tree, context in (
            (tree_zero, "projected gene tree zero"),
            (tree_one, "projected gene tree one"),
            (tree_ref, "projected-rate reference tree"),
        ):
            self.validate_tree(
                tree,
                min_tips=3,
                require_branch_lengths=True,
                context=context,
            )

        shared_taxa = self._shared_taxa(tree_zero, tree_one, tree_ref)
        tree_zero = self._prune_to_taxa(tree_zero, shared_taxa)
        tree_one = self._prune_to_taxa(tree_one, shared_taxa)
        tree_ref = self._prune_to_taxa(tree_ref, shared_taxa)

        for tree, context in (
            (tree_zero, "pruned projected gene tree zero"),
            (tree_one, "pruned projected gene tree one"),
            (tree_ref, "pruned projected-rate reference tree"),
        ):
            self.validate_tree(
                tree,
                min_tips=3,
                require_branch_lengths=True,
                context=context,
            )

        edges = self._reference_edge_basis(tree_ref, shared_taxa)
        pair_i, pair_j, selected_pairs, total_pair_count = (
            self._selected_pair_indices(
                len(shared_taxa),
                len(edges),
            )
        )
        design = self._path_design_matrix(shared_taxa, edges, pair_i, pair_j)
        self._validate_projection_design(
            design,
            edge_count=len(edges),
            was_subsampled=len(pair_i) < total_pair_count,
        )

        reference_distances = self._selected_distances(
            tree_ref,
            shared_taxa,
            pair_i,
            pair_j,
            selected_pairs,
        )
        self._validate_reference_basis(design, edges, reference_distances)
        weights = self._pair_weights(reference_distances, self.weighting)

        tree_zero_distances = self._selected_distances(
            tree_zero,
            shared_taxa,
            pair_i,
            pair_j,
            selected_pairs,
        )
        tree_one_distances = self._selected_distances(
            tree_one,
            shared_taxa,
            pair_i,
            pair_j,
            selected_pairs,
        )
        projection_zero = self._project_distances(
            design,
            tree_zero_distances,
            weights,
        )
        projection_one = self._project_distances(
            design,
            tree_one_distances,
            weights,
        )

        rows, retained_zero, retained_one = self._rate_rows(
            edges,
            projection_zero.edge_lengths,
            projection_one.edge_lengths,
        )
        self._validate_rate_vectors(retained_zero, retained_one)
        standardized_zero = _zscore(retained_zero)
        standardized_one = _zscore(retained_one)
        correlation, p_value = _pearsonr(standardized_zero, standardized_one)
        self._add_standardized_rates(rows, standardized_zero, standardized_one)
        retained_edge_indices = np.asarray(
            [
                index
                for index, row in enumerate(rows)
                if row["status"] == "retained"
            ],
            dtype=np.intp,
        )

        return ProjectedRateAnalysis(
            reference_tree=tree_ref,
            shared_taxa=shared_taxa,
            total_pair_count=total_pair_count,
            used_pair_count=len(pair_i),
            edges=edges,
            rows=rows,
            retained_edge_indices=retained_edge_indices,
            standardized_zero=np.asarray(standardized_zero, dtype=float),
            standardized_one=np.asarray(standardized_one, dtype=float),
            correlation=correlation,
            p_value=p_value,
            projection_zero=projection_zero,
            projection_one=projection_one,
        )

    def _shared_taxa(self, tree_zero, tree_one, tree_ref) -> list[str]:
        shared = (
            set(self.get_tip_names_from_tree(tree_zero))
            & set(self.get_tip_names_from_tree(tree_one))
            & set(self.get_tip_names_from_tree(tree_ref))
        )
        if len(shared) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} taxa are shared by both gene trees and "
                    "the reference tree.",
                    "At least 3 shared taxa are required for distance projection.",
                ],
                code=2,
            )
        return sorted(shared)

    def _prune_to_taxa(self, tree, shared_taxa):
        shared = set(shared_taxa)
        tips_to_prune = [
            tip for tip in self.get_tip_names_from_tree(tree) if tip not in shared
        ]
        if not tips_to_prune:
            return tree
        return self.prune_tree_using_taxa_list(
            self._fast_copy(tree),
            tips_to_prune,
        )

    @staticmethod
    def _canonical_split(side, all_taxa) -> tuple[str, ...]:
        complement = all_taxa - side
        side_key = tuple(sorted(side))
        complement_key = tuple(sorted(complement))
        if len(side_key) < len(complement_key):
            return side_key
        if len(complement_key) < len(side_key):
            return complement_key
        return min(side_key, complement_key)

    def _reference_edge_basis(
        self,
        tree,
        shared_taxa,
    ) -> list[ReferenceEdge]:
        all_taxa = frozenset(shared_taxa)
        clade_taxa = {}
        clades = []
        stack = [tree.root]
        while stack:
            clade = stack.pop()
            clades.append(clade)
            stack.extend(clade.clades)

        for clade in reversed(clades):
            if clade.clades:
                taxa = frozenset().union(
                    *(clade_taxa[id(child)] for child in clade.clades)
                )
            else:
                taxa = frozenset([clade.name])
            clade_taxa[id(clade)] = taxa

        lengths_by_split = {}
        for clade in clades:
            if clade is tree.root:
                continue
            split = self._canonical_split(clade_taxa[id(clade)], all_taxa)
            if not split or len(split) == len(all_taxa):
                continue
            lengths_by_split[split] = (
                lengths_by_split.get(split, 0.0) + float(clade.branch_length)
            )

        edges = [
            ReferenceEdge(split=split, length=length)
            for split, length in lengths_by_split.items()
        ]
        edges.sort(key=lambda edge: (len(edge.split), edge.split))
        if len(edges) < 2:
            raise PhykitUserError(
                ["The pruned reference tree has fewer than 2 identifiable edges."],
                code=2,
            )
        return edges

    def _selected_pair_indices(self, taxon_count, edge_count):
        pair_i, pair_j = np.triu_indices(taxon_count, k=1)
        total_pair_count = len(pair_i)
        cell_limited_pairs = self.MAX_DESIGN_CELLS // edge_count
        selected_count = min(
            total_pair_count,
            self.max_pairs,
            cell_limited_pairs,
        )
        if selected_count < edge_count:
            raise PhykitUserError(
                [
                    "The projected-rate design is too large for the configured "
                    "memory limit.",
                    f"At least {edge_count} distance pairs are required, but only "
                    f"{selected_count} fit. Reduce the number of shared taxa.",
                ],
                code=2,
            )
        if selected_count == total_pair_count:
            selected = np.arange(total_pair_count)
            return pair_i, pair_j, selected, total_pair_count

        rng = np.random.default_rng(self.seed)
        selected = np.sort(
            rng.choice(total_pair_count, size=selected_count, replace=False)
        )
        return pair_i[selected], pair_j[selected], selected, total_pair_count

    @staticmethod
    def _path_design_matrix(shared_taxa, edges, pair_i, pair_j):
        taxon_index = {taxon: index for index, taxon in enumerate(shared_taxa)}
        design = np.empty((len(pair_i), len(edges)), dtype=float)
        for column, edge in enumerate(edges):
            membership = np.zeros(len(shared_taxa), dtype=bool)
            membership[[taxon_index[taxon] for taxon in edge.split]] = True
            design[:, column] = membership[pair_i] != membership[pair_j]
        return design

    @staticmethod
    def _validate_projection_design(design, edge_count, was_subsampled) -> None:
        if np.any(np.count_nonzero(design, axis=0) == 0):
            raise PhykitUserError(
                [
                    "The selected taxon pairs do not observe every reference edge.",
                    "Increase --max-pairs or use a different --seed.",
                ],
                code=2,
            )
        if was_subsampled and np.linalg.matrix_rank(design) < edge_count:
            raise PhykitUserError(
                [
                    "The subsampled reference-edge design is rank deficient, so "
                    "projected edge lengths are not uniquely identifiable.",
                    "Increase --max-pairs or use a different --seed.",
                ],
                code=2,
            )

    def _selected_distances(
        self,
        tree,
        shared_taxa,
        pair_i,
        pair_j,
        selected_pairs,
    ):
        total_pair_count = len(shared_taxa) * (len(shared_taxa) - 1) // 2
        if len(selected_pairs) < total_pair_count:
            return np.asarray(
                [
                    tree.distance(shared_taxa[int(i)], shared_taxa[int(j)])
                    for i, j in zip(pair_i, pair_j)
                ],
                dtype=float,
            )

        fast_result = self.calculate_pairwise_tip_distances_fast(
            tree,
            shared_taxa,
            include_combos=False,
        )
        if fast_result is not None:
            return np.asarray(fast_result[1], dtype=float)

        return np.asarray(
            [
                tree.distance(shared_taxa[int(i)], shared_taxa[int(j)])
                for i, j in zip(pair_i, pair_j)
            ],
            dtype=float,
        )

    @staticmethod
    def _validate_reference_basis(design, edges, reference_distances) -> None:
        expected = design @ np.asarray(
            [edge.length for edge in edges],
            dtype=float,
        )
        tolerance = 1e-9 * max(1.0, float(np.max(reference_distances)))
        if not np.allclose(expected, reference_distances, atol=tolerance, rtol=1e-9):
            raise PhykitUserError(
                [
                    "The reference-tree edge basis does not reproduce its "
                    "patristic distances.",
                    "Check the reference tree for unary nodes or malformed topology.",
                ],
                code=2,
            )

    @staticmethod
    def _pair_weights(reference_distances, weighting):
        if weighting == "uniform":
            return np.ones(len(reference_distances), dtype=float)
        positive = reference_distances[reference_distances > 0.0]
        if positive.size == 0:
            raise PhykitUserError(
                ["Reference-tree pairwise distances are all zero."],
                code=2,
            )
        floor = max(float(np.min(positive)) * 1e-8, np.finfo(float).tiny)
        return 1.0 / np.maximum(reference_distances, floor)

    @staticmethod
    def _project_distances(design, distances, weights) -> ProjectionResult:
        if not np.all(np.isfinite(distances)) or np.any(distances < 0.0):
            raise PhykitUserError(
                ["Gene-tree patristic distances must be finite and nonnegative."],
                code=2,
            )
        row_scale = np.sqrt(weights)
        weighted_design = design * row_scale[:, None]
        weighted_distances = distances * row_scale
        try:
            from scipy.optimize import nnls

            edge_lengths, _ = nnls(weighted_design, weighted_distances)
        except RuntimeError:
            from scipy.optimize import lsq_linear

            fitted = lsq_linear(
                weighted_design,
                weighted_distances,
                bounds=(0.0, np.inf),
                method="trf",
                lsmr_tol="auto",
            )
            if not fitted.success:
                raise PhykitUserError(
                    [f"Nonnegative distance projection failed: {fitted.message}"],
                    code=2,
                )
            edge_lengths = fitted.x

        residuals = design @ edge_lengths - distances
        weighted_sse = float(np.dot(weights, residuals * residuals))
        denominator = float(np.dot(weights, distances * distances))
        nrmse = math.sqrt(weighted_sse / denominator) if denominator > 0.0 else 0.0
        return ProjectionResult(
            edge_lengths=np.asarray(edge_lengths, dtype=float),
            nrmse=nrmse,
            weighted_sse=weighted_sse,
            max_absolute_residual=float(np.max(np.abs(residuals))),
        )

    def _rate_rows(self, edges, projected_zero, projected_one):
        rows = []
        retained_zero = []
        retained_one = []
        for edge, length_zero, length_one in zip(
            edges,
            projected_zero,
            projected_one,
        ):
            status = "retained"
            rate_zero = None
            rate_one = None
            if edge.length <= 0.0:
                status = "zero_reference_length"
            else:
                rate_zero = float(length_zero / edge.length)
                rate_one = float(length_one / edge.length)
                if rate_zero > self.max_rate or rate_one > self.max_rate:
                    status = "rate_outlier"

            row = {
                "branch": edge.label,
                "reference_length": float(edge.length),
                "tree_zero_projected_length": float(length_zero),
                "tree_one_projected_length": float(length_one),
                "tree_zero_relative_rate": rate_zero,
                "tree_one_relative_rate": rate_one,
                "tree_zero_zscore": None,
                "tree_one_zscore": None,
                "status": status,
            }
            rows.append(row)
            if status == "retained":
                retained_zero.append(rate_zero)
                retained_one.append(rate_one)
        return rows, retained_zero, retained_one

    @staticmethod
    def _validate_rate_vectors(tree_zero_rates, tree_one_rates) -> None:
        if len(tree_zero_rates) != len(tree_one_rates):
            raise PhykitUserError(
                ["The projected branch-rate vectors have different lengths."],
                code=2,
            )
        if len(tree_zero_rates) < 3:
            raise PhykitUserError(
                [
                    "Fewer than 3 projected reference edges remain after filtering.",
                    "Check taxon overlap, reference lengths, and --max-rate.",
                ],
                code=2,
            )
        if (
            min(tree_zero_rates) == max(tree_zero_rates)
            or min(tree_one_rates) == max(tree_one_rates)
        ):
            raise PhykitUserError(
                [
                    "Pearson correlation is undefined because one gene has no "
                    "variation in projected relative rates."
                ],
                code=2,
            )

    @staticmethod
    def _add_standardized_rates(rows, standardized_zero, standardized_one) -> None:
        retained_index = 0
        for row in rows:
            if row["status"] != "retained":
                continue
            row["tree_zero_zscore"] = float(standardized_zero[retained_index])
            row["tree_one_zscore"] = float(standardized_one[retained_index])
            retained_index += 1

    def _result_payload(
        self,
        *,
        correlation,
        p_value,
        shared_taxa,
        total_pair_count,
        used_pair_count,
        edges,
        rows,
        projection_zero,
        projection_one,
    ):
        retained_edge_count = sum(row["status"] == "retained" for row in rows)
        payload = {
            "method": "projected_covarying_rates",
            "experimental": True,
            "correlation": float(correlation),
            "p_value": float(p_value),
            "shared_taxa": shared_taxa,
            "shared_taxon_count": len(shared_taxa),
            "reference_edge_count": len(edges),
            "retained_edge_count": retained_edge_count,
            "distance_pair_count": total_pair_count,
            "distance_pairs_used": used_pair_count,
            "weighting": self.weighting,
            "max_rate": self.max_rate,
            "tree_zero_projection": {
                "nrmse": projection_zero.nrmse,
                "weighted_sse": projection_zero.weighted_sse,
                "max_absolute_residual": projection_zero.max_absolute_residual,
            },
            "tree_one_projection": {
                "nrmse": projection_one.nrmse,
                "weighted_sse": projection_one.weighted_sse,
                "max_absolute_residual": projection_one.max_absolute_residual,
            },
            "branches": rows,
        }
        if self.plot:
            payload["plot_output"] = self.plot_output
        return payload

    def _print_text(self, payload) -> None:
        lines = [
            "Projected Covarying Evolutionary Rates (experimental)",
            f"Correlation: {payload['correlation']:.6f}",
            f"P-value: {payload['p_value']:.6g}",
            f"Shared taxa: {payload['shared_taxon_count']}",
            "Distance pairs: "
            f"{payload['distance_pairs_used']}/{payload['distance_pair_count']}",
            f"Reference edges: {payload['reference_edge_count']}",
            f"Retained edges: {payload['retained_edge_count']}",
            "Tree zero projection NRMSE: "
            f"{payload['tree_zero_projection']['nrmse']:.6f}",
            "Tree one projection NRMSE: "
            f"{payload['tree_one_projection']['nrmse']:.6f}",
        ]
        if self.verbose:
            lines.extend(
                [
                    "",
                    "branch\treference_length\tprojected_zero\tprojected_one\t"
                    "relative_zero\trelative_one\tz_zero\tz_one\tstatus",
                ]
            )
            for row in payload["branches"]:
                values = [
                    row["branch"],
                    self._format_optional(row["reference_length"]),
                    self._format_optional(row["tree_zero_projected_length"]),
                    self._format_optional(row["tree_one_projected_length"]),
                    self._format_optional(row["tree_zero_relative_rate"]),
                    self._format_optional(row["tree_one_relative_rate"]),
                    self._format_optional(row["tree_zero_zscore"]),
                    self._format_optional(row["tree_one_zscore"]),
                    row["status"],
                ]
                lines.append("\t".join(values))
        if self.plot:
            lines.append(f"Saved projected covarying rates plot: {self.plot_output}")
        try:
            print("\n".join(lines))
        except BrokenPipeError:
            pass

    @staticmethod
    def _format_optional(value) -> str:
        return "NA" if value is None else f"{float(value):.6f}"

    def _plot_projected_rates(
        self,
        tree_zero_rates,
        tree_one_rates,
        correlation,
        p_value,
    ) -> None:
        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            raise PhykitUserError(
                [
                    "matplotlib is required for --plot in "
                    "projected_covarying_rates."
                ],
                code=2,
            )

        x = np.asarray(tree_zero_rates, dtype=float)
        y = np.asarray(tree_one_rates, dtype=float)
        config = self.plot_config
        config.resolve(n_rows=len(x), n_cols=None)
        colors = config.merge_colors(["#2b8cbe", "#202020"])
        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        ax.scatter(x, y, s=18, alpha=0.7, color=colors[0], edgecolors="none")
        slope, intercept = np.polyfit(x, y, 1)
        x_line = np.linspace(float(np.min(x)), float(np.max(x)), 200)
        ax.plot(
            x_line,
            slope * x_line + intercept,
            color=colors[1],
            linestyle="--",
            linewidth=1.8,
        )
        if config.show_title:
            ax.set_title(
                config.title or "Projected Covarying Evolutionary Rates",
                fontsize=config.title_fontsize,
            )
        ax.set_xlabel("Gene tree zero projected rate (z-score)")
        ax.set_ylabel("Gene tree one projected rate (z-score)")
        ax.text(
            0.02,
            0.98,
            f"r={correlation:.4f}\np={p_value:.6g}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "alpha": 0.8,
                "edgecolor": "none",
            },
        )
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)
        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
