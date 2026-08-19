"""Localize episodic evolutionary-rate covariation on a reference tree."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .projected_covarying_rates import (
    ProjectedCovaryingRates,
    ProjectedRateAnalysis,
    print_json,
)
from ...errors import PhykitUserError


@dataclass(frozen=True)
class CladeScanCandidate:
    """One rooted reference clade and its local projected-edge positions."""

    label: str
    taxa: tuple[str, ...]
    edge_indices: tuple[int, ...]
    root_depth: int


@dataclass(frozen=True)
class CladeScanBasis:
    """Scannable projected edges and rooted clade candidates."""

    rate_positions: np.ndarray
    original_edge_indices: np.ndarray
    edge_depths: np.ndarray
    edge_labels: tuple[str, ...]
    candidates: tuple[CladeScanCandidate, ...]
    ambiguous_edge_count: int


class EpisodicRateCovariation(ProjectedCovaryingRates):
    """Scan a reference tree for localized projected-rate covariation."""

    PERMUTATION_BATCH_SIZE = 128

    def __init__(self, args) -> None:
        super().__init__(args)
        parsed = self._processed_args
        self.permutations = parsed["permutations"]
        self.depth_bins = parsed["depth_bins"]
        self.min_edges = parsed["min_edges"]
        self.alpha = parsed["alpha"]
        self.output = parsed["output"]
        self.annotated_tree = parsed["annotated_tree"]

    def process_args(self, args) -> dict:
        parsed = super().process_args(args)
        permutations = int(getattr(args, "permutations", 999))
        if permutations < 1:
            raise PhykitUserError(
                ["--permutations must be at least 1."],
                code=2,
            )

        depth_bins = int(getattr(args, "depth_bins", 4))
        if depth_bins < 1:
            raise PhykitUserError(
                ["--depth-bins must be at least 1."],
                code=2,
            )

        min_edges = int(getattr(args, "min_edges", 3))
        if min_edges < 2:
            raise PhykitUserError(
                ["--min-edges must be at least 2."],
                code=2,
            )

        alpha = float(getattr(args, "alpha", 0.05))
        if not math.isfinite(alpha) or not 0.0 < alpha <= 1.0:
            raise PhykitUserError(
                ["--alpha must be greater than zero and at most one."],
                code=2,
            )

        parsed.update(
            permutations=permutations,
            depth_bins=depth_bins,
            min_edges=min_edges,
            alpha=alpha,
            output=getattr(args, "output", None),
            annotated_tree=getattr(args, "annotated_tree", None),
        )
        return parsed

    def run(self) -> None:
        analysis = self._analyze_projected_rates()
        edge_influence = self._add_edge_diagnostics(analysis)
        scan_basis = self._clade_scan_basis(
            analysis.reference_tree,
            analysis.edges,
            analysis.retained_edge_indices,
            self.min_edges,
        )
        scan_zero = analysis.standardized_zero[scan_basis.rate_positions]
        scan_one = analysis.standardized_one[scan_basis.rate_positions]
        membership = self._candidate_membership_matrix(
            scan_basis.candidates,
            len(scan_zero),
        )
        strata = self._depth_strata(scan_basis.edge_depths, self.depth_bins)
        scan_rows = self._scan_clades(
            scan_zero,
            scan_one,
            membership,
            scan_basis.candidates,
            strata,
        )

        if self.output:
            self._write_scan_table(scan_rows, self.output)
        if self.annotated_tree:
            self._write_annotated_tree(
                analysis.reference_tree,
                scan_rows,
                self.annotated_tree,
            )

        payload = self._episodic_payload(
            analysis,
            scan_basis,
            strata,
            scan_rows,
            edge_influence,
        )
        if self.json_output:
            print_json(payload)
        else:
            self._print_episodic_text(payload)

    def _clade_scan_basis(
        self,
        tree,
        edges,
        retained_edge_indices,
        min_edges,
    ) -> CladeScanBasis:
        root = tree.root
        clades = []
        depths = {id(root): 0}
        stack = [root]
        while stack:
            clade = stack.pop()
            clades.append(clade)
            child_depth = depths[id(clade)] + 1
            for child in clade.clades:
                depths[id(child)] = child_depth
                stack.append(child)

        clade_taxa = {}
        for clade in reversed(clades):
            if clade.clades:
                taxa = frozenset().union(
                    *(clade_taxa[id(child)] for child in clade.clades)
                )
            else:
                taxa = frozenset((clade.name,))
            clade_taxa[id(clade)] = taxa

        all_taxa = clade_taxa[id(root)]
        occurrences = {}
        for clade in clades:
            if clade is root:
                continue
            split = self._canonical_split(clade_taxa[id(clade)], all_taxa)
            occurrences.setdefault(split, []).append(clade)

        ambiguous_splits = {
            split for split, matches in occurrences.items() if len(matches) > 1
        }
        retained_positions = {
            int(edge_index): position
            for position, edge_index in enumerate(retained_edge_indices)
        }
        scannable_original = [
            int(edge_index)
            for edge_index in retained_edge_indices
            if edges[int(edge_index)].split not in ambiguous_splits
        ]
        scan_position_by_original = {
            edge_index: position
            for position, edge_index in enumerate(scannable_original)
        }
        scan_position_by_split = {
            edges[edge_index].split: scan_position
            for edge_index, scan_position in scan_position_by_original.items()
        }

        edge_depths = np.empty(len(scannable_original), dtype=float)
        for split, scan_position in scan_position_by_split.items():
            edge_depths[scan_position] = depths[id(occurrences[split][0])]

        subtree_edges = {}
        candidates = []
        seen_edge_sets = set()
        scan_edge_count = len(scannable_original)
        for clade in reversed(clades):
            local_edges = set()
            for child in clade.clades:
                local_edges.update(subtree_edges[id(child)])
            if clade is not root:
                split = self._canonical_split(clade_taxa[id(clade)], all_taxa)
                scan_position = scan_position_by_split.get(split)
                if scan_position is not None:
                    local_edges.add(scan_position)
            subtree_edges[id(clade)] = frozenset(local_edges)

            if clade is root or not clade.clades:
                continue
            edge_indices = tuple(sorted(local_edges))
            if (
                len(edge_indices) < min_edges
                or scan_edge_count - len(edge_indices) < min_edges
                or edge_indices in seen_edge_sets
            ):
                continue
            seen_edge_sets.add(edge_indices)
            taxa = tuple(sorted(clade_taxa[id(clade)]))
            candidates.append(
                CladeScanCandidate(
                    label=";".join(taxa),
                    taxa=taxa,
                    edge_indices=edge_indices,
                    root_depth=depths[id(clade)],
                )
            )

        if scan_edge_count < min_edges * 2:
            raise PhykitUserError(
                [
                    f"Only {scan_edge_count} localizable projected edges remain.",
                    f"At least {min_edges * 2} are required with "
                    f"--min-edges {min_edges}.",
                ],
                code=2,
            )
        if not candidates:
            raise PhykitUserError(
                [
                    "No rooted reference clades satisfy the episodic scan "
                    "requirements.",
                    "Reduce --min-edges or use a reference tree with more taxa.",
                ],
                code=2,
            )

        candidates.sort(key=lambda candidate: (candidate.root_depth, candidate.taxa))
        return CladeScanBasis(
            rate_positions=np.asarray(
                [retained_positions[index] for index in scannable_original],
                dtype=np.intp,
            ),
            original_edge_indices=np.asarray(
                scannable_original,
                dtype=np.intp,
            ),
            edge_depths=edge_depths,
            edge_labels=tuple(edges[index].label for index in scannable_original),
            candidates=tuple(candidates),
            ambiguous_edge_count=sum(
                edges[int(index)].split in ambiguous_splits
                for index in retained_edge_indices
            ),
        )

    @staticmethod
    def _candidate_membership_matrix(candidates, edge_count) -> np.ndarray:
        membership = np.zeros((len(candidates), edge_count), dtype=float)
        for row, candidate in enumerate(candidates):
            membership[row, list(candidate.edge_indices)] = 1.0
        return membership

    @staticmethod
    def _depth_strata(depths, requested_bins) -> np.ndarray:
        depths = np.asarray(depths, dtype=float)
        unique_depths = np.unique(depths)
        if len(unique_depths) <= requested_bins:
            return np.searchsorted(unique_depths, depths).astype(np.intp)

        boundaries = np.unique(
            np.quantile(
                depths,
                np.linspace(0.0, 1.0, requested_bins + 1)[1:-1],
            )
        )
        return np.searchsorted(boundaries, depths, side="right").astype(
            np.intp
        )

    @staticmethod
    def _scan_scores(products, membership):
        products = np.asarray(products, dtype=float)
        edge_count = products.shape[-1]
        candidate_sizes = membership.sum(axis=1)
        denominator = np.sqrt(
            candidate_sizes * (1.0 - candidate_sizes / edge_count)
        )

        if products.ndim == 1:
            centered = products - products.mean()
            inside_sums = membership @ products
            centered_sums = membership @ centered
            total = products.sum()
        else:
            centered = products - products.mean(axis=1, keepdims=True)
            inside_sums = products @ membership.T
            centered_sums = centered @ membership.T
            total = products.sum(axis=1, keepdims=True)

        outside_sizes = edge_count - candidate_sizes
        local_means = inside_sums / candidate_sizes
        background_means = (total - inside_sums) / outside_sizes
        scores = centered_sums / denominator
        return scores, local_means, background_means

    def _scan_clades(
        self,
        tree_zero_rates,
        tree_one_rates,
        membership,
        candidates,
        strata,
    ) -> list[dict]:
        observed_products = tree_zero_rates * tree_one_rates
        observed_scores, local_means, background_means = self._scan_scores(
            observed_products,
            membership,
        )
        unadjusted_counts = np.zeros(len(candidates), dtype=np.int64)
        familywise_counts = np.zeros(len(candidates), dtype=np.int64)
        groups = [
            np.flatnonzero(strata == stratum)
            for stratum in np.unique(strata)
        ]
        if not any(len(group) > 1 for group in groups):
            raise PhykitUserError(
                [
                    "No projected edges are exchangeable within depth bins.",
                    "Decrease --depth-bins or use a larger reference tree.",
                ],
                code=2,
            )

        rng = np.random.default_rng(self.seed)
        observed_absolute = np.abs(observed_scores)
        completed = 0
        while completed < self.permutations:
            batch_size = min(
                self.PERMUTATION_BATCH_SIZE,
                self.permutations - completed,
            )
            permuted = np.tile(tree_one_rates, (batch_size, 1))
            for row in range(batch_size):
                for group in groups:
                    if len(group) > 1:
                        permuted[row, group] = tree_one_rates[
                            group[rng.permutation(len(group))]
                        ]
            permutation_products = permuted * tree_zero_rates[None, :]
            permutation_scores = self._scan_scores(
                permutation_products,
                membership,
            )[0]
            permutation_absolute = np.abs(permutation_scores)
            unadjusted_counts += np.count_nonzero(
                permutation_absolute >= observed_absolute[None, :],
                axis=0,
            )
            maximum_absolute = np.max(permutation_absolute, axis=1)
            familywise_counts += np.count_nonzero(
                maximum_absolute[:, None] >= observed_absolute[None, :],
                axis=0,
            )
            completed += batch_size

        denominator = self.permutations + 1.0
        unadjusted_p = (unadjusted_counts + 1.0) / denominator
        familywise_p = (familywise_counts + 1.0) / denominator
        rows = []
        for index, candidate in enumerate(candidates):
            local_mean = float(local_means[index])
            background_mean = float(background_means[index])
            rows.append(
                {
                    "clade": candidate.label,
                    "taxa": list(candidate.taxa),
                    "taxon_count": len(candidate.taxa),
                    "root_depth": candidate.root_depth,
                    "edge_count": len(candidate.edge_indices),
                    "score": float(observed_scores[index]),
                    "local_mean_contribution": local_mean,
                    "background_mean_contribution": background_mean,
                    "contribution_contrast": local_mean - background_mean,
                    "coupling": (
                        "concordant" if local_mean >= 0.0 else "antagonistic"
                    ),
                    "unadjusted_p_value": float(unadjusted_p[index]),
                    "fwer_p_value": float(familywise_p[index]),
                    "significant": bool(familywise_p[index] <= self.alpha),
                }
            )
        rows.sort(
            key=lambda row: (
                row["fwer_p_value"],
                -abs(row["score"]),
                row["clade"],
            )
        )
        return rows

    @staticmethod
    def _leave_one_out_correlations(tree_zero_rates, tree_one_rates):
        x = np.asarray(tree_zero_rates, dtype=float)
        y = np.asarray(tree_one_rates, dtype=float)
        retained_count = len(x) - 1
        sum_x = float(x.sum())
        sum_y = float(y.sum())
        sum_xx = float(np.dot(x, x))
        sum_yy = float(np.dot(y, y))
        sum_xy = float(np.dot(x, y))
        correlations = []
        for value_x, value_y in zip(x, y):
            reduced_x = sum_x - value_x
            reduced_y = sum_y - value_y
            covariance = (
                sum_xy
                - value_x * value_y
                - reduced_x * reduced_y / retained_count
            )
            variance_x = sum_xx - value_x * value_x - reduced_x**2 / retained_count
            variance_y = sum_yy - value_y * value_y - reduced_y**2 / retained_count
            scale = math.sqrt(max(0.0, variance_x * variance_y))
            correlations.append(None if scale == 0.0 else float(covariance / scale))
        return correlations

    def _add_edge_diagnostics(self, analysis: ProjectedRateAnalysis) -> dict:
        leave_one_out = self._leave_one_out_correlations(
            analysis.standardized_zero,
            analysis.standardized_one,
        )
        influence_rows = []
        for position, edge_index in enumerate(analysis.retained_edge_indices):
            row = analysis.rows[int(edge_index)]
            correlation_without = leave_one_out[position]
            delta = (
                None
                if correlation_without is None
                else correlation_without - analysis.correlation
            )
            row["local_contribution"] = float(
                analysis.standardized_zero[position]
                * analysis.standardized_one[position]
            )
            row["leave_one_out_correlation"] = correlation_without
            row["correlation_delta_without_edge"] = delta
            if delta is not None:
                influence_rows.append((abs(delta), row, correlation_without, delta))

        for row in analysis.rows:
            row.setdefault("local_contribution", None)
            row.setdefault("leave_one_out_correlation", None)
            row.setdefault("correlation_delta_without_edge", None)

        if not influence_rows:
            return {
                "branch": None,
                "correlation_without_edge": None,
                "delta": None,
                "absolute_delta": None,
            }
        _, row, correlation_without, delta = max(
            influence_rows,
            key=lambda item: (item[0], item[1]["branch"]),
        )
        return {
            "branch": row["branch"],
            "correlation_without_edge": correlation_without,
            "delta": delta,
            "absolute_delta": abs(delta),
        }

    def _episodic_payload(
        self,
        analysis,
        scan_basis,
        strata,
        scan_rows,
        edge_influence,
    ) -> dict:
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
        payload.update(
            method="episodic_rate_covariation",
            scan_statistic="depth_matched_max_clade_contrast",
            permutations=self.permutations,
            seed=self.seed,
            requested_depth_bins=self.depth_bins,
            effective_depth_bins=len(np.unique(strata)),
            permutable_edge_count=sum(
                np.count_nonzero(strata == value)
                for value in np.unique(strata)
                if np.count_nonzero(strata == value) > 1
            ),
            min_edges=self.min_edges,
            alpha=self.alpha,
            scan_edge_count=len(scan_basis.rate_positions),
            ambiguous_root_edge_count=scan_basis.ambiguous_edge_count,
            candidate_count=len(scan_rows),
            significant_clade_count=sum(
                row["significant"] for row in scan_rows
            ),
            most_influential_edge=edge_influence,
            clade_scans=scan_rows,
        )
        if self.output:
            payload["output"] = self.output
        if self.annotated_tree:
            payload["annotated_tree"] = self.annotated_tree
        return payload

    @staticmethod
    def _write_scan_table(rows, output_path) -> None:
        columns = (
            "clade",
            "taxon_count",
            "root_depth",
            "edge_count",
            "score",
            "local_mean_contribution",
            "background_mean_contribution",
            "contribution_contrast",
            "coupling",
            "unadjusted_p_value",
            "fwer_p_value",
            "significant",
        )
        with open(output_path, "w", encoding="utf-8") as handle:
            handle.write("\t".join(columns) + "\n")
            for row in rows:
                handle.write("\t".join(str(row[column]) for column in columns) + "\n")

    def _write_annotated_tree(self, tree, rows, output_path) -> None:
        significant_by_taxa = {
            frozenset(row["taxa"]): row
            for row in rows
            if row["significant"]
        }
        annotated = self._fast_copy(tree)
        clade_taxa = {}
        clades = list(annotated.find_clades(order="postorder"))
        for clade in clades:
            if clade.clades:
                taxa = frozenset().union(
                    *(clade_taxa[id(child)] for child in clade.clades)
                )
            else:
                taxa = frozenset((clade.name,))
            clade_taxa[id(clade)] = taxa
            row = significant_by_taxa.get(taxa)
            if row is None:
                continue
            annotation = (
                f"erc_score={row['score']:.6g},"
                f"erc_fwer_p={row['fwer_p_value']:.6g},"
                f"erc_coupling={row['coupling']}"
            )
            clade.comment = (
                f"{clade.comment},{annotation}"
                if clade.comment
                else annotation
            )
        self.write_tree_file(annotated, output_path)

    def _print_episodic_text(self, payload) -> None:
        lines = [
            "Episodic Rate Covariation (experimental)",
            f"Global correlation: {payload['correlation']:.6f}",
            f"Descriptive global p-value: {payload['p_value']:.6g}",
            f"Shared taxa: {payload['shared_taxon_count']}",
            f"Retained projected edges: {payload['retained_edge_count']}",
            f"Localizable scan edges: {payload['scan_edge_count']}",
            f"Candidate clades: {payload['candidate_count']}",
            f"Permutations: {payload['permutations']}",
            "Depth bins: "
            f"{payload['effective_depth_bins']}/"
            f"{payload['requested_depth_bins']}",
            "FWER-significant clades: "
            f"{payload['significant_clade_count']} (alpha={payload['alpha']})",
        ]
        if payload["most_influential_edge"]["branch"] is not None:
            influence = payload["most_influential_edge"]
            lines.append(
                "Most influential edge: "
                f"{influence['branch']} "
                f"(absolute delta r={influence['absolute_delta']:.6f})"
            )

        reported = [
            row for row in payload["clade_scans"] if row["significant"]
        ]
        if self.verbose or not reported:
            reported = payload["clade_scans"]
        if reported:
            lines.extend(
                [
                    "",
                    "clade\tedges\tscore\tcoupling\tunadjusted_p\tfwer_p\t"
                    "significant",
                ]
            )
            for row in reported:
                lines.append(
                    f"{row['clade']}\t{row['edge_count']}\t"
                    f"{row['score']:.6f}\t{row['coupling']}\t"
                    f"{row['unadjusted_p_value']:.6g}\t"
                    f"{row['fwer_p_value']:.6g}\t{row['significant']}"
                )
        if self.output:
            lines.append(f"Saved clade scan table: {self.output}")
        if self.annotated_tree:
            lines.append(f"Saved annotated reference tree: {self.annotated_tree}")
        try:
            print("\n".join(lines))
        except BrokenPipeError:
            pass
