from typing import Dict, List

import numpy as np

from .base import Alignment
from ...helpers.json_output import print_json


class AlignmentOutlierTaxa(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.gap_z = parsed["gap_z"]
        self.composition_z = parsed["composition_z"]
        self.distance_z = parsed["distance_z"]
        self.rcvt_z = parsed["rcvt_z"]
        self.occupancy_z = parsed["occupancy_z"]
        self.entropy_z = parsed["entropy_z"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            gap_z=args.gap_z,
            composition_z=args.composition_z,
            distance_z=args.distance_z,
            rcvt_z=args.rcvt_z,
            occupancy_z=args.occupancy_z,
            entropy_z=args.entropy_z,
            json_output=getattr(args, "json", False),
        )

    def _high_outlier_threshold(self, values: np.ndarray, z: float) -> float:
        finite_values = values[np.isfinite(values)]
        if finite_values.size < 2:
            return float("inf")

        median = float(np.median(finite_values))
        mad = float(np.median(np.abs(finite_values - median)))
        if mad > 0.0:
            robust_sigma = 1.4826 * mad
            return median + (z * robust_sigma)

        # Fallback when variation is very low: flag values above second-highest.
        unique_vals = np.unique(finite_values)
        if unique_vals.size > 1:
            return float(unique_vals[-2])
        return float("inf")

    def _low_outlier_threshold(self, values: np.ndarray, z: float) -> float:
        finite_values = values[np.isfinite(values)]
        if finite_values.size < 2:
            return float("-inf")

        median = float(np.median(finite_values))
        mad = float(np.median(np.abs(finite_values - median)))
        if mad > 0.0:
            robust_sigma = 1.4826 * mad
            return median - (z * robust_sigma)

        # Fallback when variation is very low: flag values below second-lowest.
        unique_vals = np.unique(finite_values)
        if unique_vals.size > 1:
            return float(unique_vals[1])
        return float("-inf")

    def _invalid_chars(self, is_protein: bool) -> np.ndarray:
        if is_protein:
            return np.array(["-", "?", "*", "X"], dtype="U1")
        return np.array(["-", "?", "*", "X", "N"], dtype="U1")

    def calculate_outliers(self, alignment, is_protein: bool) -> Dict[str, object]:
        taxa = [record.id for record in alignment]
        n_taxa = len(taxa)
        if n_taxa == 0:
            return dict(
                features=[
                    "gap_rate",
                    "occupancy",
                    "composition_distance",
                    "long_branch_proxy",
                    "rcvt",
                    "entropy_burden",
                ],
                thresholds=dict(
                    gap_rate=float("inf"),
                    occupancy=float("-inf"),
                    composition_distance=float("inf"),
                    long_branch_proxy=float("inf"),
                    rcvt=float("inf"),
                    entropy_burden=float("inf"),
                ),
                rows=[],
                outliers=[],
            )

        alignment_array = np.array(
            [[c.upper() for c in str(record.seq)] for record in alignment],
            dtype="U1",
        )
        aln_len = alignment_array.shape[1] if alignment_array.ndim == 2 else 0

        invalid_chars = self._invalid_chars(is_protein)
        valid_mask = ~np.isin(alignment_array, invalid_chars)
        valid_lengths = np.sum(valid_mask, axis=1).astype(np.float64)

        if aln_len > 0:
            gap_rates = 1.0 - (valid_lengths / float(aln_len))
            occupancies = valid_lengths / float(aln_len)
        else:
            gap_rates = np.ones(n_taxa, dtype=np.float64)
            occupancies = np.zeros(n_taxa, dtype=np.float64)

        valid_chars = alignment_array[valid_mask]
        if valid_chars.size > 0:
            symbols = np.unique(valid_chars)
            comp_matrix = np.zeros((n_taxa, len(symbols)), dtype=np.float64)
            for i, seq in enumerate(alignment_array):
                if valid_lengths[i] <= 0:
                    continue
                seq_valid = valid_mask[i]
                for j, symbol in enumerate(symbols):
                    comp_matrix[i, j] = np.sum((seq == symbol) & seq_valid) / valid_lengths[i]
            center = np.median(comp_matrix, axis=0)
            composition_distances = np.linalg.norm(comp_matrix - center, axis=1)
        else:
            composition_distances = np.zeros(n_taxa, dtype=np.float64)
            comp_matrix = np.zeros((n_taxa, 0), dtype=np.float64)

        long_branch_proxy = np.full(n_taxa, np.nan, dtype=np.float64)
        for i in range(n_taxa):
            pairwise_distances = []
            for j in range(n_taxa):
                if i == j:
                    continue
                overlap = valid_mask[i] & valid_mask[j]
                overlap_count = int(np.sum(overlap))
                if overlap_count == 0:
                    continue
                mismatches = int(np.sum(alignment_array[i, overlap] != alignment_array[j, overlap]))
                pairwise_distances.append(mismatches / overlap_count)
            if pairwise_distances:
                long_branch_proxy[i] = float(np.mean(pairwise_distances))

        gap_threshold = self._high_outlier_threshold(gap_rates, self.gap_z)
        composition_threshold = self._high_outlier_threshold(
            composition_distances, self.composition_z
        )
        distance_threshold = self._high_outlier_threshold(
            long_branch_proxy, self.distance_z
        )
        if comp_matrix.shape[1] > 0:
            average_counts = np.sum(comp_matrix, axis=0) / n_taxa
            deviations = np.abs(comp_matrix - average_counts)
            seq_sums = np.sum(deviations, axis=1)
            denom = n_taxa * valid_lengths
            rcvt_values = np.divide(
                seq_sums,
                denom,
                out=np.zeros_like(seq_sums, dtype=np.float64),
                where=denom > 0,
            )
        else:
            rcvt_values = np.zeros(n_taxa, dtype=np.float64)
        rcvt_threshold = self._high_outlier_threshold(rcvt_values, self.rcvt_z)
        occupancy_threshold = self._low_outlier_threshold(occupancies, self.occupancy_z)

        # Entropy burden: average site entropy over valid positions for each taxon.
        if aln_len > 0:
            site_entropies = np.zeros(aln_len, dtype=np.float64)
            for col_idx in range(aln_len):
                column = alignment_array[:, col_idx]
                col_valid = ~np.isin(column, invalid_chars)
                vals = column[col_valid]
                if vals.size == 0:
                    continue
                _, counts = np.unique(vals, return_counts=True)
                probs = counts / float(np.sum(counts))
                site_entropies[col_idx] = -float(np.sum(probs * np.log2(probs)))

            entropy_burden = np.zeros(n_taxa, dtype=np.float64)
            for i in range(n_taxa):
                if np.any(valid_mask[i]):
                    entropy_burden[i] = float(np.mean(site_entropies[valid_mask[i]]))
        else:
            entropy_burden = np.zeros(n_taxa, dtype=np.float64)
        entropy_threshold = self._high_outlier_threshold(entropy_burden, self.entropy_z)

        rows = []
        outliers = []
        for i, taxon in enumerate(taxa):
            reasons: List[Dict[str, float]] = []

            if np.isfinite(gap_rates[i]) and gap_rates[i] > gap_threshold:
                reasons.append(
                    dict(
                        feature="gap_rate",
                        value=float(gap_rates[i]),
                        threshold=float(gap_threshold),
                        explanation="High fraction of gap/ambiguous symbols compared to other taxa.",
                    )
                )
            if np.isfinite(occupancies[i]) and occupancies[i] < occupancy_threshold:
                reasons.append(
                    dict(
                        feature="occupancy",
                        value=float(occupancies[i]),
                        threshold=float(occupancy_threshold),
                        explanation="Low fraction of valid symbols compared to other taxa.",
                    )
                )
            if (
                np.isfinite(composition_distances[i])
                and composition_distances[i] > composition_threshold
            ):
                reasons.append(
                    dict(
                        feature="composition_distance",
                        value=float(composition_distances[i]),
                        threshold=float(composition_threshold),
                        explanation="Unusual sequence composition profile relative to other taxa.",
                    )
                )
            if np.isfinite(long_branch_proxy[i]) and long_branch_proxy[i] > distance_threshold:
                reasons.append(
                    dict(
                        feature="long_branch_proxy",
                        value=float(long_branch_proxy[i]),
                        threshold=float(distance_threshold),
                        explanation="High mean pairwise sequence distance to other taxa.",
                    )
                )
            if np.isfinite(rcvt_values[i]) and rcvt_values[i] > rcvt_threshold:
                reasons.append(
                    dict(
                        feature="rcvt",
                        value=float(rcvt_values[i]),
                        threshold=float(rcvt_threshold),
                        explanation="High relative composition variability for this taxon.",
                    )
                )
            if np.isfinite(entropy_burden[i]) and entropy_burden[i] > entropy_threshold:
                reasons.append(
                    dict(
                        feature="entropy_burden",
                        value=float(entropy_burden[i]),
                        threshold=float(entropy_threshold),
                        explanation="High average site entropy across this taxon's valid positions.",
                    )
                )

            row = dict(
                taxon=taxon,
                gap_rate=round(float(gap_rates[i]), 4),
                occupancy=round(float(occupancies[i]), 4),
                composition_distance=round(float(composition_distances[i]), 4),
                long_branch_proxy=(
                    None if not np.isfinite(long_branch_proxy[i]) else round(float(long_branch_proxy[i]), 4)
                ),
                rcvt=round(float(rcvt_values[i]), 4),
                entropy_burden=round(float(entropy_burden[i]), 4),
                flagged=bool(reasons),
                reasons=[
                    dict(
                        feature=r["feature"],
                        value=round(float(r["value"]), 4),
                        threshold=round(float(r["threshold"]), 4),
                        direction=("low" if r["feature"] == "occupancy" else "high"),
                        explanation=r["explanation"],
                    )
                    for r in reasons
                ],
            )
            rows.append(row)
            if reasons:
                outliers.append(row)

        return dict(
            features=[
                "gap_rate",
                "occupancy",
                "composition_distance",
                "long_branch_proxy",
                "rcvt",
                "entropy_burden",
            ],
            thresholds=dict(
                gap_rate=round(float(gap_threshold), 4)
                if np.isfinite(gap_threshold)
                else None,
                occupancy=round(float(occupancy_threshold), 4)
                if np.isfinite(occupancy_threshold)
                else None,
                composition_distance=round(float(composition_threshold), 4)
                if np.isfinite(composition_threshold)
                else None,
                long_branch_proxy=round(float(distance_threshold), 4)
                if np.isfinite(distance_threshold)
                else None,
                rcvt=round(float(rcvt_threshold), 4)
                if np.isfinite(rcvt_threshold)
                else None,
                entropy_burden=round(float(entropy_threshold), 4)
                if np.isfinite(entropy_threshold)
                else None,
            ),
            rows=rows,
            outliers=outliers,
        )

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()
        result = self.calculate_outliers(alignment, is_protein)

        if self.json_output:
            payload = dict(
                features_evaluated=result["features"],
                thresholds=result["thresholds"],
                outliers=result["outliers"],
                rows=result["rows"],
                taxa=result["rows"],
            )
            print_json(payload)
            return

        print(
            "features_evaluated\t"
            "gap_rate,occupancy,composition_distance,long_branch_proxy,rcvt,entropy_burden"
        )
        print(
            "thresholds\t"
            f"gap_rate>{result['thresholds']['gap_rate']};"
            f"occupancy<{result['thresholds']['occupancy']};"
            f"composition_distance>{result['thresholds']['composition_distance']};"
            f"long_branch_proxy>{result['thresholds']['long_branch_proxy']};"
            f"rcvt>{result['thresholds']['rcvt']};"
            f"entropy_burden>{result['thresholds']['entropy_burden']}"
        )

        if not result["outliers"]:
            print("No outlier taxa detected.")
            return

        for row in result["outliers"]:
            reason_string = ";".join(
                f"{reason['feature']}={reason['value']}"
                f"{'<' if reason.get('direction') == 'low' else '>'}{reason['threshold']}"
                for reason in row["reasons"]
            )
            explanation_string = " | ".join(
                reason["explanation"] for reason in row["reasons"]
            )
            print(f"{row['taxon']}\t{reason_string}\t{explanation_string}")
