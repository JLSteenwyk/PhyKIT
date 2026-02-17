import sys
from argparse import Namespace
from typing import Dict

import numpy as np

from .base import Alignment
from .alignment_outlier_taxa import AlignmentOutlierTaxa
from ...helpers.json_output import print_json


class PlotAlignmentQC(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.output = parsed["output"]
        self.width = parsed["width"]
        self.height = parsed["height"]
        self.dpi = parsed["dpi"]
        self.gap_z = parsed["gap_z"]
        self.composition_z = parsed["composition_z"]
        self.distance_z = parsed["distance_z"]
        self.rcvt_z = parsed["rcvt_z"]
        self.occupancy_z = parsed["occupancy_z"]
        self.entropy_z = parsed["entropy_z"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> Dict[str, object]:
        return dict(
            alignment_file_path=args.alignment,
            output=args.output,
            width=args.width,
            height=args.height,
            dpi=args.dpi,
            gap_z=args.gap_z,
            composition_z=args.composition_z,
            distance_z=args.distance_z,
            rcvt_z=args.rcvt_z,
            occupancy_z=args.occupancy_z,
            entropy_z=args.entropy_z,
            json_output=getattr(args, "json", False),
        )

    def _get_outlier_result(self, alignment, is_protein: bool) -> Dict[str, object]:
        svc = AlignmentOutlierTaxa(
            Namespace(
                alignment=self.alignment_file_path,
                gap_z=self.gap_z,
                composition_z=self.composition_z,
                distance_z=self.distance_z,
                rcvt_z=self.rcvt_z,
                occupancy_z=self.occupancy_z,
                entropy_z=self.entropy_z,
                json=False,
            )
        )
        return svc.calculate_outliers(alignment, is_protein)

    def run(self) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.lines import Line2D
        except ImportError:
            print("matplotlib is required for plot_alignment_qc. Install matplotlib and retry.")
            sys.exit(2)

        alignment, _, is_protein = self.get_alignment_and_format()
        result = self._get_outlier_result(alignment, is_protein)
        rows = result["rows"]
        thresholds = result["thresholds"]

        taxa = [row["taxon"] for row in rows]
        occupancies = np.array([row["occupancy"] for row in rows], dtype=float)
        gap_rates = np.array([row["gap_rate"] for row in rows], dtype=float)

        flagged_taxa = {row["taxon"] for row in result["outliers"]}
        flagged_mask = np.array([taxon in flagged_taxa for taxon in taxa], dtype=bool)

        # panel 1/2 ordering for easier visual scanning
        order = np.argsort(occupancies)
        taxa_ordered = [taxa[i] for i in order]
        occupancy_ordered = occupancies[order]
        gap_ordered = gap_rates[order]
        flagged_ordered = flagged_mask[order]

        fig, axes = plt.subplots(2, 2, figsize=(self.width, self.height))
        flagged_color = "#ca0020"
        normal_color = "#636363"
        legend_handles = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor=normal_color, markersize=7, label="Not flagged"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=flagged_color, markersize=7, label="Flagged"),
        ]

        ax = axes[0, 0]
        colors = [flagged_color if is_flagged else normal_color for is_flagged in flagged_ordered]
        ax.bar(np.arange(len(taxa_ordered)), occupancy_ordered, color=colors)
        occ_thr = thresholds["occupancy"]
        if occ_thr is not None:
            ax.axhline(occ_thr, color="black", linestyle="--", linewidth=1)
            ax.axhspan(0.0, occ_thr, color=flagged_color, alpha=0.08)
        ax.set_title(f"Occupancy Per Taxon (thr={occ_thr})")
        ax.set_ylabel("Occupancy")
        ax.set_xticks(np.arange(len(taxa_ordered)))
        ax.set_xticklabels(taxa_ordered, rotation=90, fontsize=7)
        ax.legend(handles=legend_handles, fontsize=8, loc="lower right")

        ax = axes[0, 1]
        ax.bar(np.arange(len(taxa_ordered)), gap_ordered, color=colors)
        gap_thr = thresholds["gap_rate"]
        if gap_thr is not None:
            ax.axhline(gap_thr, color="black", linestyle="--", linewidth=1)
            ax.axhspan(gap_thr, max(1.0, float(np.max(gap_ordered) + 0.05)), color=flagged_color, alpha=0.08)
        ax.set_title(f"Gap Rate Per Taxon (thr={gap_thr})")
        ax.set_ylabel("Gap Rate")
        ax.set_xticks(np.arange(len(taxa_ordered)))
        ax.set_xticklabels(taxa_ordered, rotation=90, fontsize=7)

        ax = axes[1, 0]
        for row in rows:
            x = row["composition_distance"]
            y = row["long_branch_proxy"] if row["long_branch_proxy"] is not None else np.nan
            color = flagged_color if row["taxon"] in flagged_taxa else normal_color
            ax.scatter(x, y, c=color, s=24, alpha=0.85)
            if row["taxon"] in flagged_taxa:
                ax.text(x, y, row["taxon"], fontsize=7, color=flagged_color, ha="left", va="bottom")
        comp_thr = thresholds["composition_distance"]
        dist_thr = thresholds["long_branch_proxy"]
        if comp_thr is not None:
            ax.axvline(comp_thr, color="black", linestyle="--", linewidth=1)
            ax.axvspan(comp_thr, max(comp_thr + 0.1, float(np.nanmax([r["composition_distance"] for r in rows]) + 0.05)), color=flagged_color, alpha=0.06)
        if dist_thr is not None:
            ax.axhline(dist_thr, color="black", linestyle="--", linewidth=1)
            long_vals = [r["long_branch_proxy"] for r in rows if r["long_branch_proxy"] is not None]
            upper = max(dist_thr + 0.1, float(np.max(long_vals) + 0.05)) if long_vals else dist_thr + 0.1
            ax.axhspan(dist_thr, upper, color=flagged_color, alpha=0.06)
        ax.set_title(f"Comp Dist vs Long-Branch (thr={comp_thr}, {dist_thr})")
        ax.set_xlabel("Composition Distance")
        ax.set_ylabel("Long-Branch Proxy")
        ax.legend(handles=legend_handles, fontsize=8, loc="upper left")

        # Panel 4: robust z-score heatmap (taxa x features)
        ax = axes[1, 1]
        features = result["features"]
        feature_vals = {
            "gap_rate": np.array([row["gap_rate"] for row in rows], dtype=float),
            "occupancy": np.array([row["occupancy"] for row in rows], dtype=float),
            "composition_distance": np.array([row["composition_distance"] for row in rows], dtype=float),
            "long_branch_proxy": np.array(
                [np.nan if row["long_branch_proxy"] is None else row["long_branch_proxy"] for row in rows], dtype=float
            ),
            "rcvt": np.array([row["rcvt"] for row in rows], dtype=float),
            "entropy_burden": np.array([row["entropy_burden"] for row in rows], dtype=float),
        }

        z_data = np.zeros((len(rows), len(features)), dtype=float)
        for col_idx, feature in enumerate(features):
            values = feature_vals[feature]
            finite_mask = np.isfinite(values)
            finite_vals = values[finite_mask]
            if finite_vals.size < 2:
                continue
            median = float(np.median(finite_vals))
            mad = float(np.median(np.abs(finite_vals - median)))
            sigma = 1.4826 * mad if mad > 0 else float(np.std(finite_vals))
            if sigma <= 0:
                continue
            if feature == "occupancy":
                z_data[:, col_idx] = (median - values) / sigma
            else:
                z_data[:, col_idx] = (values - median) / sigma

        heat_order = np.argsort(~flagged_mask)
        z_plot = z_data[heat_order, :]
        taxa_heat = [taxa[i] for i in heat_order]
        im = ax.imshow(z_plot, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
        ax.set_title("Per-Taxon Feature Robust Z-Scores")
        ax.set_yticks(np.arange(len(taxa_heat)))
        ax.set_yticklabels(taxa_heat, fontsize=7)
        ax.set_xticks(np.arange(len(features)))
        ax.set_xticklabels(features, rotation=30, ha="right", fontsize=8)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Robust z-score")

        n_taxa = len(rows)
        aln_len = alignment.get_alignment_length() if len(alignment) > 0 else 0
        fig.suptitle(
            f"Alignment QC: taxa={n_taxa}, sites={aln_len}, outliers={len(result['outliers'])}",
            fontsize=12,
        )
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        fig.savefig(self.output, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)

        if self.json_output:
            print_json(
                dict(
                    output_file=self.output,
                    n_taxa=n_taxa,
                    alignment_length=aln_len,
                    features_evaluated=result["features"],
                    thresholds=result["thresholds"],
                    outlier_count=len(result["outliers"]),
                    outliers=result["outliers"],
                )
            )
            return

        print(f"Saved alignment QC plot: {self.output}")
