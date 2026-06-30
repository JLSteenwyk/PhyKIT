from __future__ import annotations

import sys
from argparse import Namespace

from .base import Alignment


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


def AlignmentOutlierTaxa(*args, **kwargs):
    from .alignment_outlier_taxa import AlignmentOutlierTaxa as _AlignmentOutlierTaxa

    return _AlignmentOutlierTaxa(*args, **kwargs)


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


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
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> dict[str, object]:
        from ...helpers.plot_config import PlotConfig

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
            plot_config=PlotConfig.from_args(args),
        )

    def _get_outlier_result(self, alignment, is_protein: bool) -> dict[str, object]:
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

    @staticmethod
    def _prepare_plot_arrays(rows, outliers, features):
        n_rows = len(rows)
        flagged_taxa = {row["taxon"] for row in outliers}
        taxa = []
        taxa_append = taxa.append
        occupancies = np.empty(n_rows, dtype=float)
        gap_rates = np.empty(n_rows, dtype=float)
        composition_distances = np.empty(n_rows, dtype=float)
        long_branch_proxies = np.empty(n_rows, dtype=float)
        rcvt = np.empty(n_rows, dtype=float)
        entropy_burden = np.empty(n_rows, dtype=float)
        has_flagged_taxa = bool(flagged_taxa)
        flagged_mask = (
            np.empty(n_rows, dtype=bool)
            if has_flagged_taxa
            else np.zeros(n_rows, dtype=bool)
        )

        if has_flagged_taxa:
            for idx, row in enumerate(rows):
                taxon = row["taxon"]
                taxa_append(taxon)
                flagged_mask[idx] = taxon in flagged_taxa
                occupancies[idx] = row["occupancy"]
                gap_rates[idx] = row["gap_rate"]
                composition_distances[idx] = row["composition_distance"]
                long_branch_proxy = row["long_branch_proxy"]
                long_branch_proxies[idx] = (
                    np.nan if long_branch_proxy is None else long_branch_proxy
                )
                rcvt[idx] = row["rcvt"]
                entropy_burden[idx] = row["entropy_burden"]
        else:
            for idx, row in enumerate(rows):
                taxa_append(row["taxon"])
                occupancies[idx] = row["occupancy"]
                gap_rates[idx] = row["gap_rate"]
                composition_distances[idx] = row["composition_distance"]
                long_branch_proxy = row["long_branch_proxy"]
                long_branch_proxies[idx] = (
                    np.nan if long_branch_proxy is None else long_branch_proxy
                )
                rcvt[idx] = row["rcvt"]
                entropy_burden[idx] = row["entropy_burden"]

        feature_vals = {
            "gap_rate": gap_rates,
            "occupancy": occupancies,
            "composition_distance": composition_distances,
            "long_branch_proxy": long_branch_proxies,
            "rcvt": rcvt,
            "entropy_burden": entropy_burden,
        }
        return {
            "taxa": taxa,
            "occupancies": occupancies,
            "gap_rates": gap_rates,
            "composition_distances": composition_distances,
            "long_branch_proxies": long_branch_proxies,
            "flagged_taxa": flagged_taxa,
            "flagged_mask": flagged_mask,
            "feature_vals": {
                feature: feature_vals[feature]
                for feature in features
            },
        }

    @staticmethod
    def _plot_composition_distance_panel(
        ax, rows, flagged_taxa, thresholds, normal_color, flagged_color,
        legend_handles, comp_values=None, branch_values=None, flagged_mask=None,
    ) -> None:
        if comp_values is None:
            comp_values = np.array(
                [row["composition_distance"] for row in rows], dtype=float
            )
        if branch_values is None:
            branch_values = np.array(
                [
                    np.nan
                    if row["long_branch_proxy"] is None
                    else row["long_branch_proxy"]
                    for row in rows
                ],
                dtype=float,
            )
        if flagged_mask is None:
            flagged_mask = np.array(
                [row["taxon"] in flagged_taxa for row in rows], dtype=bool
            )

        normal_mask = ~flagged_mask
        if np.any(normal_mask):
            ax.scatter(
                comp_values[normal_mask],
                branch_values[normal_mask],
                c=normal_color,
                s=24,
                alpha=0.85,
            )
        if np.any(flagged_mask):
            ax.scatter(
                comp_values[flagged_mask],
                branch_values[flagged_mask],
                c=flagged_color,
                s=24,
                alpha=0.85,
            )
            for idx in np.flatnonzero(flagged_mask):
                ax.text(
                    comp_values[idx],
                    branch_values[idx],
                    rows[idx]["taxon"],
                    fontsize=7,
                    color=flagged_color,
                    ha="left",
                    va="bottom",
                )

        comp_thr = thresholds["composition_distance"]
        dist_thr = thresholds["long_branch_proxy"]
        if comp_thr is not None:
            ax.axvline(comp_thr, color="black", linestyle="--", linewidth=1)
            ax.axvspan(
                comp_thr,
                max(comp_thr + 0.1, float(np.nanmax(comp_values) + 0.05)),
                color=flagged_color,
                alpha=0.06,
            )
        if dist_thr is not None:
            ax.axhline(dist_thr, color="black", linestyle="--", linewidth=1)
            finite_branch_values = branch_values[np.isfinite(branch_values)]
            upper = (
                max(dist_thr + 0.1, float(np.max(finite_branch_values) + 0.05))
                if finite_branch_values.size
                else dist_thr + 0.1
            )
            ax.axhspan(dist_thr, upper, color=flagged_color, alpha=0.06)
        ax.set_title(f"Comp Dist vs Long-Branch (thr={comp_thr}, {dist_thr})")
        ax.set_xlabel("Composition Distance")
        ax.set_ylabel("Long-Branch Proxy")
        ax.legend(handles=legend_handles, fontsize=8, loc="upper left")

    @staticmethod
    def _flag_colors(flagged_mask, normal_color, flagged_color):
        return np.where(flagged_mask, flagged_color, normal_color)

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

        plot_arrays = self._prepare_plot_arrays(
            rows,
            result["outliers"],
            result["features"],
        )
        taxa = plot_arrays["taxa"]
        occupancies = plot_arrays["occupancies"]
        gap_rates = plot_arrays["gap_rates"]
        flagged_taxa = plot_arrays["flagged_taxa"]
        flagged_mask = plot_arrays["flagged_mask"]

        # panel 1/2 ordering for easier visual scanning
        order = np.argsort(occupancies)
        taxa_ordered = [taxa[i] for i in order]
        occupancy_ordered = occupancies[order]
        gap_ordered = gap_rates[order]
        flagged_ordered = flagged_mask[order]
        x_positions = np.arange(len(taxa_ordered))

        fig, axes = plt.subplots(2, 2, figsize=(self.width, self.height))
        flagged_color = "#ca0020"
        normal_color = "#636363"
        legend_handles = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor=normal_color, markersize=7, label="Not flagged"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=flagged_color, markersize=7, label="Flagged"),
        ]

        ax = axes[0, 0]
        colors = self._flag_colors(flagged_ordered, normal_color, flagged_color)
        ax.bar(x_positions, occupancy_ordered, color=colors)
        occ_thr = thresholds["occupancy"]
        if occ_thr is not None:
            ax.axhline(occ_thr, color="black", linestyle="--", linewidth=1)
            ax.axhspan(0.0, occ_thr, color=flagged_color, alpha=0.08)
        ax.set_title(f"Occupancy Per Taxon (thr={occ_thr})")
        ax.set_ylabel("Occupancy")
        ax.set_xticks(x_positions)
        ax.set_xticklabels(taxa_ordered, rotation=90, fontsize=7)
        ax.legend(handles=legend_handles, fontsize=8, loc="lower right")

        ax = axes[0, 1]
        ax.bar(x_positions, gap_ordered, color=colors)
        gap_thr = thresholds["gap_rate"]
        if gap_thr is not None:
            ax.axhline(gap_thr, color="black", linestyle="--", linewidth=1)
            ax.axhspan(gap_thr, max(1.0, float(np.max(gap_ordered) + 0.05)), color=flagged_color, alpha=0.08)
        ax.set_title(f"Gap Rate Per Taxon (thr={gap_thr})")
        ax.set_ylabel("Gap Rate")
        ax.set_xticks(x_positions)
        ax.set_xticklabels(taxa_ordered, rotation=90, fontsize=7)

        ax = axes[1, 0]
        self._plot_composition_distance_panel(
            ax, rows, flagged_taxa, thresholds,
            normal_color, flagged_color, legend_handles,
            plot_arrays["composition_distances"],
            plot_arrays["long_branch_proxies"],
            flagged_mask,
        )

        # Panel 4: robust z-score heatmap (taxa x features)
        ax = axes[1, 1]
        features = result["features"]
        feature_vals = plot_arrays["feature_vals"]

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
