from typing import Dict, List, Tuple, Union
from collections import Counter
import numpy as np

from scipy.stats import chisquare, false_discovery_control
from scipy.stats._stats_py import Power_divergenceResult
from Bio.Align import MultipleSeqAlignment

from .base import Alignment
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig


class CompositionalBiasPerSite(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.json_output = parsed["json_output"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()

        stat_res, p_vals_corrected = \
            self.calculate_compositional_bias_per_site(alignment, is_protein)

        rows = self._build_rows(stat_res, p_vals_corrected)

        if self.plot:
            self._plot_compositional_bias_manhattan(rows)

        if self.json_output:
            payload = dict(rows=rows, sites=rows)
            if self.plot:
                payload["plot_output"] = self.plot_output
            print_json(payload)
            return

        for row in rows:
            pval_cor_str = "nan" if row["p_value_corrected"] is None else row["p_value_corrected"]
            p_val_str = "nan" if row["p_value"] is None else row["p_value"]
            print(f"{row['site']}\t{row['chi_square']}\t{pval_cor_str}\t{p_val_str}")

        if self.plot:
            print(f"Saved compositional bias plot: {self.plot_output}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
            plot=getattr(args, "plot", False),
            plot_output=getattr(args, "plot_output", "compositional_bias_per_site_plot.png"),
            plot_config=PlotConfig.from_args(args),
        )

    def _build_rows(
        self,
        stat_res: List[Power_divergenceResult],
        p_vals_corrected: List[Union[float, str]],
    ) -> List[Dict[str, Union[int, float, None]]]:
        rows = []
        for idx, (stat_info, pval_cor) in enumerate(
            zip(stat_res, p_vals_corrected), start=1
        ):
            corrected = None if isinstance(pval_cor, str) else round(float(pval_cor), 4)
            raw_p = float(stat_info.pvalue)
            rows.append(
                dict(
                    site=idx,
                    chi_square=round(float(stat_info.statistic), 4),
                    p_value_corrected=corrected,
                    p_value=None if np.isnan(raw_p) else round(raw_p, 4),
                )
            )
        return rows

    def _plot_compositional_bias_manhattan(
        self,
        rows: List[Dict[str, Union[int, float, None]]],
        alpha: float = 0.05,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in compositional_bias_per_site. Install matplotlib and retry.")
            raise SystemExit(2)

        if not rows:
            return

        sites = np.array([int(row["site"]) for row in rows], dtype=np.int32)
        corrected_pvals = []
        for row in rows:
            value = row["p_value_corrected"]
            if value is None:
                corrected_pvals.append(np.nan)
            else:
                corrected_pvals.append(float(value))
        corrected_pvals = np.array(corrected_pvals, dtype=np.float64)

        finite_mask = np.isfinite(corrected_pvals) & (corrected_pvals > 0.0)
        y_vals = np.full_like(corrected_pvals, np.nan, dtype=np.float64)
        y_vals[finite_mask] = -np.log10(corrected_pvals[finite_mask])

        sig_threshold = -np.log10(alpha)
        sig_mask = finite_mask & (corrected_pvals < alpha)
        nonsig_mask = finite_mask & ~sig_mask

        config = self.plot_config
        config.resolve(n_rows=len(sites), n_cols=None)
        default_colors = ["#2b8cbe", "#d62728", "#000000"]
        colors = config.merge_colors(default_colors)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        if np.any(nonsig_mask):
            ax.scatter(
                sites[nonsig_mask],
                y_vals[nonsig_mask],
                s=10,
                alpha=0.75,
                color=colors[0],
                edgecolors="none",
                label="Not significant",
            )
        if np.any(sig_mask):
            ax.scatter(
                sites[sig_mask],
                y_vals[sig_mask],
                s=14,
                alpha=0.9,
                color=colors[1],
                edgecolors="none",
                label=f"FDR < {alpha}",
            )

        ax.axhline(sig_threshold, color=colors[2], linestyle="--", linewidth=1.5, label=f"-log10({alpha})")
        ax.set_xlabel("Alignment site")
        ax.set_ylabel("-log10(corrected p-value)")
        ax.set_xlim(1, int(np.max(sites)))
        if np.any(finite_mask):
            max_y = float(np.nanmax(y_vals))
            ax.set_ylim(0, max(sig_threshold * 1.1, max_y * 1.1))

        legend_loc = config.legend_position or "upper right"
        if legend_loc != "none":
            ax.legend(loc=legend_loc, frameon=False, fontsize=8)

        if config.show_title:
            ax.set_title(config.title or "Compositional Bias Per Site (Manhattan)", fontsize=config.title_fontsize)
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def get_number_of_occurrences_per_character(
        self,
        alignment: MultipleSeqAlignment,
        idx: int,
        is_protein: bool = False,
    ) -> List[int]:
        gap_chars = self.get_gap_chars(is_protein)
        seq_at_position = alignment[:, idx].upper()
        filtered_seq = "".join([char for char in seq_at_position if char not in gap_chars])

        return list(Counter(filtered_seq).values())

    def calculate_compositional_bias_per_site(
        self,
        alignment: MultipleSeqAlignment,
        is_protein: bool = False,
    ) -> Tuple[
        List[Power_divergenceResult],
        List[Union[float, str]],
    ]:
        aln_len = alignment.get_alignment_length()
        gap_chars = set(self.get_gap_chars(is_protein))

        # Convert alignment to numpy array for faster operations
        alignment_array = np.array([
            [c.upper() for c in str(record.seq)]
            for record in alignment
        ], dtype='U1')

        stat_res = []
        p_vals = []
        nan_idx = []

        # Process each column
        for col_idx in range(aln_len):
            column = alignment_array[:, col_idx]

            # Filter out gaps
            non_gap_mask = ~np.isin(column, list(gap_chars))
            filtered_column = column[non_gap_mask]

            if len(filtered_column) > 0:
                # Count occurrences using numpy
                unique_chars, counts = np.unique(filtered_column, return_counts=True)

                # Perform chi-square test
                chisquare_res = chisquare(counts)
                stat_res.append(chisquare_res)

                if not np.isnan(chisquare_res.pvalue):
                    p_vals.append(chisquare_res.pvalue)
                else:
                    nan_idx.append(col_idx)
            else:
                # Handle empty column
                dummy_res = chisquare([1])  # Create dummy result
                stat_res.append(dummy_res)
                nan_idx.append(col_idx)

        # Apply FDR correction
        if p_vals:
            p_vals_corrected = list(false_discovery_control(p_vals))
        else:
            p_vals_corrected = []

        # Insert NaNs at appropriate positions
        for idx in reversed(nan_idx):
            p_vals_corrected.insert(idx, "nan")

        return stat_res, p_vals_corrected
