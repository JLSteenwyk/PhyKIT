from __future__ import annotations

from collections import namedtuple
from math import erfc, isnan, pi, sqrt


from .base import Alignment


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)

class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()
_DNA_GAP_LOOKUP = None
_PROTEIN_GAP_LOOKUP = None
_FDR_VECTOR_MIN_LENGTH = 2048


def _column_sum_squares(counts: np.ndarray) -> np.ndarray:
    return np.einsum("ij,ij->j", counts, counts, dtype=np.float64)


def _column_totals(counts: np.ndarray) -> np.ndarray:
    if counts.shape[1] <= 20000:
        return counts.sum(axis=0)
    return np.sum(counts, axis=0)


def _get_gap_lookup(is_protein: bool):
    global _DNA_GAP_LOOKUP, _PROTEIN_GAP_LOOKUP
    if is_protein:
        if _PROTEIN_GAP_LOOKUP is None:
            lookup = np.zeros(256, dtype=np.bool_)
            lookup[np.frombuffer("-?*X".encode("ascii"), dtype=np.uint8)] = True
            _PROTEIN_GAP_LOOKUP = lookup
        return _PROTEIN_GAP_LOOKUP

    if _DNA_GAP_LOOKUP is None:
        lookup = np.zeros(256, dtype=np.bool_)
        lookup[np.frombuffer("-?*XN".encode("ascii"), dtype=np.uint8)] = True
        _DNA_GAP_LOOKUP = lookup
    return _DNA_GAP_LOOKUP


Power_divergenceResult = namedtuple(
    "Power_divergenceResult",
    ["statistic", "pvalue"],
)


def _erfc_array(values: np.ndarray) -> np.ndarray:
    sqrt_values = np.sqrt(values)
    return np.fromiter(
        (erfc(float(value)) for value in sqrt_values),
        dtype=np.float64,
        count=sqrt_values.size,
    )


def _chi2_sf_for_integer_df(statistics: np.ndarray, degrees_of_freedom: int) -> np.ndarray:
    y_values = statistics * 0.5
    if degrees_of_freedom % 2 == 0:
        term = np.exp(-y_values)
        total = term.copy()
        for idx in range(1, degrees_of_freedom // 2):
            term = term * y_values / idx
            total += term
        return np.clip(total, 0.0, 1.0)

    total = _erfc_array(y_values)
    if degrees_of_freedom == 1:
        return total

    term = np.exp(-y_values) * 2.0 * np.sqrt(y_values) / sqrt(pi)
    total += term
    current_shape = 0.5
    for _ in range(1, degrees_of_freedom // 2):
        term = term * y_values / (current_shape + 1.0)
        total += term
        current_shape += 1.0
    return np.clip(total, 0.0, 1.0)


def _chi2_sf(statistics: np.ndarray, degrees_of_freedom: np.ndarray) -> np.ndarray:
    p_values = np.empty(statistics.size, dtype=np.float64)
    for degree in np.unique(degrees_of_freedom):
        degree_mask = degrees_of_freedom == degree
        p_values[degree_mask] = _chi2_sf_for_integer_df(
            statistics[degree_mask],
            int(degree),
        )
    return p_values


def _false_discovery_control(p_values: list[float]) -> list[float]:
    p_value_count = len(p_values)
    if p_value_count == 0:
        return []

    if p_value_count < _FDR_VECTOR_MIN_LENGTH:
        ordered = sorted(enumerate(p_values), key=lambda item: item[1])
        adjusted = [0.0] * p_value_count
        previous = 1.0
        for rank_index in range(p_value_count - 1, -1, -1):
            original_index, p_value = ordered[rank_index]
            rank = rank_index + 1
            corrected = min(p_value * p_value_count / rank, previous)
            corrected = min(corrected, 1.0)
            adjusted[original_index] = corrected
            previous = corrected
        return adjusted

    p_array = np.asarray(p_values, dtype=np.float64)

    order = np.argsort(p_array)
    ordered = p_array[order]
    ranks = np.arange(1, ordered.size + 1, dtype=np.float64)
    adjusted_ordered = ordered * ordered.size / ranks
    adjusted_ordered = np.minimum.accumulate(adjusted_ordered[::-1])[::-1]
    adjusted = np.empty_like(adjusted_ordered)
    adjusted[order] = np.minimum(adjusted_ordered, 1.0)
    return adjusted.tolist()


def _restore_nan_corrected_pvalues(
    p_vals_corrected: list[float],
    nan_mask: np.ndarray,
) -> list[float | str]:
    mask_size = int(nan_mask.size)
    if len(p_vals_corrected) == mask_size:
        return p_vals_corrected

    corrected_full: list[float | str] = ["nan"] * mask_size
    if not p_vals_corrected:
        return corrected_full

    valid_indices = np.flatnonzero(~nan_mask)
    for idx, value in zip(valid_indices, p_vals_corrected):
        corrected_full[int(idx)] = value
    return corrected_full


def _power_divergence_results_from_arrays(statistics, p_values):
    result_type = Power_divergenceResult
    return [
        result_type(float(statistic), float(p_value))
        for statistic, p_value in zip(statistics, p_values)
    ]


def _pvalues_for_fdr(p_values: np.ndarray, nan_mask: np.ndarray) -> list[float]:
    if not nan_mask.any():
        return p_values.tolist()
    return p_values[~nan_mask].tolist()


def _prepare_compositional_bias_plot_series(p_vals_corrected: list[float | str]):
    count = len(p_vals_corrected)
    sites = np.arange(1, count + 1, dtype=np.int32)
    corrected_pvals = np.fromiter(
        (
            np.nan if isinstance(value, str) else round(float(value), 4)
            for value in p_vals_corrected
        ),
        dtype=np.float64,
        count=count,
    )
    return sites, corrected_pvals


def _prepare_compositional_bias_row_plot_series(
    rows: list[dict[str, int | float | None]],
):
    sites = np.array([int(row["site"]) for row in rows], dtype=np.int32)
    corrected_pvals = np.fromiter(
        (
            np.nan if row["p_value_corrected"] is None else float(row["p_value_corrected"])
            for row in rows
        ),
        dtype=np.float64,
        count=len(rows),
    )
    return sites, corrected_pvals


def _column_count_stats_from_ascii_codes(
    alignment_array: np.ndarray,
    valid_mask: np.ndarray | None,
    valid_symbols: np.ndarray,
    block_size: int = 8192,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    aln_len = alignment_array.shape[1]
    category_counts = np.zeros(aln_len, dtype=np.int16)
    totals = np.zeros(aln_len, dtype=np.float64)
    sum_squares = np.zeros(aln_len, dtype=np.float64)
    symbol_codes = valid_symbols.astype(np.intp, copy=False)

    for start in range(0, aln_len, block_size):
        stop = min(start + block_size, aln_len)
        width = stop - start
        block = alignment_array[:, start:stop]

        block_by_site = block.T.astype(np.intp, copy=False)
        site_offsets = np.arange(width, dtype=np.intp)[:, None] * 256
        if valid_mask is None:
            count_index = (block_by_site + site_offsets).ravel()
        else:
            block_valid = valid_mask[:, start:stop]
            if not block_valid.any():
                continue
            count_index = (block_by_site + site_offsets)[block_valid.T]
        counts = np.bincount(
            count_index,
            minlength=width * 256,
        ).reshape(width, 256)[:, symbol_codes].T

        slc = slice(start, stop)
        category_counts[slc] = np.count_nonzero(counts, axis=0)
        totals[slc] = counts.sum(axis=0, dtype=np.float64)
        sum_squares[slc] = _column_sum_squares(counts)

    return category_counts, totals, sum_squares


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

        rows = None

        if self.plot and self.json_output:
            rows = self._build_rows(stat_res, p_vals_corrected)
            self._plot_compositional_bias_manhattan(rows)
        elif self.plot:
            self._plot_compositional_bias_corrected_pvalues(p_vals_corrected)

        if self.json_output:
            if rows is None:
                rows = self._build_rows(stat_res, p_vals_corrected)
            payload = dict(rows=rows, sites=rows)
            if self.plot:
                payload["plot_output"] = self.plot_output
            print_json(payload)
            return

        lines = []
        for idx, (stat_info, pval_cor) in enumerate(
            zip(stat_res, p_vals_corrected), start=1
        ):
            pval_cor_str = "nan" if isinstance(pval_cor, str) else round(float(pval_cor), 4)
            raw_p = float(stat_info.pvalue)
            p_val_str = "nan" if isnan(raw_p) else round(raw_p, 4)
            lines.append(
                f"{idx}\t{round(float(stat_info.statistic), 4)}\t{pval_cor_str}\t{p_val_str}"
            )
        if lines:
            print("\n".join(lines))

        if self.plot:
            print(f"Saved compositional bias plot: {self.plot_output}")

    def process_args(self, args) -> dict[str, str]:
        from ...helpers.plot_config import PlotConfig

        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
            plot=getattr(args, "plot", False),
            plot_output=getattr(args, "plot_output", "compositional_bias_per_site_plot.png"),
            plot_config=PlotConfig.from_args(args),
        )

    def _build_rows(
        self,
        stat_res: list[Power_divergenceResult],
        p_vals_corrected: list[float | str],
    ) -> list[dict[str, int | float | None]]:
        rows = []
        append_row = rows.append
        for idx, (stat_info, pval_cor) in enumerate(
            zip(stat_res, p_vals_corrected), start=1
        ):
            corrected = None if isinstance(pval_cor, str) else round(float(pval_cor), 4)
            raw_p = float(stat_info.pvalue)
            append_row(
                {
                    "site": idx,
                    "chi_square": round(float(stat_info.statistic), 4),
                    "p_value_corrected": corrected,
                    "p_value": None if isnan(raw_p) else round(raw_p, 4),
                }
            )
        return rows

    def _plot_compositional_bias_manhattan(
        self,
        rows: list[dict[str, int | float | None]],
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

        sites, corrected_pvals = _prepare_compositional_bias_row_plot_series(rows)
        self._render_compositional_bias_manhattan(sites, corrected_pvals, plt, alpha)

    def _plot_compositional_bias_corrected_pvalues(
        self,
        p_vals_corrected: list[float | str],
        alpha: float = 0.05,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in compositional_bias_per_site. Install matplotlib and retry.")
            raise SystemExit(2)

        if not p_vals_corrected:
            return

        sites, corrected_pvals = _prepare_compositional_bias_plot_series(p_vals_corrected)
        self._render_compositional_bias_manhattan(sites, corrected_pvals, plt, alpha)

    def _render_compositional_bias_manhattan(
        self,
        sites,
        corrected_pvals,
        plt,
        alpha: float,
    ) -> None:
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
        if nonsig_mask.any():
            ax.scatter(
                sites[nonsig_mask],
                y_vals[nonsig_mask],
                s=10,
                alpha=0.75,
                color=colors[0],
                edgecolors="none",
                label="Not significant",
            )
        if sig_mask.any():
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
        if finite_mask.any():
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
        alignment: "MultipleSeqAlignment",
        idx: int,
        is_protein: bool = False,
    ) -> list[int]:
        gap_chars = set(self.get_gap_chars(is_protein))
        counts = {}
        for record in alignment:
            char = record.seq[idx].upper()
            if char not in gap_chars:
                counts[char] = counts.get(char, 0) + 1

        return list(counts.values())

    def calculate_compositional_bias_per_site(
        self,
        alignment: "MultipleSeqAlignment",
        is_protein: bool = False,
    ) -> tuple[
        list[Power_divergenceResult],
        list[float | str],
    ]:
        aln_len = alignment.get_alignment_length()
        num_records = len(alignment)

        if num_records == 0:
            return [], []

        raw_sequences = [str(record.seq) for record in alignment]
        first_raw_sequence = raw_sequences[0]
        first_sequence = first_raw_sequence.upper()
        for idx in range(1, num_records):
            sequence = raw_sequences[idx]
            if sequence != first_raw_sequence and sequence.upper() != first_sequence:
                break
        else:
            stat_res = [
                Power_divergenceResult(0.0, float("nan"))
                for _ in range(aln_len)
            ]
            return stat_res, ["nan"] * aln_len
        sequences = [sequence.upper() for sequence in raw_sequences]

        gap_chars = {char.upper() for char in self.get_gap_chars(is_protein)}
        valid_mask = None
        try:
            alignment_bytes = "".join(sequences).encode("ascii")
            alignment_array = np.frombuffer(
                alignment_bytes,
                dtype=np.uint8,
            ).reshape(len(sequences), aln_len)
            gap_lookup = _get_gap_lookup(is_protein)
            if is_protein:
                if any(alignment_bytes.find(gap_code) != -1 for gap_code in b"-?*X"):
                    valid_mask = ~gap_lookup[alignment_array]
                    valid_symbols = np.unique(alignment_array[valid_mask])
                else:
                    valid_symbols = np.unique(alignment_array)
            else:
                observed_symbols = np.unique(alignment_array)
                valid_symbols = observed_symbols[~gap_lookup[observed_symbols]]
                if valid_symbols.size > 8 and valid_symbols.size != observed_symbols.size:
                    valid_mask = ~gap_lookup[alignment_array]
        except UnicodeEncodeError:
            alignment_array = np.array([list(seq) for seq in sequences], dtype="U1")
            gap_chars_array = np.array(list(gap_chars), dtype="U1")
            valid_mask = ~np.isin(alignment_array, gap_chars_array)
            valid_symbols = np.unique(alignment_array[valid_mask])

        if valid_symbols.size == 1:
            stat_res = [
                Power_divergenceResult(0.0, float("nan"))
                for _ in range(aln_len)
            ]
            return stat_res, ["nan"] * aln_len

        if valid_symbols.size:
            if alignment_array.dtype == np.uint8 and valid_symbols.size > 8:
                (
                    category_counts,
                    totals,
                    sum_squares,
                ) = _column_count_stats_from_ascii_codes(
                    alignment_array,
                    valid_mask,
                    valid_symbols,
                )
            else:
                counts = np.array(
                    [
                        np.sum(alignment_array == symbol, axis=0)
                        for symbol in valid_symbols
                    ],
                    dtype=np.float64,
                )
                category_counts = np.count_nonzero(counts, axis=0)
                totals = _column_totals(counts)
                sum_squares = _column_sum_squares(counts)
            statistics = np.divide(
                category_counts * sum_squares,
                totals,
                out=np.zeros(aln_len, dtype=np.float64),
                where=totals > 0,
            ) - totals
            p_values = np.full(aln_len, np.nan, dtype=np.float64)
            valid_pvalue_mask = category_counts > 1
            p_values[valid_pvalue_mask] = _chi2_sf(
                statistics[valid_pvalue_mask],
                category_counts[valid_pvalue_mask] - 1,
            )
        else:
            statistics = np.zeros(aln_len, dtype=np.float64)
            p_values = np.full(aln_len, np.nan, dtype=np.float64)

        stat_res = _power_divergence_results_from_arrays(statistics, p_values)

        # Apply FDR correction
        nan_mask = np.isnan(p_values)
        p_vals = _pvalues_for_fdr(p_values, nan_mask)
        if p_vals:
            p_vals_corrected = _false_discovery_control(p_vals)
        else:
            p_vals_corrected = []

        return stat_res, _restore_nan_corrected_pvalues(p_vals_corrected, nan_mask)
