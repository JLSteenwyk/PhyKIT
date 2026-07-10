from __future__ import annotations

from collections import Counter

from .base import Alignment, _all_sequences_identical


class _LazyNumpy:
    _module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = _np
            self._module = module
        value = getattr(module, name)
        setattr(self, name, value)
        return value


np = _LazyNumpy()
_DNA_GAP_LOOKUP = None
_PROTEIN_GAP_LOOKUP = None
_GAP_DELETE_TABLES = {}
_PLOT_DIRECT_MAX_LIMIT = 100_000
_ASCII_COUNT_BLOCK_SIZE = 512
_ERPS_SCALAR_MAX_CELLS = 8192
_NUCLEOTIDE_ALPHABET_BYTES = b"ACGTU"


def _column_sum_squares(counts: np.ndarray) -> np.ndarray:
    return np.einsum("ij,ij->j", counts, counts, dtype=np.float64)


def _column_totals(counts: np.ndarray) -> np.ndarray:
    if counts.shape[1] <= 20000:
        return counts.sum(axis=0)
    return np.sum(counts, axis=0)


def _bounded_ascii_symbol_counts(alignment_array, valid_symbols):
    num_records = alignment_array.shape[0]
    if num_records <= 0xFFFF:
        count_dtype = np.uint16
    elif num_records <= 0xFFFFFFFF:
        count_dtype = np.uint32
    else:
        count_dtype = np.uint64

    return np.array(
        [
            (alignment_array == symbol).sum(axis=0, dtype=count_dtype)
            for symbol in valid_symbols
        ],
        dtype=np.float64,
    )


def _plot_max(values):
    if values.size <= _PLOT_DIRECT_MAX_LIMIT:
        return values.max()
    return np.max(values)


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


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


def _column_totals_and_sum_squares_from_ascii_codes(
    alignment_array,
    valid_mask,
    valid_symbols,
    block_size: int = _ASCII_COUNT_BLOCK_SIZE,
):
    aln_len = alignment_array.shape[1]
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
        totals[slc] = counts.sum(axis=0, dtype=np.float64)
        sum_squares[slc] = _column_sum_squares(counts)

    return totals, sum_squares


def _prepare_evolutionary_rate_plot_series(pic_values):
    count = len(pic_values)
    sites = np.arange(1, count + 1, dtype=np.int32)
    rates = np.fromiter(
        (round(value, 4) for value in pic_values),
        dtype=np.float64,
        count=count,
    )
    return sites, rates


def _prepare_evolutionary_rate_row_plot_series(rows):
    sites = np.array([row["site"] for row in rows], dtype=np.int32)
    rates = np.array([row["evolutionary_rate"] for row in rows], dtype=np.float64)
    return sites, rates


def _evolutionary_rate_per_site_scalar(
    raw_sequences: list[str],
    aln_len: int,
    is_protein: bool,
):
    gap_chars = frozenset("-?*X" if is_protein else "-?*XN")
    counts_by_site = [dict() for _ in range(aln_len)]
    observed_valid = False

    for sequence in raw_sequences:
        if len(sequence) != aln_len:
            return None
        for idx, char in enumerate(sequence.upper()):
            if char in gap_chars:
                continue
            observed_valid = True
            counts = counts_by_site[idx]
            counts[char] = counts.get(char, 0) + 1

    if not observed_valid:
        return [0] * aln_len

    pic_values = []
    for counts in counts_by_site:
        if not counts:
            pic_values.append(0.0)
            continue
        total = 0
        sum_squares = 0
        for count in counts.values():
            total += count
            sum_squares += count * count
        pic_values.append(1 - (sum_squares / (total * total)))
    return pic_values


class EvolutionaryRatePerSite(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.json_output = parsed["json_output"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def run(self):
        alignment, _, is_protein = self.get_alignment_and_format()
        pic_values = self.calculate_evolutionary_rate_per_site(alignment, is_protein)
        rows = None

        if self.plot and self.json_output:
            rows = [
                {"site": idx, "evolutionary_rate": round(value, 4)}
                for idx, value in enumerate(pic_values, start=1)
            ]
            self._plot_evolutionary_rate_per_site(rows)
        elif self.plot:
            self._plot_evolutionary_rate_values(pic_values)

        if self.json_output:
            if rows is None:
                rows = [
                    {"site": idx, "evolutionary_rate": round(value, 4)}
                    for idx, value in enumerate(pic_values, start=1)
                ]
            payload = dict(rows=rows, sites=rows)
            if self.plot:
                payload["plot_output"] = self.plot_output
            print_json(payload)
            return

        lines = [
            f"{idx}\t{round(value, 4)}"
            for idx, value in enumerate(pic_values, start=1)
        ]
        if lines:
            print("\n".join(lines))

        if self.plot:
            print(f"Saved evolutionary-rate plot: {self.plot_output}")

    def process_args(self, args):
        plot = getattr(args, "plot", False)
        plot_config = None
        if plot:
            from ...helpers.plot_config import PlotConfig

            plot_config = PlotConfig.from_args(args)

        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
            plot=plot,
            plot_output=getattr(args, "plot_output", "evolutionary_rate_per_site_plot.png"),
            plot_config=plot_config,
        )

    def _plot_evolutionary_rate_per_site(self, rows):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in evolutionary_rate_per_site. Install matplotlib and retry.")
            raise SystemExit(2)

        if not rows:
            return

        sites, rates = _prepare_evolutionary_rate_row_plot_series(rows)
        self._render_evolutionary_rate_plot(sites, rates, plt)

    def _plot_evolutionary_rate_values(self, pic_values):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in evolutionary_rate_per_site. Install matplotlib and retry.")
            raise SystemExit(2)

        if not pic_values:
            return

        sites, rates = _prepare_evolutionary_rate_plot_series(pic_values)
        self._render_evolutionary_rate_plot(sites, rates, plt)

    def _render_evolutionary_rate_plot(self, sites, rates, plt):
        config = self.plot_config
        config.resolve(n_rows=len(sites), n_cols=None)
        default_colors = ["#2b8cbe"]
        colors = config.merge_colors(default_colors)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        ax.plot(sites, rates, color=colors[0], linewidth=1.2, alpha=0.9)
        ax.scatter(sites, rates, s=8, color=colors[0], alpha=0.75, edgecolors="none")
        ax.set_xlabel("Alignment site")
        ax.set_ylabel("Evolutionary rate")
        ax.set_xlim(1, int(_plot_max(sites)))
        ax.set_ylim(0, max(1.0, float(_plot_max(rates) * 1.1)))

        if config.show_title:
            ax.set_title(config.title or "Evolutionary Rate Per Site", fontsize=config.title_fontsize)
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def remove_gap_characters(self, seq: str, gap_chars: list[str]) -> str:
        key = tuple(gap_chars)
        table = _GAP_DELETE_TABLES.get(key)
        if table is None:
            if not all(len(char) == 1 for char in gap_chars):
                return ''.join([char for char in seq if char not in gap_chars]).upper()
            table = str.maketrans("", "", "".join(gap_chars))
            _GAP_DELETE_TABLES[key] = table
        return seq.translate(table).upper()

    def get_number_of_occurrences_per_character(
        self,
        alignment: MultipleSeqAlignment,
        idx: int,
        gap_chars: list[str]
    ) -> dict[str, int]:
        gap_chars_set = set(gap_chars)
        counts = {}
        for record in alignment:
            char = record.seq[idx]
            if char not in gap_chars_set:
                char = char.upper()
                counts[char] = counts.get(char, 0) + 1
        return Counter(counts)

    def calculate_pic(
        self,
        num_occurrences: dict[str, int],
    ) -> float:
        if not num_occurrences:
            return 1

        total_frequencies = 0
        sum_squares = 0
        for frequency in num_occurrences.values():
            total_frequencies += frequency
            sum_squares += frequency * frequency

        return 1 - (sum_squares / (total_frequencies * total_frequencies))

    def calculate_evolutionary_rate_per_site(
        self,
        alignment: MultipleSeqAlignment,
        is_protein: bool = False,
    ) -> list[float]:
        aln_len = alignment.get_alignment_length()
        num_records = len(alignment)

        if num_records == 0:
            return []

        raw_sequences = [str(record.seq) for record in alignment]
        first_raw_sequence = raw_sequences[0]
        first_sequence = first_raw_sequence.upper()
        all_identical = True
        for idx in range(1, num_records):
            sequence = raw_sequences[idx]
            if sequence != first_raw_sequence and sequence.upper() != first_sequence:
                all_identical = False
                break

        if all_identical:
            return [0.0] * aln_len

        if num_records * aln_len <= _ERPS_SCALAR_MAX_CELLS:
            scalar_values = _evolutionary_rate_per_site_scalar(
                raw_sequences,
                aln_len,
                is_protein,
            )
            if scalar_values is not None:
                return scalar_values

        sequences = [sequence.upper() for sequence in raw_sequences]

        gap_chars = {char.upper() for char in self.get_gap_chars(is_protein)}
        ascii_matrix = False
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
                if not alignment_bytes.translate(
                    None,
                    _NUCLEOTIDE_ALPHABET_BYTES,
                ):
                    observed_symbols = np.fromiter(
                        (
                            code
                            for code in _NUCLEOTIDE_ALPHABET_BYTES
                            if code in alignment_bytes
                        ),
                        dtype=np.uint8,
                    )
                else:
                    observed_symbols = np.unique(alignment_array)
                valid_symbols = observed_symbols[~gap_lookup[observed_symbols]]
                if valid_symbols.size > 8 and valid_symbols.size != observed_symbols.size:
                    valid_mask = ~gap_lookup[alignment_array]
            ascii_matrix = True
        except UnicodeEncodeError:
            alignment_array = np.array([list(seq) for seq in sequences], dtype="U1")
            gap_chars_array = np.array(list(gap_chars), dtype="U1")
            valid_mask = ~np.isin(alignment_array, gap_chars_array)
            valid_symbols = np.unique(alignment_array[valid_mask])

        if valid_symbols.size == 0:
            return [0] * aln_len
        if valid_symbols.size == 1:
            return [0.0] * aln_len

        if ascii_matrix and valid_symbols.size > 8:
            totals, sum_squares = _column_totals_and_sum_squares_from_ascii_codes(
                alignment_array,
                valid_mask,
                valid_symbols,
            )
        elif ascii_matrix:
            counts = _bounded_ascii_symbol_counts(
                alignment_array,
                valid_symbols,
            )
            totals = _column_totals(counts)
            sum_squares = _column_sum_squares(counts)
        else:
            counts = np.array(
                [
                    (alignment_array == symbol).sum(axis=0)
                    for symbol in valid_symbols
                ],
                dtype=np.float64,
            )
            totals = _column_totals(counts)
            sum_squares = _column_sum_squares(counts)
        squared_frequency_sums = np.divide(
            sum_squares,
            totals * totals,
            out=np.ones(aln_len, dtype=np.float64),
            where=totals > 0,
        )
        pic_values = 1.0 - squared_frequency_sums
        return pic_values.tolist()
