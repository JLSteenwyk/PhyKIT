from .base import Alignment, _all_sequences_identical


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()
_DNA_INVALID_LOOKUP = None
_PROTEIN_INVALID_LOOKUP = None


def _get_invalid_lookup(is_protein: bool):
    global _DNA_INVALID_LOOKUP, _PROTEIN_INVALID_LOOKUP
    if is_protein:
        if _PROTEIN_INVALID_LOOKUP is None:
            lookup = np.zeros(256, dtype=np.bool_)
            lookup[np.frombuffer("-?*X".encode("ascii"), dtype=np.uint8)] = True
            _PROTEIN_INVALID_LOOKUP = lookup
        return _PROTEIN_INVALID_LOOKUP

    if _DNA_INVALID_LOOKUP is None:
        lookup = np.zeros(256, dtype=np.bool_)
        lookup[np.frombuffer("-?*XN".encode("ascii"), dtype=np.uint8)] = True
        _DNA_INVALID_LOOKUP = lookup
    return _DNA_INVALID_LOOKUP


class RelativeCompositionVariabilityTaxon(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.json_output = parsed["json_output"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def run(self):
        alignment, _, is_protein = self.get_alignment_and_format()
        rows = self.calculate_rows(alignment, is_protein)

        if self.json_output:
            payload = dict(rows=rows, taxa=rows)
            if self.plot:
                self._plot_rcvt(rows)
                payload["plot_output"] = self.plot_output
            print_json(payload)
            return

        lines = [f"{row['taxon']}\t{row['rcvt']}" for row in rows]
        if lines:
            print("\n".join(lines))

        if self.plot:
            self._plot_rcvt(rows)
            print(f"Saved RCVT plot: {self.plot_output}")

    def _plot_rcvt(self, rows):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in rcvt. Install matplotlib and retry.")
            raise SystemExit(2)

        taxa, rcvt_values = self._prepare_rcvt_plot_series(rows)

        config = self.plot_config
        config.resolve(n_rows=len(taxa), n_cols=None)
        default_colors = ["#2b8cbe"]
        colors = config.merge_colors(default_colors)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        ax.bar(np.arange(len(taxa)), rcvt_values, color=colors[0])
        ax.set_ylabel("RCVT")
        ax.set_xticks(np.arange(len(taxa)))
        if config.xlabel_fontsize and config.xlabel_fontsize > 0:
            ax.set_xticklabels(taxa, rotation=90, fontsize=config.xlabel_fontsize)
        else:
            ax.set_xticklabels([])

        if config.show_title:
            ax.set_title(config.title or "RCVT Per Taxon", fontsize=config.title_fontsize)
        if config.axis_fontsize:
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    @staticmethod
    def _prepare_rcvt_plot_series(rows):
        taxa = [row["taxon"] for row in rows]
        rcvt_values = np.fromiter(
            (row["rcvt"] for row in rows),
            dtype=np.float64,
            count=len(rows),
        )
        order = np.argsort(-rcvt_values, kind="mergesort")
        return [taxa[idx] for idx in order], rcvt_values[order]

    def calculate_rows(self, alignment, is_protein: bool):
        num_records = len(alignment)
        if num_records == 0:
            return []
        if num_records == 1:
            return [dict(taxon=alignment[0].id, rcvt=0.0)]

        raw_sequences = [str(record.seq) for record in alignment]
        aln_len = alignment.get_alignment_length()
        first_raw_sequence = raw_sequences[0]
        first_sequence = first_raw_sequence.upper()
        all_identical = True
        for idx in range(1, num_records):
            sequence = raw_sequences[idx]
            if sequence != first_raw_sequence and sequence.upper() != first_sequence:
                all_identical = False
                break

        if all_identical:
            return [dict(taxon=record.id, rcvt=0.0) for record in alignment]
        sequences = [sequence.upper() for sequence in raw_sequences]

        if is_protein:
            invalid_chars = ["-", "?", "*", "X"]
        else:
            invalid_chars = ["-", "?", "*", "X", "N"]

        try:
            alignment_array = np.frombuffer(
                "".join(sequences).encode("ascii"),
                dtype=np.uint8,
            ).reshape(num_records, aln_len)
            invalid_lookup = _get_invalid_lookup(is_protein)
            observed_chars = np.unique(alignment_array)
            unique_chars = observed_chars[~invalid_lookup[observed_chars]]
            if unique_chars.size == 0:
                return [dict(taxon=record.id, rcvt=0.0) for record in alignment]
            if unique_chars.size == observed_chars.size:
                valid_mask = None
                valid_lengths = np.full(num_records, aln_len, dtype=np.float64)
            else:
                valid_mask = ~invalid_lookup[alignment_array]
                valid_lengths = np.count_nonzero(valid_mask, axis=1).astype(np.float64)
        except UnicodeEncodeError:
            alignment_array = np.array([list(seq) for seq in sequences], dtype="U1")
            invalid_chars_array = np.array(invalid_chars, dtype="U1")
            valid_mask = ~np.isin(alignment_array, invalid_chars_array)
            valid_lengths = np.count_nonzero(valid_mask, axis=1).astype(np.float64)

            valid_chars = alignment_array[valid_mask]
            if valid_chars.size == 0:
                return [dict(taxon=record.id, rcvt=0.0) for record in alignment]

            unique_chars = np.unique(valid_chars)
        if alignment_array.dtype == np.uint8 and len(unique_chars) > 8:
            count_matrix = self._ascii_count_matrix(
                alignment_array,
                unique_chars,
                valid_mask,
            )
        else:
            count_matrix = np.array(
                [
                    np.sum(alignment_array == char, axis=1)
                    for char in unique_chars
                ],
                dtype=np.float32,
            ).T

        average_counts = count_matrix.sum(axis=0) / num_records
        deviations = np.abs(count_matrix - average_counts)
        seq_sums = np.sum(deviations, axis=1)
        denom = num_records * valid_lengths
        rcv_values = np.divide(
            seq_sums,
            denom,
            out=np.zeros_like(seq_sums, dtype=np.float64),
            where=denom > 0,
        )

        return [
            dict(taxon=record.id, rcvt=round(float(value), 4))
            for record, value in zip(alignment, rcv_values)
        ]

    @staticmethod
    def _ascii_count_matrix(alignment_array, unique_chars, valid_mask):
        num_records, aln_len = alignment_array.shape
        if valid_mask is None and num_records >= 10_000 and aln_len <= 256:
            encoded = alignment_array.astype(np.int64)
            encoded += (np.arange(num_records, dtype=np.int64) * 256)[:, None]
            counts = np.bincount(
                encoded.ravel(),
                minlength=num_records * 256,
            ).reshape(num_records, 256)
            return counts[:, unique_chars].astype(np.float32, copy=False)

        count_matrix = np.zeros(
            (num_records, len(unique_chars)),
            dtype=np.float32,
        )
        for row_idx, row in enumerate(alignment_array):
            row_values = row if valid_mask is None else row[valid_mask[row_idx]]
            count_matrix[row_idx] = np.bincount(
                row_values,
                minlength=256,
            )[unique_chars]
        return count_matrix

    def process_args(self, args):
        from ...helpers.plot_config import PlotConfig

        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
            plot=getattr(args, "plot", False),
            plot_output=getattr(args, "plot_output", "rcvt_plot.png"),
            plot_config=PlotConfig.from_args(args),
        )
