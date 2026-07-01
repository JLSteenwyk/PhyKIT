from .base import Alignment


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()
_DNA_GAP_CODES = None
_PROTEIN_GAP_CODES = None


def _get_gap_codes(is_protein: bool):
    global _DNA_GAP_CODES, _PROTEIN_GAP_CODES
    if is_protein:
        if _PROTEIN_GAP_CODES is None:
            _PROTEIN_GAP_CODES = np.frombuffer(b"-?*X", dtype=np.uint8)
        return _PROTEIN_GAP_CODES

    if _DNA_GAP_CODES is None:
        _DNA_GAP_CODES = np.frombuffer(b"-?*XN", dtype=np.uint8)
    return _DNA_GAP_CODES


def _entropy_from_ascii_codes(
    alignment_array,
    valid_mask,
    valid_chars,
    block_size: int = 8192,
    all_valid: bool = False,
) -> list[float]:
    aln_len = alignment_array.shape[1]
    entropies = np.zeros(aln_len, dtype=np.float64)
    symbol_codes = valid_chars.astype(np.intp, copy=False)

    for start in range(0, aln_len, block_size):
        stop = min(start + block_size, aln_len)
        width = stop - start
        block = alignment_array[:, start:stop]

        block_by_site = block.T.astype(np.intp, copy=False)
        site_offsets = np.arange(width, dtype=np.intp)[:, None] * 256
        if all_valid:
            count_index = (block_by_site + site_offsets).ravel()
        else:
            block_valid = valid_mask[:, start:stop]
            if not block_valid.any():
                continue
            count_index = (block_by_site + site_offsets)[block_valid.T]
        counts = np.bincount(
            count_index,
            minlength=width * 256,
        ).reshape(width, 256)[:, symbol_codes].T.astype(np.float64, copy=False)
        totals = counts.sum(axis=0)

        entropies[start:stop] = _entropy_from_counts(counts, totals)

    return _entropy_values_to_list(entropies)


def _entropy_values_to_list(entropies):
    return entropies.tolist()


def _entropy_columns_from_probabilities(probs, log_probs):
    return -np.einsum("ij,ij->j", probs, log_probs)


def _entropy_from_counts(counts, totals):
    with np.errstate(divide="ignore", invalid="ignore"):
        probs = np.divide(
            counts,
            totals,
            out=np.zeros_like(counts, dtype=np.float64),
            where=totals > 0,
        )
        positive = probs > 0
        log_probs = np.zeros_like(probs, dtype=np.float64)
        np.log2(probs, out=log_probs, where=positive)

    entropies = _entropy_columns_from_probabilities(probs, log_probs)
    entropies[totals == 0] = 0.0
    return entropies


def _prepare_entropy_plot_series(entropies):
    count = len(entropies)
    sites = np.arange(1, count + 1, dtype=np.int32)
    entropy_values = np.fromiter(
        (round(entropy, 4) for entropy in entropies),
        dtype=np.float64,
        count=count,
    )
    return sites, entropy_values


def _prepare_entropy_row_plot_series(rows):
    sites = np.array([int(row["site"]) for row in rows], dtype=np.int32)
    entropies = np.array([float(row["entropy"]) for row in rows], dtype=np.float64)
    return sites, entropies


class AlignmentEntropy(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"], verbose=parsed["verbose"])
        self.json_output = parsed["json_output"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()
        entropies = self.calculate_site_entropies(alignment, is_protein)
        rows = None

        if self.plot and self.verbose:
            rows = [
                dict(site=idx, entropy=round(entropy, 4))
                for idx, entropy in enumerate(entropies, start=1)
            ]
            self._plot_alignment_entropy(rows)
        elif self.plot:
            self._plot_alignment_entropy_values(entropies)

        if self.json_output:
            if self.verbose:
                if rows is None:
                    rows = [
                        dict(site=idx, entropy=round(entropy, 4))
                        for idx, entropy in enumerate(entropies, start=1)
                    ]
                payload = dict(
                    verbose=True,
                    rows=rows,
                    sites=rows,
                )
            else:
                payload = dict(
                    verbose=False,
                    mean_entropy=round(sum(entropies) / len(entropies), 4) if entropies else 0.0,
                )
            if self.plot:
                payload["plot_output"] = self.plot_output
            print_json(payload)
            return

        if self.verbose:
            if rows is None:
                lines = [
                    f"{idx}\t{round(entropy, 4)}"
                    for idx, entropy in enumerate(entropies, start=1)
                ]
            else:
                lines = [
                    f"{row['site']}\t{row['entropy']}"
                    for row in rows
                ]
            if lines:
                print("\n".join(lines))
        else:
            if entropies:
                print(round(sum(entropies) / len(entropies), 4))
            else:
                print(0.0)

        if self.plot:
            print(f"Saved alignment entropy plot: {self.plot_output}")

    def process_args(self, args) -> dict[str, str]:
        from ...helpers.plot_config import PlotConfig

        return dict(
            alignment_file_path=args.alignment,
            verbose=args.verbose,
            json_output=getattr(args, "json", False),
            plot=getattr(args, "plot", False),
            plot_output=getattr(args, "plot_output", "alignment_entropy_plot.png"),
            plot_config=PlotConfig.from_args(args),
        )

    def _plot_alignment_entropy(self, rows: list[dict[str, float]]) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in alignment_entropy. Install matplotlib and retry.")
            raise SystemExit(2)

        if not rows:
            return

        sites, entropies = _prepare_entropy_row_plot_series(rows)
        self._render_alignment_entropy_plot(sites, entropies, plt)

    def _plot_alignment_entropy_values(self, entropies) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in alignment_entropy. Install matplotlib and retry.")
            raise SystemExit(2)

        if not entropies:
            return

        sites, entropies = _prepare_entropy_plot_series(entropies)
        self._render_alignment_entropy_plot(sites, entropies, plt)

    def _render_alignment_entropy_plot(self, sites, entropies, plt) -> None:
        config = self.plot_config
        config.resolve(n_rows=len(sites), n_cols=None)
        default_colors = ["#2b8cbe"]
        colors = config.merge_colors(default_colors)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        ax.plot(sites, entropies, color=colors[0], linewidth=1.2, alpha=0.9)
        ax.scatter(sites, entropies, s=8, color=colors[0], alpha=0.75, edgecolors="none")
        ax.set_xlabel("Alignment site")
        ax.set_ylabel("Shannon entropy")
        ax.set_xlim(1, int(np.max(sites)))
        max_entropy = float(np.max(entropies)) if entropies.size else 0.0
        ax.set_ylim(0, max(1.0, max_entropy * 1.1))

        if config.show_title:
            ax.set_title(config.title or "Alignment Entropy Per Site", fontsize=config.title_fontsize)
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.tight_layout()
        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def calculate_site_entropies(self, alignment, is_protein: bool) -> list[float]:
        sequences = [str(record.seq).upper() for record in alignment]
        if not sequences:
            return []

        aln_len = len(sequences[0])
        first_sequence = sequences[0]
        for sequence in sequences:
            if sequence != first_sequence:
                break
        else:
            return [0.0] * aln_len

        try:
            alignment_array = np.frombuffer(
                "".join(sequences).encode("ascii"),
                dtype=np.uint8,
            ).reshape(len(sequences), aln_len)
            invalid_mask = np.zeros(alignment_array.shape, dtype=np.bool_)
            for gap_code in _get_gap_codes(is_protein):
                invalid_mask |= alignment_array == gap_code
            valid_mask = ~invalid_mask
        except UnicodeEncodeError:
            invalid_chars = {char.upper() for char in self.get_gap_chars(is_protein)}
            alignment_array = np.array([list(seq) for seq in sequences], dtype="U1")
            invalid_chars_array = np.array(list(invalid_chars), dtype="U1")
            valid_mask = ~np.isin(alignment_array, invalid_chars_array)

        if alignment_array.size == 0:
            return []

        valid_chars = np.unique(alignment_array[valid_mask])
        if valid_chars.size == 0:
            return [0.0] * alignment_array.shape[1]
        if valid_chars.size == 1:
            return [0.0] * alignment_array.shape[1]

        if alignment_array.dtype == np.uint8 and valid_chars.size > 8:
            return _entropy_from_ascii_codes(
                alignment_array,
                valid_mask,
                valid_chars,
                all_valid=bool(valid_mask.all()),
            )

        count_reducer = (
            np.count_nonzero if alignment_array.dtype == np.uint8 else np.sum
        )
        counts = np.vstack(
            [
                count_reducer(alignment_array == char, axis=0)
                for char in valid_chars
            ]
        ).astype(np.float64)
        totals = np.sum(counts, axis=0)

        entropies = _entropy_from_counts(counts, totals)
        return [float(entropy) for entropy in entropies]
