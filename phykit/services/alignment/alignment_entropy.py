from .base import Alignment


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyNumpy:
    def __init__(self):
        self._module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = _np
            self._module = module

        attr = getattr(module, name)
        setattr(self, name, attr)
        return attr


np = _LazyNumpy()
_DNA_GAP_CODES = None
_PROTEIN_GAP_CODES = None
_PROTEIN_GAP_BYTES = b"-?*X"
_PLOT_DIRECT_MAX_LIMIT = 100_000
_ASCII_ENTROPY_BLOCK_SIZE = 512


def _get_gap_codes(is_protein: bool):
    global _DNA_GAP_CODES, _PROTEIN_GAP_CODES
    if is_protein:
        if _PROTEIN_GAP_CODES is None:
            _PROTEIN_GAP_CODES = np.frombuffer(_PROTEIN_GAP_BYTES, dtype=np.uint8)
        return _PROTEIN_GAP_CODES

    if _DNA_GAP_CODES is None:
        _DNA_GAP_CODES = np.frombuffer(b"-?*XN", dtype=np.uint8)
    return _DNA_GAP_CODES


def _entropy_from_ascii_codes(
    alignment_array,
    valid_mask,
    valid_chars,
    block_size: int = _ASCII_ENTROPY_BLOCK_SIZE,
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
        totals = _entropy_count_totals(counts)

        entropies[start:stop] = _entropy_from_counts(counts, totals)

    return _entropy_values_to_list(entropies)


def _entropy_values_to_list(entropies):
    return entropies.tolist()


def _entropy_count_totals(counts):
    if counts.shape[1] <= 20000 or counts.shape[0] >= 8:
        return counts.sum(axis=0)
    return np.sum(counts, axis=0)


def _entropy_columns_from_probabilities(probs, log_probs):
    return -np.einsum("ij,ij->j", probs, log_probs)


def _plot_max(values):
    if values.size <= _PLOT_DIRECT_MAX_LIMIT:
        return values.max()
    return np.max(values)


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
    count = len(rows)
    sites = np.fromiter(
        (int(row["site"]) for row in rows),
        dtype=np.int32,
        count=count,
    )
    entropies = np.fromiter(
        (float(row["entropy"]) for row in rows),
        dtype=np.float64,
        count=count,
    )
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
        mean_entropy = None
        if not self.verbose and not self.plot:
            mean_entropy = self._mean_entropy_if_identical(alignment)

        if mean_entropy is None:
            entropies = self.calculate_site_entropies(alignment, is_protein)
        else:
            entropies = None
        rows = None

        if self.plot and self.verbose:
            rows = [
                {"site": idx, "entropy": round(entropy, 4)}
                for idx, entropy in enumerate(entropies, start=1)
            ]
            self._plot_alignment_entropy(rows)
        elif self.plot:
            self._plot_alignment_entropy_values(entropies)

        if self.json_output:
            if self.verbose:
                if rows is None:
                    rows = [
                        {"site": idx, "entropy": round(entropy, 4)}
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
                    mean_entropy=round(
                        mean_entropy
                        if mean_entropy is not None
                        else (sum(entropies) / len(entropies) if entropies else 0.0),
                        4,
                    ),
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
            if mean_entropy is not None:
                print(round(mean_entropy, 4))
            elif entropies:
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

    @staticmethod
    def _mean_entropy_if_identical(alignment) -> float | None:
        try:
            record_count = len(alignment)
        except TypeError:
            return None

        if record_count == 0:
            return 0.0

        try:
            first_record = alignment[0]
            first_raw_sequence = str(first_record.seq)
        except (AttributeError, IndexError, TypeError):
            return None

        first_sequence = first_raw_sequence.upper()
        for sample_idx in {record_count // 2, record_count - 1}:
            if sample_idx == 0:
                continue
            try:
                sample_sequence = str(alignment[sample_idx].seq)
            except (AttributeError, IndexError, TypeError):
                return None
            if (
                sample_sequence != first_raw_sequence
                and sample_sequence.upper() != first_sequence
            ):
                return None

        for idx in range(1, record_count):
            sequence = str(alignment[idx].seq)
            if (
                sequence != first_raw_sequence
                and sequence.upper() != first_sequence
            ):
                return None
        return 0.0

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
        ax.set_xlim(1, int(_plot_max(sites)))
        max_entropy = float(_plot_max(entropies)) if entropies.size else 0.0
        ax.set_ylim(0, max(1.0, max_entropy * 1.1))

        if config.show_title:
            ax.set_title(config.title or "Alignment Entropy Per Site", fontsize=config.title_fontsize)
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def calculate_site_entropies(self, alignment, is_protein: bool) -> list[float]:
        raw_sequences = [str(record.seq) for record in alignment]
        if not raw_sequences:
            return []

        first_raw_sequence = raw_sequences[0]
        aln_len = len(first_raw_sequence)
        first_sequence = first_raw_sequence.upper()
        for idx in range(1, len(raw_sequences)):
            sequence = raw_sequences[idx]
            if sequence != first_raw_sequence and sequence.upper() != first_sequence:
                break
        else:
            return [0.0] * aln_len
        sequences = [sequence.upper() for sequence in raw_sequences]

        try:
            alignment_bytes = "".join(sequences).encode("ascii")
            alignment_array = np.frombuffer(
                alignment_bytes,
                dtype=np.uint8,
            ).reshape(len(sequences), aln_len)
            if is_protein and not any(
                gap_code in alignment_bytes for gap_code in _PROTEIN_GAP_BYTES
            ):
                valid_mask = None
                valid_chars = np.unique(alignment_array)
                all_valid = True
            else:
                invalid_mask = np.zeros(alignment_array.shape, dtype=np.bool_)
                for gap_code in _get_gap_codes(is_protein):
                    invalid_mask |= alignment_array == gap_code
                valid_mask = ~invalid_mask
                valid_chars = np.unique(alignment_array[valid_mask])
                all_valid = bool(valid_mask.all())
        except UnicodeEncodeError:
            invalid_chars = {char.upper() for char in self.get_gap_chars(is_protein)}
            alignment_array = np.array([list(seq) for seq in sequences], dtype="U1")
            invalid_chars_array = np.array(list(invalid_chars), dtype="U1")
            valid_mask = ~np.isin(alignment_array, invalid_chars_array)
            valid_chars = np.unique(alignment_array[valid_mask])
            all_valid = False

        if alignment_array.size == 0:
            return []

        if valid_chars.size == 0:
            return [0.0] * alignment_array.shape[1]
        if valid_chars.size == 1:
            return [0.0] * alignment_array.shape[1]

        if alignment_array.dtype == np.uint8 and valid_chars.size > 8:
            return _entropy_from_ascii_codes(
                alignment_array,
                valid_mask,
                valid_chars,
                all_valid=all_valid,
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
        totals = _entropy_count_totals(counts)

        entropies = _entropy_from_counts(counts, totals)
        return _entropy_values_to_list(entropies)
