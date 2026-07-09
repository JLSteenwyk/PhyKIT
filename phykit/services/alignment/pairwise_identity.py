from __future__ import annotations

import itertools
import sys
import os
from math import sqrt

from .base import Alignment


class _LazyMultiprocessing:
    _module = None

    def _load(self):
        module = self._module
        if module is None:
            import multiprocessing as module

            self._module = module
        return module

    def cpu_count(self):
        return self._load().cpu_count()

    def Pool(self, *args, **kwargs):
        return self._load().Pool(*args, **kwargs)


mp = _LazyMultiprocessing()


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
_SQUAREFORM = None
_PAIRWISE_IDENTITY_SCALAR_STATS_MAX_CELLS = 8192
_PAIRWISE_IDENTITY_SCALAR_RESULT_MAX_CELLS = 8192


def squareform(*args, **kwargs):
    global _SQUAREFORM
    if _SQUAREFORM is None:
        from scipy.spatial.distance import squareform as _squareform

        _SQUAREFORM = _squareform

    return _SQUAREFORM(*args, **kwargs)


def calculate_summary_statistics_from_dict(*args, **kwargs):
    from ...helpers.stats_summary import (
        calculate_summary_statistics_from_dict as _calculate_summary_statistics_from_dict,
    )

    return _calculate_summary_statistics_from_dict(*args, **kwargs)


def calculate_summary_statistics_from_arr(*args, **kwargs):
    from ...helpers.stats_summary import (
        calculate_summary_statistics_from_arr as _calculate_summary_statistics_from_arr,
    )

    return _calculate_summary_statistics_from_arr(*args, **kwargs)


def print_summary_statistics(stats):
    try:
        print(
            (
                "mean: %s\n"
                "median: %s\n"
                "25th percentile: %s\n"
                "75th percentile: %s\n"
                "minimum: %s\n"
                "maximum: %s\n"
                "standard deviation: %s\n"
                "variance: %s"
            )
            % (
                round(stats["mean"], 4),
                round(stats["median"], 4),
                round(stats["twenty_fifth"], 4),
                round(stats["seventy_fifth"], 4),
                round(stats["minimum"], 4),
                round(stats["maximum"], 4),
                round(stats["standard_deviation"], 4),
                round(stats["variance"], 4),
            )
        )
    except BrokenPipeError:
        pass


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def tqdm(iterable, *args, **kwargs):
    try:
        from tqdm import tqdm as _tqdm
    except ImportError:
        return iterable
    return _tqdm(iterable, *args, **kwargs)


def _pair_ids_follow_taxa_upper_triangle(taxa, pair_ids) -> bool:
    n_taxa = len(taxa)
    expected_pairs = n_taxa * (n_taxa - 1) // 2
    if len(pair_ids) != expected_pairs:
        return False

    cursor = 0
    for idx1 in range(n_taxa - 1):
        left_taxon = taxa[idx1]
        for idx2 in range(idx1 + 1, n_taxa):
            pair = pair_ids[cursor]
            if pair[0] != left_taxon or pair[1] != taxa[idx2]:
                return False
            cursor += 1
    return True


def _pairwise_identity_matrix_from_pairs(taxa, pair_ids, pairwise_identities):
    n_taxa = len(taxa)

    if _pair_ids_follow_taxa_upper_triangle(taxa, pair_ids):
        values = np.fromiter(
            pairwise_identities.values(),
            dtype=np.float32,
            count=len(pair_ids),
        )
        matrix = squareform(values, checks=False)
        np.fill_diagonal(matrix, 1.0)
        return matrix

    matrix = np.ones((n_taxa, n_taxa), dtype=np.float32)
    taxon_to_index = {taxon: idx for idx, taxon in enumerate(taxa)}
    for pair, identity in zip(pair_ids, pairwise_identities.values()):
        idx1 = taxon_to_index[pair[0]]
        idx2 = taxon_to_index[pair[1]]
        matrix[idx1, idx2] = identity
        matrix[idx2, idx1] = identity
    return matrix


def _constant_identity_stats(identity: float, n_pairs: int):
    if n_pairs < 2:
        return calculate_summary_statistics_from_arr(
            np.full(n_pairs, identity, dtype=np.float64)
        )
    return dict(
        mean=identity,
        median=identity,
        twenty_fifth=identity,
        seventy_fifth=identity,
        minimum=identity,
        maximum=identity,
        standard_deviation=0.0,
        variance=0.0,
    )


def _empty_pairwise_identity_stats():
    return calculate_summary_statistics_from_dict({})


def _alignment_size(alignment):
    try:
        return len(alignment)
    except TypeError:
        return None


def _identity_for_identical_sequence(
    first_sequence: str,
    is_protein: bool,
    exclude_gaps: bool,
) -> float:
    aln_len = len(first_sequence)
    if aln_len == 0:
        return 0.0
    if not exclude_gaps:
        return 1.0

    if not is_protein:
        try:
            return (
                len(first_sequence.encode("ascii").translate(None, b"-?*XN"))
                / float(aln_len)
            )
        except UnicodeEncodeError:
            pass

    gap_chars = "-?*X" if is_protein else "-?*XN"
    nongap_count = aln_len
    for gap_char in gap_chars:
        nongap_count -= first_sequence.count(gap_char)
    return nongap_count / float(aln_len)


def _all_sequences_identical(sequences: list[str]) -> bool:
    empty = object()
    iterator = iter(sequences)
    first_sequence = next(iterator, empty)
    if first_sequence is empty:
        return True
    for sequence in iterator:
        if sequence != first_sequence:
            return False
    return True


def _linear_percentile(sorted_values, position):
    lower = int(position)
    upper = lower + 1
    if upper >= len(sorted_values):
        return sorted_values[lower]
    fraction = position - lower
    if fraction == 0:
        return sorted_values[lower]
    return sorted_values[lower] * (1.0 - fraction) + sorted_values[upper] * fraction


def _summary_statistics_from_small_values(values: list[float]) -> dict[str, float] | None:
    count = len(values)
    if count < 2:
        return calculate_summary_statistics_from_arr(values)
    sorted_values = sorted(values)
    minimum = sorted_values[0]
    maximum = sorted_values[-1]
    if minimum == maximum:
        return dict(
            mean=minimum,
            median=minimum,
            twenty_fifth=minimum,
            seventy_fifth=minimum,
            minimum=minimum,
            maximum=maximum,
            standard_deviation=0.0,
            variance=0.0,
        )

    mean = sum(sorted_values) / count
    median = _linear_percentile(sorted_values, (count - 1) * 0.5)
    twenty_fifth = _linear_percentile(sorted_values, (count - 1) * 0.25)
    seventy_fifth = _linear_percentile(sorted_values, (count - 1) * 0.75)
    sum_squared_deviations = 0.0
    for value in sorted_values:
        delta = value - mean
        sum_squared_deviations += delta * delta
    variance = sum_squared_deviations / (count - 1)
    return dict(
        mean=mean,
        median=median,
        twenty_fifth=twenty_fifth,
        seventy_fifth=seventy_fifth,
        minimum=minimum,
        maximum=maximum,
        standard_deviation=sqrt(variance),
        variance=variance,
    )


def _pairwise_identity_stats_scalar(records, is_protein: bool, exclude_gaps: bool):
    raw_sequences = []
    total_cells = 0
    aln_len = None
    for record in records:
        sequence = str(record.seq)
        if aln_len is None:
            aln_len = len(sequence)
        elif len(sequence) != aln_len:
            return None
        total_cells += len(sequence)
        if total_cells > _PAIRWISE_IDENTITY_SCALAR_STATS_MAX_CELLS:
            return None
        raw_sequences.append(sequence.upper())

    num_records = len(raw_sequences)
    if num_records < 2:
        return None
    if aln_len is None:
        return None

    gap_chars = frozenset("-?*X" if is_protein else "-?*XN")
    denominator = float(aln_len)
    identities = []
    for idx1 in range(num_records - 1):
        seq_one = raw_sequences[idx1]
        for idx2 in range(idx1 + 1, num_records):
            seq_two = raw_sequences[idx2]
            matches = 0
            if exclude_gaps:
                for char_one, char_two in zip(seq_one, seq_two):
                    if (
                        char_one == char_two
                        and (char_one not in gap_chars or char_two not in gap_chars)
                    ):
                        matches += 1
            else:
                for char_one, char_two in zip(seq_one, seq_two):
                    if char_one == char_two:
                        matches += 1
            identities.append(matches / denominator if aln_len > 0 else 0.0)

    return _summary_statistics_from_small_values(identities)


def _pairwise_identities_scalar(records, is_protein: bool, exclude_gaps: bool):
    raw_records = []
    total_cells = 0
    aln_len = None
    for record in records:
        sequence = str(record.seq)
        if aln_len is None:
            aln_len = len(sequence)
        elif len(sequence) != aln_len:
            return None
        total_cells += len(sequence)
        if total_cells > _PAIRWISE_IDENTITY_SCALAR_RESULT_MAX_CELLS:
            return None
        raw_records.append((record.id, sequence.upper()))

    num_records = len(raw_records)
    if num_records < 2 or aln_len is None:
        return None

    gap_chars = frozenset("-?*X" if is_protein else "-?*XN")
    denominator = float(aln_len)
    pair_ids = []
    pairwise_identities = {}
    identities = []

    for idx1 in range(num_records - 1):
        left_id, seq_one = raw_records[idx1]
        for idx2 in range(idx1 + 1, num_records):
            right_id, seq_two = raw_records[idx2]
            matches = 0
            if exclude_gaps:
                for char_one, char_two in zip(seq_one, seq_two):
                    if (
                        char_one == char_two
                        and (char_one not in gap_chars or char_two not in gap_chars)
                    ):
                        matches += 1
            else:
                for char_one, char_two in zip(seq_one, seq_two):
                    if char_one == char_two:
                        matches += 1

            identity = matches / denominator if aln_len > 0 else 0.0
            pair_id = [left_id, right_id]
            pair_ids.append(pair_id)
            pairwise_identities[f"{left_id}-{right_id}"] = identity
            identities.append(identity)

    stats = _summary_statistics_from_small_values(identities)
    return pair_ids, pairwise_identities, stats


class PairwiseIdentity(Alignment):
    MP_MIN_PAIRS = 500_000
    MAX_MP_WORKERS = 8
    _DNA_GAP_LOOKUP = None
    _PROTEIN_GAP_LOOKUP = None

    @classmethod
    def _get_gap_lookup(cls, is_protein: bool):
        if is_protein:
            if cls._PROTEIN_GAP_LOOKUP is None:
                lookup = np.zeros(256, dtype=np.bool_)
                lookup[np.frombuffer("-?*X".encode("ascii"), dtype=np.uint8)] = True
                cls._PROTEIN_GAP_LOOKUP = lookup
            return cls._PROTEIN_GAP_LOOKUP

        if cls._DNA_GAP_LOOKUP is None:
            lookup = np.zeros(256, dtype=np.bool_)
            lookup[np.frombuffer("-?*XN".encode("ascii"), dtype=np.uint8)] = True
            cls._DNA_GAP_LOOKUP = lookup
        return cls._DNA_GAP_LOOKUP

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            alignment_file_path=parsed["alignment_file_path"],
            verbose=parsed["verbose"],
            exclude_gaps=parsed["exclude_gaps"],
        )
        self.json_output = parsed["json_output"]
        self.plot = parsed["plot"]
        self.plot_output = parsed["plot_output"]
        self.plot_config = parsed["plot_config"]

    def _should_use_multiprocessing(self, n_pairs: int) -> bool:
        if os.environ.get("PHYKIT_DISABLE_MP", "0") == "1":
            return False
        if os.environ.get("PHYKIT_FORCE_MP", "0") == "1":
            return True
        return n_pairs >= self.MP_MIN_PAIRS

    def run(self):
        alignment, _, is_protein = self.get_alignment_and_format()
        taxa = [record.id for record in alignment] if self.plot else []

        summary_only = not self.verbose and not self.plot
        if summary_only:
            stats = self.calculate_pairwise_identity_stats(
                alignment, self.exclude_gaps, is_protein
            )
            pair_ids = []
            pairwise_identities = {}
        else:
            pair_ids, pairwise_identities, stats = \
                self.calculate_pairwise_identities(
                    alignment, self.exclude_gaps, is_protein
                )

        if self.plot:
            self._plot_pairwise_identity_heatmap(taxa, pair_ids, pairwise_identities)

        if self.json_output:
            self._print_json_output(pair_ids, pairwise_identities, stats)
            return

        if self.verbose:
            try:
                lines = [
                    f"{pair[0]}\t{pair[1]}\t{round(identity, 4)}"
                    for pair, identity in zip(
                        pair_ids, pairwise_identities.values()
                    )
                ]
                if lines:
                    print("\n".join(lines))
            except BrokenPipeError:
                pass
        else:
            print_summary_statistics(stats)

        if self.plot:
            print(f"Saved pairwise identity heatmap: {self.plot_output}")

    def process_args(self, args) -> dict[str, str]:
        plot = getattr(args, "plot", False)
        plot_config = None
        if plot:
            from ...helpers.plot_config import PlotConfig

            plot_config = PlotConfig.from_args(args)

        return dict(
            alignment_file_path=args.alignment,
            verbose=args.verbose,
            exclude_gaps=args.exclude_gaps,
            json_output=getattr(args, "json", False),
            plot=plot,
            plot_output=getattr(args, "plot_output", "pairwise_identity_heatmap.png"),
            plot_config=plot_config,
        )

    def _print_json_output(
        self,
        pair_ids: list[list[str]],
        pairwise_identities: dict[str, float],
        stats: dict[str, float],
    ) -> None:
        payload = dict(verbose=self.verbose, exclude_gaps=self.exclude_gaps)
        if self.verbose:
            rows = [
                {
                    "taxon_a": pair[0],
                    "taxon_b": pair[1],
                    "identity": round(identity, 4),
                }
                for pair, identity in zip(pair_ids, pairwise_identities.values())
            ]
            payload["rows"] = rows
            payload["pairs"] = rows
        else:
            payload["summary"] = stats
        if self.plot:
            payload["plot_output"] = self.plot_output
        print_json(payload)

    def _plot_pairwise_identity_heatmap(
        self,
        taxa: list[str],
        pair_ids: list[list[str]],
        pairwise_identities: dict[str, float],
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for --plot in pairwise_identity. Install matplotlib and retry.")
            raise SystemExit(2)

        n_taxa = len(taxa)
        if n_taxa == 0:
            return

        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform

        matrix = _pairwise_identity_matrix_from_pairs(
            taxa,
            pair_ids,
            pairwise_identities,
        )

        if n_taxa >= 3:
            distance_matrix = np.clip(1.0 - matrix, 0.0, 1.0)
            np.fill_diagonal(distance_matrix, 0.0)
            condensed = squareform(distance_matrix, checks=False)
            order = leaves_list(linkage(condensed, method="average"))
        else:
            order = np.arange(n_taxa)

        ordered_matrix = matrix[np.ix_(order, order)]
        ordered_taxa = [taxa[idx] for idx in order]

        config = self.plot_config
        config.resolve(n_rows=n_taxa, n_cols=n_taxa)

        fig_w = config.fig_width or max(6, min(20, n_taxa * 0.35))
        fig_h = config.fig_height or fig_w
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        image = ax.imshow(ordered_matrix, cmap="viridis", vmin=0, vmax=1, interpolation="nearest")

        if config.ylabel_fontsize and config.ylabel_fontsize > 0:
            ax.set_xticks(np.arange(n_taxa))
            ax.set_yticks(np.arange(n_taxa))
            ax.set_xticklabels(ordered_taxa, rotation=90, fontsize=config.xlabel_fontsize or config.ylabel_fontsize)
            ax.set_yticklabels(ordered_taxa, fontsize=config.ylabel_fontsize)
        else:
            ax.set_xticks([])
            ax.set_yticks([])

        ax.set_xlabel("Taxa (clustered)")
        ax.set_ylabel("Taxa (clustered)")
        colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        colorbar.set_label("Identity")

        if config.show_title:
            ax.set_title(config.title or "Pairwise Identity Heatmap", fontsize=config.title_fontsize)
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.savefig(self.plot_output, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def _calculate_identity_vectorized(self, seq_arr1, seq_arr2, gap_mask=None, exclude_gaps=False):
        """Vectorized calculation of sequence identity."""
        matches = (seq_arr1 == seq_arr2)

        if exclude_gaps and gap_mask is not None:
            # Match original behavior: count identities when at least one doesn't have a gap
            # This matches the original "res_one not in gap_chars or res_two not in gap_chars"
            valid_for_identity = ~gap_mask[0] | ~gap_mask[1]
            identities = np.count_nonzero(matches & valid_for_identity)
        else:
            identities = np.count_nonzero(matches)

        # Total compared is always the full length (matching original behavior)
        total_compared = len(seq_arr1)

        return identities / total_compared if total_compared > 0 else 0

    @staticmethod
    def _sequence_to_array(sequence):
        upper_sequence = str(sequence).upper()
        try:
            return np.frombuffer(upper_sequence.encode("ascii"), dtype="S1")
        except UnicodeEncodeError:
            return np.array(list(upper_sequence), dtype="U1")

    @staticmethod
    def _gap_chars_array_for_sequence(seq_array, gap_chars):
        gap_chars_upper = [char.upper() for char in gap_chars]
        if getattr(seq_array.dtype, "kind", None) == "S":
            try:
                return np.array(
                    [char.encode("ascii") for char in gap_chars_upper],
                    dtype="S1",
                )
            except UnicodeEncodeError:
                pass
        return np.array(gap_chars_upper, dtype="U1")

    def _process_pair_batch(self, alignment_data, pair_indices, exclude_gaps, gap_chars):
        """Process a batch of sequence pairs."""
        results = []
        append_result = results.append
        calculate_identity = self._calculate_identity_vectorized
        gap_array_for_sequence = self._gap_chars_array_for_sequence
        for idx1, idx2 in pair_indices:
            record_one = alignment_data[idx1]
            record_two = alignment_data[idx2]
            seq_one = record_one['seq']
            seq_two = record_two['seq']

            if exclude_gaps:
                gap_mask1 = record_one.get('gap_mask')
                gap_mask2 = record_two.get('gap_mask')
                if gap_mask1 is None:
                    gap_chars_array = gap_array_for_sequence(seq_one, gap_chars)
                    gap_mask1 = np.isin(seq_one, gap_chars_array)
                if gap_mask2 is None:
                    gap_chars_array = gap_array_for_sequence(seq_two, gap_chars)
                    gap_mask2 = np.isin(seq_two, gap_chars_array)
                identity = calculate_identity(
                    seq_one, seq_two, (gap_mask1, gap_mask2), exclude_gaps
                )
            else:
                identity = calculate_identity(seq_one, seq_two)

            append_result({
                'pair_id': [record_one['id'], record_two['id']],
                'identity': identity
            })
        return results

    @staticmethod
    def _batched_pair_indices(num_records: int, chunk_size: int):
        pair_indices = itertools.combinations(range(num_records), 2)
        while True:
            chunk = list(itertools.islice(pair_indices, chunk_size))
            if not chunk:
                break
            yield chunk

    def _calculate_pairwise_identities_matrix(
        self,
        alignment: "MultipleSeqAlignment",
        is_protein: bool,
        exclude_gaps: bool,
    ):
        records = list(alignment)
        num_records = len(records)
        if num_records < 2:
            return None

        raw_sequences = [str(record.seq) for record in records]
        lengths = {len(sequence) for sequence in raw_sequences}
        if len(lengths) != 1:
            return None
        aln_len = lengths.pop()

        first_raw_sequence = raw_sequences[0]
        first_sequence = first_raw_sequence.upper()
        all_identical = True
        for idx in range(1, num_records):
            sequence = raw_sequences[idx]
            if sequence != first_raw_sequence and sequence.upper() != first_sequence:
                all_identical = False
                break

        if all_identical:
            identity = _identity_for_identical_sequence(
                first_sequence,
                is_protein,
                exclude_gaps,
            )
            pair_ids = []
            pairwise_identities = {}
            for record_one, record_two in itertools.combinations(records, 2):
                pair_id = [record_one.id, record_two.id]
                pair_ids.append(pair_id)
                pairwise_identities["-".join(pair_id)] = identity
            stats = _constant_identity_stats(identity, len(pair_ids))
            return pair_ids, pairwise_identities, stats
        sequences = [sequence.upper() for sequence in raw_sequences]

        try:
            joined_bytes = "".join(sequences).encode("ascii")
            sequence_matrix = np.frombuffer(
                joined_bytes,
                dtype=np.uint8,
            ).reshape(num_records, aln_len)
        except UnicodeEncodeError:
            return None

        if exclude_gaps:
            gap_bytes = b"-?*X" if is_protein else b"-?*XN"
            if any(code in joined_bytes for code in gap_bytes):
                gap_lookup = self._get_gap_lookup(is_protein)
                identity_counts = np.zeros((num_records, num_records), dtype=np.int32)
                for symbol in np.unique(sequence_matrix):
                    if gap_lookup[symbol]:
                        continue
                    symbol_mask = (sequence_matrix == symbol).astype(
                        np.int32,
                        copy=False,
                    )
                    identity_counts += symbol_mask @ symbol_mask.T
            else:
                exclude_gaps = False
        else:
            exclude_gaps = False

        if not exclude_gaps:
            target_bytes = 16 * 1024 * 1024
            block_size = max(
                1,
                min(64, target_bytes // max(1, num_records * max(1, aln_len))),
            )
            identity_counts = np.empty((num_records, num_records), dtype=np.int32)
            for start in range(0, num_records, block_size):
                stop = min(num_records, start + block_size)
                identity_counts[start:stop] = (
                    sequence_matrix[start:stop, None, :] == sequence_matrix[None, :, :]
                ).sum(axis=2, dtype=np.int32)

        condensed_counts = squareform(identity_counts, checks=False)
        if aln_len > 0:
            identities = condensed_counts.astype(np.float64, copy=False)
            identities /= float(aln_len)
        else:
            identities = np.zeros(len(condensed_counts), dtype=np.float64)

        pairwise_identities = {}
        pair_ids = []
        cursor = 0
        for idx1 in range(num_records - 1):
            left_id = records[idx1].id
            for idx2 in range(idx1 + 1, num_records):
                right_id = records[idx2].id
                pair_id = [left_id, right_id]
                pair_ids.append(pair_id)
                pairwise_identities[f"{left_id}-{right_id}"] = float(identities[cursor])
                cursor += 1

        stats = calculate_summary_statistics_from_arr(identities)
        return pair_ids, pairwise_identities, stats

    def _calculate_pairwise_identity_stats_matrix(
        self,
        alignment: "MultipleSeqAlignment",
        is_protein: bool,
        exclude_gaps: bool,
    ):
        records = list(alignment)
        num_records = len(records)
        if num_records < 2:
            return None

        raw_sequences = [str(record.seq) for record in records]
        lengths = {len(sequence) for sequence in raw_sequences}
        if len(lengths) != 1:
            return None
        aln_len = lengths.pop()

        first_raw_sequence = raw_sequences[0]
        first_sequence = first_raw_sequence.upper()
        all_identical = True
        for idx in range(1, num_records):
            sequence = raw_sequences[idx]
            if sequence != first_raw_sequence and sequence.upper() != first_sequence:
                all_identical = False
                break

        if all_identical:
            identity = _identity_for_identical_sequence(
                first_sequence,
                is_protein,
                exclude_gaps,
            )
            n_pairs = num_records * (num_records - 1) // 2
            return _constant_identity_stats(identity, n_pairs)
        sequences = [sequence.upper() for sequence in raw_sequences]

        try:
            joined_bytes = "".join(sequences).encode("ascii")
            sequence_matrix = np.frombuffer(
                joined_bytes,
                dtype=np.uint8,
            ).reshape(num_records, aln_len)
        except UnicodeEncodeError:
            return None

        if exclude_gaps:
            gap_bytes = b"-?*X" if is_protein else b"-?*XN"
            if any(code in joined_bytes for code in gap_bytes):
                gap_lookup = self._get_gap_lookup(is_protein)
                identity_counts = np.zeros((num_records, num_records), dtype=np.int32)
                for symbol in np.unique(sequence_matrix):
                    if gap_lookup[symbol]:
                        continue
                    symbol_mask = (sequence_matrix == symbol).astype(
                        np.int32,
                        copy=False,
                    )
                    identity_counts += symbol_mask @ symbol_mask.T
                condensed_counts = squareform(identity_counts, checks=False)
                if aln_len > 0:
                    identities = (
                        condensed_counts.astype(np.float64, copy=False)
                        / float(aln_len)
                    )
                else:
                    identities = np.zeros(len(condensed_counts), dtype=np.float64)
                return calculate_summary_statistics_from_arr(identities)

        target_bytes = 16 * 1024 * 1024
        block_size = max(
            1,
            min(64, target_bytes // max(1, num_records * max(1, aln_len))),
        )
        values = np.empty(num_records * (num_records - 1) // 2, dtype=np.float64)
        cursor = 0
        denominator = float(aln_len)
        col_indices = np.arange(num_records)[None, :]
        for start in range(0, num_records, block_size):
            stop = min(num_records, start + block_size)
            identity_counts = (
                sequence_matrix[start:stop, None, :] == sequence_matrix[None, :, :]
            ).sum(axis=2, dtype=np.int32)
            upper_mask = col_indices > np.arange(start, stop)[:, None]
            block_values = identity_counts[upper_mask]
            next_cursor = cursor + block_values.size
            if aln_len > 0:
                values[cursor:next_cursor] = block_values / denominator
            else:
                values[cursor:next_cursor] = 0.0
            cursor = next_cursor

        return calculate_summary_statistics_from_arr(values)

    def calculate_pairwise_identity_stats(
        self,
        alignment: "MultipleSeqAlignment",
        exclude_gaps: bool,
        is_protein: bool = False,
    ) -> dict[str, float]:
        alignment_size = _alignment_size(alignment)
        if alignment_size is not None and alignment_size < 2:
            return _empty_pairwise_identity_stats()

        if alignment_size is not None:
            scalar_stats = _pairwise_identity_stats_scalar(
                alignment,
                is_protein,
                exclude_gaps,
            )
            if scalar_stats is not None:
                return scalar_stats

        matrix_stats = self._calculate_pairwise_identity_stats_matrix(
            alignment,
            is_protein,
            exclude_gaps,
        )
        if matrix_stats is not None:
            return matrix_stats

        _, _, stats = self.calculate_pairwise_identities(
            alignment,
            exclude_gaps,
            is_protein,
        )
        return stats

    def calculate_pairwise_identities(
        self,
        alignment: "MultipleSeqAlignment",
        exclude_gaps: bool,
        is_protein: bool = False,
    ) -> tuple[list[list[str]], dict[str, float], dict[str, float]]:
        alignment_size = _alignment_size(alignment)
        if alignment_size is not None and alignment_size < 2:
            return [], {}, _empty_pairwise_identity_stats()

        if alignment_size is not None:
            scalar_result = _pairwise_identities_scalar(
                alignment,
                is_protein,
                exclude_gaps,
            )
            if scalar_result is not None:
                return scalar_result

        matrix_result = self._calculate_pairwise_identities_matrix(
            alignment,
            is_protein,
            exclude_gaps,
        )
        if matrix_result is not None:
            return matrix_result

        gap_chars = self.get_gap_chars(is_protein)

        # Convert sequences to numpy arrays for faster comparison
        alignment_data = []
        for record in alignment:
            seq_array = self._sequence_to_array(record.seq)
            data = {
                'id': record.id,
                'seq': seq_array
            }
            if exclude_gaps:
                gap_chars_array = self._gap_chars_array_for_sequence(seq_array, gap_chars)
                data['gap_mask'] = np.isin(seq_array, gap_chars_array)
            alignment_data.append(data)

        num_records = len(alignment_data)
        n_pairs = num_records * (num_records - 1) // 2

        pairwise_identities = {}
        pair_ids = []

        # For small/medium workloads, multiprocessing overhead dominates.
        if not self._should_use_multiprocessing(n_pairs):
            # Process all pairs without multiprocessing
            results = self._process_pair_batch(
                alignment_data,
                itertools.combinations(range(num_records), 2),
                exclude_gaps,
                gap_chars,
            )
            for result in results:
                pair_id = result['pair_id']
                pair_ids.append(pair_id)
                pairwise_identities["-".join(pair_id)] = result['identity']
        else:
            from functools import partial

            # Use multiprocessing for larger datasets
            num_workers = min(mp.cpu_count(), self.MAX_MP_WORKERS)
            chunk_size = max(1, n_pairs // (num_workers * 4))
            pair_chunks = self._batched_pair_indices(num_records, chunk_size)

            # Create partial function
            process_func = partial(
                self._process_pair_batch,
                alignment_data,
                exclude_gaps=exclude_gaps,
                gap_chars=gap_chars
            )

            # Process in parallel with progress bar
            with mp.Pool(processes=num_workers) as pool:
                # Only show progress bar if stderr is a tty (not redirected)
                if sys.stderr.isatty():
                    chunk_results = list(tqdm(
                        pool.imap(process_func, pair_chunks),
                        total=len(pair_chunks),
                        desc="Calculating pairwise identities",
                        unit="batch"
                    ))
                else:
                    chunk_results = pool.map(process_func, pair_chunks)

            # Combine results
            for chunk_result in chunk_results:
                for result in chunk_result:
                    pair_id = result['pair_id']
                    pair_ids.append(pair_id)
                    pairwise_identities["-".join(pair_id)] = result['identity']

        stats = calculate_summary_statistics_from_dict(pairwise_identities)

        return pair_ids, pairwise_identities, stats
