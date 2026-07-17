from __future__ import annotations

import json
from pathlib import Path

from ...errors import PhykitUserError
from .base import Alignment


TAPER_DOI = "10.1111/2041-210X.13696"
TAPER_CITATION = (
    "Zhang C, Zhao Y, Braun EL, and Mirarab S. 2021. TAPER: Pinpointing "
    "errors in multiple sequence alignments despite varying rates of evolution. "
    "Methods in Ecology and Evolution 12:2145-2158. "
    f"doi:{TAPER_DOI}"
)

# Published TAPER defaults: window size, across-taxon tail, within-sequence
# tail, and the maximum retained region length for that scale.
_DEFAULT_SCALES = (
    (5, 0.25, 0.10, 30),
    (9, 0.25, 0.25, 54),
    (17, 0.10, 0.50, None),
)
_SCORE_BLOCK_SIZE = 2_048


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


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def _ascii_valid_mask(alignment_array, is_protein: bool):
    invalid = b"-?*X" if is_protein else b"-?*XN"
    lookup = np.zeros(256, dtype=np.bool_)
    lookup[np.frombuffer(invalid, dtype=np.uint8)] = True
    return ~lookup[alignment_array]


def _column_divergence_scores(alignment_array, valid_mask):
    """Score valid cells as 1 / (unique symbols * symbol frequency)."""

    n_taxa, n_sites = alignment_array.shape
    scores = np.zeros((n_taxa, n_sites), dtype=np.float64)
    if n_taxa == 0 or n_sites == 0:
        return scores

    if alignment_array.dtype == np.uint8:
        for start in range(0, n_sites, _SCORE_BLOCK_SIZE):
            stop = min(start + _SCORE_BLOCK_SIZE, n_sites)
            width = stop - start
            block = alignment_array[:, start:stop].T
            block_valid = valid_mask[:, start:stop].T
            if not block_valid.any():
                continue

            offsets = np.arange(width, dtype=np.intp)[:, None] * 256
            encoded = block.astype(np.intp, copy=False) + offsets
            counts = np.bincount(
                encoded[block_valid],
                minlength=width * 256,
            ).reshape(width, 256)
            totals = counts.sum(axis=1, dtype=np.int64)
            unique_counts = np.count_nonzero(counts, axis=1)
            observed_counts = counts[
                np.arange(width, dtype=np.intp)[:, None],
                block,
            ]
            denominator = unique_counts[:, None] * observed_counts
            block_scores = np.divide(
                totals[:, None],
                denominator,
                out=np.zeros_like(observed_counts, dtype=np.float64),
                where=block_valid & (denominator > 0),
            )
            scores[:, start:stop] = block_scores.T
        return scores

    for site in range(n_sites):
        column_valid = valid_mask[:, site]
        symbols = alignment_array[column_valid, site]
        if symbols.size == 0:
            continue
        unique_symbols, counts = np.unique(symbols, return_counts=True)
        count_by_symbol = dict(zip(unique_symbols.tolist(), counts.tolist()))
        total = float(symbols.size)
        unique_count = float(unique_symbols.size)
        for row in np.flatnonzero(column_valid):
            scores[row, site] = total / (
                unique_count * count_by_symbol[alignment_array[row, site]]
            )
    return scores


def _jenks_two_class_cutoff(values) -> float:
    """Return the low-class ceiling for a two-class natural break."""

    ordered = np.sort(values)
    size = ordered.size
    if size == 0:
        return 0.0
    if size == 1 or ordered[0] == ordered[-1]:
        return float(ordered[-1])

    prefix = np.cumsum(ordered, dtype=np.float64)
    low_sizes = np.arange(1, size, dtype=np.float64)
    high_sizes = float(size) - low_sizes
    objective = (
        (prefix[:-1] * prefix[:-1]) / low_sizes
        + ((prefix[-1] - prefix[:-1]) ** 2) / high_sizes
    )
    split = int(np.argmax(objective))
    return float(ordered[split])


def _upper_tail_boundary(ordered_values, tail_fraction: float) -> float:
    if len(ordered_values) == 0:
        return 0.0
    index = len(ordered_values) - int(len(ordered_values) * tail_fraction) - 1
    index = max(0, min(index, len(ordered_values) - 1))
    return float(ordered_values[index])


def _segment_outlier_windows(outlier_windows, window_size: int):
    """Smooth window calls into residue regions with a two-state DP."""

    n_windows = len(outlier_windows)
    if n_windows == 0:
        return np.zeros(0, dtype=np.bool_)

    normal_scores = np.zeros(n_windows, dtype=np.int32)
    outlier_scores = np.zeros(n_windows, dtype=np.int32)
    normal_previous = np.zeros(n_windows, dtype=np.int8)
    outlier_previous = np.zeros(n_windows, dtype=np.int8)

    normal_scores[0] = int(not outlier_windows[0])
    outlier_scores[0] = int(outlier_windows[0])

    for index in range(1, n_windows):
        normal_match = int(not outlier_windows[index])
        outlier_match = int(outlier_windows[index])

        normal_scores[index] = normal_scores[index - 1] + normal_match
        outlier_scores[index] = outlier_scores[index - 1] + outlier_match

        if index >= window_size:
            switch_to_normal = (
                outlier_scores[index - window_size] + normal_match
            )
            if switch_to_normal > normal_scores[index]:
                normal_scores[index] = switch_to_normal
                normal_previous[index] = 1

            switch_to_outlier = (
                normal_scores[index - window_size] + outlier_match
            )
            if switch_to_outlier > outlier_scores[index]:
                outlier_scores[index] = switch_to_outlier
                outlier_previous[index] = 1

    residue_mask = np.zeros(n_windows + window_size - 1, dtype=np.bool_)
    state = 0 if normal_scores[-1] >= outlier_scores[-1] else 1
    index = n_windows - 1
    if state == 1:
        residue_mask[index:index + window_size] = True

    while index > 0:
        if state == 0 and normal_previous[index] == 0:
            index -= 1
        elif state == 0:
            index -= window_size
            state = 1
            residue_mask[index:index + window_size] = True
        elif outlier_previous[index] == 1:
            index -= window_size
            state = 0
        else:
            index -= 1
            residue_mask[index] = True

    if index == 0 and state == 1:
        residue_mask[0] = True

    return residue_mask


def _remove_long_regions(mask, maximum_length: int | None):
    if maximum_length is None or not mask.any():
        return mask

    filtered = mask.copy()
    padded = np.pad(mask.astype(np.int8), (1, 1))
    boundaries = np.flatnonzero(np.diff(padded))
    for start, stop in boundaries.reshape(-1, 2):
        if stop - start > maximum_length:
            filtered[start:stop] = False
    return filtered


def _detect_at_scale(
    residue_scores_by_taxon,
    window_size: int,
    taxon_tail: float,
    sequence_tail: float,
    score_cutoff: float,
    maximum_length: int | None,
):
    window_scores = []
    initial_cutoffs = []
    eligible_cutoffs = []

    for residue_scores in residue_scores_by_taxon:
        if residue_scores.size < (2 * window_size):
            window_scores.append(None)
            initial_cutoffs.append(None)
            continue
        windows = np.lib.stride_tricks.sliding_window_view(
            residue_scores,
            window_size,
        )
        medians = np.median(windows, axis=1)
        cutoff = _jenks_two_class_cutoff(medians)
        window_scores.append(medians)
        initial_cutoffs.append(cutoff)
        eligible_cutoffs.append(cutoff)

    if eligible_cutoffs:
        ordered_cutoffs = np.sort(np.asarray(eligible_cutoffs, dtype=np.float64))
        across_taxa_cutoff = _upper_tail_boundary(
            ordered_cutoffs,
            taxon_tail,
        )
    else:
        across_taxa_cutoff = score_cutoff

    masks = []
    final_cutoffs = []
    for residue_scores, medians, initial_cutoff in zip(
        residue_scores_by_taxon,
        window_scores,
        initial_cutoffs,
    ):
        if medians is None:
            masks.append(np.zeros(residue_scores.size, dtype=np.bool_))
            final_cutoffs.append(None)
            continue

        ordered = np.sort(medians)
        within_sequence_cutoff = _upper_tail_boundary(
            ordered,
            sequence_tail,
        )
        final_cutoff = max(
            float(initial_cutoff),
            across_taxa_cutoff,
            within_sequence_cutoff,
            score_cutoff,
        )
        outlier_windows = medians > final_cutoff
        residue_mask = _segment_outlier_windows(
            outlier_windows,
            window_size,
        )
        masks.append(_remove_long_regions(residue_mask, maximum_length))
        final_cutoffs.append(final_cutoff)

    return masks, final_cutoffs


class AlignmentOutlierRegions(Alignment):
    """Detect localized sequence outliers with TAPER-inspired 2D scoring."""

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.cutoff = parsed["cutoff"]
        self.mask_output = parsed["mask_output"]
        self.mask_character = parsed["mask_character"]
        self.report = parsed["report"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> dict[str, object]:
        return dict(
            alignment_file_path=args.alignment,
            cutoff=args.cutoff,
            mask_output=getattr(args, "mask_output", None),
            mask_character=getattr(args, "mask_character", "?"),
            report=getattr(args, "report", None),
            json_output=getattr(args, "json", False),
        )

    def _validate_options(self) -> None:
        messages = []
        if self.cutoff <= 1.0:
            messages.append("cutoff must be greater than 1.0.")
        if len(self.mask_character) != 1 or self.mask_character.isspace():
            messages.append("mask_character must be exactly one non-whitespace character.")
        elif self.mask_character == "-":
            messages.append("mask_character cannot be '-'; use an ambiguity symbol such as '?'.")

        input_path = Path(self.alignment_file_path).resolve()
        output_paths = [
            Path(path).resolve()
            for path in (self.mask_output, self.report)
            if path is not None
        ]
        if input_path in output_paths:
            messages.append("Output paths cannot overwrite the input alignment.")
        if len(output_paths) == 2 and output_paths[0] == output_paths[1]:
            messages.append("report and mask_output must use different paths.")
        if messages:
            raise PhykitUserError(messages, code=2)

    @staticmethod
    def _alignment_arrays(alignment, is_protein: bool):
        sequences = [str(record.seq).upper() for record in alignment]
        n_taxa = len(sequences)
        n_sites = alignment.get_alignment_length()
        if n_taxa == 0:
            empty = np.zeros((0, n_sites), dtype=np.uint8)
            return sequences, empty, empty.astype(np.bool_)

        try:
            alignment_array = np.frombuffer(
                "".join(sequences).encode("ascii"),
                dtype=np.uint8,
            ).reshape(n_taxa, n_sites)
            valid_mask = _ascii_valid_mask(alignment_array, is_protein)
        except UnicodeEncodeError:
            alignment_array = np.asarray(
                [list(sequence) for sequence in sequences],
                dtype="U1",
            )
            invalid = ["-", "?", "*", "X"]
            if not is_protein:
                invalid.append("N")
            valid_mask = ~np.isin(alignment_array, invalid)
        return sequences, alignment_array, valid_mask

    def detect(self, alignment, is_protein: bool) -> dict[str, object]:
        sequences, alignment_array, valid_mask = self._alignment_arrays(
            alignment,
            is_protein,
        )
        scores = _column_divergence_scores(alignment_array, valid_mask)
        valid_positions_by_taxon = [
            np.flatnonzero(valid_mask[row])
            for row in range(len(alignment))
        ]
        residue_scores_by_taxon = [
            scores[row, positions]
            for row, positions in enumerate(valid_positions_by_taxon)
        ]

        masks_by_scale = {}
        cutoffs_by_scale = {}
        for window_size, taxon_tail, sequence_tail, maximum_length in _DEFAULT_SCALES:
            scale_masks, scale_cutoffs = _detect_at_scale(
                residue_scores_by_taxon,
                window_size,
                taxon_tail,
                sequence_tail,
                self.cutoff,
                maximum_length,
            )
            masks_by_scale[window_size] = scale_masks
            cutoffs_by_scale[window_size] = scale_cutoffs

        combined_masks = []
        for row, residue_scores in enumerate(residue_scores_by_taxon):
            combined = np.zeros(residue_scores.size, dtype=np.bool_)
            for window_size, *_ in _DEFAULT_SCALES:
                combined |= masks_by_scale[window_size][row]
            combined_masks.append(combined)

        regions = self._region_rows(
            alignment,
            residue_scores_by_taxon,
            valid_positions_by_taxon,
            combined_masks,
            masks_by_scale,
        )
        masked_sequences = self._masked_sequences(
            sequences,
            valid_positions_by_taxon,
            combined_masks,
        )
        return dict(
            regions=regions,
            masked_sequences=masked_sequences,
            masked_residues=sum(int(mask.sum()) for mask in combined_masks),
            affected_taxa=sum(bool(mask.any()) for mask in combined_masks),
            cutoffs=cutoffs_by_scale,
        )

    @staticmethod
    def _region_rows(
        alignment,
        residue_scores_by_taxon,
        valid_positions_by_taxon,
        combined_masks,
        masks_by_scale,
    ) -> list[dict[str, object]]:
        regions = []
        for row, record in enumerate(alignment):
            mask = combined_masks[row]
            if not mask.any():
                continue
            padded = np.pad(mask.astype(np.int8), (1, 1))
            boundaries = np.flatnonzero(np.diff(padded)).reshape(-1, 2)
            positions = valid_positions_by_taxon[row]
            scores = residue_scores_by_taxon[row]
            for start, stop in boundaries:
                region_scores = scores[start:stop]
                scales = [
                    window_size
                    for window_size, scale_masks in masks_by_scale.items()
                    if scale_masks[row][start:stop].any()
                ]
                regions.append(
                    {
                        "taxon": record.id,
                        "alignment_start": int(positions[start]) + 1,
                        "alignment_end": int(positions[stop - 1]) + 1,
                        "sequence_start": int(start) + 1,
                        "sequence_end": int(stop),
                        "length": int(stop - start),
                        "mean_divergence_score": round(
                            float(np.mean(region_scores)),
                            4,
                        ),
                        "max_divergence_score": round(
                            float(np.max(region_scores)),
                            4,
                        ),
                        "window_sizes": scales,
                    }
                )
        return regions

    def _masked_sequences(
        self,
        sequences,
        valid_positions_by_taxon,
        combined_masks,
    ) -> list[str]:
        masked_sequences = []
        for sequence, positions, mask in zip(
            sequences,
            valid_positions_by_taxon,
            combined_masks,
        ):
            characters = list(sequence)
            for alignment_position in positions[mask]:
                characters[int(alignment_position)] = self.mask_character
            masked_sequences.append("".join(characters))
        return masked_sequences

    def _write_masked_alignment(self, alignment, sequences) -> None:
        try:
            with Path(self.mask_output).open("w") as handle:
                for record, sequence in zip(alignment, sequences):
                    header = record.description or record.id
                    handle.write(f">{header}\n")
                    for start in range(0, len(sequence), 60):
                        handle.write(sequence[start:start + 60] + "\n")
        except OSError as error:
            raise PhykitUserError(
                [f"Could not write masked alignment: {error}"],
                code=2,
            ) from error

    @staticmethod
    def _render_tsv(regions) -> str:
        header = (
            "taxon\talignment_start\talignment_end\tsequence_start\t"
            "sequence_end\tlength\tmean_divergence_score\t"
            "max_divergence_score\twindow_sizes"
        )
        lines = [header]
        lines.extend(
            "\t".join(
                [
                    str(region["taxon"]),
                    str(region["alignment_start"]),
                    str(region["alignment_end"]),
                    str(region["sequence_start"]),
                    str(region["sequence_end"]),
                    str(region["length"]),
                    str(region["mean_divergence_score"]),
                    str(region["max_divergence_score"]),
                    ",".join(str(value) for value in region["window_sizes"]),
                ]
            )
            for region in regions
        )
        return "\n".join(lines)

    def _payload(self, alignment, is_protein: bool, result) -> dict[str, object]:
        scales = [
            {
                "window_size": window_size,
                "taxon_tail": taxon_tail,
                "sequence_tail": sequence_tail,
                "maximum_region_length": maximum_length,
            }
            for window_size, taxon_tail, sequence_tail, maximum_length in _DEFAULT_SCALES
        ]
        return {
            "method": "TAPER-inspired two-dimensional alignment outlier detection",
            "citation": TAPER_CITATION,
            "doi": TAPER_DOI,
            "alignment": {
                "path": self.alignment_file_path,
                "taxa": len(alignment),
                "sites": alignment.get_alignment_length(),
                "data_type": "protein" if is_protein else "nucleotide",
            },
            "parameters": {
                "cutoff": self.cutoff,
                "mask_character": self.mask_character,
                "scales": scales,
            },
            "summary": {
                "regions": len(result["regions"]),
                "masked_residues": result["masked_residues"],
                "affected_taxa": result["affected_taxa"],
            },
            "mask_output": self.mask_output,
            "regions": result["regions"],
        }

    def _write_report(self, content: str) -> None:
        try:
            Path(self.report).write_text(content.rstrip() + "\n")
        except OSError as error:
            raise PhykitUserError(
                [f"Could not write report: {error}"],
                code=2,
            ) from error

    def run(self) -> None:
        self._validate_options()
        alignment, _, is_protein = self.get_alignment_and_format()
        result = self.detect(alignment, is_protein)

        if self.mask_output is not None:
            self._write_masked_alignment(
                alignment,
                result["masked_sequences"],
            )

        if self.json_output:
            payload = self._payload(alignment, is_protein, result)
            if self.report is not None:
                self._write_report(json.dumps(payload, indent=2, sort_keys=True))
            else:
                print_json(payload)
        else:
            report = self._render_tsv(result["regions"])
            if self.report is not None:
                self._write_report(report)
            else:
                print(report)

        if self.report is not None:
            print(f"report\t{self.report}")
        if self.mask_output is not None and self.report is not None:
            print(f"masked_alignment\t{self.mask_output}")
