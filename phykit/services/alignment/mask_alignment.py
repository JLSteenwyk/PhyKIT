from __future__ import annotations

import sys

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


def _column_entropies_from_ascii_codes(
    alignment_array,
    valid_mask,
    valid_symbols,
    block_size: int = 8192,
):
    aln_len = alignment_array.shape[1]
    entropies = np.zeros(aln_len, dtype=np.float64)
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
            if not np.any(block_valid):
                continue
            count_index = (block_by_site + site_offsets)[block_valid.T]
        counts = np.bincount(
            count_index,
            minlength=width * 256,
        ).reshape(width, 256)[:, symbol_codes].T.astype(np.float64, copy=False)
        totals = counts.sum(axis=0)
        probs = np.divide(
            counts,
            totals,
            out=np.zeros_like(counts, dtype=np.float64),
            where=totals > 0,
        )
        log_probs = np.zeros_like(probs, dtype=np.float64)
        positive_probs = probs > 0
        np.log2(probs, out=log_probs, where=positive_probs)
        probs *= log_probs
        entropies[start:stop] = -np.sum(probs, axis=0)

    return entropies


class MaskAlignment(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.max_gap = parsed["max_gap"]
        self.min_occupancy = parsed["min_occupancy"]
        self.max_entropy = parsed["max_entropy"]
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        self._validate_thresholds()
        alignment, _, is_protein = self.get_alignment_and_format()
        keep_mask = self.calculate_keep_mask(alignment, is_protein)
        masked = self.apply_mask(alignment, keep_mask)

        if self.json_output:
            rows = [dict(taxon=taxon, sequence=seq) for taxon, seq in masked.items()]
            print_json(
                dict(
                    thresholds=dict(
                        max_gap=self.max_gap,
                        min_occupancy=self.min_occupancy,
                        max_entropy=self.max_entropy,
                    ),
                    kept_sites=int(np.count_nonzero(keep_mask)),
                    total_sites=int(len(keep_mask)),
                    rows=rows,
                    taxa=rows,
                )
            )
            return

        if masked:
            print("\n".join(
                f">{taxon}\n{seq}"
                for taxon, seq in masked.items()
            ))

    def process_args(self, args) -> dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            max_gap=args.max_gap,
            min_occupancy=args.min_occupancy,
            max_entropy=args.max_entropy,
            json_output=getattr(args, "json", False),
        )

    def _validate_thresholds(self) -> None:
        if self.max_gap < 0.0 or self.max_gap > 1.0:
            print("max_gap must be between 0 and 1.")
            sys.exit(2)
        if self.min_occupancy < 0.0 or self.min_occupancy > 1.0:
            print("min_occupancy must be between 0 and 1.")
            sys.exit(2)
        if self.max_entropy is not None and self.max_entropy < 0.0:
            print("max_entropy must be >= 0.")
            sys.exit(2)

    @staticmethod
    def _valid_sites_for_identical_sequence(
        sequence: str,
        aln_len: int,
        is_protein: bool,
        invalid_chars,
    ) -> np.ndarray:
        try:
            sequence_array = np.frombuffer(
                sequence[:aln_len].encode("ascii"),
                dtype=np.uint8,
            )
            invalid_lookup = _get_invalid_lookup(is_protein)
            return ~invalid_lookup[sequence_array]
        except UnicodeEncodeError:
            return np.fromiter(
                (char not in invalid_chars for char in sequence),
                dtype=np.bool_,
                count=aln_len,
            )

    def calculate_keep_mask(self, alignment, is_protein: bool) -> np.ndarray:
        sequences = [str(record.seq).upper() for record in alignment]
        aln_len = alignment.get_alignment_length()
        if (
            self.max_entropy is None
            and self.max_gap >= 1.0
            and self.min_occupancy <= 0.0
        ):
            return np.ones(aln_len, dtype=np.bool_)

        if sequences:
            if _all_sequences_identical(sequences):
                first_sequence = sequences[0]
                if is_protein:
                    invalid_chars = {"-", "?", "*", "X"}
                else:
                    invalid_chars = {"-", "?", "*", "X", "N"}
                if self.max_gap >= 1.0 and self.min_occupancy <= 0.0:
                    keep_mask = np.ones(aln_len, dtype=np.bool_)
                else:
                    keep_mask = self._valid_sites_for_identical_sequence(
                        first_sequence,
                        aln_len,
                        is_protein,
                        invalid_chars,
                    )
                if self.max_entropy is not None:
                    keep_mask &= 0.0 <= self.max_entropy
                return keep_mask

        if is_protein:
            invalid_chars = ["-", "?", "*", "X"]
        else:
            invalid_chars = ["-", "?", "*", "X", "N"]

        ascii_matrix = False
        try:
            alignment_bytes = "".join(sequences).encode("ascii")
            invalid_codes = b"-?*X" if is_protein else b"-?*XN"
            clean_ascii = not any(
                alignment_bytes.find(invalid_code) != -1
                for invalid_code in invalid_codes
            )
            if (
                self.max_entropy is None
                and self.max_gap >= 0.0
                and self.min_occupancy <= 1.0
                and clean_ascii
            ):
                return np.ones(aln_len, dtype=np.bool_)

            alignment_array = np.frombuffer(
                alignment_bytes,
                dtype=np.uint8,
            ).reshape(len(sequences), aln_len)
            if clean_ascii:
                valid_mask = None
                keep_mask = np.full(
                    aln_len,
                    self.max_gap >= 0.0 and self.min_occupancy <= 1.0,
                    dtype=np.bool_,
                )
            else:
                invalid_lookup = _get_invalid_lookup(is_protein)
                valid_mask = ~invalid_lookup[alignment_array]
                occupancy = np.mean(valid_mask, axis=0).astype(np.float64)
                gap_fraction = 1.0 - occupancy
                keep_mask = (
                    (gap_fraction <= self.max_gap)
                    & (occupancy >= self.min_occupancy)
                )
            ascii_matrix = True
        except UnicodeEncodeError:
            alignment_array = np.array([list(seq) for seq in sequences], dtype="U1")
            invalid_chars_array = np.array(invalid_chars, dtype="U1")
            valid_mask = ~np.isin(alignment_array, invalid_chars_array)
            occupancy = np.mean(valid_mask, axis=0).astype(np.float64)
            gap_fraction = 1.0 - occupancy
            keep_mask = (
                (gap_fraction <= self.max_gap)
                & (occupancy >= self.min_occupancy)
            )

        if self.max_entropy is not None:
            entropies = np.zeros(alignment_array.shape[1], dtype=np.float64)
            if valid_mask is None:
                valid_symbols = np.unique(alignment_array)
            else:
                valid_symbols = np.unique(alignment_array[valid_mask])
            if valid_symbols.size <= 1:
                keep_mask &= 0.0 <= self.max_entropy
            else:
                if ascii_matrix and valid_symbols.size > 8:
                    entropies = _column_entropies_from_ascii_codes(
                        alignment_array,
                        valid_mask,
                        valid_symbols,
                    )
                else:
                    count_reducer = np.count_nonzero if ascii_matrix else np.sum
                    counts = np.array(
                        [
                            count_reducer(alignment_array == symbol, axis=0)
                            for symbol in valid_symbols
                        ],
                        dtype=np.float64,
                    )
                    totals = np.sum(counts, axis=0)
                    probs = np.divide(
                        counts,
                        totals,
                        out=np.zeros_like(counts, dtype=np.float64),
                        where=totals > 0,
                    )
                    log_probs = np.zeros_like(probs, dtype=np.float64)
                    positive_probs = probs > 0
                    np.log2(probs, out=log_probs, where=positive_probs)
                    probs *= log_probs
                    entropies = -np.sum(probs, axis=0)
            keep_mask &= entropies <= self.max_entropy

        return keep_mask

    def apply_mask(self, alignment, keep_mask: np.ndarray) -> dict[str, str]:
        records = list(alignment)
        if not records:
            return {}

        mask_len = len(keep_mask)
        if (
            mask_len
            and not keep_mask[0]
            and not np.any(keep_mask)
            and all(len(str(record.seq)) == mask_len for record in records)
        ):
            return {
                record.id: ""
                for record in records
            }

        if len(keep_mask) and np.all(keep_mask):
            return {
                record.id: str(record.seq).upper()
                for record in records
            }

        sequences = [str(record.seq).upper() for record in records]
        try:
            alignment_array = np.frombuffer(
                "".join(sequences).encode("ascii"),
                dtype="S1",
            ).reshape(len(sequences), len(keep_mask))
            masked_array = alignment_array[:, keep_mask]
            return {
                record.id: row.tobytes().decode("ascii")
                for record, row in zip(records, masked_array)
            }
        except UnicodeEncodeError:
            alignment_array = np.array([list(seq) for seq in sequences], dtype="U1")
            masked_array = alignment_array[:, keep_mask]
            return {
                record.id: "".join(row.tolist())
                for record, row in zip(records, masked_array)
            }
