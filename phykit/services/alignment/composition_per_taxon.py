from __future__ import annotations

from .base import Alignment


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()
_COMPOSITION_ROW_ZIP_MIN_COUNT = 50_000


def _composition_output_rows(record_ids, freqs):
    if len(record_ids) >= _COMPOSITION_ROW_ZIP_MIN_COUNT:
        return list(zip(record_ids, freqs))
    return [
        (record_id, freqs[row_idx])
        for row_idx, record_id in enumerate(record_ids)
    ]


class CompositionPerTaxon(Alignment):
    _INVALID_LOOKUP_CACHE = {}

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()
        symbols, rows = self.calculate_composition_per_taxon(alignment, is_protein)
        if not symbols:
            return

        if self.json_output:
            payload_rows = [
                dict(
                    taxon=taxon,
                    composition={
                        symbol: round(float(value), 4)
                        for symbol, value in zip(symbols, comps)
                    },
                )
                for taxon, comps in rows
            ]
            print_json(dict(symbols=symbols, rows=payload_rows, taxa=payload_rows))
            return

        lines = []
        for taxon, comps in rows:
            comp_str = ";".join(
                f"{symbol}:{round(float(value), 4)}"
                for symbol, value in zip(symbols, comps)
            )
            lines.append(f"{taxon}\t{comp_str}")
        if lines:
            print("\n".join(lines))

    def process_args(self, args) -> dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
        )

    @classmethod
    def _invalid_lookup_for_chars(cls, invalid_chars):
        invalid_chars = frozenset(char.upper() for char in invalid_chars)
        lookup = cls._INVALID_LOOKUP_CACHE.get(invalid_chars)
        if lookup is None:
            lookup = np.zeros(256, dtype=bool)
            lookup[
                np.fromiter((ord(char) for char in invalid_chars), dtype=np.uint8)
            ] = True
            cls._INVALID_LOOKUP_CACHE[invalid_chars] = lookup
        return lookup

    def calculate_composition_per_taxon(
        self, alignment, is_protein: bool
    ) -> tuple[list[str], list[tuple[str, np.ndarray]]]:
        invalid_chars = {char.upper() for char in self.get_gap_chars(is_protein)}
        raw_records = [(record.id, str(record.seq)) for record in alignment]
        if not raw_records:
            return [], []

        record_ids = [record_id for record_id, _ in raw_records]
        first_raw_sequence = raw_records[0][1]
        first_sequence = first_raw_sequence.upper()
        all_identical = True
        for idx in range(1, len(raw_records)):
            sequence = raw_records[idx][1]
            if sequence != first_raw_sequence and sequence.upper() != first_sequence:
                all_identical = False
                break

        if all_identical:
            try:
                sequence_array = np.frombuffer(
                    first_sequence.encode("ascii"),
                    dtype=np.uint8,
                )
                invalid_lookup = self._invalid_lookup_for_chars(invalid_chars)
                valid_symbols_raw, counts = np.unique(
                    sequence_array[~invalid_lookup[sequence_array]],
                    return_counts=True,
                )
                symbols = [chr(int(symbol)) for symbol in valid_symbols_raw]
            except UnicodeEncodeError:
                sequence_array = np.array(list(first_sequence), dtype="U1")
                invalid_chars_array = np.array(list(invalid_chars), dtype="U1")
                valid_symbols_raw, counts = np.unique(
                    sequence_array[~np.isin(sequence_array, invalid_chars_array)],
                    return_counts=True,
                )
                symbols = valid_symbols_raw.tolist()

            if not symbols:
                return [], []

            freqs = counts.astype(np.float64) / float(counts.sum())
            return symbols, [
                (record_id, freqs.copy())
                for record_id in record_ids
            ]

        records = [
            (record_id, sequence.upper())
            for record_id, sequence in raw_records
        ]
        sequences = [seq for _, seq in records]
        aln_len = len(sequences[0])
        try:
            alignment_bytes = "".join(sequences).encode("ascii")
            alignment_array = np.frombuffer(
                alignment_bytes,
                dtype=np.uint8,
            ).reshape(len(sequences), aln_len)
            invalid_lookup = self._invalid_lookup_for_chars(invalid_chars)
            invalid_codes = tuple(ord(char) for char in invalid_chars)
            if any(code in alignment_bytes for code in invalid_codes):
                valid_mask = ~invalid_lookup[alignment_array]
                valid_lengths = np.count_nonzero(valid_mask, axis=1)
                valid_symbols_raw = np.unique(
                    alignment_array[valid_mask]
                )
            else:
                valid_mask = None
                valid_lengths = np.full(len(sequences), aln_len, dtype=np.intp)
                valid_symbols_raw = np.unique(alignment_array)
            symbols = [chr(int(symbol)) for symbol in valid_symbols_raw]
            symbol_values = valid_symbols_raw
        except UnicodeEncodeError:
            alignment_array = np.array([list(seq) for seq in sequences], dtype="U1")
            invalid_chars_array = np.array(list(invalid_chars), dtype="U1")
            valid_mask = ~np.isin(alignment_array, invalid_chars_array)
            valid_lengths = np.count_nonzero(valid_mask, axis=1)
            symbol_values = np.unique(
                alignment_array[valid_mask]
            )
            symbols = symbol_values.tolist()

        if not symbols:
            return [], []

        if len(symbol_values) == 1:
            freqs = np.zeros((len(sequences), 1), dtype=np.float64)
            freqs[valid_lengths > 0, 0] = 1.0
            return symbols, _composition_output_rows(record_ids, freqs)

        if alignment_array.dtype == np.uint8:
            if len(symbol_values) <= 8 or (
                len(sequences) >= 1024 and aln_len <= 512
            ):
                counts = np.array(
                    [
                        np.count_nonzero(alignment_array == symbol, axis=1)
                        for symbol in symbol_values
                    ],
                    dtype=np.float64,
                ).T
            else:
                counts = np.zeros((len(sequences), len(symbol_values)), dtype=np.float64)
                for row_idx, row in enumerate(alignment_array):
                    row_values = row if valid_mask is None else row[valid_mask[row_idx]]
                    counts[row_idx] = np.bincount(
                        row_values,
                        minlength=256,
                    )[symbol_values]
        else:
            counts = np.array(
                [
                    np.count_nonzero(alignment_array == symbol, axis=1)
                    for symbol in symbol_values
                ],
                dtype=np.float64,
            ).T
        freqs = np.divide(
            counts,
            valid_lengths[:, None],
            out=np.zeros_like(counts, dtype=np.float64),
            where=valid_lengths[:, None] > 0,
        )

        output = _composition_output_rows(record_ids, freqs)

        return symbols, output
