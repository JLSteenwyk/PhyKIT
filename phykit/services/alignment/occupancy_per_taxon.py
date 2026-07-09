from __future__ import annotations

from .base import Alignment


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
_DNA_VALID_LOOKUP = None
_PROTEIN_VALID_LOOKUP = None
_DNA_INVALID_BYTES = b"-?*XNxn"
_PROTEIN_INVALID_BYTES = b"-?*Xx"
_DNA_INVALID_COUNT_CHARS = "-?*XNxn"
_PROTEIN_INVALID_COUNT_CHARS = "-?*Xx"
_INVALID_SCAN_BYTES = 4096


class _IdenticalOccupancyRows(list):
    __slots__ = ("shared_occupancy",)

    def __init__(self, rows, shared_occupancy):
        super().__init__(rows)
        self.shared_occupancy = shared_occupancy


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def _get_valid_lookup(is_protein: bool):
    global _DNA_VALID_LOOKUP, _PROTEIN_VALID_LOOKUP

    if is_protein:
        if _PROTEIN_VALID_LOOKUP is None:
            lookup = np.ones(256, dtype=np.bool_)
            lookup[np.frombuffer("-?*Xx".encode("ascii"), dtype=np.uint8)] = False
            _PROTEIN_VALID_LOOKUP = lookup
        return _PROTEIN_VALID_LOOKUP

    if _DNA_VALID_LOOKUP is None:
        lookup = np.ones(256, dtype=np.bool_)
        lookup[np.frombuffer("-?*XNxn".encode("ascii"), dtype=np.uint8)] = False
        _DNA_VALID_LOOKUP = lookup
    return _DNA_VALID_LOOKUP


def _has_invalid_bytes(sequence_bytes: bytes, invalid_bytes: bytes) -> bool:
    if len(sequence_bytes) <= _INVALID_SCAN_BYTES:
        return any(code in sequence_bytes for code in invalid_bytes)

    has_invalid = any(
        code in sequence_bytes[:_INVALID_SCAN_BYTES]
        for code in invalid_bytes
    )
    if not has_invalid:
        has_invalid = any(
            code in sequence_bytes[-_INVALID_SCAN_BYTES:]
            for code in invalid_bytes
        )
    if not has_invalid:
        has_invalid = any(code in sequence_bytes for code in invalid_bytes)
    return has_invalid


def _identical_occupancy_rows(record_data, occupancy):
    return _IdenticalOccupancyRows(
        ((record_id, occupancy) for record_id, _ in record_data),
        occupancy,
    )


def _occupancy_from_ascii_matrix(record_data, is_protein: bool):
    if not record_data:
        return []

    record_iter = iter(record_data)
    _, first_sequence = next(record_iter)
    seq_len = len(first_sequence)
    all_identical = True
    for _, sequence in record_iter:
        if len(sequence) != seq_len:
            return None
        if sequence != first_sequence:
            all_identical = False

    if all_identical:
        try:
            sequence_bytes = first_sequence.encode("ascii")
        except UnicodeEncodeError:
            return None
        invalid_bytes = _PROTEIN_INVALID_BYTES if is_protein else _DNA_INVALID_BYTES
        if not _has_invalid_bytes(sequence_bytes, invalid_bytes):
            return _identical_occupancy_rows(record_data, 1.0)
        occupancy = (
            0.0
            if seq_len == 0
            else len(sequence_bytes.translate(None, invalid_bytes)) / seq_len
        )
        return _identical_occupancy_rows(record_data, occupancy)

    sequences = [sequence for _, sequence in record_data]
    try:
        alignment_bytes = "".join(sequences).encode("ascii")
    except UnicodeEncodeError:
        return None

    invalid_bytes = _PROTEIN_INVALID_BYTES if is_protein else _DNA_INVALID_BYTES
    if not any(code in alignment_bytes for code in invalid_bytes):
        occupancy = 0.0 if seq_len == 0 else 1.0
        return _identical_occupancy_rows(record_data, occupancy)

    try:
        alignment_array = np.frombuffer(
            alignment_bytes,
            dtype=np.uint8,
        ).reshape(len(sequences), seq_len)
    except ValueError:
        return None

    valid_lookup = _get_valid_lookup(is_protein)
    valid_counts = np.count_nonzero(valid_lookup[alignment_array], axis=1)
    if seq_len == 0:
        occupancies = np.zeros(len(sequences), dtype=np.float64)
    else:
        occupancies = valid_counts / seq_len
    return [
        (record_id, float(occupancy))
        for (record_id, _), occupancy in zip(record_data, occupancies)
    ]


def _occupancy_for_sequence(sequence: str, is_protein: bool) -> float:
    if is_protein:
        invalid_bytes = _PROTEIN_INVALID_BYTES
    else:
        invalid_bytes = _DNA_INVALID_BYTES

    try:
        seq_bytes = sequence.encode("ascii")
        if _has_invalid_bytes(seq_bytes, invalid_bytes):
            valid_count = len(seq_bytes.translate(None, invalid_bytes))
        else:
            valid_count = len(seq_bytes)
    except UnicodeEncodeError:
        invalid_count_chars = (
            _PROTEIN_INVALID_COUNT_CHARS
            if is_protein
            else _DNA_INVALID_COUNT_CHARS
        )
        valid_count = len(sequence) - sum(
            sequence.count(char) for char in invalid_count_chars
        )
    return (valid_count / len(sequence)) if len(sequence) > 0 else 0.0


def _occupancy_json_rows(occupancies):
    shared_occupancy = getattr(occupancies, "shared_occupancy", None)
    if shared_occupancy is not None:
        occupancy = 1.0 if shared_occupancy == 1.0 else round(shared_occupancy, 4)
        return [
            {"taxon": taxon, "occupancy": occupancy}
            for taxon, _ in occupancies
        ]

    if occupancies and occupancies[0][1] == 1.0:
        return [
            {
                "taxon": taxon,
                "occupancy": 1.0 if occupancy == 1.0 else round(occupancy, 4),
            }
            for taxon, occupancy in occupancies
        ]
    return [
        {"taxon": taxon, "occupancy": round(occupancy, 4)}
        for taxon, occupancy in occupancies
    ]


def _alignment_size(alignment):
    try:
        return len(alignment)
    except TypeError:
        return None


class OccupancyPerTaxon(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()
        occupancies = self.calculate_occupancy_per_taxon(alignment, is_protein)

        if self.json_output:
            rows = _occupancy_json_rows(occupancies)
            print_json(
                dict(
                    rows=rows,
                    taxa=rows,
                )
            )
            return

        shared_occupancy = getattr(occupancies, "shared_occupancy", None)
        if shared_occupancy is None:
            lines = [
                f"{taxon}\t{round(occupancy, 4)}"
                for taxon, occupancy in occupancies
            ]
        else:
            occupancy = round(shared_occupancy, 4)
            lines = [f"{taxon}\t{occupancy}" for taxon, _ in occupancies]
        if lines:
            print("\n".join(lines))

    def process_args(self, args) -> dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
        )

    def calculate_occupancy_per_taxon(self, alignment, is_protein: bool):
        alignment_size = _alignment_size(alignment)
        if alignment_size == 1:
            record = alignment[0]
            return [
                (
                    record.id,
                    _occupancy_for_sequence(str(record.seq), is_protein),
                )
            ]

        record_data = [(record.id, str(record.seq)) for record in alignment]
        matrix_result = _occupancy_from_ascii_matrix(record_data, is_protein)
        if matrix_result is not None:
            return matrix_result

        output = []
        for record_id, seq in record_data:
            output.append((record_id, _occupancy_for_sequence(seq, is_protein)))
        return output
