from __future__ import annotations

import sys


from .base import Alignment


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def get_alignment_and_format(*args, **kwargs):
    from ...helpers.files import get_alignment_and_format as _get_alignment_and_format

    return _get_alignment_and_format(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()
_DNA_VALID_LOOKUP = None
_PROTEIN_VALID_LOOKUP = None
_GC_LOOKUP = None
_DNA_INVALID_BYTES = b"-?*XNxn"
_PROTEIN_INVALID_BYTES = b"-?*Xx"
_DNA_INVALID_CHARS = "-?*XN"
_PROTEIN_INVALID_CHARS = "-?*X"


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


def _get_gc_lookup():
    global _GC_LOOKUP
    if _GC_LOOKUP is None:
        lookup = np.zeros(256, dtype=np.bool_)
        lookup[np.frombuffer("GCgc".encode("ascii"), dtype=np.uint8)] = True
        _GC_LOOKUP = lookup
    return _GC_LOOKUP


def _gc_counts_from_upper_sequence(sequence: str, is_protein: bool) -> tuple[int, int]:
    invalid_chars = _PROTEIN_INVALID_CHARS if is_protein else _DNA_INVALID_CHARS
    invalid_bytes = _PROTEIN_INVALID_BYTES if is_protein else _DNA_INVALID_BYTES
    try:
        sequence_bytes = sequence.encode("ascii")
    except UnicodeEncodeError:
        valid_count = len(sequence) - sum(sequence.count(char) for char in invalid_chars)
        gc_count = sequence.count("G") + sequence.count("C")
        return valid_count, gc_count

    if any(code in sequence_bytes for code in invalid_bytes):
        valid_count = len(sequence_bytes.translate(None, invalid_bytes))
    else:
        valid_count = len(sequence_bytes)
    gc_count = sequence_bytes.count(b"G") + sequence_bytes.count(b"C")
    return valid_count, gc_count


def _common_upper_sequence(sequences: list[str]) -> str | None:
    first_raw = sequences[0]
    first_upper = first_raw.upper()
    for sequence in sequences:
        if sequence != first_raw and sequence.upper() != first_upper:
            return None
    return first_upper


def _gc_counts_from_ascii_matrix(records, is_protein: bool):
    record_data = [(record.id, str(record.seq)) for record in records]
    if not record_data:
        return [], None, None

    sequences = [sequence for _, sequence in record_data]
    first_sequence = _common_upper_sequence(sequences)
    if first_sequence is not None:
        valid_count, gc_count = _gc_counts_from_upper_sequence(
            first_sequence,
            is_protein,
        )
        return (
            record_data,
            np.full(len(sequences), valid_count, dtype=np.intp),
            np.full(len(sequences), gc_count, dtype=np.intp),
        )

    seq_len = len(sequences[0])
    if any(len(sequence) != seq_len for sequence in sequences):
        return record_data, None, None

    try:
        alignment_bytes = "".join(sequences).encode("ascii")
        alignment_array = np.frombuffer(
            alignment_bytes,
            dtype=np.uint8,
        ).reshape(len(sequences), seq_len)
    except UnicodeEncodeError:
        return record_data, None, None

    gc_lookup = _get_gc_lookup()
    gc_counts = np.count_nonzero(gc_lookup[alignment_array], axis=1)
    invalid_bytes = _PROTEIN_INVALID_BYTES if is_protein else _DNA_INVALID_BYTES
    if any(code in alignment_bytes for code in invalid_bytes):
        valid_lookup = _get_valid_lookup(is_protein)
        valid_counts = np.count_nonzero(valid_lookup[alignment_array], axis=1)
    else:
        valid_counts = np.full(len(sequences), seq_len, dtype=np.intp)
    return record_data, valid_counts, gc_counts


def _gc_total_from_ascii(records, is_protein: bool):
    sequences = [str(record.seq) for record in records]
    if not sequences:
        return None

    first_sequence = _common_upper_sequence(sequences)
    if first_sequence is not None:
        valid_count, gc_count = _gc_counts_from_upper_sequence(
            first_sequence,
            is_protein,
        )
        n_sequences = len(sequences)
        return valid_count * n_sequences, gc_count * n_sequences

    try:
        alignment_bytes = "".join(sequences).encode("ascii")
    except UnicodeEncodeError:
        return None

    seq_array = np.frombuffer(alignment_bytes, dtype=np.uint8)
    invalid_bytes = _PROTEIN_INVALID_BYTES if is_protein else _DNA_INVALID_BYTES
    has_invalid = any(code in alignment_bytes[:4096] for code in invalid_bytes)
    if not has_invalid:
        has_invalid = any(code in alignment_bytes[-4096:] for code in invalid_bytes)
    if not has_invalid:
        has_invalid = any(code in alignment_bytes for code in invalid_bytes)
    if has_invalid:
        valid_lookup = _get_valid_lookup(is_protein)
        valid_count = int(np.count_nonzero(valid_lookup[seq_array]))
    else:
        valid_count = len(seq_array)
    gc_lookup = _get_gc_lookup()
    gc_count = int(np.count_nonzero(gc_lookup[seq_array]))
    return valid_count, gc_count


class GCContent(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(fasta=parsed["fasta"], verbose=parsed["verbose"])
        self.json_output = parsed["json_output"]

    def run(self):
        records, _, is_protein = get_alignment_and_format(self.fasta)

        if is_protein:
            print("GC content can't be calculated for protein sequences")
            sys.exit(2)

        if self.json_output:
            if self.verbose:
                rows = self.calculate_gc_per_sequence_data(records, is_protein)
                row_payload = [
                    {"taxon": taxon, "gc_content": round(gc_content, 4)}
                    for taxon, gc_content in rows
                ]
                print_json(
                    {
                        "verbose": True,
                        "rows": row_payload,
                        "sequences": row_payload,
                    }
                )
            else:
                print_json(
                    dict(
                        verbose=False,
                        gc_content=self.calculate_gc_total_value(records, is_protein),
                    )
                )
            return

        if self.verbose:
            self.calculate_gc_per_sequence(records, is_protein)
        else:
            self.calculate_gc_total(records, is_protein)

    def process_args(self, args) -> dict[str, str]:
        return dict(
            fasta=args.fasta,
            verbose=args.verbose,
            json_output=getattr(args, "json", False),
        )

    def calculate_gc_per_sequence_data(
        self, records: MultipleSeqAlignment, is_protein: bool = False
    ) -> list[tuple[str, float]]:
        record_data, valid_counts, gc_counts = _gc_counts_from_ascii_matrix(
            records,
            is_protein,
        )
        if valid_counts is not None and gc_counts is not None:
            gc_contents = np.divide(
                gc_counts,
                valid_counts,
                out=np.zeros_like(gc_counts, dtype=np.float64),
                where=valid_counts > 0,
            )
            return [
                (record_id, float(gc_content))
                for (record_id, _), gc_content in zip(record_data, gc_contents)
            ]

        gap_chars = {char.upper() for char in self.get_gap_chars(is_protein)}
        valid_lookup = _get_valid_lookup(is_protein)
        gc_lookup = _get_gc_lookup()
        output = []
        for record_id, seq in record_data:
            try:
                seq_array = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
                valid_len = int(np.count_nonzero(valid_lookup[seq_array]))
                gc_count = int(np.count_nonzero(gc_lookup[seq_array]))
            except UnicodeEncodeError:
                seq = seq.upper()
                valid_len = len(seq) - sum(seq.count(char) for char in gap_chars)
                gc_count = seq.count("G") + seq.count("C")
            if valid_len > 0:
                gc_content = gc_count / valid_len
            else:
                gc_content = 0.0
            output.append((record_id, float(gc_content)))
        return output

    def calculate_gc_per_sequence(
        self, records: MultipleSeqAlignment, is_protein: bool = False
    ) -> None:
        rows = self.calculate_gc_per_sequence_data(records, is_protein)
        lines = [
            f"{record_id}\t{round(gc_content, 4)}"
            for record_id, gc_content in rows
        ]
        if not lines:
            return

        try:
            print("\n".join(lines))
        except BrokenPipeError:
            pass

    def calculate_gc_total_value(
        self, records: MultipleSeqAlignment, is_protein: bool = False
    ) -> float:
        records = list(records)
        total_counts = _gc_total_from_ascii(records, is_protein)
        if total_counts is not None:
            valid_count, gc_count = total_counts
            if valid_count > 0:
                return round(gc_count / valid_count, 4)

            print(
                "Input file has an unacceptable format. Please check input file argument."
            )
            sys.exit(2)

        gap_chars = {char.upper() for char in self.get_gap_chars(is_protein)}
        valid_lookup = _get_valid_lookup(is_protein)
        gc_lookup = _get_gc_lookup()
        valid_count = 0
        gc_count = 0
        ascii_chunks = []
        for record in records:
            seq = str(record.seq)
            try:
                ascii_chunks.append(seq.encode("ascii"))
            except UnicodeEncodeError:
                seq = seq.upper()
                valid_count += len(seq) - sum(seq.count(char) for char in gap_chars)
                gc_count += seq.count("G") + seq.count("C")
        if ascii_chunks:
            seq_array = np.frombuffer(b"".join(ascii_chunks), dtype=np.uint8)
            valid_count += int(np.count_nonzero(valid_lookup[seq_array]))
            gc_count += int(np.count_nonzero(gc_lookup[seq_array]))
        if valid_count > 0:
            return round(gc_count / valid_count, 4)

        print(
            "Input file has an unacceptable format. Please check input file argument."
        )
        sys.exit(2)

    def calculate_gc_total(
        self, records: MultipleSeqAlignment, is_protein: bool = False
    ) -> None:
        print(self.calculate_gc_total_value(records, is_protein))
