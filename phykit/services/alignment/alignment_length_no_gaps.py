from __future__ import annotations

from .base import Alignment


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)

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
_DNA_GAP_BYTES = b"-?*XNxn"
_PROTEIN_GAP_BYTES = b"-?*Xx"
_DNA_GAP_CHARS = {"-", "?", "*", "X", "N"}
_PROTEIN_GAP_CHARS = {"-", "?", "*", "X"}
_DNA_GAP_COUNT_CHARS = "-?*XNxn"
_PROTEIN_GAP_COUNT_CHARS = "-?*Xx"
_DNA_GAP_COUNT_CHAR_SET = frozenset(_DNA_GAP_COUNT_CHARS)
_PROTEIN_GAP_COUNT_CHAR_SET = frozenset(_PROTEIN_GAP_COUNT_CHARS)
_NO_GAP_SCALAR_MAX_CELLS = 8192
_DENSE_GAP_SCAN_MIN_BYTES = 1_000_000
_DENSE_GAP_SAMPLE_BYTES = 8192
_DENSE_GAP_SAMPLE_MIN_FRACTION = 0.03


def _count_no_gap_sites_in_identical_sequence(sequence: str, is_protein: bool) -> int:
    try:
        sequence_bytes = sequence.encode("ascii")
        gap_bytes = _PROTEIN_GAP_BYTES if is_protein else _DNA_GAP_BYTES
        if not any(code in sequence_bytes for code in gap_bytes):
            return len(sequence)
        return len(sequence_bytes.translate(None, gap_bytes))
    except UnicodeEncodeError:
        gap_chars = _PROTEIN_GAP_COUNT_CHARS if is_protein else _DNA_GAP_COUNT_CHARS
        return len(sequence) - sum(sequence.count(char) for char in gap_chars)


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


def _count_columns_without_gap_bytes(
    alignment_bytes: bytes,
    aln_len: int,
    gap_bytes: bytes,
) -> int:
    columns_with_gaps = bytearray(aln_len)
    gap_column_count = 0
    for gap_byte in gap_bytes:
        gap_position = alignment_bytes.find(gap_byte)
        while gap_position != -1:
            column = gap_position % aln_len
            if not columns_with_gaps[column]:
                columns_with_gaps[column] = 1
                gap_column_count += 1
                if gap_column_count == aln_len:
                    return 0
            gap_position = alignment_bytes.find(gap_byte, gap_position + 1)
    return aln_len - gap_column_count


def _should_use_dense_gap_matrix_scan(
    alignment_bytes: bytes,
    gap_bytes: bytes,
) -> bool:
    if len(alignment_bytes) < _DENSE_GAP_SCAN_MIN_BYTES:
        return False

    sample_len = min(_DENSE_GAP_SAMPLE_BYTES, len(alignment_bytes))
    sample = alignment_bytes[:sample_len]
    gap_count = sum(sample.count(bytes((gap_byte,))) for gap_byte in gap_bytes)
    return (gap_count / sample_len) >= _DENSE_GAP_SAMPLE_MIN_FRACTION


def _count_columns_without_gap_bytes_matrix(
    alignment_bytes: bytes,
    aln_len: int,
    gap_bytes: bytes,
) -> int:
    alignment_array = np.frombuffer(
        alignment_bytes,
        dtype=np.uint8,
    ).reshape(len(alignment_bytes) // aln_len, aln_len)
    gap_lookup = np.zeros(256, dtype=np.bool_)
    gap_lookup[np.frombuffer(gap_bytes, dtype=np.uint8)] = True
    columns_with_gaps = np.any(gap_lookup[alignment_array], axis=0)
    return int(np.count_nonzero(~columns_with_gaps))


def _count_no_gap_sites_scalar(
    sequences: list[str],
    aln_len: int,
    is_protein: bool,
) -> int | None:
    if not sequences or (len(sequences) * aln_len) > _NO_GAP_SCALAR_MAX_CELLS:
        return None
    if any(len(sequence) != aln_len for sequence in sequences):
        return None

    gap_chars = _PROTEIN_GAP_COUNT_CHAR_SET if is_protein else _DNA_GAP_COUNT_CHAR_SET
    no_gap_sites = 0
    for column_idx in range(aln_len):
        for sequence in sequences:
            if sequence[column_idx] in gap_chars:
                break
        else:
            no_gap_sites += 1
    return no_gap_sites


class AlignmentLengthNoGaps(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()
        (
            aln_len_no_gaps,
            aln_len,
            aln_len_no_gaps_per,
        ) = self.calculate_alignment_length_no_gaps(alignment, is_protein)
        if self.json_output:
            print_json(
                dict(
                    alignment_length_no_gaps=aln_len_no_gaps,
                    alignment_length=aln_len,
                    percent_no_gaps=round(aln_len_no_gaps_per, 4),
                )
            )
            return
        print(f"{aln_len_no_gaps}\t{aln_len}\t{round(aln_len_no_gaps_per, 4)}")

    def process_args(
        self,
        args,
    ) -> dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
        )

    def calculate_alignment_length_no_gaps(
        self,
        alignment: MultipleSeqAlignment,
        is_protein: bool,
    ) -> tuple[int, int, float]:
        aln_len = alignment.get_alignment_length()
        aln_len_no_gaps = self.get_sites_no_gaps_count(
            alignment,
            aln_len,
            is_protein
        )

        aln_len_no_gaps_per = (aln_len_no_gaps / aln_len) * 100

        return aln_len_no_gaps, aln_len, aln_len_no_gaps_per

    def get_sites_no_gaps_count(
        self,
        alignment: MultipleSeqAlignment,
        aln_len: int,
        is_protein: bool,
    ) -> int:
        """
        Count sites in the alignment with no gaps
        """
        sequences = [str(record.seq) for record in alignment]

        if not sequences:
            return aln_len

        first_sequence = sequences[0]
        if _all_sequences_identical(sequences):
            return _count_no_gap_sites_in_identical_sequence(
                first_sequence,
                is_protein,
            )

        scalar_no_gap_sites = _count_no_gap_sites_scalar(
            sequences,
            aln_len,
            is_protein,
        )
        if scalar_no_gap_sites is not None:
            return scalar_no_gap_sites

        try:
            alignment_bytes = "".join(sequences).encode("ascii")
            gap_bytes = _PROTEIN_GAP_BYTES if is_protein else _DNA_GAP_BYTES
            if not any(gap_code in alignment_bytes for gap_code in gap_bytes):
                return aln_len
            if _should_use_dense_gap_matrix_scan(alignment_bytes, gap_bytes):
                return _count_columns_without_gap_bytes_matrix(
                    alignment_bytes,
                    aln_len,
                    gap_bytes,
                )
            return _count_columns_without_gap_bytes(
                alignment_bytes,
                aln_len,
                gap_bytes,
            )
        except UnicodeEncodeError:
            gap_chars = {char.upper() for char in self.get_gap_chars(is_protein)}
            sequences = [seq.upper() for seq in sequences]
            alignment_array = np.array([list(seq) for seq in sequences], dtype="U1")
            gap_chars_array = np.array(list(gap_chars), dtype="U1")
            columns_with_gaps = np.any(
                np.isin(alignment_array, gap_chars_array),
                axis=0,
            )
        return int(np.count_nonzero(~columns_with_gaps))
