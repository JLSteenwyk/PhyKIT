from __future__ import annotations

from argparse import Namespace

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
_DNA_GAP_BYTES = b"-?*XNxn"
_PROTEIN_GAP_BYTES = b"-?*Xx"
_DNA_GAP_CHARS = {"-", "?", "*", "X", "N"}
_PROTEIN_GAP_CHARS = {"-", "?", "*", "X"}
_DNA_GAP_COUNT_CHARS = "-?*XNxn"
_PROTEIN_GAP_COUNT_CHARS = "-?*Xx"


def _get_gap_codes(is_protein: bool):
    global _DNA_GAP_CODES, _PROTEIN_GAP_CODES
    if is_protein:
        if _PROTEIN_GAP_CODES is None:
            _PROTEIN_GAP_CODES = np.frombuffer(b"-?*Xx", dtype=np.uint8)
        return _PROTEIN_GAP_CODES

    if _DNA_GAP_CODES is None:
        _DNA_GAP_CODES = np.frombuffer(b"-?*XNxn", dtype=np.uint8)
    return _DNA_GAP_CODES


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
    if not sequences:
        return True

    first_sequence = sequences[0]
    for idx in range(1, len(sequences)):
        if sequences[idx] != first_sequence:
            return False
    return True


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
        args: Namespace,
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
        gap_chars = {char.upper() for char in self.get_gap_chars(is_protein)}
        sequences = [str(record.seq) for record in alignment]

        if not sequences:
            return aln_len

        first_sequence = sequences[0]
        if _all_sequences_identical(sequences):
            return _count_no_gap_sites_in_identical_sequence(
                first_sequence,
                is_protein,
            )

        try:
            alignment_bytes = "".join(sequences).encode("ascii")
            gap_bytes = _PROTEIN_GAP_BYTES if is_protein else _DNA_GAP_BYTES
            if not any(gap_code in alignment_bytes for gap_code in gap_bytes):
                return aln_len
            gap_codes = _get_gap_codes(is_protein)
            alignment_array = np.frombuffer(
                alignment_bytes,
                dtype=np.uint8,
            ).reshape(len(sequences), aln_len)
            columns_with_gaps = np.zeros(aln_len, dtype=np.bool_)
            for gap_code in gap_codes:
                columns_with_gaps |= np.any(alignment_array == gap_code, axis=0)
        except UnicodeEncodeError:
            sequences = [seq.upper() for seq in sequences]
            alignment_array = np.array([list(seq) for seq in sequences], dtype="U1")
            gap_chars_array = np.array(list(gap_chars), dtype="U1")
            columns_with_gaps = np.any(
                np.isin(alignment_array, gap_chars_array),
                axis=0,
            )
        return int(np.count_nonzero(~columns_with_gaps))
