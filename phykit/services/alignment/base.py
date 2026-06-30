from __future__ import annotations

from ..base import BaseService


def get_alignment_and_format_helper(*args, **kwargs):
    from ...helpers.files import get_alignment_and_format

    return get_alignment_and_format(*args, **kwargs)


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()
_DNA_INVALID_LOOKUP = None
_PROTEIN_INVALID_LOOKUP = None


def _all_sequences_identical(sequences) -> bool:
    if len(sequences) < 2:
        return True
    first_sequence = sequences[0]
    for idx in range(1, len(sequences)):
        if sequences[idx] != first_sequence:
            return False
    return True


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


class Alignment(BaseService):
    def __init__(
        self,
        *args,
        alignment_file_path=None,
        code=None,
        fasta=None,
        output_file_path=None,
        protein_file_path=None,
        nucleotide_file_path=None,
        alignment_list_path=None,
        prefix=None,
        idmap=None,
        reference=None,
        verbose=None,
        entry=None,
        exclude_gaps=None,
    ):
        self.alignment_file_path = alignment_file_path
        self.code = code
        self.output_file_path = output_file_path
        self.protein_file_path = protein_file_path
        self.nucleotide_file_path = nucleotide_file_path
        self.alignment_list_path = alignment_list_path
        self.prefix = prefix
        self.fasta = fasta
        self.idmap = idmap
        self.reference = reference
        self.verbose = verbose
        self.entry = entry
        self.exclude_gaps = exclude_gaps

    def get_alignment_and_format(self):
        """
        automatic file type determination
        """
        return get_alignment_and_format_helper(self.alignment_file_path)

    def calculate_rcv(self) -> float:
        alignment, _, is_protein = self.get_alignment_and_format()
        num_records = len(alignment)

        if num_records == 0:
            return 0.0

        sequences = [str(record.seq).upper() for record in alignment]
        if num_records == 1:
            return 0.0
        if _all_sequences_identical(sequences):
            return 0.0

        aln_len = alignment.get_alignment_length()

        if is_protein:
            invalid_chars = ["-", "?", "*", "X"]
        else:
            invalid_chars = ["-", "?", "*", "X", "N"]

        try:
            alignment_bytes = "".join(sequences).encode("ascii")
            alignment_array = np.frombuffer(
                alignment_bytes,
                dtype=np.uint8,
            ).reshape(num_records, aln_len)
            invalid_lookup = _get_invalid_lookup(is_protein)
            observed_chars = np.unique(alignment_array)
            unique_chars = observed_chars[~invalid_lookup[observed_chars]]
            if unique_chars.size == 0:
                return 0.0
            if unique_chars.size == observed_chars.size:
                valid_mask = None
                valid_lengths = np.full(num_records, aln_len, dtype=np.float64)
            else:
                valid_mask = ~invalid_lookup[alignment_array]
                valid_lengths = np.count_nonzero(valid_mask, axis=1).astype(np.float64)
        except UnicodeEncodeError:
            alignment_array = np.array([list(seq) for seq in sequences], dtype="U1")
            invalid_chars_array = np.array(invalid_chars, dtype="U1")
            valid_mask = ~np.isin(alignment_array, invalid_chars_array)
            valid_lengths = np.count_nonzero(valid_mask, axis=1).astype(np.float64)

            # Get all unique valid characters in the alignment
            valid_chars = alignment_array[valid_mask]
            if valid_chars.size == 0:
                return 0.0

            unique_chars = np.unique(valid_chars)
        if alignment_array.dtype == np.uint8 and len(unique_chars) > 8:
            count_matrix = np.zeros(
                (num_records, len(unique_chars)),
                dtype=np.float64,
            )
            for row_idx, row in enumerate(alignment_array):
                row_values = row if valid_mask is None else row[valid_mask[row_idx]]
                count_matrix[row_idx] = np.bincount(
                    row_values,
                    minlength=256,
                )[unique_chars]
        else:
            count_matrix = np.array(
                [
                    np.sum(alignment_array == char, axis=1)
                    for char in unique_chars
                ],
                dtype=np.float64,
            ).T

        # Calculate total counts and averages using matrix operations
        total_counts = np.sum(count_matrix, axis=0)
        average_counts = total_counts / num_records

        # Calculate RCV values using vectorized operations
        # Compute absolute differences from average for all sequences at once
        abs_diffs = np.abs(count_matrix - average_counts)

        # Sum across characters for each sequence
        seq_rcv_sums = np.sum(abs_diffs, axis=1)

        # Normalize each sequence by its valid (non-gap/non-ambiguous) length.
        # Sequences with no valid symbols contribute 0.
        denom = num_records * valid_lengths
        indiv_rcv_values = np.divide(
            seq_rcv_sums,
            denom,
            out=np.zeros_like(seq_rcv_sums, dtype=np.float64),
            where=denom > 0,
        )

        return float(np.sum(indiv_rcv_values))

    @staticmethod
    def get_gap_chars(is_protein: bool = False) -> list[str]:
        if is_protein:
            return ["-", "?", "*", "X", "x"]
        else:
            return ["-", "?", "*", "X", "x", "N", "n"]
