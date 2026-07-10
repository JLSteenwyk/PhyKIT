from __future__ import annotations

from ..base import BaseService


def get_alignment_and_format_helper(*args, **kwargs):
    from ...helpers.files import get_alignment_and_format

    return get_alignment_and_format(*args, **kwargs)


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
_DNA_INVALID_LOOKUP = None
_PROTEIN_INVALID_LOOKUP = None
_ASCII_RCV_GLOBAL_BINCOUNT_MIN_RECORDS = 2_048
_ASCII_RCV_GLOBAL_BINCOUNT_MAX_LENGTH = 512
_NUCLEOTIDE_ALPHABET_BYTES = b"ACGTU"


def _all_sequences_identical(sequences) -> bool:
    empty = object()
    iterator = iter(sequences)
    first_sequence = next(iterator, empty)
    if first_sequence is empty:
        return True
    for sequence in iterator:
        if sequence != first_sequence:
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


def _ascii_rcv_count_matrix(alignment_array, unique_chars, valid_mask):
    num_records, aln_len = alignment_array.shape
    if (
        num_records >= _ASCII_RCV_GLOBAL_BINCOUNT_MIN_RECORDS
        and aln_len <= _ASCII_RCV_GLOBAL_BINCOUNT_MAX_LENGTH
    ):
        encoded = alignment_array.astype(np.int64)
        encoded += (np.arange(num_records, dtype=np.int64) * 256)[:, None]
        counts = np.bincount(
            encoded.ravel(),
            minlength=num_records * 256,
        ).reshape(num_records, 256)
        return counts[:, unique_chars].astype(np.float64, copy=False)

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
    return count_matrix


def _rcv_row_sums(abs_diffs):
    if abs_diffs.shape[1] <= 4:
        return abs_diffs.sum(axis=1)
    return np.sum(abs_diffs, axis=1)


def _rcv_column_totals(count_matrix):
    if count_matrix.shape[0] <= 1200:
        return count_matrix.sum(axis=0)
    return np.sum(count_matrix, axis=0)


def _bounded_ascii_count_dtype(max_count: int):
    if max_count <= 0xFFFF:
        return np.uint16
    if max_count <= 0xFFFFFFFF:
        return np.uint32
    return np.uint64


def _bounded_ascii_valid_lengths(valid_mask):
    count_dtype = _bounded_ascii_count_dtype(valid_mask.shape[1])
    return valid_mask.sum(axis=1, dtype=count_dtype).astype(np.float64)


def _bounded_ascii_row_symbol_counts(alignment_array, unique_chars):
    count_dtype = _bounded_ascii_count_dtype(alignment_array.shape[1])
    return np.array(
        [
            (alignment_array == char).sum(axis=1, dtype=count_dtype)
            for char in unique_chars
        ],
        dtype=np.float64,
    ).T


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
        if num_records == 1:
            return 0.0

        raw_sequences = [str(record.seq) for record in alignment]
        first_raw_sequence = raw_sequences[0]
        first_sequence = first_raw_sequence.upper()
        all_identical = True
        for idx in range(1, num_records):
            sequence = raw_sequences[idx]
            if sequence != first_raw_sequence and sequence.upper() != first_sequence:
                all_identical = False
                break

        if all_identical:
            return 0.0
        sequences = [sequence.upper() for sequence in raw_sequences]

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
            if (
                not is_protein
                and not alignment_bytes.translate(
                    None,
                    _NUCLEOTIDE_ALPHABET_BYTES,
                )
            ):
                observed_chars = np.fromiter(
                    (
                        code
                        for code in _NUCLEOTIDE_ALPHABET_BYTES
                        if code in alignment_bytes
                    ),
                    dtype=np.uint8,
                )
            else:
                observed_chars = np.unique(alignment_array)
            unique_chars = observed_chars[~invalid_lookup[observed_chars]]
            if unique_chars.size == 0:
                return 0.0
            if unique_chars.size == observed_chars.size:
                valid_mask = None
                valid_lengths = np.full(num_records, aln_len, dtype=np.float64)
            else:
                valid_mask = ~invalid_lookup[alignment_array]
                valid_lengths = _bounded_ascii_valid_lengths(valid_mask)
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
            count_matrix = _ascii_rcv_count_matrix(
                alignment_array,
                unique_chars,
                valid_mask,
            )
        elif alignment_array.dtype == np.uint8:
            count_matrix = _bounded_ascii_row_symbol_counts(
                alignment_array,
                unique_chars,
            )
        else:
            count_matrix = np.array(
                [
                    np.sum(alignment_array == char, axis=1)
                    for char in unique_chars
                ],
                dtype=np.float64,
            ).T

        # Calculate total counts and averages using matrix operations
        total_counts = _rcv_column_totals(count_matrix)
        average_counts = total_counts / num_records

        # Calculate RCV values using vectorized operations
        # Compute absolute differences from average for all sequences at once
        abs_diffs = np.abs(count_matrix - average_counts)

        # Sum across characters for each sequence
        seq_rcv_sums = _rcv_row_sums(abs_diffs)

        # Normalize each sequence by its valid (non-gap/non-ambiguous) length.
        # Sequences with no valid symbols contribute 0.
        denom = num_records * valid_lengths
        indiv_rcv_values = np.divide(
            seq_rcv_sums,
            denom,
            out=np.zeros_like(seq_rcv_sums, dtype=np.float64),
            where=denom > 0,
        )

        return float(indiv_rcv_values.sum())

    @staticmethod
    def get_gap_chars(is_protein: bool = False) -> list[str]:
        if is_protein:
            return ["-", "?", "*", "X", "x"]
        else:
            return ["-", "?", "*", "X", "x", "N", "n"]
