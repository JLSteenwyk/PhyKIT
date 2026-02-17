import sys

from typing import List

from ..base import BaseService
from ...helpers.files import (
    get_alignment_and_format as get_alignment_and_format_helper
)


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
        try:
            return get_alignment_and_format_helper(self.alignment_file_path)
        except FileNotFoundError:
            print("Input corresponds to no such file or directory.")
            print("Please double check pathing and filenames")
            sys.exit(2)

    def calculate_rcv(self) -> float:
        import numpy as np

        alignment, _, is_protein = self.get_alignment_and_format()
        num_records = len(alignment)

        if num_records == 0:
            return 0.0

        # Convert alignment to numpy array for faster operations
        alignment_array = np.array([
            list(str(record.seq).upper()) for record in alignment
        ], dtype='U1')

        if is_protein:
            invalid_chars = np.array(["-", "?", "*", "X"], dtype='U1')
        else:
            invalid_chars = np.array(["-", "?", "*", "X", "N"], dtype='U1')
        valid_mask = ~np.isin(alignment_array, invalid_chars)
        valid_lengths = np.sum(valid_mask, axis=1).astype(np.float64)

        # Get all unique valid characters in the alignment
        valid_chars = alignment_array[valid_mask]
        if valid_chars.size == 0:
            return 0.0

        unique_chars = np.unique(valid_chars)

        # Shape: (num_records, num_unique_chars)
        count_matrix = np.zeros((num_records, len(unique_chars)), dtype=np.int32)

        # Count characters for each sequence using vectorized operations
        for seq_idx, seq in enumerate(alignment_array):
            seq_valid = valid_mask[seq_idx]
            for char_idx, char in enumerate(unique_chars):
                count_matrix[seq_idx, char_idx] = np.sum((seq == char) & seq_valid)

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

    def get_gap_chars(is_protein: bool) -> List[str]:
        if is_protein:
            return ["-", "?", "*", "X", "x"]
        else:
            return ["-", "?", "*", "X", "x", "N", "n"]
