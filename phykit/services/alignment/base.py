from collections import Counter
import sys
import numpy as np

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
        self.code = code,
        self.output_file_path = output_file_path
        self.protein_file_path = (protein_file_path,)
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
        alignment, _, _ = self.get_alignment_and_format()
        aln_len = alignment.get_alignment_length()
        num_records = len(alignment)

        # Convert alignment to numpy array for faster operations
        alignment_array = np.array([
            list(str(record.seq)) for record in alignment
        ], dtype='U1')

        # Get all unique characters in the alignment
        unique_chars = np.unique(alignment_array)

        # Vectorized approach: create a count matrix for all sequences and characters at once
        # Shape: (num_records, num_unique_chars)
        count_matrix = np.zeros((num_records, len(unique_chars)), dtype=np.int32)

        # Build character index mapping for fast lookup
        char_to_idx = {char: idx for idx, char in enumerate(unique_chars)}

        # Count characters for each sequence using vectorized operations
        for seq_idx in range(num_records):
            seq = alignment_array[seq_idx]
            for char_idx, char in enumerate(unique_chars):
                count_matrix[seq_idx, char_idx] = np.sum(seq == char)

        # Calculate total counts and averages using matrix operations
        total_counts = np.sum(count_matrix, axis=0)
        average_counts = total_counts / num_records

        # Calculate RCV values using vectorized operations
        # Compute absolute differences from average for all sequences at once
        abs_diffs = np.abs(count_matrix - average_counts)

        # Sum across characters for each sequence
        seq_rcv_sums = np.sum(abs_diffs, axis=1)

        # Normalize and sum
        indiv_rcv_values = seq_rcv_sums / (num_records * aln_len)
        relative_composition_variability = np.sum(indiv_rcv_values)

        return float(relative_composition_variability)

    def get_gap_chars(is_protein: bool) -> List[str]:
        if is_protein:
            return ["-", "?", "*", "X", "x"]
        else:
            return ["-", "?", "*", "X", "x", "N", "n"]
