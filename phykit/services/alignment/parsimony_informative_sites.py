from collections import Counter
from typing import Dict, Tuple
import numpy as np

from Bio.Align import MultipleSeqAlignment

from .base import Alignment


class ParsimonyInformative(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _, is_protein = self.get_alignment_and_format()
        pi_sites, aln_len, pi_sites_per = self.calculate_parsimony_informative_sites(
            alignment
        )

        print(f"{pi_sites}\t{aln_len}\t{round(pi_sites_per, 4)}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(alignment_file_path=args.alignment)

    def get_number_of_occurrences_per_character(
        self,
        alignment: MultipleSeqAlignment,
        idx: int
    ) -> Counter:
        gap_chars = self.get_gap_chars()
        seq_at_position = alignment[:, idx].upper()
        filtered_seq = filter(lambda c: c not in gap_chars, seq_at_position)

        return Counter(filtered_seq)

    def is_parsimony_informative(
        self,
        num_occurrences: Counter,
    ) -> bool:
        """
        Check if a site is parsimony informative.
        That is, the site has two characters that appear at least twice.
        """
        informative_char_count = sum(1 for count in num_occurrences.values() if count >= 2)
        return informative_char_count >= 2

    def calculate_parsimony_informative_sites(
        self,
        alignment: MultipleSeqAlignment,
    ) -> Tuple[int, int, float]:
        aln_len = alignment.get_alignment_length()
        gap_chars = self.get_gap_chars()

        # Convert alignment to numpy array for vectorized operations
        alignment_array = np.array([
            [c.upper() for c in str(record.seq)]
            for record in alignment
        ], dtype='U1')

        pi_sites = 0

        # Process each column more efficiently
        for col_idx in range(aln_len):
            column = alignment_array[:, col_idx]

            # Filter out gaps
            non_gap_mask = ~np.isin(column, list(gap_chars))
            filtered_column = column[non_gap_mask]

            if len(filtered_column) > 0:
                # Count occurrences of each character
                unique_chars, counts = np.unique(filtered_column, return_counts=True)

                # Check if parsimony informative (at least 2 chars appearing >= 2 times)
                chars_appearing_twice = np.sum(counts >= 2)
                if chars_appearing_twice >= 2:
                    pi_sites += 1

        pi_sites_per = (pi_sites / aln_len) * 100

        return pi_sites, aln_len, pi_sites_per
