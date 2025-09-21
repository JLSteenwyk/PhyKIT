from typing import Dict, Tuple
import numpy as np

from Bio.Align import MultipleSeqAlignment

from .base import Alignment


class VariableSites(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _, is_protein = self.get_alignment_and_format()
        var_sites, aln_len, var_sites_per = \
            self.calculate_variable_sites(alignment)

        print(f"{var_sites}\t{aln_len}\t{round(var_sites_per, 4)}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(alignment_file_path=args.alignment)

    def calculate_variable_sites(
        self,
        alignment: MultipleSeqAlignment
    ) -> Tuple[int, int, float]:
        aln_len = alignment.get_alignment_length()
        gap_chars = self.get_gap_chars()

        # Convert alignment to numpy array for vectorized operations
        alignment_array = np.array([
            [c.upper() for c in str(record.seq)]
            for record in alignment
        ], dtype='U1')

        var_sites = 0

        # Process each column
        for col in range(aln_len):
            column = alignment_array[:, col]

            # Filter out gap characters
            non_gap_mask = ~np.isin(column, list(gap_chars))
            filtered_column = column[non_gap_mask]

            # Check if variable (more than one unique character)
            if len(filtered_column) > 0:
                unique_chars = np.unique(filtered_column)
                if len(unique_chars) > 1:
                    var_sites += 1

        var_sites_per = (var_sites / aln_len) * 100

        return var_sites, aln_len, var_sites_per
