from typing import Dict

import numpy as np

from .base import Alignment


class OccupancyPerTaxon(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()
        occupancies = self.calculate_occupancy_per_taxon(alignment, is_protein)

        for taxon, occupancy in occupancies:
            print(f"{taxon}\t{round(occupancy, 4)}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(alignment_file_path=args.alignment)

    def calculate_occupancy_per_taxon(self, alignment, is_protein: bool):
        if is_protein:
            invalid_chars = np.array(["-", "?", "*", "X"], dtype="U1")
        else:
            invalid_chars = np.array(["-", "?", "*", "X", "N"], dtype="U1")

        output = []
        for record in alignment:
            seq_arr = np.array(list(str(record.seq).upper()), dtype="U1")
            valid_count = int(np.sum(~np.isin(seq_arr, invalid_chars)))
            occupancy = (valid_count / len(seq_arr)) if len(seq_arr) > 0 else 0.0
            output.append((record.id, occupancy))
        return output
