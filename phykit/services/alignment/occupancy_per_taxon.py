from typing import Dict

import numpy as np

from .base import Alignment
from ...helpers.json_output import print_json


class OccupancyPerTaxon(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()
        occupancies = self.calculate_occupancy_per_taxon(alignment, is_protein)

        if self.json_output:
            rows = [
                dict(taxon=taxon, occupancy=round(occupancy, 4))
                for taxon, occupancy in occupancies
            ]
            print_json(
                dict(
                    rows=rows,
                    taxa=rows,
                )
            )
            return

        for taxon, occupancy in occupancies:
            print(f"{taxon}\t{round(occupancy, 4)}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
        )

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
