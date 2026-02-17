import numpy as np

from .base import Alignment
from ...helpers.json_output import print_json


class RelativeCompositionVariabilityTaxon(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.json_output = parsed["json_output"]

    def run(self):
        alignment, _, is_protein = self.get_alignment_and_format()
        num_records = len(alignment)

        if num_records == 0:
            return

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

        # Get all unique valid symbols
        valid_chars = alignment_array[valid_mask]
        if valid_chars.size == 0:
            if self.json_output:
                rows = [dict(taxon=record.id, rcvt=0.0) for record in alignment]
                print_json(dict(rows=rows, taxa=rows))
                return
            for record in alignment:
                print(f"{record.id}\t0.0")
            return

        unique_chars = np.unique(valid_chars)

        # Count valid symbols for each sequence and character
        count_matrix = np.zeros((num_records, len(unique_chars)), dtype=np.float32)
        for seq_idx, seq in enumerate(alignment_array):
            seq_valid = valid_mask[seq_idx]
            for char_idx, char in enumerate(unique_chars):
                count_matrix[seq_idx, char_idx] = np.sum((seq == char) & seq_valid)

        # Calculate average counts per sequence (total counts / num_records)
        average_counts = np.sum(count_matrix, axis=0) / num_records

        # Vectorized RCV calculation for all sequences at once
        deviations = np.abs(count_matrix - average_counts)
        seq_sums = np.sum(deviations, axis=1)
        denom = num_records * valid_lengths
        rcv_values = np.divide(
            seq_sums,
            denom,
            out=np.zeros_like(seq_sums, dtype=np.float64),
            where=denom > 0,
        )

        if self.json_output:
            rows = [
                dict(taxon=record.id, rcvt=round(float(rcv_values[i]), 4))
                for i, record in enumerate(alignment)
            ]
            print_json(
                dict(
                    rows=rows,
                    taxa=rows,
                )
            )
            return

        # Print results - convert to float64 for consistent rounding
        for i, record in enumerate(alignment):
            rcv_val = float(rcv_values[i])
            print(f"{record.id}\t{round(rcv_val, 4)}")

    def process_args(self, args):
        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
        )
