import numpy as np

from .base import Alignment


class RelativeCompositionVariabilityTaxon(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _, _ = self.get_alignment_and_format()
        aln_len = alignment.get_alignment_length()
        num_records = len(alignment)

        # Convert alignment to numpy array for faster operations
        alignment_array = np.array([
            list(str(record.seq)) for record in alignment
        ], dtype='U1')

        # Get all unique characters and create mapping
        unique_chars = np.unique(alignment_array)
        char_to_idx = {char: i for i, char in enumerate(unique_chars)}

        # Create integer representation for faster counting
        alignment_int = np.zeros_like(alignment_array, dtype=np.int8)
        for char, idx in char_to_idx.items():
            alignment_int[alignment_array == char] = idx

        # Vectorized counting for all sequences and characters
        count_matrix = np.zeros((num_records, len(unique_chars)), dtype=np.float32)
        for i in range(len(unique_chars)):
            count_matrix[:, i] = np.sum(alignment_int == i, axis=1)

        # Calculate average counts per sequence (total counts / num_records)
        average_counts = np.sum(count_matrix, axis=0) / num_records

        # Vectorized RCV calculation for all sequences at once
        deviations = np.abs(count_matrix - average_counts)
        rcv_values = np.sum(deviations, axis=1) / (num_records * aln_len)

        # Print results - convert to float64 for consistent rounding
        for i, record in enumerate(alignment):
            rcv_val = float(rcv_values[i])
            print(f"{record.id}\t{round(rcv_val, 4)}")

    def process_args(self, args):
        return dict(alignment_file_path=args.alignment)
