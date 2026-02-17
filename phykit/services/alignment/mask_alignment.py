import sys
from typing import Dict

import numpy as np

from .base import Alignment


class MaskAlignment(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.max_gap = parsed["max_gap"]
        self.min_occupancy = parsed["min_occupancy"]
        self.max_entropy = parsed["max_entropy"]

    def run(self) -> None:
        self._validate_thresholds()
        alignment, _, is_protein = self.get_alignment_and_format()
        keep_mask = self.calculate_keep_mask(alignment, is_protein)
        masked = self.apply_mask(alignment, keep_mask)

        for taxon, seq in masked.items():
            print(f">{taxon}\n{seq}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            max_gap=args.max_gap,
            min_occupancy=args.min_occupancy,
            max_entropy=args.max_entropy,
        )

    def _validate_thresholds(self) -> None:
        if self.max_gap < 0.0 or self.max_gap > 1.0:
            print("max_gap must be between 0 and 1.")
            sys.exit(2)
        if self.min_occupancy < 0.0 or self.min_occupancy > 1.0:
            print("min_occupancy must be between 0 and 1.")
            sys.exit(2)
        if self.max_entropy is not None and self.max_entropy < 0.0:
            print("max_entropy must be >= 0.")
            sys.exit(2)

    def calculate_keep_mask(self, alignment, is_protein: bool) -> np.ndarray:
        if is_protein:
            invalid_chars = np.array(["-", "?", "*", "X"], dtype="U1")
        else:
            invalid_chars = np.array(["-", "?", "*", "X", "N"], dtype="U1")

        alignment_array = np.array(
            [[c.upper() for c in str(record.seq)] for record in alignment], dtype="U1"
        )
        valid_mask = ~np.isin(alignment_array, invalid_chars)

        occupancy = np.mean(valid_mask, axis=0).astype(np.float64)
        gap_fraction = 1.0 - occupancy

        keep_mask = (gap_fraction <= self.max_gap) & (occupancy >= self.min_occupancy)

        if self.max_entropy is not None:
            entropies = np.zeros(alignment_array.shape[1], dtype=np.float64)
            for col_idx in range(alignment_array.shape[1]):
                column = alignment_array[:, col_idx]
                valid = column[valid_mask[:, col_idx]]
                if valid.size == 0:
                    entropies[col_idx] = 0.0
                    continue
                _, counts = np.unique(valid, return_counts=True)
                probs = counts / counts.sum()
                entropies[col_idx] = float(-np.sum(probs * np.log2(probs)))
            keep_mask &= entropies <= self.max_entropy

        return keep_mask

    def apply_mask(self, alignment, keep_mask: np.ndarray) -> Dict[str, str]:
        output = {}
        for record in alignment:
            seq_arr = np.array(list(str(record.seq).upper()), dtype="U1")
            output[record.id] = "".join(seq_arr[keep_mask].tolist())
        return output
