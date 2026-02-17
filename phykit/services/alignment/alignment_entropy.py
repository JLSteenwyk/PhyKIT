from typing import Dict, List

import numpy as np

from .base import Alignment


class AlignmentEntropy(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()
        entropies = self.calculate_site_entropies(alignment, is_protein)

        if self.verbose:
            for idx, entropy in enumerate(entropies, start=1):
                print(f"{idx}\t{round(entropy, 4)}")
            return

        if entropies:
            print(round(float(np.mean(entropies)), 4))
        else:
            print(0.0)

    def process_args(self, args) -> Dict[str, str]:
        return dict(alignment_file_path=args.alignment, verbose=args.verbose)

    def calculate_site_entropies(self, alignment, is_protein: bool) -> List[float]:
        if is_protein:
            invalid_chars = np.array(["-", "?", "*", "X"], dtype="U1")
        else:
            invalid_chars = np.array(["-", "?", "*", "X", "N"], dtype="U1")
        alignment_array = np.array(
            [[c.upper() for c in str(record.seq)] for record in alignment],
            dtype="U1",
        )

        site_entropies: List[float] = []
        for col_idx in range(alignment_array.shape[1]):
            column = alignment_array[:, col_idx]
            valid = column[~np.isin(column, invalid_chars)]
            if valid.size == 0:
                site_entropies.append(0.0)
                continue

            _, counts = np.unique(valid, return_counts=True)
            probs = counts / counts.sum()
            entropy = float(-np.sum(probs * np.log2(probs)))
            site_entropies.append(entropy)

        return site_entropies
