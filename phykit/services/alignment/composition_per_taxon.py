from typing import Dict, List, Tuple

import numpy as np

from .base import Alignment


class CompositionPerTaxon(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()
        symbols, rows = self.calculate_composition_per_taxon(alignment, is_protein)
        if not symbols:
            return

        for taxon, comps in rows:
            comp_str = ";".join(
                f"{symbol}:{round(comps[idx], 4)}"
                for idx, symbol in enumerate(symbols)
            )
            print(f"{taxon}\t{comp_str}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(alignment_file_path=args.alignment)

    def calculate_composition_per_taxon(
        self, alignment, is_protein: bool
    ) -> Tuple[List[str], List[Tuple[str, np.ndarray]]]:
        if is_protein:
            invalid_chars = np.array(["-", "?", "*", "X"], dtype="U1")
        else:
            invalid_chars = np.array(["-", "?", "*", "X", "N"], dtype="U1")

        alignment_array = np.array(
            [[c.upper() for c in str(record.seq)] for record in alignment], dtype="U1"
        )
        valid_mask = ~np.isin(alignment_array, invalid_chars)
        valid_chars = alignment_array[valid_mask]

        if valid_chars.size == 0:
            return [], []

        symbols = sorted(np.unique(valid_chars).tolist())
        output: List[Tuple[str, np.ndarray]] = []
        for record_idx, record in enumerate(alignment):
            seq = alignment_array[record_idx]
            seq_valid = valid_mask[record_idx]
            valid_len = int(np.sum(seq_valid))
            freqs = np.zeros(len(symbols), dtype=np.float64)
            if valid_len > 0:
                for sym_idx, symbol in enumerate(symbols):
                    freqs[sym_idx] = np.sum((seq == symbol) & seq_valid) / valid_len
            output.append((record.id, freqs))

        return symbols, output
