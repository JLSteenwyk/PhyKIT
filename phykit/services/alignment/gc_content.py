from enum import Enum
import re
import sys
from typing import Dict, Tuple
from collections import Counter
import numpy as np

from Bio.Align import MultipleSeqAlignment

from .base import Alignment
from ...helpers.files import get_alignment_and_format


class FileFormat(Enum):
    fasta = "fasta"
    clustal = "clustal"
    maf = "maf"
    mauve = "mauve"
    phylip = "phylip"
    phylip_seq = "phylip-sequential"
    stockholm = "stockholm"


class GCContent(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        records, _, is_protein = get_alignment_and_format(self.fasta)

        if is_protein:
            print("GC content can't be calculated for protein sequences")
            sys.exit(2)

        if self.verbose:
            self.calculate_gc_per_sequence(records)
        else:
            self.calculate_gc_total(records)

    def process_args(self, args) -> Dict[str, str]:
        return dict(fasta=args.fasta, verbose=args.verbose)

    def calculate_gc_per_sequence(self, records: MultipleSeqAlignment) -> None:
        gap_chars = set(self.get_gap_chars())

        for record in records:
            # Convert to numpy array for faster operations
            seq_arr = np.array(list(str(record.seq).upper()), dtype='U1')

            # Filter out gaps
            non_gap_mask = ~np.isin(seq_arr, list(gap_chars))
            cleaned_seq = seq_arr[non_gap_mask]

            if len(cleaned_seq) > 0:
                # Count G and C
                gc_count = np.sum((cleaned_seq == 'G') | (cleaned_seq == 'C'))
                gc_content = gc_count / len(cleaned_seq)
            else:
                gc_content = 0

            try:
                print(f"{record.id}\t{round(gc_content, 4)}")
            except BrokenPipeError:
                pass

    def calculate_gc_total(self, records: MultipleSeqAlignment) -> None:
        gap_chars = set(self.get_gap_chars())

        # Combine all sequences into one array
        all_seqs = [list(str(record.seq).upper()) for record in records]
        combined_arr = np.concatenate([np.array(seq, dtype='U1') for seq in all_seqs])

        # Filter out gaps
        non_gap_mask = ~np.isin(combined_arr, list(gap_chars))
        cleaned_seq = combined_arr[non_gap_mask]

        if len(cleaned_seq) > 0:
            # Count G and C
            gc_count = np.sum((cleaned_seq == 'G') | (cleaned_seq == 'C'))
            gc_content = round(gc_count / len(cleaned_seq), 4)
            print(gc_content)
        else:
            print(
                "Input file has an unacceptable format. Please check input file argument."
            )
            sys.exit(2)

    def remove_gaps_and_count_gc(self, seq: str) -> Tuple[str, float]:
        gap_chars = self.get_gap_chars()
        pattern = "[" + "".join(re.escape(char) for char in gap_chars) + "]"
        cleaned_seq = re.sub(pattern, "", seq)
        gc_count = Counter(cleaned_seq.upper())["G"] + Counter(cleaned_seq.upper())["C"]

        return cleaned_seq, gc_count
