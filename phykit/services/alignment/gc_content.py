from enum import Enum
import re
import sys
from typing import Dict, Tuple
from collections import Counter

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
            sys.exit()

        if self.verbose:
            self.calculate_gc_per_sequence(records)
        else:
            self.calculate_gc_total(records)

    def process_args(self, args) -> Dict[str, str]:
        return dict(fasta=args.fasta, verbose=args.verbose)

    def calculate_gc_per_sequence(self, records: MultipleSeqAlignment) -> None:
        for record in records:
            seq = str(record.seq)
            seq, gc_count = self.remove_gaps_and_count_gc(seq)
            try:
                print(f"{record.id}\t{round(gc_count / len(seq), 4) if seq else 0}")
            except BrokenPipeError:
                pass

    def calculate_gc_total(self, records: MultipleSeqAlignment) -> None:
        combined_seq = "".join(str(record.seq) for record in records)
        combined_seq, gc_count = self.remove_gaps_and_count_gc(combined_seq)
        try:
            gc_content = round(gc_count / len(combined_seq), 4)
            print(gc_content)
        except (BrokenPipeError, ZeroDivisionError):
            print(
                "Input file has an unacceptable format. Please check input file argument."
            )
            sys.exit()

    def remove_gaps_and_count_gc(self, seq: str) -> Tuple[str, float]:
        gap_chars = self.get_gap_chars()
        pattern = "[" + "".join(re.escape(char) for char in gap_chars) + "]"
        cleaned_seq = re.sub(pattern, "", seq)
        gc_count = Counter(cleaned_seq.upper())["G"] + Counter(cleaned_seq.upper())["C"]

        return cleaned_seq, gc_count
