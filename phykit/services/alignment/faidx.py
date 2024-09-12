from typing import Dict

from Bio import SeqIO

from .base import Alignment


class Faidx(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        record_dict = SeqIO.index(self.fasta, "fasta")

        # Split entries and iterate
        for e in map(str.strip, self.entry.split(",")):
            record = record_dict[e]
            print(f">{record.name}\n{record.seq}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(fasta=args.fasta, entry=args.entry)
