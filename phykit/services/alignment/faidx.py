from Bio import SeqIO

from .base import Alignment

class Faidx(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        record_dict = SeqIO.index(self.fasta, "fasta")
        print(f">{record_dict[self.entry].name}\n{record_dict[self.entry].seq}")

    def process_args(self, args):
        return dict(fasta=args.fasta, entry=args.entry)
