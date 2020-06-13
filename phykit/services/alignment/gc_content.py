import re

from Bio import SeqIO

from .base import Alignment

class GCContent(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # create biopython object of sequences
        records = SeqIO.parse(self.fasta, "fasta")
        
        # initialize and populate dict for
        # holding the entry sequences
        entry_and_seq = {}
        for record in records:
            entry_and_seq[record.id] = str(record.seq)

        regex_pattern = re.compile('[GgCc]')
        if self.verbose:
            for entry, seq in entry_and_seq.items():
                seq = seq.replace('-', '')
                matches = regex_pattern.findall(seq)
                print(f"{entry}\t{len(matches)/len(seq)}")
        else:
            all_seqs = ''
            for entry, seq in entry_and_seq.items():
                all_seqs += seq
            all_seqs = all_seqs.replace('-', '')
            matches = regex_pattern.findall(all_seqs)
            print(matches)
            print(f"{len(matches)/len(all_seqs)}")
            


    def process_args(self, args):
        return dict(
            fasta=args.fasta, 
            verbose=args.verbose
        )