from enum import Enum
import re
import sys

from Bio import SeqIO

from .base import Alignment

class FileFormat(Enum):
    fasta = "fasta"
    clustal = "clustal"
    maf = "maf"
    mauve = "mauve"
    phylip = "phylip"
    phylip_seq = "phylip-sequential"
    phylip_rel = "phylip-relaxed"
    stockholm = "stockholm"

class GCContent(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # create biopython object of sequences
        records, file_format = self.determine_file_type_fasta_file(self.fasta)
        
        # initialize and populate dict for
        # holding the entry sequences
        entry_and_seq = {}
        for record in records:
            entry_and_seq[record.id] = str(record.seq)

        if self.verbose:
            for entry, seq in entry_and_seq.items():
                seq, matches = self.find_matches_and_remove_gaps(seq)
                print(f"{entry}\t{round(len(matches)/len(seq), 4)}")
        else:
            all_seqs = []
            for entry, seq in entry_and_seq.items():
                all_seqs.append(seq)
            all_seqs = ''.join(all_seqs)
            all_seqs, matches = self.find_matches_and_remove_gaps(all_seqs)
            try:
                gc_content = round(len(matches)/len(all_seqs), 4)
            except ZeroDivisionError:
                print("Input file has an unacceptable format. Please check input file argument.")
                sys.exit()
            print(gc_content)
            

    def process_args(self, args):
        return dict(
            fasta=args.fasta, 
            verbose=args.verbose
        )

    def find_matches_and_remove_gaps(self, seq: str):
        regex_pattern = re.compile('[GgCc]')
        seq = seq.replace('-', '')
        matches = regex_pattern.findall(seq)
        return seq, matches


    def determine_file_type_fasta_file(self, fasta):
        # automatically determine file type and read in the fasta file
        for fileFormat in FileFormat:
            try:
                records = SeqIO.parse(fasta, fileFormat.value)
                return records, fileFormat.value
            # the following exceptions refer to skipping over errors
            # associated with reading the wrong input file
            except ValueError:
                continue
            except AssertionError:
                continue

        raise Exception("Input file could not be read")