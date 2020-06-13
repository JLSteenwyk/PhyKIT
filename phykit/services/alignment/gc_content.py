from enum import Enum
import re

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
            print(f"{len(matches)/len(all_seqs)}")
            

    def process_args(self, args):
        return dict(
            fasta=args.fasta, 
            verbose=args.verbose
        )

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