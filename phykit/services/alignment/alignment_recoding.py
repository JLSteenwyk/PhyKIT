from os import path
import sys
from typing import Dict, List

from Bio.Align import MultipleSeqAlignment

from .base import Alignment

here = path.dirname(__file__)


class AlignmentRecoding(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()

        recoding_table = self.read_recoding_table(self.code[0])

        recoded_alignment = self.recode_alignment(
            alignment, recoding_table, is_protein
        )

        for k, v in recoded_alignment.items():
            print(f">{k}\n{''.join(v)}")

    def recode_alignment(
        self,
        alignment: MultipleSeqAlignment,
        recoding_table: Dict[str, str],
        is_protein: bool,
    ) -> Dict[str, List[str]]:

        gap_chars = self.get_gap_chars()
        recoded_alignment = dict()

        for record in alignment:
            recoded_sequence = [
                recoding_table.get(base.upper(), base)
                if base not in gap_chars else base
                for base in record.seq
            ]
            recoded_alignment[record.id] = recoded_sequence

        return recoded_alignment

    def read_recoding_table(
        self,
        recoding: str
    ) -> Dict[str, str]:
        """
        return translation table with codons as keys and amino acids as values
        """

        recoding_table = dict()

        if recoding is None:
            print("Please specify a recoding table")
            sys.exit()

        recoding_paths = {
            "RY-nucleotide": "../../recoding_tables/RY-nucleotide.txt",
            "SandR-6": "../../recoding_tables/S_and_R-6.txt",
            "KGB-6": "../../recoding_tables/KGB-6.txt",
            "Dayhoff-6": "../../recoding_tables/Dayhoff-6.txt",
            "Dayhoff-9": "../../recoding_tables/Dayhoff-9.txt",
            "Dayhoff-12": "../../recoding_tables/Dayhoff-12.txt",
            "Dayhoff-15": "../../recoding_tables/Dayhoff-15.txt",
            "Dayhoff-18": "../../recoding_tables/Dayhoff-18.txt",
        }
        pathing = recoding_paths.get(recoding, str(recoding))

        try:
            with open(path.join(here, pathing)) as code:
                for line in code:
                    parts = line.split()
                    recoding_table[parts[1].upper()] = parts[0].upper()
        except FileNotFoundError:
            print(f"Recoding table file '{pathing}' not found.")
            sys.exit()

        return recoding_table

    def process_args(self, args):
        return dict(
            alignment_file_path=args.alignment,
            code=args.code
        )
