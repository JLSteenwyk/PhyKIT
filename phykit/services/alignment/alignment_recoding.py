from __future__ import annotations

from os import path
import sys


from .base import Alignment


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)

here = path.dirname(__file__)


class AlignmentRecoding(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"], code=parsed["code"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()

        recoding_table = self.read_recoding_table(self.code)

        recoded_alignment = self._recode_alignment_strings(
            alignment, recoding_table, is_protein
        )

        if self.json_output:
            print_json(
                dict(
                    code=self.code,
                    taxa=[
                        dict(taxon=taxon, sequence=seq)
                        for taxon, seq in recoded_alignment.items()
                    ],
                )
            )
            return

        if recoded_alignment:
            print("\n".join(
                f">{taxon}\n{sequence}"
                for taxon, sequence in recoded_alignment.items()
            ))

    def recode_alignment(
        self,
        alignment: MultipleSeqAlignment,
        recoding_table: dict[str, str],
        is_protein: bool,
    ) -> dict[str, list[str]]:

        recoded_alignment = self._recode_alignment_strings(
            alignment,
            recoding_table,
            is_protein,
        )

        return {
            taxon: list(sequence)
            for taxon, sequence in recoded_alignment.items()
        }

    def _recode_alignment_strings(
        self,
        alignment: MultipleSeqAlignment,
        recoding_table: dict[str, str],
        is_protein: bool,
    ) -> dict[str, str]:
        translate_table = self._build_translation_table(
            recoding_table,
            self.get_gap_chars(is_protein),
        )
        recoded_alignment = dict()

        for record in alignment:
            recoded_alignment[record.id] = str(record.seq).translate(translate_table)

        return recoded_alignment

    @staticmethod
    def _build_translation_table(
        recoding_table: dict[str, str],
        gap_chars: list[str],
    ) -> dict[int, str]:
        gap_chars_set = set(gap_chars)
        translate_table: dict[int, str] = {}
        for original, recoded in recoding_table.items():
            upper_original = original.upper()
            if upper_original not in gap_chars_set:
                translate_table[ord(upper_original)] = recoded
            lower_original = upper_original.lower()
            if lower_original not in gap_chars_set:
                translate_table[ord(lower_original)] = recoded
        return translate_table

    def read_recoding_table(
        self,
        recoding: str
    ) -> dict[str, str]:
        """
        return translation table with codons as keys and amino acids as values
        """

        recoding_table = dict()

        if recoding is None:
            print("Please specify a recoding table")
            sys.exit(2)

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
                    parts = line.split(None, 2)
                    recoding_table[parts[1].upper()] = parts[0].upper()
        except FileNotFoundError:
            print(f"Recoding table file '{pathing}' not found.")
            sys.exit(2)

        return recoding_table

    def process_args(self, args):
        return dict(
            alignment_file_path=args.alignment,
            code=args.code,
            json_output=getattr(args, "json", False),
        )
