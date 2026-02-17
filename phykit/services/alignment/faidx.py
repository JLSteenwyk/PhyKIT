from typing import Dict

from Bio import SeqIO

from .base import Alignment
from ...helpers.json_output import print_json


class Faidx(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(fasta=parsed["fasta"], entry=parsed["entry"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        record_dict = SeqIO.index(self.fasta, "fasta")
        entries = [e for e in map(str.strip, self.entry.split(",")) if e]

        if self.json_output:
            rows = [
                dict(entry=e, name=record_dict[e].name, sequence=str(record_dict[e].seq))
                for e in entries
            ]
            print_json(
                dict(
                    rows=rows,
                    entries=rows,
                )
            )
            return

        # Split entries and iterate
        for e in entries:
            record = record_dict[e]
            print(f">{record.name}\n{record.seq}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            fasta=args.fasta,
            entry=args.entry,
            json_output=getattr(args, "json", False),
        )
