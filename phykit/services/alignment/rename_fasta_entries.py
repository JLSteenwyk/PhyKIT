import sys
from typing import Dict

from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaIterator

from .base import Alignment
from ...helpers.json_output import print_json


class RenameFastaEntries(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            fasta=parsed["fasta"],
            idmap=parsed["idmap"],
            output_file_path=parsed["output_file_path"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        try:
            records = SeqIO.parse(self.fasta, "fasta")
        except FileNotFoundError:
            print("FASTA file path corresponds to no such file. Please check the path.")
            sys.exit(2)

        idmap = self.load_idmap(self.idmap)

        renamed_count, total_records = self.replace_ids_and_write(
            self.output_file_path, records, idmap
        )

        if self.json_output:
            print_json(
                dict(
                    input_fasta=self.fasta,
                    idmap=self.idmap,
                    output_file=self.output_file_path,
                    total_records=total_records,
                    renamed_records=renamed_count,
                )
            )

    def process_args(self, args) -> Dict[str, str]:
        output_file_path = f"{args.output or args.fasta}.renamed.fa"
        return dict(
            fasta=args.fasta,
            idmap=args.idmap,
            output_file_path=output_file_path,
            json_output=getattr(args, "json", False),
        )

    def load_idmap(self, idmap_file: str) -> Dict[str, str]:
        try:
            with open(idmap_file) as f:
                return dict(line.split() for line in f)
        except FileNotFoundError:
            print("Idmap path corresponds to no such file. Please check the path.")
            sys.exit(2)

    def replace_ids_and_write(
        self,
        output_file_path: str,
        records: FastaIterator,
        idmap: Dict[str, str]
    ) -> tuple[int, int]:
        renamed_count = 0
        total_records = 0
        with open(output_file_path, "w") as output_file:
            for record in records:
                total_records += 1
                if record.id in idmap:
                    record.id = idmap[record.id]
                    record.description = ""
                    renamed_count += 1
                SeqIO.write(record, output_file, "fasta")
        return renamed_count, total_records
