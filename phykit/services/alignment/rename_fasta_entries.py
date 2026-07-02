from __future__ import annotations

import os
import sys

from .base import Alignment

_path_exists = os.path.exists


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


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
        if not _path_exists(self.fasta):
            print("FASTA file path corresponds to no such file. Please check the path.")
            sys.exit(2)

        idmap = self.load_idmap(self.idmap)
        renamed_count, total_records = self.replace_ids_in_file_and_write(
            self.output_file_path,
            self.fasta,
            idmap,
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

    def process_args(self, args) -> dict[str, str]:
        output_file_path = f"{args.output or args.fasta}.renamed.fa"
        return dict(
            fasta=args.fasta,
            idmap=args.idmap,
            output_file_path=output_file_path,
            json_output=getattr(args, "json", False),
        )

    def load_idmap(self, idmap_file: str) -> dict[str, str]:
        try:
            with open(idmap_file) as f:
                idmap = {}
                split = str.split
                for line in f:
                    key, val = split(line)
                    idmap[key] = val
                return idmap
        except FileNotFoundError:
            print("Idmap path corresponds to no such file. Please check the path.")
            sys.exit(2)

    def replace_ids_and_write(
        self,
        output_file_path: str,
        records: FastaIterator,
        idmap: dict[str, str]
    ) -> tuple[int, int]:
        from Bio import SeqIO

        renamed_count = 0

        def renamed_records():
            nonlocal renamed_count
            for record in records:
                if record.id in idmap:
                    record.id = idmap[record.id]
                    record.description = ""
                    renamed_count += 1
                yield record

        with open(output_file_path, "w") as output_file:
            total_records = SeqIO.write(renamed_records(), output_file, "fasta")
        return renamed_count, total_records

    def replace_ids_in_file_and_write(
        self,
        output_file_path: str,
        fasta_path: str,
        idmap: dict[str, str],
    ) -> tuple[int, int]:
        from Bio.SeqIO.FastaIO import SimpleFastaParser

        renamed_count = 0
        total_records = 0
        missing = object()

        with open(fasta_path) as input_file, open(output_file_path, "w") as output_file:
            write = output_file.write
            width = 60
            for title, sequence in SimpleFastaParser(input_file):
                total_records += 1
                record_id = title.split(None, 1)[0]
                mapped_id = idmap.get(record_id, missing)
                if mapped_id is not missing:
                    header = mapped_id
                    renamed_count += 1
                else:
                    header = title
                if not sequence:
                    write(f">{header}\n")
                elif len(sequence) <= width:
                    write(f">{header}\n{sequence}\n")
                else:
                    wrapped_sequence = "\n".join(
                        [
                            sequence[idx:idx + width]
                            for idx in range(0, len(sequence), width)
                        ]
                    )
                    write(f">{header}\n{wrapped_sequence}\n")

        return renamed_count, total_records

    @staticmethod
    def _write_wrapped_fasta_sequence(handle, sequence: str, width: int = 60) -> None:
        if not sequence:
            return
        handle.write(
            "\n".join(
                [
                    sequence[idx:idx + width]
                    for idx in range(0, len(sequence), width)
                ]
            )
        )
        handle.write("\n")
