from __future__ import annotations

import os
import sys

from .base import Alignment

_path_exists = os.path.exists
_FASTA_WRAP_WIDTH = 60
_WRAPPED_FASTA_BATCH_CHUNKS = 8_192
_WRAPPED_FASTA_BATCH_MIN_LENGTH = 1_000_000


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
            width = _FASTA_WRAP_WIDTH
            two_line_limit = width * 2
            batch_min_length = _WRAPPED_FASTA_BATCH_MIN_LENGTH
            for title, sequence in SimpleFastaParser(input_file):
                total_records += 1
                record_id = title.split(None, 1)[0]
                mapped_id = idmap.get(record_id, missing)
                if mapped_id is not missing:
                    header = mapped_id
                    renamed_count += 1
                else:
                    header = title
                sequence_length = len(sequence)
                if sequence_length == 0:
                    write(f">{header}\n")
                elif sequence_length <= width:
                    write(f">{header}\n{sequence}\n")
                elif sequence_length <= two_line_limit:
                    write(f">{header}\n{sequence[:width]}\n{sequence[width:]}\n")
                elif sequence_length >= batch_min_length:
                    write(f">{header}\n")
                    self._write_wrapped_fasta_sequence(output_file, sequence, width)
                else:
                    wrapped_sequence = "\n".join(
                        [
                            sequence[idx:idx + width]
                            for idx in range(0, sequence_length, width)
                        ]
                    )
                    write(f">{header}\n{wrapped_sequence}\n")

        return renamed_count, total_records

    @staticmethod
    def _write_wrapped_fasta_sequence(handle, sequence: str, width: int = 60) -> None:
        sequence_length = len(sequence)
        if sequence_length == 0:
            return
        if sequence_length < _WRAPPED_FASTA_BATCH_MIN_LENGTH:
            handle.write(
                "\n".join(
                    [
                        sequence[idx:idx + width]
                        for idx in range(0, sequence_length, width)
                    ]
                )
            )
            handle.write("\n")
            return

        write = handle.write
        batch_size = width * _WRAPPED_FASTA_BATCH_CHUNKS
        for start in range(0, sequence_length, batch_size):
            stop = min(start + batch_size, sequence_length)
            write(
                "\n".join(
                    [
                        sequence[idx:idx + width]
                        for idx in range(start, stop, width)
                    ]
                )
            )
            write("\n")
