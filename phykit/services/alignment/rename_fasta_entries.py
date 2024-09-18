from multiprocessing import Pool
import sys
from typing import Dict, List

from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.SeqRecord import SeqRecord

from .base import Alignment


class RenameFastaEntries(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        try:
            records = SeqIO.parse(self.fasta, "fasta")
        except FileNotFoundError:
            print("FASTA file path corresponds to no such file. Please check the path.")
            sys.exit()

        idmap = self.load_idmap(self.idmap)

        self.replace_ids_and_write(self.output_file_path, records, idmap)

    def process_args(self, args) -> Dict[str, str]:
        output_file_path = f"{args.output or args.fasta}.renamed.fa"
        return dict(
            fasta=args.fasta,
            idmap=args.idmap,
            output_file_path=output_file_path,
        )

    def load_idmap(self, idmap_file: str) -> Dict[str, str]:
        try:
            with open(idmap_file) as f:
                return dict(line.split() for line in f)
        except FileNotFoundError:
            print("Idmap path corresponds to no such file. Please check the path.")
            sys.exit()

    def replace_ids_and_write(
        self,
        output_file_path: str,
        records: List[FastaIterator],
        idmap: Dict[str, str]
    ) -> None:
        cpu = self.set_cpu()
        with Pool(cpu) as pool:
            updated_records = pool.starmap(
                self.replace_record_id,
                [
                    (record, idmap)
                    for record in records
                ]
            )

        with open(output_file_path, "w") as output_file:
            SeqIO.write(updated_records, output_file, "fasta")

    def replace_record_id(self, record, idmap: Dict[str, str]) -> SeqRecord:
        if record.id in idmap:
            record.id = idmap[record.id]
            record.description = ""
        return record
