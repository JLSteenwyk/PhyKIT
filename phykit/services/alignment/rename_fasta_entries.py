# import sys

# from Bio import SeqIO

# from .base import Alignment


# class RenameFastaEntries(Alignment):
#     def __init__(self, args) -> None:
#         super().__init__(**self.process_args(args))

#     def run(self):
#         try:
#             records = SeqIO.parse(self.fasta, "fasta")
#         except FileNotFoundError:
#             print(f"{self.fasta} corresponds to no such file or directory.")
#             print("Please double check pathing and filenames")
#             sys.exit()

#         # save idmap to a dictionary
#         idmap = self.idmap_to_dictionary(self.idmap)

#         # replace and write out
#         self.replace_ids_and_write(self.output_file_path, records, idmap)

#     def process_args(self, args):
#         if args.output is None:
#             output_file_path = f"{args.fasta}.renamed.fa"
#         else:
#             output_file_path = f"{args.output}"

#         return dict(
#             fasta=args.fasta, idmap=args.idmap, output_file_path=output_file_path
#         )

#     def replace_ids_and_write(self, output_file_path, records, idmap):
#         """
#         for tips with a name as a key in the idmap,
#         replace that tip name with the value in the
#         idmap and write to output file
#         """
#         with open(output_file_path, "w") as output_file_path:
#             for record in records:
#                 if record.id in idmap:
#                     # replace ID
#                     record.id = idmap[record.id]
#                     # remove description
#                     record.description = ""
#                 SeqIO.write(record, output_file_path, "fasta")

#     def idmap_to_dictionary(self, idmap: str) -> dict:
#         """
#         read idmap into a dictionary
#         """
#         idmap = {}
#         try:
#             with open(self.idmap) as identifiers:
#                 for line in identifiers:
#                     (key, val) = line.split()
#                     idmap[key] = val
#             return idmap
#         except FileNotFoundError:
#             print(f"{self.idmap} corresponds to no such file or directory.")
#             print("Please double check pathing and filenames")
#             sys.exit()

import sys
from Bio import SeqIO
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

    def process_args(self, args):
        output_file_path = f"{args.output or args.fasta}.renamed.fa"
        return dict(fasta=args.fasta, idmap=args.idmap, output_file_path=output_file_path)

    def load_idmap(self, idmap_file: str) -> dict:
        try:
            with open(idmap_file) as f:
                return dict(line.split() for line in f)
        except FileNotFoundError:
            print("Idmap path corresponds to no such file. Please check the path.")
            sys.exit()

    def replace_ids_and_write(self, output_file_path: str, records, idmap: dict):
        with open(output_file_path, "w") as output_file:
            for record in records:
                if record.id in idmap:
                    record.id = idmap[record.id]  # Replace ID
                    record.description = ""       # Remove description
                SeqIO.write(record, output_file, "fasta")
