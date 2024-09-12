# import sys
# from textwrap import dedent
# from typing import Dict, List, Tuple, Union

# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from collections import defaultdict

# from .base import Alignment

# from ...helpers.files import read_single_column_file_to_list


# class CreateConcatenationMatrix(Alignment):
#     def __init__(self, args) -> None:
#         super().__init__(**self.process_args(args))

#     def run(self) -> None:
#         self.create_concatenation_matrix(
#             self.alignment_list_path,
#             self.prefix
#         )

#     def process_args(self, args) -> Dict[str, str]:
#         return dict(alignment_list_path=args.alignment_list, prefix=args.prefix)

#     def read_alignment_paths(self, alignment_list_path: str) -> list:
#         try:
#             return read_single_column_file_to_list(alignment_list_path)
#         except FileNotFoundError:
#             print("Alignment list file (-a) is not found. Please check pathing.")
#             sys.exit()

#     def get_taxa_names(self, alignment_paths: List[str]) -> List[str]:
#         taxa = set()
#         for alignment_path in alignment_paths:
#             alignment = SeqIO.parse(alignment_path, "fasta")
#             for seq_record in alignment:
#                 taxa.add(seq_record.name)

#         return sorted(taxa)

#     def print_start_message(
#         self,
#         taxa: List[str],
#         alignment_paths: List[str],
#         file_partition: str,
#         fasta_output: str,
#         file_occupancy: str,
#     ) -> None:
#         start_message = dedent(
#             f"""
#             --------------------
#             | General features |
#             --------------------
#             Total number of taxa: {len(taxa)}
#             Total number of alignments: {len(alignment_paths)}


#             ----------------
#             | Output files |
#             ----------------
#             Partition file output: {file_partition}
#             Concatenated fasta output: {fasta_output}
#             Occupancy report: {file_occupancy}
#         """
#         )

#         print(start_message)

#     def get_list_of_taxa_and_records(
#         self,
#         alignment_path: str
#     ) -> Tuple[List[str], List[SeqRecord]]:
#         og_taxa = []
#         records = list(SeqIO.parse(alignment_path, "fasta"))

#         for record in records:
#             og_taxa.append(record.id)

#         return og_taxa, records

#     def create_missing_seq_str(
#         self,
#         records: List[SeqRecord],
#         alignment_path: str
#     ) -> Tuple[str, int]:
#         """
#         create a placeholder string for sequences with missing taxa
#         """
#         og_len = ""
#         missing_seq = ""
#         try:
#             og_len = len(records[0].seq)
#         except IndexError:
#             print(f"There are no sequence records in {alignment_path}.\nExiting now...")
#             sys.exit()
#         missing_seq = og_len * "?"

#         return missing_seq, og_len

#     def create_missing_taxa_sequence(
#         self,
#         missing_taxa: List[str],
#         missing_seq: str,
#         concat: List[Union[Seq, str]],
#     ) -> List[Union[Seq, str]]:
#         for indiv in missing_taxa:
#             concat[indiv].append(missing_seq)
#         return concat

#     def create_record_for_present_taxa(
#         self,
#         records: List[SeqRecord],
#         og_taxa: List[str],
#         concat: Dict[str, List[Union[Seq, str]]],
#     ) -> Dict[str, List[Seq]]:
#         for record in records:
#             if record.id in og_taxa:
#                 concat[record.id].append(record.seq)

#         return concat

#     def add_to_partition_info(
#         self,
#         partition_info: List[str],
#         og_len: int,
#         field_one: str,
#         fasta: str,
#         first_len: int,
#         second_len: int,
#     ) -> Tuple[List[str], int, int]:
#         second_len += og_len
#         entry = f"{field_one}, {str(fasta)}={str(first_len)}-{str(second_len)}\n"
#         partition_info.append(entry)

#         first_len += og_len

#         return partition_info, first_len, second_len

#     def add_to_occupancy_info(
#         self,
#         occupancy_info: List[str],
#         og_taxa: List[str],
#         missing_taxa: List[str],
#         fasta: str,
#     ) -> List[str]:
#         num_present = len(og_taxa)
#         num_missing = len(missing_taxa)
#         percent_occupancy = num_present / (num_present + num_missing)

#         if num_missing == 0:
#             missing_taxa = ["None"]
#         else:
#             missing_taxa.sort()
#         og_taxa.sort()
#         entry = f"{str(fasta)}\t{str(num_present)}\t{str(num_missing)}\t{str(percent_occupancy)}\t{';'.join(missing_taxa)}\t{';'.join(og_taxa)}\n"
#         occupancy_info.append(entry)

#         return occupancy_info

#     def fasta_file_write(
#         self,
#         fasta_output: str,
#         concat: List[Union[Seq, str]]
#     ) -> None:
#         with open(fasta_output, "w") as final_fasta_file:
#             for x in concat:
#                 concatenated = []
#                 for s in concat[x]:
#                     try:
#                         # if a seq object
#                         concatenated.append(str(s._data.decode("utf-8")))
#                     except AttributeError:
#                         # if a string
#                         concatenated.append(str(s))
#                 concat[x] = "".join(concatenated)
#                 entry = f">{x}\n{''.join(str(concat[x]))}\n"
#                 final_fasta_file.write(str(entry))

#     def write_occupancy_or_partition_file(
#         self,
#         occu_info: List[str],
#         output_file_name: str,
#     ) -> None:
#         with open(output_file_name, "w") as f:
#             for entry in occu_info:
#                 f.write(entry)

#     def create_concatenation_matrix(
#         self,
#         alignment_list_path: str,
#         prefix: str,
#     ) -> None:
#         alignment_paths = self.read_alignment_paths(alignment_list_path)

#         taxa = self.get_taxa_names(alignment_paths)

#         # assign the output file names
#         file_partition = f"{prefix}.partition"
#         fasta_output = f"{prefix}.fa"
#         file_occupancy = f"{prefix}.occupancy"

#         self.print_start_message(
#             taxa, alignment_paths, file_partition, fasta_output, file_occupancy
#         )

#         # assign placeholders for lengths in the partition file
#         first_len = 1
#         second_len = 0

#         # create a dictionary with keys as taxa and empty list values
#         concatenated_seqs = {key: [] for key in taxa}

#         # string for first field of partition file
#         field_one = "AUTO"

#         # create empty lists to contain information for the
#         # partition and occupancy files
#         partition_info = []
#         occupancy_info = []

#         # loop through alignment files
#         for alignment_path in alignment_paths:
#             # create a list of taxa and records
#             og_taxa, records = self.get_list_of_taxa_and_records(alignment_path)

#             # initialize list and determine missing taxa
#             missing_taxa = list(set(taxa) - set(og_taxa))

#             # create string for missing seq
#             missing_seq, og_len = self.create_missing_seq_str(records, alignment_path)

#             # create missing taxa sequence
#             concatenated_seqs = self.create_missing_taxa_sequence(
#                 missing_taxa, missing_seq, concatenated_seqs
#             )

#             # create record for present data
#             concatenated_seqs = self.create_record_for_present_taxa(
#                 records, og_taxa, concatenated_seqs
#             )

#             # append partition to to partition info
#             partition_info, first_len, second_len = self.add_to_partition_info(
#                 partition_info, og_len, field_one, alignment_path, first_len, second_len
#             )

#             # append occupancy data to occupancy info
#             occupancy_info = self.add_to_occupancy_info(
#                 occupancy_info, og_taxa, missing_taxa, alignment_path
#             )

#         # write output fasta file
#         self.fasta_file_write(fasta_output, concatenated_seqs)

#         # write out occupancy and partition file
#         self.write_occupancy_or_partition_file(occupancy_info, file_occupancy)
#         self.write_occupancy_or_partition_file(partition_info, file_partition)

#         print("Complete!\n")

import sys
from textwrap import dedent
from typing import Dict, List, Tuple, Union
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

from .base import Alignment
from ...helpers.files import read_single_column_file_to_list


class CreateConcatenationMatrix(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        self.create_concatenation_matrix(
            self.alignment_list_path,
            self.prefix
        )

    def process_args(self, args) -> Dict[str, str]:
        return dict(alignment_list_path=args.alignment_list, prefix=args.prefix)

    def read_alignment_paths(self, alignment_list_path: str) -> List[str]:
        try:
            return read_single_column_file_to_list(alignment_list_path)
        except FileNotFoundError:
            print("Alignment list file (-a) is not found. Please check pathing.")
            sys.exit()

    def get_taxa_names(self, alignment_paths: List[str]) -> List[str]:
        taxa = set()
        for alignment_path in alignment_paths:
            for seq_record in SeqIO.parse(alignment_path, "fasta"):
                taxa.add(seq_record.id)
        return sorted(taxa)

    def print_start_message(
        self,
        taxa: List[str],
        alignment_paths: List[str],
        file_partition: str,
        fasta_output: str,
        file_occupancy: str,
    ) -> None:
        start_message = dedent(f"""
            --------------------
            | General features |
            --------------------
            Total number of taxa: {len(taxa)}
            Total number of alignments: {len(alignment_paths)}


            ----------------
            | Output files |
            ----------------
            Partition file output: {file_partition}
            Concatenated fasta output: {fasta_output}
            Occupancy report: {file_occupancy}
        """)
        print(start_message)

    def get_list_of_taxa_and_records(
        self, alignment_path: str
    ) -> Tuple[set, List[SeqRecord]]:
        records = list(SeqIO.parse(alignment_path, "fasta"))
        og_taxa = {record.id for record in records}
        return og_taxa, records

    def create_missing_seq_str(self, records: List[SeqRecord]) -> Tuple[str, int]:
        """Create a placeholder string for sequences with missing taxa."""
        if not records:
            print(f"No sequence records found. Exiting...")
            sys.exit()

        og_len = len(records[0].seq)
        missing_seq = '?' * og_len
        return missing_seq, og_len

    def process_taxa_sequences(
        self,
        records: List[SeqRecord],
        taxa: List[str],
        concatenated_seqs: Dict[str, List[str]],
        missing_seq: str,
    ) -> None:
        present_taxa = {record.id for record in records}
        missing_taxa = set(taxa) - present_taxa

        # Add sequences for present taxa
        for record in records:
            concatenated_seqs[record.id].append(str(record.seq))

        # Add missing sequences for missing taxa
        for taxon in missing_taxa:
            concatenated_seqs[taxon].append(missing_seq)

    def add_to_partition_info(
        self,
        partition_info: List[str],
        og_len: int,
        field_one: str,
        fasta: str,
        first_len: int,
        second_len: int,
    ) -> Tuple[List[str], int, int]:
        second_len += og_len
        partition_info.append(f"{field_one}, {fasta}={first_len}-{second_len}\n")
        return partition_info, second_len + 1, second_len

    def add_to_occupancy_info(
        self,
        occupancy_info: List[str],
        present_taxa: set,
        taxa: List[str],
        fasta: str,
    ) -> List[str]:
        missing_taxa = sorted(set(taxa) - present_taxa)
        num_present = len(present_taxa)
        num_missing = len(missing_taxa)
        percent_occupancy = num_present / len(taxa)
        occupancy_info.append(f"{fasta}\t{num_present}\t{num_missing}\t{percent_occupancy:.4f}\t{';'.join(missing_taxa)}\n")
        return occupancy_info

    def fasta_file_write(self, fasta_output: str, concatenated_seqs: Dict[str, List[str]]) -> None:
        with open(fasta_output, "w") as final_fasta_file:
            for taxon, sequences in concatenated_seqs.items():
                final_fasta_file.write(f">{taxon}\n{''.join(sequences)}\n")

    def write_occupancy_or_partition_file(self, info: List[str], output_file_name: str) -> None:
        with open(output_file_name, "w") as f:
            f.writelines(info)

    def create_concatenation_matrix(self, alignment_list_path: str, prefix: str) -> None:
        alignment_paths = self.read_alignment_paths(alignment_list_path)
        taxa = self.get_taxa_names(alignment_paths)

        # Assign output file names
        file_partition = f"{prefix}.partition"
        fasta_output = f"{prefix}.fa"
        file_occupancy = f"{prefix}.occupancy"

        self.print_start_message(taxa, alignment_paths, file_partition, fasta_output, file_occupancy)

        # Initialize placeholders for partition info
        first_len, second_len = 1, 0
        partition_info, occupancy_info = [], []
        concatenated_seqs = {taxon: [] for taxon in taxa}

        # Process each alignment file
        for alignment_path in alignment_paths:
            present_taxa, records = self.get_list_of_taxa_and_records(alignment_path)
            missing_seq, og_len = self.create_missing_seq_str(records)

            # Process taxa sequences and add to the concatenated sequences
            self.process_taxa_sequences(records, taxa, concatenated_seqs, missing_seq)

            # Add to partition and occupancy info
            partition_info, first_len, second_len = self.add_to_partition_info(
                partition_info, og_len, "AUTO", alignment_path, first_len, second_len
            )
            occupancy_info = self.add_to_occupancy_info(occupancy_info, present_taxa, taxa, alignment_path)

        # Write output files
        self.fasta_file_write(fasta_output, concatenated_seqs)
        self.write_occupancy_or_partition_file(occupancy_info, file_occupancy)
        self.write_occupancy_or_partition_file(partition_info, file_partition)

        print("Complete!\n")
