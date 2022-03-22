import statistics as stat
import sys
from textwrap import dedent
from typing import Tuple

from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq

from .base import Alignment

from ...helpers.files import read_single_column_file_to_list

class CreateConcatenationMatrix(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        self.create_concatenation_matrix(self.alignment_list_path, self.prefix)

    def process_args(self, args) -> dict:
        return dict(
            alignment_list_path=args.alignment_list,
            prefix=args.prefix
        )

    def read_alignment_paths(
        self,
        alignment_list_path: str
    ) -> list:
        """
        Read alignment paths into a list 
        """
        try:
            return read_single_column_file_to_list(alignment_list_path)
        except FileNotFoundError:
            print("Alignment list file (-a) is not found. Please check pathing.")
            sys.exit()

    def get_taxa_names(
        self,
        alignment_paths: list
    ) -> list:
        # get taxa names
        taxa = set()
        for alignment_path in alignment_paths:
            alignment = SeqIO.parse(alignment_path, "fasta")
            for seq_record in alignment:
                taxa.add(seq_record.name)

        return list(sorted(taxa))

    def print_start_message(
        self,
        taxa: list,
        alignment_paths: list,
        file_partition: str,
        fasta_output: str,
        file_occupancy: str
    ) -> None:
        """
        print a log message to the user about general features of
        input and output file names
        """

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
        self,
        alignment_path: str
    ) -> Tuple[list, list]:
        """
        get list of a taxa and their records
        """
        # list for keeping track of taxa per USCO
        og_taxa = []
        # open alignment file
        records = list(SeqIO.parse(alignment_path, "fasta"))

        # create lists of USCOtaxa and missing_taxa
        for record in records:
            # append taxa per USCO to USCOtaxa
            og_taxa.append(record.id)

        return og_taxa, records

    def determine_missing_sequences(
        self,
        taxa: list,
        og_taxa: list
    ) -> list:
        """
        determine taxa not present in an alignment file
        """
        missing_taxa = list(set(taxa)-set(og_taxa))

        return missing_taxa

    def create_missing_seq_str(
        self,
        records: list,
        alignment_path: str
    ) -> Tuple[str, int]:
        """
        create a placeholder string for sequences with missing taxa
        """
        og_len = ''
        missing_seq = ''
        try:
            og_len = len(records[0].seq)
        except IndexError:
            print(f"There are no sequence records in {alignment_path}.\nExiting now...")
            sys.exit()
        missing_seq = og_len * '?'
    
        return missing_seq, og_len

    def create_missing_taxa_sequence(
        self,
        missing_taxa: list,
        missing_seq: str,
        concat: dict
    ) -> dict:

        for indiv in missing_taxa:
            concat[indiv].append(missing_seq)
            
        return concat 

    def create_record_for_present_taxa(
        self,
        records: list,
        og_taxa: list,
        concat: dict
    ) -> dict:
        """
        create a record for taxa that are present in the concatenation dict
        """
        for record in records:
            if record.id in og_taxa:
                concat[record.id].append(record.seq)

        return concat

    def add_to_partition_info(
        self,
        partition_info: list,
        og_len: int,
        field_one: str,
        fasta: str,
        first_len: int,
        second_len: int
    ) -> Tuple[list, int, int]:
        """
        append to partition file
        """
        # second value in partition file
        second_len += og_len
        # add partition file entry to partition_info
        entry = f"{field_one}, {str(fasta)}={str(first_len)}-{str(second_len)}\n"
        partition_info.append(entry)
        # add to first value for partition file
        first_len += og_len

        return partition_info, first_len, second_len

    def add_to_occupancy_info(
        self,
        occupancy_info: list,
        og_taxa: list,
        missing_taxa: list,
        fasta: str
    ) -> list:
        """
        write partition file to detail boundaries in the concat alignment
        """

        # determine number of present taxa
        num_present = len(og_taxa)
        # determine number of missing taxa
        num_missing = len(missing_taxa)
        # determine percent occupied
        percent_occupancy = num_present/(num_present+num_missing)
        # handle missing taxa
        if num_missing == 0:
            missing_taxa = ["None"]
        else:
            missing_taxa.sort()
        og_taxa.sort()
        entry = f"{str(fasta)}\t{str(num_present)}\t{str(num_missing)}\t{str(percent_occupancy)}\t{';'.join(missing_taxa)}\t{';'.join(og_taxa)}\n"
        occupancy_info.append(entry)
        
        return occupancy_info

    def fasta_file_write(
        self,
        fasta_output: str,
        concat: dict
    ) -> None:
        """
        join seqs of genes in value (list) and write to concat_fa 
        """
        with open(fasta_output, "w") as final_fasta_file: 
            for x in concat:
                concatenated = []
                for s in concat[x]:
                    try:
                        # if a seq object
                        concatenated.append(str(s._data.decode("utf-8")))
                    except AttributeError:
                        # if a string
                        concatenated.append(str(s))
                concat[x] = ''.join(concatenated)
                entry = f">{x}\n{''.join(str(concat[x]))}\n"
                final_fasta_file.write(str(entry))   

    def write_occupancy_or_partition_file(
        self,
        list_of_lists: list,
        output_file_name: str
    ) -> None:
        """
        loops through occupancy or partition file information and writes
        them to the corresponding output file
        """

        with open(output_file_name, "w") as f:
            for entry in list_of_lists:
                f.write(entry)

    def create_concatenation_matrix(
        self,
        alignment_list_path: str,
        prefix: str
    ) -> None:
        alignment_paths = self.read_alignment_paths(alignment_list_path)

        taxa = self.get_taxa_names(alignment_paths)

        # assign the output file names
        file_partition = f"{prefix}.partition"
        fasta_output = f"{prefix}.fa"
        file_occupancy = f"{prefix}.occupancy"

        # print log message
        self.print_start_message(
            taxa,
            alignment_paths,
            file_partition,
            fasta_output,
            file_occupancy
        )

        # assign placeholders for lengths in the partition file
        first_len    = 1
        second_len   = 0

        # create a dictionary with keys as taxa and empty list values
        concatenated_seqs = {key: [] for key in taxa}

        # string for first field of partition file
        field_one = 'AUTO'

        # create empty lists to contain information for the 
        # partition and occupancy files
        partition_info = []
        occupancy_info = []

        # loop through alignment files
        for alignment_path in alignment_paths:

            # create a list of taxa and records
            og_taxa, records = self.get_list_of_taxa_and_records(
                alignment_path
            )
                
            # initialize list and determine missing taxa
            missing_taxa = self.determine_missing_sequences(
                taxa,
                og_taxa
            )

            # create string for missing seq
            missing_seq, og_len = self.create_missing_seq_str(
                records,
                alignment_path
            )

            # create missing taxa sequence
            concatenated_seqs = self.create_missing_taxa_sequence(
                missing_taxa,
                missing_seq,
                concatenated_seqs
            )      

            # create record for present data
            concatenated_seqs = self.create_record_for_present_taxa(
                records,
                og_taxa,
                concatenated_seqs
            )

            # append partition to to partition info
            partition_info, first_len, second_len = self.add_to_partition_info(
                partition_info,
                og_len,
                field_one,
                alignment_path,
                first_len,
                second_len
            )

            # append occupancy data to occupancy info
            occupancy_info = self.add_to_occupancy_info(
                occupancy_info,
                og_taxa,
                missing_taxa,
                alignment_path
            )

        # write output fasta file
        self.fasta_file_write(
            fasta_output,
            concatenated_seqs
        )

        # write out occupancy and partition file
        self.write_occupancy_or_partition_file(occupancy_info, file_occupancy)
        self.write_occupancy_or_partition_file(partition_info, file_partition)

        print("Complete!\n")
