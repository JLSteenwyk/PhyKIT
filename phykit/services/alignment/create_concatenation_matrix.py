import statistics as stat
from textwrap import dedent
from typing import Tuple

from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq

from .base import Alignment


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

        return list(taxa)
    
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
                    concatenated.append(s._data)
                concat[x] = concatenated
                entry = f">{x}\n{''.join(concat[x])}\n"
                final_fasta_file.write(str(entry))   

    def occupancy_file_write(
        self,
        file_occupancy: str,
        og_taxa: list,
        missing_taxa: list,
        fasta: str
    ) -> None:
        """
        write partition file to detail boundaries in the concat alignment
        """
        with open(file_occupancy, "a") as f:
            # determine number of present taxa
            num_present = len(og_taxa)
            # determine number of missing taxa
            num_missing = len(missing_taxa)
            # determine percent occupied
            percent_occupancy = num_present/(num_present+num_missing)
            if num_missing == 0:
                missing_taxa = ["None"]
            entry = f"{str(fasta)}\t{str(num_present)}\t{str(num_missing)}\t{str(percent_occupancy)}\t{';'.join(missing_taxa)}\t{';'.join(og_taxa)}\n"
            f.write(str(entry))

    def partition_file_write(
        self,
        file_partition: str,
        og_len: int,
        field_one: str,
        fasta: str,
        first_len: int,
        second_len: int
    ) -> Tuple[int, int]:
        """
        append to partition file
        """
        with open(file_partition, "a") as f:
            # second value in partition file
            second_len += og_len
            entry = f"{field_one}, {str(fasta)}={str(first_len)}-{str(second_len)}\n"
            f.write(str(entry))
            # add to first value for partition file
            first_len += og_len

        return first_len, second_len

    def create_missing_taxa_sequence(
        self,
        missing_taxa: list,
        missing_seq: str,
        concat: dict
    ) -> dict:
        for indiv in missing_taxa:
            indiv_record = SeqRecord(Seq(missing_seq), 
                id = indiv, name = indiv, description = indiv)
            concat[indiv].append(indiv_record.seq)
            
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

    def create_missing_seq_str(
        self,
        records: list
    ) -> Tuple[str, int]:
        """
        create a placeholder string for sequences with missing taxa
        """
        og_len = ''
        missing_seq = ''
        og_len = len(records[0].seq)
        missing_seq = og_len * '?'
    
        return missing_seq, og_len

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

    def read_alignment_paths(self, alignment_list_path: str) -> list:
        """
        Read alignment paths into a list 
        """
        # TODO: handle I/O exception here (try & except)
        with open(alignment_list_path) as f:
            return [line.rstrip('\n') for line in f]

    def create_output_files(self, file_partition, file_occupancy):
        # create an empty file_partition
        open(file_partition, "w")
        # create an empty occupancy per OG file
        open(file_occupancy, "w")

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

        self.create_output_files(file_partition, file_occupancy)

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
                records
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

            # append to partition file
            # TODO: refactor partition file and occupancy file writing
            # write data structures to list of lists
            first_len, second_len = self.partition_file_write(
                file_partition,
                og_len,
                field_one,
                alignment_path,
                first_len,
                second_len
            )

            # append to occupancy file
            self.occupancy_file_write(
                file_occupancy,
                og_taxa,
                missing_taxa,
                alignment_path
            )

        # write output fasta file
        self.fasta_file_write(
            fasta_output,
            concatenated_seqs
        )

        print("Complete!\n")
