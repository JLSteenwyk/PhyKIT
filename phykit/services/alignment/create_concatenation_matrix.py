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
        alignment_list = self.alignment_list
        prefix = self.prefix
        
        self.create_concatenation_matrix(alignment_list, prefix)

    def process_args(self, args) -> dict:
        return dict(
            alignment_list=args.alignment_list,
            prefix=args.prefix
        )

    def get_taxa_names(
        self,
        alignments: list
    ) -> list:
        # get taxa names
        taxa = []
        for alignment in alignments:
            alignment = SeqIO.parse(alignment, "fasta")
            for indiv in alignment:
                if indiv.name not in taxa:
                    taxa.append(indiv.name)

        return taxa
    
    def start_message_print(
        self,
        taxa: list,
        alignments: list,
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
            Total number of alignments: {len(alignments)}


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
            num_present  = len(og_taxa)
            # determine number of missing taxa
            num_missing  = len(missing_taxa)
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
            entry = field_one+", "+str(fasta)+"="+str(first_len)+"-"+str(second_len)+"\n"
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

        og_len     = ''
        missing_seq = ''
        og_len   = len(records[0].seq)
        missing_seq = og_len*'?'
    
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
        fasta: str
    ) -> Tuple[list, list]:
        """
        get list of a taxa and their records
        """
        # list for keeping track of taxa per USCO
        og_taxa    = []
        # open alignment file
        records = list(SeqIO.parse(fasta, "fasta"))

        # create lists of USCOtaxa and missing_taxa
        for record in records:
            # append taxa per USCO to USCOtaxa
            og_taxa.append(record.id)

        return og_taxa, records

    def create_concatenation_matrix(
        self,
        alignment_list: str,
        prefix: str
    ):
        # read lists into python lists
        alignments  = [line.rstrip('\n') for line in open(alignment_list)]

        taxa = self.get_taxa_names(alignments)

        # assignment the output file names
        file_partition = prefix + ".partition"
        fasta_output   = prefix + ".fa"
        file_occupancy = prefix + ".occupancy"

        # print log message
        self.start_message_print(
            taxa,
            alignments,
            file_partition,
            fasta_output,
            file_occupancy
        )

        # assign placeholders for lengths in the partition file
        first_len    = 1
        second_len   = 0

        # create a dictionary with keys as taxa and empty list values
        concat = {key: [] for key in taxa}

        # string for first field of partition file
        field_one = 'AUTO'

        # create an empty file_partition
        open(file_partition, "w")
        # create an empty occupancy per OG file
        open(file_occupancy, "w")

        # loop through USCOs files
        for fasta in alignments:

            # create a list of taxa and records
            og_taxa, records = self.get_list_of_taxa_and_records(
                fasta
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
            concat = self.create_missing_taxa_sequence(
                missing_taxa,
                missing_seq,
                concat
            )      

            # create record for present data
            concat = self.create_record_for_present_taxa(
                records,
                og_taxa,
                concat
            )

            # append to partition file
            first_len, second_len = self.partition_file_write(
                file_partition,
                og_len,
                field_one,
                fasta,
                first_len,
                second_len
            )

            # append to occupancy file
            self.occupancy_file_write(
                file_occupancy,
                og_taxa,
                missing_taxa,
                fasta
            )

        # write output fasta file
        self.fasta_file_write(
            fasta_output,
            concat
        )

        # print log message
        print("Complete!\n")
            
            

