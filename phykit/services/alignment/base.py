import sys

from ..base import BaseService
from ...helpers.files import get_alignment_and_format as get_alignment_and_format_helper

class Alignment(BaseService):
    def __init__(
        self,
        *args,
        alignment_file_path=None,
        fasta=None,
        output_file_path=None,
        protein_file_path=None,
        nucleotide_file_path=None,
        alignment_list_path=None,
        prefix=None,
        idmap=None,
        reference=None,
        verbose=None,
        entry=None,
        clipkit_log_file=None,
    ):
        self.alignment_file_path = alignment_file_path
        self.output_file_path = output_file_path
        self.protein_file_path = protein_file_path,
        self.nucleotide_file_path = nucleotide_file_path 
        self.alignment_list_path = alignment_list_path
        self.prefix = prefix
        self.fasta = fasta
        self.idmap = idmap
        self.reference = reference
        self.verbose = verbose
        self.entry = entry
        self.clipkit_log_file = clipkit_log_file

    def get_alignment_and_format(self):
        """
        automatic file type determination
        """
        try:
            return get_alignment_and_format_helper(self.alignment_file_path)
        except FileNotFoundError:
            print("Input corresponds to no such file or directory.")
            print("Please double check pathing and filenames")
            sys.exit()
            
    def calculate_rcv(self):
        alignment, _ = self.get_alignment_and_format()
        aln_len = alignment.get_alignment_length()

        # string to hold all sequences
        concat_seq = ''
        # initialize a counter for the number of sequences in the input fasta file
        num_records = 0

        # for each record join concatSeq string and sequence as well as keeping track 
        # of the number of records
        for record in alignment:
            concat_seq  += record.seq
            num_records += 1

        # dictionary to hold the average occurence of each sequence letter
        average_d = {}
        # loop through the different sequences that appear in the fasta file
        # population dictionary with how many times that sequence appears
        for seq in set(concat_seq):
            average_d[seq] = (concat_seq.count(seq)/num_records)

        # intiailize list to hold the RCV values per ith taxa 
        # that will later be summed
        indiv_rcv_values = []

        # loop through records again and calculate RCV for 
        # each taxa and append to indivRCVvalues
        for record in alignment:
            # temp holds a temporary value of the numerator before appending
            # to numeratorRCVvalues and then is reassigned to 0 when it goes
            # through the loop again
            temp = 0
            # calculates the absolute value of the ith sequence letter minus the average
            for seq_letter in set(concat_seq):
                temp += abs(record.seq.count(seq_letter)-average_d[seq_letter])
            indiv_rcv_values.append(temp/(num_records*aln_len))

        relative_composition_variability = sum(indiv_rcv_values)

        # print the sum of all RCV values
        return relative_composition_variability

