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
        alignment_list=None,
        prefix=None,
        idmap=None,
        verbose=None
    ):
        self.alignment_file_path = alignment_file_path
        self.output_file_path = output_file_path
        self.protein_file_path = protein_file_path,
        self.nucleotide_file_path = nucleotide_file_path 
        self.alignment_list = alignment_list
        self.prefix = prefix
        self.fasta = fasta
        self.idmap = idmap
        self.verbose = verbose

    def get_alignment_and_format(self):
        """
        automatic file type determination
        """
        return get_alignment_and_format_helper(self.alignment_file_path)

    def calculate_alignment_length(self, alignment):
        """
        calculate alignment length
        """
        aln_len = alignment.get_alignment_length()
        return aln_len

    def calculate_rcv(self, alignment, aln_len):
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

