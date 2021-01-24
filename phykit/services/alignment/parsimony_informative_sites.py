import sys

from Bio.Align import MultipleSeqAlignment

from .base import Alignment

class ParsimonyInformative(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _ = self.get_alignment_and_format()
        pi_sites, aln_len, pi_sites_per = self.calculate_parsimony_informative_sites(alignment)
        if (pi_sites, aln_len, pi_sites_per):
            print(f"{pi_sites}\t{aln_len}\t{round(pi_sites_per, 4)}")


    def process_args(self, args):
        return dict(alignment_file_path=args.alignment)

    def get_number_of_occurences_per_character(
        self,
        alignment,
        idx: int
    ) -> dict:

        # obtain sequence at position, remove gaps, and make
        # all characters uppercase
        seq_at_position = ""
        seq_at_position += alignment[:, idx]
        seq_at_position = seq_at_position.upper().replace("-", "")
        num_occurences = {}
        for char in set(seq_at_position.replace("-", "")):
            num_occurences[char] = seq_at_position.count(char)

        return num_occurences

    def count_if_parsimony_informative(
        self,
        num_occurences: dict,
        pi_sites: int
    ):
        # create a dictionary of characters that occur at least twice
        d = dict((k, v) for k, v in num_occurences.items() if v >= 2)

        # determine number of characters that occur at least twice
        if len(d) >= 2:
            pi_sites += 1 
        
        return pi_sites

    def calculate_parsimony_informative_sites(
        self,
        alignment
    ):
        # get aln length
        aln_len = alignment.get_alignment_length()
        
        pi_sites = 0

        # count number of parsimony informative sites
        for idx in range(0, aln_len, int(1)):
            # count occurneces of each character at site idx
            num_occurences = self.get_number_of_occurences_per_character(alignment, idx) 

            # add one to pi_sites if site idx is parsimony informative
            pi_sites = self.count_if_parsimony_informative(num_occurences, pi_sites)

        # calculate percent of variable sites
        pi_sites_per = (pi_sites / aln_len)*100
        
        return pi_sites, aln_len, pi_sites_per