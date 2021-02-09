from collections import Counter
import itertools
import sys

from Bio import SeqIO

from .base import Alignment

class SumOfPairsScore(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # create biopython object of sequences for query and 
        # reference alignments
        query_records = SeqIO.to_dict(SeqIO.parse(self.fasta, "fasta"))
        reference_records = SeqIO.to_dict(SeqIO.parse(self.reference, "fasta"))

        # get all record pairs
        record_id_pairs = self.get_record_ids(reference_records)

        # calculate how many matches there are and how many total pairs there are
        number_of_matches, number_of_total_pairs = self.determine_number_of_matches_and_total_pairs(
            record_id_pairs,
            reference_records,
            query_records
        )
        
        # print res
        print(round(number_of_matches/number_of_total_pairs, 4))

    def process_args(self, args) -> dict:
        return dict(
            fasta=args.fasta,
            reference=args.reference
        )

    def get_record_ids(
        self,
        reference_records: dict
    ) -> list:
        # loop through record names and save each to 
        record_ids = []
        for entry_name in reference_records:
            record_ids.append(entry_name)
        # create all pairwise combinations
        record_id_pairs = list(itertools.combinations(record_ids, 2))
        return record_id_pairs

    def determine_number_of_matches_and_total_pairs(
        self,
        record_id_pairs: list,
        reference_records: dict,
        query_records: dict
    ):
        number_of_matches = 0
        number_of_total_pairs = 0
        # loop through each pair
        for record_pair in record_id_pairs:
            first_in_pair = record_pair[0]
            second_in_pair = record_pair[1]
            
            pairs_in_reference = []
            pairs_in_query = []
            # for each pair, loop through the length of the alignment and get sequence at each site
            for i in range(0, len(reference_records[first_in_pair].seq)):
                pairs_in_reference.append(reference_records[first_in_pair].seq[i] + reference_records[second_in_pair].seq[i])
            for i in range(0, len(query_records[first_in_pair].seq)):
                pairs_in_query.append(query_records[first_in_pair].seq[i] + query_records[second_in_pair].seq[i])

            # count the number of matches and total pairs
            matches = list((Counter(pairs_in_reference) & Counter(pairs_in_query)).elements())
            number_of_matches+=len(matches)
            number_of_total_pairs+=len(pairs_in_reference)

        return number_of_matches, number_of_total_pairs





