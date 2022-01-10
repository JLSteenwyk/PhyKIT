from collections import Counter
import itertools
import sys

from Bio import SeqIO, AlignIO

from .base import Alignment

class ColumnScore(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # create biopython object of sequences for query and 
        # reference alignments
        query_records = AlignIO.read(self.fasta, "fasta")
        reference_records = AlignIO.read(self.reference, "fasta")

        # create lists with strings of every columns
        ref_columns, query_columns = self.get_columns_from_alignments(reference_records, query_records)
        
        # count the number of matches and total pairs
        number_of_matches, number_of_total_columns = self.calculate_matches_between_ref_and_query_columns(ref_columns, query_columns)

        # print res
        print(round(number_of_matches/number_of_total_columns, 4))

    def process_args(self, args) -> dict:
        return dict(
            fasta=args.fasta,
            reference=args.reference
        )

    def get_columns_from_alignments(
        self,
        reference_records,
        query_records
    ):
        # loop through reference alignment and get each column
        ref_columns = []
        for i in range(reference_records.get_alignment_length()):
            ref_seq_at_position = ''
            ref_seq_at_position += reference_records[:,i]
            ref_seq_at_position = ref_seq_at_position.upper()
            ref_columns.append(ref_seq_at_position)

        # loop through query alignment and get each column
        query_columns = []
        for i in range(query_records.get_alignment_length()):
            query_seq_at_position = ''
            query_seq_at_position += query_records[:,i]
            query_seq_at_position = query_seq_at_position.upper()
            query_columns.append(query_seq_at_position)

        return ref_columns, query_columns

    def calculate_matches_between_ref_and_query_columns(
        self,
        ref_columns: list,
        query_columns :list
    ):
        # determine the number of matches
        matches = list((Counter(ref_columns) & Counter(query_columns)).elements())
        return len(matches), len(query_columns)






