import itertools
from typing import Dict, List, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .base import Alignment


class SumOfPairsScore(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        query_records = SeqIO.to_dict(SeqIO.parse(self.fasta, "fasta"))
        reference_records = SeqIO.to_dict(SeqIO.parse(self.reference, "fasta"))

        record_id_pairs = list(
            itertools.combinations(reference_records.keys(), 2)
        )

        number_of_matches, number_of_total_pairs = \
            self.determine_number_of_matches_and_total_pairs(
                record_id_pairs, reference_records, query_records
            )

        print(round(number_of_matches / number_of_total_pairs, 4))

    def process_args(self, args) -> Dict[str, str]:
        return dict(fasta=args.fasta, reference=args.reference)

    def determine_number_of_matches_and_total_pairs(
        self,
        record_id_pairs: List[Tuple[str, str]],
        reference_records: Dict[str, SeqRecord],
        query_records: Dict[str, SeqRecord],
    ) -> Tuple[int, int]:
        number_of_matches = 0
        number_of_total_pairs = 0

        for first_in_pair, second_in_pair in record_id_pairs:
            ref_seq1 = reference_records[first_in_pair].seq
            ref_seq2 = reference_records[second_in_pair].seq
            query_seq1 = query_records[first_in_pair].seq
            query_seq2 = query_records[second_in_pair].seq

            for ref_res1, ref_res2, query_res1, query_res2 in zip(
                ref_seq1, ref_seq2, query_seq1, query_seq2
            ):
                number_of_total_pairs += 1
                if ref_res1 == query_res1 and ref_res2 == query_res2:
                    number_of_matches += 1

        return number_of_matches, number_of_total_pairs
