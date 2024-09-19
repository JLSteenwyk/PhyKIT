from multiprocessing import Pool
from typing import Dict, List, Tuple

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from .base import Alignment


class ColumnScore(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        query_records = AlignIO.read(self.fasta, "fasta")
        reference_records = AlignIO.read(self.reference, "fasta")

        ref_columns, query_columns = self.get_columns_from_alignments(
            reference_records, query_records
        )

        number_of_matches, number_of_total_columns = \
            self.calculate_matches_between_ref_and_query_columns(
                ref_columns, query_columns
            )

        print(round(number_of_matches / number_of_total_columns, 4))

    def process_args(self, args) -> Dict[str, str]:
        return dict(fasta=args.fasta, reference=args.reference)

    def get_columns_from_alignments(
        self,
        reference_records: MultipleSeqAlignment,
        query_records: MultipleSeqAlignment,
    ) -> Tuple[List[str], List[str]]:
        ref_aln_len = reference_records.get_alignment_length()
        query_aln_len = query_records.get_alignment_length()

        cpu = self.set_cpu()

        with Pool(cpu) as pool:
            ref_columns = pool.starmap(
                self.get_column, [
                    (reference_records, i) for i in range(ref_aln_len)
                ]
            )
            query_columns = pool.starmap(
                self.get_column, [
                    (query_records, i) for i in range(query_aln_len)
                ]
            )

        return ref_columns, query_columns

    def get_column(self, alignment: MultipleSeqAlignment, index: int) -> str:
        return alignment[:, index].upper()

    def calculate_matches_between_ref_and_query_columns(
        self,
        ref_columns: List[str],
        query_columns: List[str],
    ) -> Tuple[int, int]:
        set1 = set(ref_columns)
        set2 = set(query_columns)

        matches = set1.intersection(set2)

        return len(matches), len(query_columns)
