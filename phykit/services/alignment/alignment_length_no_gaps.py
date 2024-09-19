from argparse import Namespace
from typing import Dict, Tuple

from Bio.Align import MultipleSeqAlignment

from .base import Alignment


class AlignmentLengthNoGaps(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()
        (
            aln_len_no_gaps,
            aln_len,
            aln_len_no_gaps_per,
        ) = self.calculate_alignment_length_no_gaps(alignment, is_protein)
        print(f"{aln_len_no_gaps}\t{aln_len}\t{round(aln_len_no_gaps_per, 4)}")

    def process_args(
        self,
        args: Namespace,
    ) -> Dict[str, str]:
        return dict(alignment_file_path=args.alignment)

    def calculate_alignment_length_no_gaps(
        self,
        alignment: MultipleSeqAlignment,
        is_protein: bool,
    ) -> Tuple[int, int, float]:
        aln_len = alignment.get_alignment_length()
        aln_len_no_gaps = self.get_sites_no_gaps_count(
            alignment,
            aln_len,
            is_protein
        )

        aln_len_no_gaps_per = (aln_len_no_gaps / aln_len) * 100

        return aln_len_no_gaps, aln_len, aln_len_no_gaps_per

    def get_sites_no_gaps_count(
        self,
        alignment: MultipleSeqAlignment,
        aln_len: int,
        is_protein: bool,
    ) -> int:
        """
        Count sites in the alignment with no gaps
        """
        aln_len_no_gaps = 0

        gap_chars = self.get_gap_chars()

        for i in range(aln_len):
            column = set(alignment[:, i])
            if column.isdisjoint(gap_chars):
                aln_len_no_gaps += 1

        return aln_len_no_gaps
