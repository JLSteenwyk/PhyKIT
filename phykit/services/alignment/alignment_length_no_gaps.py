from argparse import Namespace
from multiprocessing import Pool, cpu_count
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
        return dict(alignment_file_path=args.alignment, cpu=args.cpu)

    def calculate_alignment_length_no_gaps(
        self,
        alignment: MultipleSeqAlignment,
        is_protein: bool,
    ) -> Tuple[int, int, float]:
        aln_len = alignment.get_alignment_length()
        aln_len_no_gaps = self.get_sites_no_gaps_count(
            alignment,
            aln_len,
            is_protein,
        )

        aln_len_no_gaps_per = (aln_len_no_gaps / aln_len) * 100

        return aln_len_no_gaps, aln_len, aln_len_no_gaps_per

    def get_sites_no_gaps_count(
        self,
        alignment: MultipleSeqAlignment,
        aln_len: int,
        is_protein: bool,
    ) -> int:
        gap_chars = self.get_gap_chars()

        cpu = self.set_cpu()

        with Pool(cpu) as pool:
            aln_len_no_gaps = pool.starmap(
                self.is_column_no_gap,
                [(alignment[:, i], gap_chars) for i in range(aln_len)]
            )

        return sum(aln_len_no_gaps)

    def is_column_no_gap(self, column: str, gap_chars: set) -> int:
        return 1 if set(column).isdisjoint(gap_chars) else 0
