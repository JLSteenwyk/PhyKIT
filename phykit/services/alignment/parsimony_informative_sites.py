from collections import Counter
from multiprocessing import Pool
from typing import Dict, Tuple

from Bio.Align import MultipleSeqAlignment

from .base import Alignment


class ParsimonyInformative(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _, _ = self.get_alignment_and_format()
        pi_sites, aln_len, pi_sites_per = \
            self.calculate_parsimony_informative_sites(alignment)

        print(f"{pi_sites}\t{aln_len}\t{round(pi_sites_per, 4)}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(alignment_file_path=args.alignment, cpu=args.cpu)

    def get_number_of_occurrences_per_character(
        self,
        alignment: MultipleSeqAlignment,
        idx: int
    ) -> Counter:
        gap_chars = self.get_gap_chars()
        seq_at_position = alignment[:, idx].upper()
        filtered_seq = filter(lambda c: c not in gap_chars, seq_at_position)

        return Counter(filtered_seq)

    def is_parsimony_informative(
        self,
        num_occurrences: Counter,
    ) -> bool:
        informative_char_count = sum(
            1 for count in num_occurrences.values() if count >= 2
        )
        return informative_char_count >= 2

    def calculate_parsimony_informative_sites(
        self,
        alignment: MultipleSeqAlignment,
    ) -> Tuple[int, int, float]:
        cpu = self.set_cpu()

        aln_len = alignment.get_alignment_length()

        with Pool(cpu) as pool:
            results = pool.map(
                self.check_site_parsimony_informative,
                [(alignment, idx) for idx in range(aln_len)]
            )

        pi_sites = sum(results)
        pi_sites_per = (pi_sites / aln_len) * 100

        return pi_sites, aln_len, pi_sites_per

    def check_site_parsimony_informative(
        self,
        args: Tuple[MultipleSeqAlignment, int]
    ) -> int:
        alignment, idx = args
        num_occurrences = self.get_number_of_occurrences_per_character(alignment, idx)
        return 1 if self.is_parsimony_informative(num_occurrences) else 0
