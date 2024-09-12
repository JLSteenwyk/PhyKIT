from collections import Counter
from typing import List, Dict

from Bio.Align import MultipleSeqAlignment

from .base import Alignment


class EvolutionaryRatePerSite(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _, is_protein = self.get_alignment_and_format()
        pic_values = self.calculate_evolutionary_rate_per_site(alignment)

        for idx, value in enumerate(pic_values):
            print(f"{idx + 1}\t{round(value, 4)}")

    def process_args(self, args):
        return dict(alignment_file_path=args.alignment)

    def remove_gap_characters(self, seq: str, gap_chars: List[str]) -> str:
        return ''.join([char for char in seq if char not in gap_chars]).upper()

    def get_number_of_occurrences_per_character(
        self,
        alignment: MultipleSeqAlignment,
        idx: int,
        gap_chars: List[str]
    ) -> Dict[str, int]:
        seq_at_position = alignment[:, idx]
        clean_seq = self.remove_gap_characters(seq_at_position, gap_chars)

        return Counter(clean_seq)

    def calculate_pic(
        self,
        num_occurrences: Dict[str, int],
    ) -> float:
        total_frequencies = sum(num_occurrences.values())
        sum_of_frequencies = sum(
            (frequency / total_frequencies) ** 2
            for frequency in num_occurrences.values()
        )
        return 1 - sum_of_frequencies

    def calculate_evolutionary_rate_per_site(
        self,
        alignment: MultipleSeqAlignment,
    ) -> List[float]:
        aln_len = alignment.get_alignment_length()
        pic_values = []

        gap_chars = self.get_gap_chars()

        # calculate PIC for each site in the alignment
        for idx in range(aln_len):
            num_occurrences = self.get_number_of_occurrences_per_character(
                alignment, idx, gap_chars
            )
            pic_values.append(self.calculate_pic(num_occurrences))

        return pic_values
