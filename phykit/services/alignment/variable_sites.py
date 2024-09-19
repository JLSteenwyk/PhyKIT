from typing import Dict, Tuple

from Bio.Align import MultipleSeqAlignment

from .base import Alignment


class VariableSites(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _, is_protein = self.get_alignment_and_format()
        var_sites, aln_len, var_sites_per = \
            self.calculate_variable_sites(alignment)

        print(f"{var_sites}\t{aln_len}\t{round(var_sites_per, 4)}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(alignment_file_path=args.alignment)

    def calculate_variable_sites(
        self,
        alignment: MultipleSeqAlignment
    ) -> Tuple[int, int, float]:
        aln_len = alignment.get_alignment_length()

        gap_chars = self.get_gap_chars()

        var_sites = 0

        for i in range(aln_len):
            seq_at_position = [
                residue.upper()
                for residue in alignment[:, i]
                if residue not in gap_chars
            ]

            if len(set(seq_at_position)) > 1:
                var_sites += 1

        var_sites_per = (var_sites / aln_len) * 100

        return var_sites, aln_len, var_sites_per
