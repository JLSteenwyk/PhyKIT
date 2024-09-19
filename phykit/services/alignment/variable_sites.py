from multiprocessing import Pool
from typing import Dict, Tuple

from Bio.Align import MultipleSeqAlignment

from .base import Alignment


class VariableSites(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _, _ = self.get_alignment_and_format()
        var_sites, aln_len, var_sites_per = \
            self.calculate_variable_sites(alignment)

        print(f"{var_sites}\t{aln_len}\t{round(var_sites_per, 4)}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(alignment_file_path=args.alignment, cpu=args.cpu)

    def calculate_variable_sites(
        self,
        alignment: MultipleSeqAlignment
    ) -> Tuple[int, int, float]:
        aln_len = alignment.get_alignment_length()
        gap_chars = self.get_gap_chars()

        cpu = self.set_cpu()

        with Pool(cpu) as pool:
            results = pool.map(
                self.check_site_variability,
                [(alignment[:, i], gap_chars) for i in range(aln_len)]
            )

        var_sites = sum(results)
        var_sites_per = (var_sites / aln_len) * 100

        return var_sites, aln_len, var_sites_per

    def check_site_variability(self, args: Tuple[str, set]) -> int:
        seq_at_position, gap_chars = args
        seq_at_position = [
            residue.upper()
            for residue in seq_at_position
            if residue not in gap_chars
        ]

        return 1 if len(set(seq_at_position)) > 1 else 0
