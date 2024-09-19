from typing import Dict, List, Tuple, Union
from collections import Counter

from scipy.stats import chisquare, false_discovery_control
from scipy.stats._stats_py import Power_divergenceResult
from Bio.Align import MultipleSeqAlignment

from .base import Alignment


class CompositionalBiasPerSite(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()

        stat_res, p_vals_corrected = \
            self.calculate_compositional_bias_per_site(alignment)

        for idx, (stat_info, pval_cor) in enumerate(
            zip(stat_res, p_vals_corrected), start=1
        ):
            pval_cor_str = "nan" if isinstance(pval_cor, str) else round(pval_cor, 4)
            print(f"{idx}\t{round(stat_info.statistic, 4)}\t{pval_cor_str}\t{round(stat_info.pvalue, 4)}")

    def process_args(self, args) -> Dict[str, str]:
        return dict(alignment_file_path=args.alignment)

    def get_number_of_occurrences_per_character(
        self,
        alignment: MultipleSeqAlignment,
        idx: int,
    ) -> List[int]:
        gap_chars = self.get_gap_chars()
        seq_at_position = alignment[:, idx].upper()
        filtered_seq = "".join([char for char in seq_at_position if char not in gap_chars])

        return list(Counter(filtered_seq).values())

    def calculate_compositional_bias_per_site(
        self,
        alignment: MultipleSeqAlignment,
    ) -> Tuple[
        List[Power_divergenceResult],
        List[Union[float, str]],
    ]:
        aln_len = alignment.get_alignment_length()

        stat_res = []
        p_vals = []
        nan_idx = []

        for idx in range(aln_len):
            num_occurrences = \
                self.get_number_of_occurrences_per_character(alignment, idx)

            chisquare_res = chisquare(num_occurrences)
            stat_res.append(chisquare_res)

            if str(chisquare_res.pvalue) != "nan":
                p_vals.append(chisquare_res.pvalue)
            else:
                nan_idx.append(idx)

        p_vals_corrected = list(false_discovery_control(p_vals))

        for idx in reversed(nan_idx):
            p_vals_corrected.insert(idx, "nan")

        return stat_res, p_vals_corrected
