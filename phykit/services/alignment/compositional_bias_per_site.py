from scipy.stats import chisquare, false_discovery_control

from .base import Alignment


class CompositionalBiasPerSite(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _ = self.get_alignment_and_format()
        stat_res, p_vals_corrected = self.calculate_compositional_bias_per_site(
            alignment
        )
        idx = 0
        for stat_info, pval_cor in zip(stat_res, p_vals_corrected):
            if isinstance(pval_cor, str):
                print(f"{idx+1}\t{round(stat_info.statistic,4)}\tnan\t{round(stat_info.pvalue,4)}")
            else:
                print(f"{idx+1}\t{round(stat_info.statistic,4)}\t{round(pval_cor,4)}\t{round(stat_info.pvalue,4)}")
            idx += 1

    def process_args(self, args):
        return dict(alignment_file_path=args.alignment)

    def get_number_of_occurences_per_character(self, alignment, idx: int) -> list:
        # obtain sequence at position, remove gaps, and make
        # all characters uppercase
        seq_at_position = ""
        seq_at_position += alignment[:, idx]
        seq_at_position = seq_at_position.upper().replace("-", "")
        num_occurences = []
        for char in set(seq_at_position.replace("-", "")):
            num_occurences.append(seq_at_position.count(char))

        return num_occurences

    def calculate_compositional_bias_per_site(self, alignment):
        # get aln length
        aln_len = alignment.get_alignment_length()

        stat_res = []
        p_vals = []
        nan_idx = []

        for idx in range(0, aln_len):
            # count occurneces of each character at site idx
            num_occurences = self.get_number_of_occurences_per_character(
                alignment, idx
            )
            
            chisquare_res = chisquare(num_occurences)

            stat_res.append(chisquare_res)

            if str(chisquare_res.pvalue) != "nan":
                p_vals.append(chisquare_res.pvalue)
            else:
                nan_idx.append(idx)

        p_vals_corrected = list(false_discovery_control(p_vals))

        for idx in reversed(nan_idx):
            p_vals_corrected.insert(idx, "nan")

        return stat_res, p_vals_corrected
