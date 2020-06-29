from Bio.Align import MultipleSeqAlignment

from .base import Alignment

class AlignmentLengthNoGaps(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, alignment_format = self.get_alignment_and_format()
        aln_len_no_gaps, aln_len, aln_len_no_gaps_per = self.calculate_alignment_length_no_gaps(alignment)
        print(f"{aln_len_no_gaps}\t{aln_len}\t{round(aln_len_no_gaps_per, 4)}")

    def process_args(self, args):
        return dict(alignment_file_path=args.alignment)

    def calculate_alignment_length_no_gaps(self, alignment):
        aln_len = alignment.get_alignment_length()
        aln_len_no_gaps = self.get_sites_no_gaps_count(alignment, aln_len)

        # calculate percent of variable sites
        aln_len_no_gaps_per = (aln_len_no_gaps / aln_len) * 100
        
        return aln_len_no_gaps, aln_len, aln_len_no_gaps_per

    def get_sites_no_gaps_count(self, alignment, aln_len):
        """
        Count sites in the alignment with no gaps
        """
        aln_len_no_gaps = 0
        for i in range(0, aln_len):
            seq_at_position = ""
            seq_at_position += alignment[:, i]
            if "-" not in seq_at_position:
                aln_len_no_gaps += 1

        return aln_len_no_gaps