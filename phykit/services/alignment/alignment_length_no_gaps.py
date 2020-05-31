from Bio.Align import MultipleSeqAlignment

from .base import Alignment

class AlignmentLengthNoGaps(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, alignment_format = self.get_alignment_and_format()
        aln_len_no_gaps, aln_len, aln_len_no_gaps_per = self.calculate_alignment_length_no_gaps(alignment)
        if (aln_len_no_gaps, aln_len, aln_len_no_gaps_per):
            print(f"{aln_len_no_gaps}\t{aln_len}\t{aln_len_no_gaps_per}")

    def process_args(self, args):
        return dict(alignment_file_path=args.alignment)

    def calculate_alignment_length_no_gaps(self, alignment):
        aln_len = alignment.get_alignment_length()
        
        aln_len_no_gaps = 0
        # count sites with no gaps
        for i in range(0, aln_len, int(1)):
            seq_at_position = ""
            seq_at_position += alignment[:, i]
            if "-" in seq_at_position:
                aln_len_no_gaps += 1
        # calculate percent of variable sites
        aln_len_no_gaps_per = (aln_len_no_gaps / aln_len)*100
        
        return aln_len_no_gaps, aln_len, aln_len_no_gaps_per