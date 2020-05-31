from Bio.Align import MultipleSeqAlignment

from .base import Alignment

class AlignmentLength(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, file_format = self.get_alignment_and_format()
        aln_len = self.calculate_alignment_length(alignment)
        if aln_len:
            print(f"{aln_len}")

    def process_args(self, args):
        return dict(alignment_file_path=args.alignment)
