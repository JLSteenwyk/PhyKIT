from Bio.Align import MultipleSeqAlignment

from .base import Alignment

class AlignmentLength(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        alignment, _ = self.get_alignment_and_format()
        aln_len = alignment.get_alignment_length()
        print(aln_len)

    def process_args(self, args):
        return dict(alignment_file_path=args.alignment)
