from Bio import SeqIO
from Bio.Seq import Seq

from .base import Alignment

class RelativeCompositionVariability(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # get alignment and alignment format
        alignment, alignment_format = self.get_alignment_and_format()

        # get alignment length
        aln_len = self.calculate_alignment_length(alignment)

        # calc rcv and print val
        relative_composition_variability = self.calculate_rcv(alignment, aln_len)
        print(f"{relative_composition_variability}")

    def process_args(self, args):
        return dict(alignment_file_path=args.alignment)