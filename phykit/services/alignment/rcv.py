from Bio import SeqIO
from Bio.Seq import Seq

from .base import Alignment

class RelativeCompositionVariability(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # calc rcv and print val
        relative_composition_variability = self.calculate_rcv()
        try:
            print(round(relative_composition_variability, 4))
        except BrokenPipeError:
            pass

    def process_args(self, args):
        return dict(alignment_file_path=args.alignment)