from .base import Alignment

from typing import Dict


class AlignmentLength(Alignment):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self) -> None:
        alignment, _, _ = self.get_alignment_and_format()
        aln_len = alignment.get_alignment_length()
        print(aln_len)

    def process_args(self, args) -> Dict[str, str]:
        return dict(alignment_file_path=args.alignment)
