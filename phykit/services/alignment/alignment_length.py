from .base import Alignment
from ...helpers.json_output import print_json

from typing import Dict


class AlignmentLength(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        alignment, _, _ = self.get_alignment_and_format()
        aln_len = alignment.get_alignment_length()
        if self.json_output:
            print_json(dict(alignment_length=aln_len))
            return
        print(aln_len)

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
        )
