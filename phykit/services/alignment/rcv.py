from .base import Alignment
from ...helpers.json_output import print_json


class RelativeCompositionVariability(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.json_output = parsed["json_output"]

    def run(self):
        # calc rcv and print val
        relative_composition_variability = self.calculate_rcv()
        rcv = round(relative_composition_variability, 4)
        if self.json_output:
            print_json(dict(rcv=rcv))
            return
        print(rcv)

    def process_args(self, args):
        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
        )
