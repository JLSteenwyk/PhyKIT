from typing import Dict

from .base import Tree
from ..alignment.base import Alignment
from ...helpers.json_output import print_json

class TreenessOverRCV(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            alignment_file_path=parsed["alignment_file_path"],
        )
        self.json_output = parsed["json_output"]

    def run(self):
        treeness = self.calculate_treeness()

        aln = Alignment(alignment_file_path=self.alignment_file_path)
        relative_composition_variability = aln.calculate_rcv()

        treeness_over_rcv = treeness / relative_composition_variability

        if self.json_output:
            print_json(
                dict(
                    treeness_over_rcv=round(treeness_over_rcv, 4),
                    treeness=round(treeness, 4),
                    rcv=round(relative_composition_variability, 4),
                )
            )
            return

        print(
            f"{round(treeness_over_rcv, 4)}\t{round(treeness, 4)}\t{round(relative_composition_variability, 4)}"
        )

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
        )
