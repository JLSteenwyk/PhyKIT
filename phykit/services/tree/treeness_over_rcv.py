from enum import Enum
from typing import Dict

from .base import Tree
from ..alignment.base import Alignment


# TODO: Here next
class FileFormat(Enum):
    fasta = "fasta"
    clustal = "clustal"
    maf = "maf"
    mauve = "mauve"
    phylip = "phylip"
    phylip_seq = "phylip-sequential"
    phylip_rel = "phylip-relaxed"
    stockholm = "stockholm"


class TreenessOverRCV(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        treeness = self.calculate_treeness()

        aln = Alignment(alignment_file_path=self.alignment_file_path)
        relative_composition_variability = aln.calculate_rcv()

        treeness_over_rcv = treeness / relative_composition_variability

        print(
            f"{round(treeness_over_rcv, 4)}\t{round(treeness, 4)}\t{round(relative_composition_variability, 4)}"
        )

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            alignment_file_path=args.alignment,
        )
