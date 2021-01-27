from enum import Enum

from Bio import AlignIO

from .base import Tree
from ..alignment.base import Alignment
from ...helpers.files import get_alignment_and_format

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
        # calculate treeness
        treeness = self.calculate_treeness()
        
        # calculate rcv
        aln = Alignment(alignment_file_path=self.alignment_file_path)
        relative_composition_variability = aln.calculate_rcv()

        # calculate treeness/rcv
        treeness_over_rcv = self.calculate_treeness_over_rcv(treeness, relative_composition_variability)
        
        # print results
        print(f"{round(treeness_over_rcv, 4)}\t{round(treeness, 4)}\t{round(relative_composition_variability, 4)}")


    def process_args(self, args):
        return dict(tree_file_path=args.tree, alignment_file_path=args.alignment)
    
    def calculate_treeness_over_rcv(self, treeness, relative_composition_variability):
        treeness_over_rcv = (treeness / relative_composition_variability)
        return treeness_over_rcv

