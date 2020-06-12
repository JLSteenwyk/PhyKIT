from enum import Enum

from Bio import AlignIO

from .base import Tree
from ..alignment.base import Alignment
from ...helpers.files import get_alignment_and_format as get_alignment_and_format_helper

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
        alignment, alignment_format = get_alignment_and_format_helper(self.alignment_file_path)
        
        # calculate treeness
        treeness = self.calculate_treeness()
        
        # calculate rcv
        aln_len = alignment.get_alignment_length()
        # TODO: check this with Thomas...still using self as arg...
        relative_composition_variability = Alignment.calculate_rcv(self, alignment, aln_len)
        
        # calculate treeness/rcv
        treeness_over_rcv = self.calculate_treeness_over_rcv(treeness, relative_composition_variability)
        
        # print results
        if (treeness, relative_composition_variability, treeness_over_rcv):
            print(f"{treeness_over_rcv}\t{treeness}\t{relative_composition_variability}")

    def process_args(self, args):
        return dict(tree_file_path=args.tree, alignment_file_path=args.alignment)
    
    def calculate_treeness_over_rcv(self, treeness, relative_composition_variability):
        treeness_over_rcv = (treeness / relative_composition_variability)
        return treeness_over_rcv

