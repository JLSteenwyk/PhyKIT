from enum import Enum

from Bio import AlignIO

from .base import Tree
from .treeness import Treeness
from ..alignment.rcv import RelativeCompositionVariability

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
        tree = self.read_tree_file()
        alignment, alignment_format = self.get_alignment_and_format()
        
        # calculate treeness
        treeness = Treeness.calculate_treeness(self, tree)
        
        # calculate rcv
        aln_len = alignment.get_alignment_length()
        relative_composition_variability = RelativeCompositionVariability.calculate_rcv(self, alignment, aln_len)
        
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

    # TODO: import this function from alignment base
    def get_alignment_and_format(self):
        """
        automatic file type determination
        """

        # if file format is provided, read the file according to the user's file format
        for fileFormat in FileFormat:
            try:
                alignment = AlignIO.read(open(self.alignment_file_path), fileFormat.value)
                return alignment, fileFormat.value
            # the following exceptions refer to skipping over errors
            # associated with reading the wrong input file
            except ValueError:
                continue
            except AssertionError:
                continue

        raise Exception("Input file could not be read")

