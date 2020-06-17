from enum import Enum
import itertools

from Bio import AlignIO
import scipy

from .base import Tree
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

class Saturation(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # read in alignment
        alignment, alignment_format = get_alignment_and_format_helper(self.alignment_file_path)
        
        # read in tree
        tree = self.read_tree_file()

        # combinations of tip names
        tips = []
        for term in tree.get_terminals():
            tips.append(term.name)
        combos = list(itertools.combinations(tips, 2))

        # get pds and pairwise identifies per pair
        patristic_distances = []
        pairwise_identities = []
        aln_len = alignment.get_alignment_length()

        for combo in combos:
            # calculate pd
            patristic_distances.append(tree.distance(combo[0], combo[1]))
            # calculate pairwise identity
            identities = 0
            for record in alignment:
                if record.name == combo[0]:
                    seq_one = record.seq
                elif record.name == combo[1]:
                    seq_two = record.seq
            for idx in range(0, aln_len):
                if seq_one[idx] == seq_two[idx]:
                    identities += 1
            pairwise_identities.append(identities / aln_len)

        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(pairwise_identities, patristic_distances)
    
        if self.verbose:
            for combo, pairwise_identity, patristic_distance in zip(combos, pairwise_identities, patristic_distances):
                print(f"{combo[0]}-{combo[1]}\t{pairwise_identity}\t{patristic_distance}")
        else:
            print(r_value**2)

    def process_args(self, args):
        return dict(tree_file_path=args.tree, alignment_file_path=args.alignment, verbose=args.verbose)

