from enum import Enum
import itertools
from typing import Tuple

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

        # get tip and tip combinations
        tips = self.get_tip_names_from_tree(tree)
        combos = list(itertools.combinations(tips, 2))

        # for pairwise combinations, calculate patristic
        # distances and pairwise identities
        patristic_distances, pairwise_identities = self.loop_through_combos_and_calculate_pds_and_pis(
            combos,
            alignment,
            tree
        )
            
        # calculate linear regression
        _, _, r_value, _, _ = scipy.stats.linregress(pairwise_identities, patristic_distances)
    
        # report res
        self.print_res(
            self.verbose,
            combos,
            pairwise_identities,
            patristic_distances,
            r_value
        )

    def process_args(self, args):
        return dict(
            tree_file_path=args.tree,
            alignment_file_path=args.alignment,
            verbose=args.verbose
        )

    def loop_through_combos_and_calculate_pds_and_pis(
        self,
        combos: list, 
        alignment,
        tree
    ) -> Tuple[list, list]:
        """
        loop through all taxon combinations and determine
        their patristic distance and pairwise identity
        """
        patristic_distances = []
        pairwise_identities = []
        aln_len = alignment.get_alignment_length()
        for combo in combos:
            # calculate pd
            patristic_distances.append(tree.distance(combo[0], combo[1]))
            # calculate pairwise identity
            pairwise_identities = self.calculate_pairwise_identities(
                alignment, pairwise_identities, aln_len, combo
            )
        return patristic_distances, pairwise_identities

    def calculate_pairwise_identities(
        self,
        alignment,
        pairwise_identities: list,
        aln_len: int,
        combo: tuple
    ) -> list:
        """
        calculate pairwise identities for a given combo
        """
        identities = 0
        seq_one = ''
        seq_two = ''
        for record in alignment:
            if record.name == combo[0]:
                seq_one = record.seq
            elif record.name == combo[1]:
                seq_two = record.seq
        for idx in range(0, aln_len):
            if seq_one[idx] == seq_two[idx]:
                identities += 1
        pairwise_identities.append(identities / aln_len)

        return pairwise_identities

    def print_res(
        self,
        verbose: bool,
        combos: list,
        pairwise_identities: list,
        patristic_distances: list,
        r_value: float
    ) -> None:
        """
        print results to stdout
        """
        try:
            if verbose:
                for combo, pairwise_identity, patristic_distance in zip(combos, pairwise_identities, patristic_distances):
                    print(f"{combo[0]}-{combo[1]}\t{round(pairwise_identity,4)}\t{round(patristic_distance, 4)}")
            else:
                print(round(r_value**2, 4))
        except BrokenPipeError:
            pass
