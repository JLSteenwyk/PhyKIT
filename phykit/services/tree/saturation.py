from enum import Enum
import itertools
import sys
from typing import Tuple

import numpy as np
from sklearn.linear_model import LinearRegression

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
        alignment, alignment_format = get_alignment_and_format_helper(
            self.alignment_file_path
        )

        # read in tree
        tree = self.read_tree_file()

        # get tip and tip combinations
        tips = self.get_tip_names_from_tree(tree)
        combos = list(itertools.combinations(tips, 2))

        # for pairwise combinations, calculate patristic
        # distances and pairwise identities
        (
            patristic_distances,
            uncorrected_distances,
        ) = self.loop_through_combos_and_calculate_pds_and_pis(combos, alignment, tree)



        # calculate slope and fit the y-intercept to zero
        # Fitting the y-intercept to zero follows Jeffroy et al.
        # See fig 2 https://www.cell.com/trends/genetics/fulltext/S0168-9525(06)00051-5
        model = LinearRegression(fit_intercept=False)
        model.fit(
            np.array(patristic_distances).reshape(-1, 1),
            np.array(uncorrected_distances)
        )
        # report res
        self.print_res(
            self.verbose, combos, uncorrected_distances, patristic_distances, model.coef_[0]
        )

    def process_args(self, args):
        return dict(
            tree_file_path=args.tree,
            alignment_file_path=args.alignment,
            verbose=args.verbose,
        )

    def loop_through_combos_and_calculate_pds_and_pis(
        self, combos: list, alignment, tree
    ) -> Tuple[list, list]:
        """
        loop through all taxon combinations and determine
        their patristic distance and pairwise identity
        """
        patristic_distances = []
        uncorrected_distances = []
        aln_len = alignment.get_alignment_length()
        for combo in combos:
            # calculate pd
            patristic_distances.append(tree.distance(combo[0], combo[1]))
            # calculate pairwise identity
            uncorrected_distances = self.calculate_uncorrected_distances(
                alignment, uncorrected_distances, aln_len, combo
            )

        return patristic_distances, uncorrected_distances

    def calculate_uncorrected_distances(
        self, alignment, uncorrected_distances: list, aln_len: int, combo: tuple
    ) -> list:
        """
        calculate uncorrected distances for a given combo
        """
        identities = 0
        seq_one = ""
        seq_two = ""
        for record in alignment:
            if record.name == combo[0]:
                seq_one = record.seq
            elif record.name == combo[1]:
                seq_two = record.seq
        for idx in range(0, aln_len):
            try:
                if seq_one[idx] == seq_two[idx]:
                    identities += 1
            except IndexError:
                print("Error: the alignment FASTA headers and tree tip labels are different.\nDouble check that the sequence alignment and phylogenetic tree have the same labels.")
                sys.exit()
        uncorrected_distances.append(1-(identities / aln_len))

        return uncorrected_distances

    def print_res(
        self,
        verbose: bool,
        combos: list,
        uncorrected_distances: list,
        patristic_distances: list,
        slope: float,
    ) -> None:
        """
        print results to stdout
        """
        try:
            if verbose:
                for combo, dist, patristic_distance in zip(
                    combos, uncorrected_distances, patristic_distances
                ):
                    print(
                        f"{combo[0]}\t{combo[1]}\t{round(dist,4)}\t{round(patristic_distance, 4)}"
                    )
            else:
                print(f"{round(slope, 4)}\t{abs(round(1-slope, 4))}")
        except BrokenPipeError:
            pass
