from enum import Enum
import itertools
import sys
from typing import Dict, List, Tuple

from Bio import Align
from Bio.Phylo import Newick
import numpy as np
from sklearn.linear_model import LinearRegression

from .base import Tree
from ...helpers.files import (
    get_alignment_and_format as get_alignment_and_format_helper
)


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

    def run(self) -> None:
        alignment, _, is_protein = get_alignment_and_format_helper(
            self.alignment_file_path
        )

        tree = self.read_tree_file()

        tips = self.get_tip_names_from_tree(tree)
        combos = list(itertools.combinations(tips, 2))

        (
            patristic_distances,
            uncorrected_distances,
        ) = self.loop_through_combos_and_calculate_pds_and_pis(
            combos, alignment, tree, self.exclude_gaps
        )

        # calculate slope and fit the y-intercept to zero
        # Fitting the y-intercept to zero follows Jeffroy et al.
        # See fig 2 https://www.cell.com/trends/genetics/fulltext/S0168-9525(06)00051-5
        model = LinearRegression(fit_intercept=False)
        model.fit(
            np.array(patristic_distances).reshape(-1, 1),
            np.array(uncorrected_distances)
        )

        self.print_res(
            self.verbose, combos, uncorrected_distances, patristic_distances, model.coef_[0]
        )

    def process_args(self, args) -> Dict[str, str]:
        return dict(
            tree_file_path=args.tree,
            alignment_file_path=args.alignment,
            exclude_gaps=args.exclude_gaps,
            verbose=args.verbose,
        )

    def loop_through_combos_and_calculate_pds_and_pis(
        self,
        combos: List[Tuple[str, str]],
        alignment: Align.MultipleSeqAlignment,
        tree: Newick.Tree,
        exclude_gaps: bool,
    ) -> Tuple[
        List[float],
        List[float]
    ]:
        """
        loop through all taxon combinations and determine
        their patristic distance and pairwise identity
        """
        patristic_distances = []
        uncorrected_distances = []
        gap_chars = self.get_gap_chars()
        aln_len = alignment.get_alignment_length()
        seq_dict = {record.name: record.seq for record in alignment}
        for combo in combos:
            # calculate pds
            patristic_distances.append(tree.distance(combo[0], combo[1]))

            # calcualte uncorrected distances
            seq_one = seq_dict[combo[0]]
            seq_two = seq_dict[combo[1]]
            if exclude_gaps:
                valid_positions = [
                    idx for idx in range(aln_len)
                    if seq_one[idx] not in gap_chars and seq_two[idx] not in gap_chars
                ]
                adjusted_len = len(valid_positions)
                identities = sum(
                    1 for idx in valid_positions if seq_one[idx].upper() == seq_two[idx].upper()
                )
            else:
                adjusted_len = aln_len
                identities = sum(
                    1 for idx in range(aln_len) if seq_one[idx].upper() == seq_two[idx].upper()
                )

            if adjusted_len > 0:
                uncorrected_distances.append(1 - (identities / adjusted_len))
            else:
                uncorrected_distances.append(float('nan'))
                uncorrected_distances.append(1 - (identities / aln_len))

        return patristic_distances, uncorrected_distances

    def print_res(
        self,
        verbose: bool,
        combos: List[str],
        uncorrected_distances: List[float],
        patristic_distances: List[float],
        slope: float,
    ) -> None:
        """
        print results to stdout
        """
        try:
            if verbose:
                for cbo, dist, pd in zip(
                    combos, uncorrected_distances, patristic_distances
                ):
                    print(
                        f"{cbo[0]}\t{cbo[1]}\t{round(dist,4)}\t{round(pd, 4)}"
                    )
            else:
                print(f"{round(slope, 4)}\t{abs(round(1-slope, 4))}")
        except BrokenPipeError:
            pass
