import math
from typing import Dict

from Bio.Phylo import Newick

from .base import Tree


class DVMC(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        dvmc = self.determine_dvmc(tree)
        print(round(dvmc, 4))

    def process_args(self, args) -> Dict[str, str]:
        return dict(tree_file_path=args.tree)

    def determine_dvmc(self, tree: Newick.Tree) -> float:
        num_spp = tree.count_terminals()
        sum_dist = float()
        sumi2N = float()

        for term in tree.get_terminals():
            dist = tree.distance(term)
            sum_dist += dist
            sumi2N += dist ** 2

        avg_dist = sum_dist / num_spp

        squared_diff_sum = sumi2N - num_spp * (avg_dist ** 2)

        return math.sqrt(squared_diff_sum / (num_spp - 1))
