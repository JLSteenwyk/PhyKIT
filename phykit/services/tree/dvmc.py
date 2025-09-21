import math
from typing import Dict
import numpy as np

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

        # Collect all distances at once for vectorized operations
        distances = np.array([tree.distance(term) for term in tree.get_terminals()])

        # Calculate statistics using numpy
        sum_dist = np.sum(distances)
        sumi2N = np.sum(distances ** 2)
        avg_dist = np.mean(distances)

        # Calculate variance more efficiently
        squared_diff_sum = sumi2N - num_spp * (avg_dist ** 2)

        # Return standard deviation
        return np.sqrt(squared_diff_sum / (num_spp - 1))
