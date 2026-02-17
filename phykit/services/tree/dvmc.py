import math
from typing import Dict
import numpy as np

from Bio.Phylo import Newick

from .base import Tree
from ...helpers.json_output import print_json


class DVMC(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file()
        dvmc = self.determine_dvmc(tree)
        if self.json_output:
            print_json(dict(dvmc=round(dvmc, 4)))
            return
        print(round(dvmc, 4))

    def process_args(self, args) -> Dict[str, str]:
        return dict(tree_file_path=args.tree, json_output=getattr(args, "json", False))

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
