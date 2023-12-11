import math
from typing import Tuple

from .base import Tree

from ...helpers.files import read_single_column_file_to_list


class DVMC(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        tree = self.read_tree_file()
        dvmc = self.determine_dvmc(tree)

        print(round(dvmc, 4))

    def process_args(self, args):
        return dict(tree_file_path=args.tree)

    def get_term_to_root_dist_and_sum_of_distances(self, tree) -> Tuple[float, float]:
        """
        calculate root to tip distances and
        the sum of root to tip distances
        """
        dist = []
        sum_dist = 0
        for term in tree.get_terminals():
            # append root to tip distance to dist
            dist.append((tree.distance(term)))
            # keep running sum of tree distances
            sum_dist += tree.distance(term)

        return dist, sum_dist

    def calculate_dvmc(self, dist: float, sum_dist: float, num_spp: int) -> float:
        """
        calculate dvmc from tip to root distances
        """
        # determine dvmc
        # calculate average tree distance
        avg_dist = sum_dist / num_spp

        # determine the sum of i=1 to N for (x_i-x_bar)^2
        sumi2N = float(0.0)
        for x_i in dist:
            sumi2N += (x_i - avg_dist) ** 2

        # multiply sumi2N by 1/(N-1) where N is the number of spp
        # and then take the square root
        dvmc = float(0.0)
        dvmc = math.sqrt((1 / (num_spp - 1)) * sumi2N)

        return dvmc

    def determine_dvmc(self, tree: Tree):
        # loop through terminal branches and store
        # distances from the root to the tip in a list.
        # Also, calc the sum of all tip to root distances
        dist, sum_dist = self.get_term_to_root_dist_and_sum_of_distances(tree)

        # determine number of taxa in the tree
        num_spp = tree.count_terminals()

        # calculate and return dvmc
        return self.calculate_dvmc(dist, sum_dist, num_spp)
