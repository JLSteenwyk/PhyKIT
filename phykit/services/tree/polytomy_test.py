import sys
import getopt
import os.path
import itertools
from scipy.stats import chisquare

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

from .base import Tree

from ...helpers.files import read_single_column_file_to_list

class PolytomyTest(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # read in groups
        groups_arr = self.read_in_groups(self.groups)

        # determine groups of groups
        groups_of_groups = self.determine_groups_of_groups(groups_arr)

        # read trees into list
        trees_file_path = read_single_column_file_to_list(self.trees)

        # loop through trees
        summary = {}
        for tree_file in trees_file_path:
            tree = Phylo.read(tree_file, 'newick')
            # get tip names
            tips = self.get_tip_names_from_tree(tree)

            # get all combinations of three tips
            triplet_tips = itertools.combinations(tips, 3)
            for triplet in triplet_tips:
                # determine tips that are not in the triplet of interest
                tips_to_prune = list(set(tips)- set(list(triplet)))
                tree = Phylo.read(tree_file, 'newick')
                
                # prune to a triplet
                tree = self.prune_tree_using_taxa_list(tree, tips_to_prune)
                
                for _, groups in groups_of_groups.items():
                    # see if there any intersctions between
                    # the triplet and the the group
                    num_groups_represented = self.count_number_of_groups_in_triplet(triplet, groups)
                    
                    # if one taxa is represented from each group, 
                    # use the triplet
                    tip_names = []
                    if num_groups_represented == 3:
                        # get names in triplet and set tree branch lengths to 1
                        tip_names = self.get_tip_names_from_tree(tree)
                        self.set_branch_lengths_in_tree_to_one(tree)
                        # get pairs from tip names
                        pairs = list(itertools.combinations(tip_names, 2))
                        for pair in pairs:
                            is_polytomy = self.check_if_triplet_is_a_polytomy(tree)
                            # if distance between pair is 2 and the triplet is 
                            # not a polytomy (i.e., having only 1 internal branch)
                            # then report the sisters in the triplet
                            if tree.distance(pair[0], pair[1]) == 2 and not is_polytomy:
                                sisters = self.determine_sisters_from_triplet(groups, pair)
                                
                                if tree_file not in summary.keys():
                                    summary[str(tree_file)]={}
                                # if the sister relationship is not in the tree file dict, create a key for it
                                if sisters not in summary[str(tree_file)].keys():
                                    summary[str(tree_file)][sisters] = 1
                                else:
                                    summary[str(tree_file)][sisters] += 1

        # Also, count the total number of sister pairings
        # for the three possible pairs.    
        group0_and_group1_count = 0
        group0_and_group2_count = 0
        group1_and_group2_count = 0

        # Also, keep track of which and how many genes 
        # support each sister pairing
        gene_support_freq = {
            '0-1' : 0,
            '1-2' : 0,
            '0-2' : 0
        }

        for tree in summary:
            # create empty key value pairs in case sister 
            # pairing was never observed
            if '0-1' not in summary[tree].keys():
                summary[tree]['0-1'] = 0 
            if '0-2' not in summary[tree].keys():
                summary[tree]['0-2'] = 0 
            if '1-2' not in summary[tree].keys():
                summary[tree]['1-2'] = 0
            # create a running value of triplets that support each sister pair
            group0_and_group1_count += summary[tree]['0-1']
            group0_and_group2_count += summary[tree]['0-2']
            group1_and_group2_count += summary[tree]['1-2']
            # determine which sister pairing is best supported and add one to
            # the corresponding gene support frequency count
            gene_support_freq[max(summary[tree], key=summary[tree].get)] += 1

        triplet_res = chisquare(
            [group0_and_group1_count, group0_and_group2_count, group1_and_group2_count]
        )
        
        gene_support_freq_res = chisquare(
            [gene_support_freq['0-1'], gene_support_freq['0-2'], gene_support_freq['1-2']]
        )

        self.print_gene_support_freq_res(gene_support_freq_res, gene_support_freq, trees_file_path)
        self.print_triplet_based_res(
            triplet_res,
            group0_and_group1_count,
            group0_and_group2_count,
            group1_and_group2_count
        )
    
    def process_args(self, args):
        return dict(trees=args.trees, groups=args.groups)

    def read_in_groups(self, groups) -> list:
        groups_arr = []
        for line in open(self.groups):
            line = line.strip()
            if not line.startswith("#"):
                line = line.split("\t")
                temp = []
                temp.append(line[0])
                temp.append(line[1].split(";"))
                temp.append(line[2].split(";"))
                temp.append(line[3].split(";"))
                groups_arr.append(temp)
        return groups_arr

    def determine_groups_of_groups(self, groups_arr: list) -> dict:
        groups_of_groups = {}

        for group in groups_arr:
            temp = []
            for i in range(1, 4):
                temp.append([taxon_name for taxon_name in group[i]])
            groups_of_groups[group[0]] = (temp)
        
        return groups_of_groups

    def count_number_of_groups_in_triplet(self, triplet: list, groups: tuple) -> int:
        """
        determine how many groups are represented in a triplet
        """
        num_groups_represented = 0
        for group in groups:
            temp = [value for value in list(triplet) if value in group]
            if len(temp) >= 1:
                num_groups_represented +=1
        return num_groups_represented

    def set_branch_lengths_in_tree_to_one(self, tree):
        for term in tree.get_terminals():
            term.branch_length = 1
        for internode in tree.get_nonterminals():
            internode.branch_length = 1

    def check_if_triplet_is_a_polytomy(self, tree):
        """
        count the number of internal branches. If 1, then the triplet is a polytomy
        """
        num_int=0
        # check if the triplet is a polytomy
        for internal_branch in tree.get_nonterminals():
            num_int+=1

        if num_int == 1:
            return True
        else:
            return False

        #return num_int

    def determine_sisters_from_triplet(self, groups: list, pair: tuple) -> str:
        """
        determine sister taxa from a triplet
        """
        sisters = sorted([[i for i, lst in enumerate(groups) if pair[0] in lst][0], [i for i, lst in enumerate(groups) if pair[1] in lst][0]])
        sisters = [str(i) for i in sisters] 
        sisters = str("-".join(sisters))

        return sisters

    def print_triplet_based_res(
        self,
        triplet_res,
        group0_and_group1_count: int,
        group0_and_group2_count: int,
        group1_and_group2_count: int
    ) -> None:
        """
        print results to stdout for user
        """
        print(f"\nTriplet Results")
        print(f"===============")
        print(f"chi-squared: {round(triplet_res.statistic, 5)}")
        print(f"p-value: {round(triplet_res.pvalue, 5)}")
        print(f"total triplets: {group0_and_group1_count+group0_and_group2_count+group1_and_group2_count}")
        print(f"0-1: {group0_and_group1_count}")
        print(f"0-2: {group0_and_group2_count}")
        print(f"1-2: {group1_and_group2_count}")

    def print_gene_support_freq_res(
        self,
        gene_support_freq_res,
        gene_support_freq: dict,
        trees_file_path: list
    ) -> None:
        """
        print results to stdout for user
        """
        print(f"Gene Support Frequency Results")
        print(f"==============================")
        print(f"chi-squared: {round(gene_support_freq_res.statistic, 5)}")
        print(f"p-value: {round(gene_support_freq_res.pvalue, 5)}")
        print(f"total genes: {len(trees_file_path)}")
        print(f"0-1: {gene_support_freq['0-1']}")
        print(f"0-2: {gene_support_freq['0-2']}")
        print(f"1-2: {gene_support_freq['1-2']}")

    def read_groups(
        groups: str
    ):
        """
        read groups into an array
        """

        groups_arr = []
        for line in open(groups):
            line = line.strip()
            if not line.startswith("#"):
                line = line.split("\t")
                temp = []
                temp.append(line[0])
                temp.append(line[1].split(";"))
                temp.append(line[2].split(";"))
                temp.append(line[3].split(";"))
                groups_arr.append(temp)

        groups_fin = {}

        for group in groups_arr:
            temp = []
            for i in range(1, 4):
                temp.append([taxon_name for taxon_name in group[i]])
            groups_fin[group[0]] = (temp)

        return groups_fin