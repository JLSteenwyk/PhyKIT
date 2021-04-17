import sys
import getopt
import os.path
import itertools
from scipy.stats import chisquare
from typing import Tuple

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
        groups_of_groups, outgroup_taxa = self.determine_groups_of_groups(groups_arr)

        # read trees into list
        trees_file_path = read_single_column_file_to_list(self.trees)

        # go through all triplets of all trees and  
        # examine sister relationships among all triplets
        summary = self.loop_through_trees_and_examine_sister_support_among_triplets(
            trees_file_path,
            groups_of_groups,
            outgroup_taxa
        )

        # count triplet and gene support frequencies for different sister relationships
        triplet_group_counts, gene_support_freq = self.get_triplet_and_gene_support_freq_counts(summary)

        # conduct chisquare tests
        triplet_res, gene_support_freq_res = self.chisquare_tests(triplet_group_counts, gene_support_freq)

        # print results
        self.print_gene_support_freq_res(gene_support_freq_res, gene_support_freq, trees_file_path)
        # self.print_triplet_based_res(triplet_res, triplet_group_counts)
    
    def process_args(self, args):
        return dict(trees=args.trees, groups=args.groups)

    def read_in_groups(self, groups) -> list:
        groups_arr = []
        try:
            for line in open(self.groups):
                line = line.strip()
                if not line.startswith("#"):
                    try:
                        line = line.split("\t")
                        temp = []
                        temp.append(line[0])
                        temp.append(line[1].split(";"))
                        temp.append(line[2].split(";"))
                        temp.append(line[3].split(";"))
                        temp.append(line[4].split(";"))
                        groups_arr.append(temp)
                    except IndexError:
                        try:
                            print(f"{self.groups} contains an indexing error.")
                            print("Please format the groups file (-g) as a four column tab-delimited file with column 1 being the name of the test")
                            print("col2: the tip names of one group (; separated)")
                            print("col3: the tip names of a second group (; separated)")
                            print("col4: the tip names of a third group (; separated)")
                            print("col5: the tip names of the outgroup taxa (; separated)")
                            sys.exit()
                        except BrokenPipeError:
                            pass
                        
        except FileNotFoundError:
            try:
                print(f"{self.groups} corresponds to no such file.")
                print("Please check filename and pathing again.")
                sys.exit()
            except BrokenPipeError:
                pass
        return groups_arr

    def loop_through_trees_and_examine_sister_support_among_triplets(
        self,
        trees_file_path: str,
        groups_of_groups: dict,
        outgroup_taxa: list
    ) -> dict:
        """
        go through all trees and all triplets of all trees. For each triplet,
        determine which two taxa are sister to one another
        """
        summary = {}
        # loop through trees
        try:
            #cnt = 0
            for tree_file in trees_file_path:
                #cnt+=1
                #print(f"processing tree {cnt} of {len(trees_file_path)}")
                tree = Phylo.read(tree_file, 'newick')
                # get tip names
                tips = self.get_tip_names_from_tree(tree)

                # examine all triplets and their support for 
                # any sister pairing
                summary = self.examine_all_triplets_and_sister_pairing(
                    tips,
                    tree_file,
                    summary,
                    groups_of_groups,
                    outgroup_taxa
                )
        except FileNotFoundError:
            try:
                print(f"{tree_file} corresponds to no such file.")
                print("Please check file name and pathing")
                sys.exit()
            except BrokenPipeError:
                pass
        
        return summary

    def determine_groups_of_groups(
        self,
        groups_arr: list
    ) -> dict:
        groups_of_groups = {}

        for group in groups_arr:
            temp = []
            for i in range(1, 4):
                temp.append([taxon_name for taxon_name in group[i]])
            groups_of_groups[group[0]] = (temp)

        outgroup_taxa = [taxon_name for taxon_name in group[4]]
        
        return groups_of_groups, outgroup_taxa

    def examine_all_triplets_and_sister_pairing(
        self,
        tips: list,
        tree_file: str,
        summary: dict,
        groups_of_groups: dict,
        outgroup_taxa: list
    ) -> dict:
        """
        evaluate all triplets for sister relationships. Polytomies
        in input trees are accounted for
        """
        # get all combinations of three tips
        identifier = list(groups_of_groups.keys())[0]
        triplet_tips = (list(itertools.product(*groups_of_groups[identifier])))
        for triplet in triplet_tips:
            # obtain tree of the triplet
            tree = self.get_triplet_tree(tips, triplet, tree_file, outgroup_taxa)
            cnt=0
            try:
                for leaf in tree.get_terminals():
                    cnt+=1
            except AttributeError:
                continue
            if tree and cnt==3:
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
                        # determine sisters and add to sister pair counter
                        summary = self.determine_sisters_and_add_to_counter(
                            tip_names,
                            tree,
                            tree_file,
                            groups,
                            summary
                        )
            else:
                continue

        return summary

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

    def set_branch_lengths_in_tree_to_one(self, tree: Tree) -> None:
        for term in tree.get_terminals():
            term.branch_length = 1
        for internode in tree.get_nonterminals():
            internode.branch_length = 1

    def check_if_triplet_is_a_polytomy(self, tree: Tree) -> Tree:
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

    def sister_relationship_counter(
        self,
        tree_file: str,
        summary: dict,
        sisters: str
    ) -> dict:
        """
        counter for how many times a particular sister relationship is observed
        """
        # if tree is not in summary, create a key for it
        if tree_file not in summary.keys():
            summary[str(tree_file)]={}
        # if the sister relationship is not in the tree file dict, create a key for it
        if sisters not in summary[str(tree_file)].keys():
            summary[str(tree_file)][sisters] = 1
        else:
            summary[str(tree_file)][sisters] += 1

        return summary

    def get_triplet_tree(
        self,
        tips: list,
        triplet: tuple,
        tree_file: str, 
        outgroup_taxa: list
    ) -> Tree:
        """
        get a tree object of only the triplet of interest
        """
        # determine tips that are not in the triplet of interest
        tips_to_prune = list(set(tips)- set(list(triplet)))
        # determine tips that in the outgroup
        outgroup_present = [value for value in tips if value in outgroup_taxa]
        tree = Phylo.read(tree_file, 'newick')
        
        # root tree on outgroup taxa
        try:
            tree.root_with_outgroup(outgroup_present)

            # prune to a triplet
            tree = self.prune_tree_using_taxa_list(tree, tips_to_prune)

            return tree
        except ValueError:
            tree = False
            return tree

    def determine_sisters_from_triplet(
        self,
        groups: list,
        pair: tuple
    ) -> str:
        """
        determine sister taxa from a triplet
        """
        sisters = sorted([[i for i, lst in enumerate(groups) if pair[0] in lst][0], [i for i, lst in enumerate(groups) if pair[1] in lst][0]])
        sisters = [str(i) for i in sisters] 
        sisters = str("-".join(sisters))

        return sisters

    def determine_sisters_and_add_to_counter(
        self,
        tip_names: list,
        tree: Tree,
        tree_file: str,
        groups: list,
        summary: dict
    ) -> dict:
        """
        determine which pair of taxa are sister to one another
        and add 1 to the counter for the sister pair
        """
        # get pairs from tip names
        pairs = list(itertools.combinations(tip_names, 2))
        for pair in pairs:
            is_polytomy = self.check_if_triplet_is_a_polytomy(tree)
            # if distance between pair is 2 and the triplet is 
            # not a polytomy (i.e., having only 1 internal branch)
            # then report the sisters in the triplet
            if tree.distance(pair[0], pair[1]) == 2 and not is_polytomy:
                # determine which two tips are sisters
                sisters = self.determine_sisters_from_triplet(groups, pair)
                # add to summary dictionary of how many times that sister
                # relationship is observed
                summary = self.sister_relationship_counter(tree_file, summary, sisters)
        return summary

    def get_triplet_and_gene_support_freq_counts(
        self,
        summary: dict
    ) -> Tuple[dict, dict]:
        """
        count how many triplets and genes support the various sister relationships
        """
        # Count the total number of sister pairings
        # for the three possible pairs for triplets
        triplet_group_counts = {
            'g0g1_count' : 0,
            'g0g2_count' : 0, 
            'g1g2_count' : 0
        }

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
            triplet_group_counts['g0g1_count'] += summary[tree]['0-1']
            triplet_group_counts['g0g2_count'] += summary[tree]['0-2']
            triplet_group_counts['g1g2_count'] += summary[tree]['1-2']
            # determine which sister pairing is best supported in a single gene
            # and add one to the corresponding gene support frequency count
            gene_support_freq[max(summary[tree], key=summary[tree].get)] += 1

        return triplet_group_counts, gene_support_freq

    def chisquare_tests(
        self,
        triplet_group_counts: dict,
        gene_support_freq: dict
    ): 
        
        triplet_res = chisquare(
            [
                triplet_group_counts['g0g1_count'],
                triplet_group_counts['g0g2_count'],
                triplet_group_counts['g1g2_count']
            ]
        )
        
        gene_support_freq_res = chisquare(
            [
                gene_support_freq['0-1'],
                gene_support_freq['0-2'],gene_support_freq['1-2']
            ]
        )

        return triplet_res, gene_support_freq_res

    # def print_triplet_based_res(
    #     self,
    #     triplet_res,
    #     triplet_group_counts: dict
    # ) -> None:
    #     """
    #     print results to stdout for user
    #     """
    #     try:
    #         print(f"\nTriplet Results")
    #         print(f"===============")
    #         print(f"chi-squared: {round(triplet_res.statistic, 4)}")
    #         print(f"p-value: {round(triplet_res.pvalue, 6)}")
    #         print(f"total triplets: {sum(triplet_group_counts.values())}")
    #         print(f"0-1: {triplet_group_counts['g0g1_count']}")
    #         print(f"0-2: {triplet_group_counts['g0g2_count']}")
    #         print(f"1-2: {triplet_group_counts['g1g2_count']}")
    #     except BrokenPipeError:
    #         pass

    def print_gene_support_freq_res(
        self,
        gene_support_freq_res,
        gene_support_freq: dict,
        trees_file_path: list
    ) -> None:
        """
        print results to stdout for user
        """
        try:
            print(f"Gene Support Frequency Results")
            print(f"==============================")
            print(f"chi-squared: {round(gene_support_freq_res.statistic, 4)}")
            print(f"p-value: {round(gene_support_freq_res.pvalue, 6)}")
            print(f"total genes: {(gene_support_freq['0-1'] + gene_support_freq['0-2'] + gene_support_freq['1-2'])}")
            print(f"0-1: {gene_support_freq['0-1']}")
            print(f"0-2: {gene_support_freq['0-2']}")
            print(f"1-2: {gene_support_freq['1-2']}")
        except BrokenPipeError:
            pass
