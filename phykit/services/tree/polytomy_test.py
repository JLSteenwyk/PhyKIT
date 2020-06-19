import sys
import getopt
import os.path
import itertools
from scipy.stats import chisquare

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin


from .base import Tree

class PolytomyTest(Tree):
    def __init__(self, args) -> None:
        super().__init__(**self.process_args(args))

    def run(self):
        # read in groups
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

        groups_of_groups = {}

        for group in groups_arr:
            temp = []
            for i in range(1, 4):
                temp.append([taxon_name for taxon_name in group[i]])
            groups_of_groups[group[0]] = (temp)
        
        # read trees into list
        trees_file_path = [line.rstrip('\n') for line in open(self.trees)]

        # loop through trees
        summary = {}
        for tree_file in trees_file_path:
            tips = []
            tree = Phylo.read(tree_file, 'newick')
            for term in tree.get_terminals():
                term.branch_length = None
                tips.append(term.name)
            for internode in tree.get_nonterminals():
                internode.branch_length = None

            # get all combinations of three tips
            triplet_tips = itertools.combinations(tips, 3)
            for triplet in triplet_tips:
                # determine tips that are not in the triplet of interest
                tips_to_prune = list(set(tips)- set(list(triplet)))
                tree = Phylo.read(tree_file, 'newick')
                
                # prune to a triplet
                for tip in tips_to_prune:
                    tree.prune(tip)
                
                for test, groups in groups_of_groups.items():
                    # see if there any intersctions between
                    # the triplet and the the group
                    num_groups_represented = 0
                    for group in groups:
                        temp = [value for value in list(triplet) if value in group]
                        if len(temp) >= 1:
                            num_groups_represented +=1
                    
                    # if one taxa is represented from each group, 
                    # use the triplet
                    tip_names = []
                    if num_groups_represented == 3:
                        for term in tree.get_terminals():
                            term.branch_length = 1
                            tip_names.append(term.name)
                        for internode in tree.get_nonterminals():
                            internode.branch_length = 1
                        pairs = list(itertools.combinations(tip_names, 2))
                        for pair in pairs:
                            num_int=0
                            # check if the triplet is a polytomy
                            for internal_branch in tree.get_nonterminals():
                                num_int+=1
                            # if distance between pair is 2 and the triplet is 
                            # not a polytomy (i.e., having only 1 internal branch)
                            # then report the sisters in the triplet
                            if tree.distance(pair[0], pair[1]) == 2 and num_int != 1:
                                sisters = sorted([[i for i, lst in enumerate(groups) if pair[0] in lst][0], [i for i, lst in enumerate(groups) if pair[1] in lst][0]])
                                sisters = [str(i) for i in sisters] 
                                sisters = str("-".join(sisters))
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

        print(f"Gene Support Frequency Results")
        print(f"==============================")
        print(f"chi-squared: {round(gene_support_freq_res.statistic, 5)}")
        print(f"p-value: {round(gene_support_freq_res.pvalue, 5)}")
        print(f"total genes: {len(trees_file_path)}")
        print(f"0-1: {gene_support_freq['0-1']}")
        print(f"0-2: {gene_support_freq['0-2']}")
        print(f"1-2: {gene_support_freq['1-2']}")

        print(f"\nTriplet Results")
        print(f"===============")
        print(f"chi-squared: {round(triplet_res.statistic, 5)}")
        print(f"p-value: {round(triplet_res.pvalue, 5)}")
        print(f"total triplets: {group0_and_group1_count+group0_and_group2_count+group1_and_group2_count}")
        print(f"0-1: {group0_and_group1_count}")
        print(f"0-2: {group0_and_group2_count}")
        print(f"1-2: {group1_and_group2_count}")
    
    
    def process_args(self, args):
        return dict(trees=args.trees, groups=args.groups)

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