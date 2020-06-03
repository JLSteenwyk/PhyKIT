#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
import itertools
import statistics as stat
import numpy as np

def calculate_LBscore(
    inTree
    ):
    """
    calculates LB score of a newick tree file
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3972080/pdf/ebo-10-2014-051.pdf
    
    Parameters
    ----------
    argv: inTree
        newick tree file
    """

    # read in tree
    tree = Phylo.read(inTree, 'newick')
        
    # get tip labels
    tips = []
    for tip in tree.get_terminals():
        tips.append(tip.name)
    
    # determine pairwise combinations
    combos = list(itertools.combinations(tips, 2))

    # determine average distance between tips
    # avgDist is PDa
    avgDist = float()
    for combo in combos:
        avgDist+=tree.distance(combo[0], combo[1])
    avgDist = avgDist/len(combos)
    
    # calculate average distance of taxon i to all other taxa
    # or PDi and save each result to LBi
    avgPDis = []
    for tip in tips:
        tipsMinusI = list(set(tips) - set(tip))
        PDi = []
        for tipMinus in tipsMinusI:    
            PDi.append(tree.distance(tip, tipMinus))
        PDi = sum(PDi) / len(PDi)
        avgPDis.append(PDi)
    
    # use PDis and avgDist to calculate LB values for each taxon
    LBis = []
    for PDi in avgPDis:
        LBis.append((((PDi/avgDist)-1)*100))

    with open(inTree + ".LBi-scores", 'w') as f:
        for tip, LBi in zip(tips, LBis):
            f.write(str(tip) + "\t" + str(LBi) + "\n")

    mean          = stat.mean(LBis)
    median        = stat.median(LBis)
    twenty_fifth  = np.percentile(LBis, 25)
    seventy_fifth = np.percentile(LBis, 75)

    print("{}\t{}\t{}\t{}".format(mean, median, twenty_fifth, seventy_fifth))



def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    tree  = ''

    try:
        opts, args = getopt.getopt(argv, "hi:")
    except getopt.GetoptError:
        # error message
        print("Error\nFor help use -h argument\n")
        sys.exit(2)
    # if no arguments are used print a help message
    if len(opts) == 0:
        # error message
        print("\nNo arguments provided...")
        print("For help use -h argument\n")
        sys.exit(2)
    
    # test for arguments
    for opt, arg in opts:
        if opt == '-h':
            # general script explanation
            print("\nusage: python LB_score.py -i newick.tre")
            print("\nCalculates LB score of a phylogenetic tree in newick format")
            print("Stdout is the mean, median, twenty fith and seventy fifth percentiles of LB scores")
            print("across all taxa in a phylogeny. Additional output includes a file with the suffix")
            print("'.LBi-scores' which, summarizes LB scores for each taxa individually.")
            print("LB score is calcualted according to Struck. TreSpEx—Detection of Misleading Signal ")
            print("in Phylogenetic Reconstructions Based on Tree Information. Evolutionary Bioinformatics 2014:10")
            print("51–67 doi: 10.4137/EBO.S14239.\n")
            print("Author: Jacob L. Steenwyk")
            print("Citation: NA\n")
            sys.exit()
        elif opt == '-i':
            if os.path.isfile(arg):
                inTree = arg
            else:
                # error message
                print("\n\nThe specified tree file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to calculate_LBscore function
    calculate_LBscore(
        inTree
        )

if __name__ == '__main__':
    main(sys.argv[1:])