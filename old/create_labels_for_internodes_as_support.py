#!/usr/bin/env python

import sys
import getopt
import os.path
#from Bio.Phylo.BaseTree import TreeMixin
from Bio.Phylo.BaseTree import Clade
from Bio import Phylo

def modTree(
    tree
    ):
    """
    creates numerical identifiers for internodes
    
    Parameters
    ----------
    argv: tree
        newick tree file
    """

    # create file name 
    filename=tree+".internodeLabels.tree"

    # internode label starter
    label=1

    # read in tree
    tree = Phylo.read(tree, 'newick')
    for i in tree.get_nonterminals():
        i.confidence=label
        label+=1

    Phylo.write(tree, filename, 'newick')

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    tree  = ''

    try:
        opts, args = getopt.getopt(argv, "ht:")
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
            print("\nCreates labels for internodes in a tree. Non-terminal branches will be looped through")
            print("Using .get_nonterminals in BioPython and each internode will get a numerical label")
            print("starting from 1 to n where n is the number of internodes in the tree. Each internode")
            print("label is reported as a confidence score that can be easily viewed in FigTree. The only")
            print("option is -t for a newick tree file.")
            sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                tree = arg
            else:
                # error message
                print("\n\nThe specified tree file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to calculate_treeness function
    modTree(
        tree
        )

if __name__ == '__main__':
    main(sys.argv[1:])