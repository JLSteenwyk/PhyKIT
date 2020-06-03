#!/usr/bin/env python

import sys
import getopt
import os.path
#from Bio.Phylo.BaseTree import TreeMixin
from Bio import Phylo

def calculate_tree_length(
    tree
    ):
    """
    calculates tree length of a newick tree file
    
    Parameters
    ----------
    argv: tree
        newick tree file
    """

    # read in tree
    tree = Phylo.read(tree, 'newick')

    # initialize variable for total branch length
    totalLen = float(0.0)

    # determine total branch length
    totalLen = tree.total_branch_length()

    # print total tree length
    print('{:.15f}'.format(totalLen))


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
            print("\nCalculates tree length of a phylogenetic tree in newick format")
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
    calculate_tree_length(
        tree
        )

if __name__ == '__main__':
    main(sys.argv[1:])