#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
import math

def calculate_DVMC(
	tree, numSpp
	):
    """ 
    calculates DVMC. DVMC is calculated using the
    following formula: The square root of (1/N-1) times
    the sum of i=1 to the N of (x_i - x_bar)^2 where x_i
    is the distance between the root of the tree and 
    species i and N is the number of species

    Parameters
    ----------
    argv: tree
        newick tree file
    argv: numSpp
        number of species in the tree
    """
    
    # initialize list of root to tip distances, the sum
    # of distances in the tree, the average of distances
    # in the tree
    distL = []
    sumDist = float(0.0)
    avgDist = float(0.0)

    # loop through terminal branches and store distances from
    # the root to the tip in a list
    for term in tree.get_terminals():
        # append root to tip distance to distL
        distL.append((tree.distance(term)))
        # keep running sum of tree distances
        sumDist += tree.distance(term)

    # calculate average tree distance
    avgDist = sumDist/numSpp

    # determine the sum of i=1 to N for (x_i-x_bar)^2
    sumi2N = float(0.0)
    for x_i in distL:
    	sumi2N += ((x_i-avgDist)**2)

    # multiple sumi2N by 1/(N-1) where N is the number of spp
    # and take the square root
    DVMC = float(0.0)
    DVMC = math.sqrt((1/(numSpp-1))*sumi2N) 

    # print result
    print(DVMC)




def root_and_prune(
    tree, outgroup
    ):
    """
    roots tree and prunes the outgroup taxa from
    the phylogenetic tree in newick format. Lastly,
    after pruning, the number of taxa present will
    be determined.
    
    Parameters
    ----------
    argv: tree
        newick tree file
    argv: outgroup
        outgroup taxa in tree
    """

    # read in tree
    tree = Phylo.read(tree, 'newick')

    # read in outgroup taxa
    cladeOutL = [line.rstrip('\n') for line in open(outgroup)]
    # initialize list for outgroup taxa present in the tree
    outPres   = []

    # initialize value to hold the number of species
    numSpp = float(0.0)

    # loop through terminal branch
    for term in tree.get_terminals():
        # if outgroup taxa is present, append it to outPres
        if term.name in cladeOutL:
            outPres.append(term.name)

    # root tree on outgroup
    tree.root_with_outgroup(outPres)

    # prune outgroiup taxa from the tree
    for taxon in outPres:
        tree.prune(taxon)

    # determine number of taxa in the tree
    for term in tree.get_terminals():
        # add 1 for each tip that is iterated through
        numSpp+=1

    #Phylo.draw_ascii(tree)

    # pass to calculate DVMC function
    calculate_DVMC(
        tree, numSpp
        )

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    tree  = ''

    try:
        opts, args = getopt.getopt(argv, "ht:o:")
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
            print("\nCalculates the degree of violation of a molecular clock (DVMC) for a tree.")
            print("The tree is first rooted on the outgroup taxa and then those outgroup taxa are")
            print("pruned from the tree. After that, DVMC is calculated. The outgroup taxa should be")
            print("the taxa that are not used in the timetree but present in the phylogenetic tree.")
            print("This is calculated according to the following paper: http://www.pnas.org/content/114/35/E7282")
            sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                tree = arg
            else:
                # error message
                print("\n\nThe specified tree file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-o':
            if os.path.isfile(arg):
                outgroup = arg
            else:
                # error message
                print("\n\nThe specified outgroup file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to calculate_treeness function
    root_and_prune(
        tree, outgroup
        )

if __name__ == '__main__':
    main(sys.argv[1:])
