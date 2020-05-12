#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import SeqIO
from Bio.Seq import Seq

def RCV(
    alignment
    ):
    """
    determines if input fasta is nucleotide or protein

    Parameters
    ----------
    argv: alignment
        alignment fasta file
    """

    # string to hold all sequences
    concatSeq = ''
    # initialize a counter for the number of sequences in the input fasta file
    numRecords = 0

    # read in records of a fasta file
    records = list(SeqIO.parse(alignment, "fasta"))
    
    # for each record join concatSeq string and sequence as well as keeping track 
    # of the number of records
    for record in records:
        concatSeq  += record.seq
        numRecords += 1

    # dictionary to hold the average occurence of each sequence letter
    averageD = {}
    # loop through the different sequences that appear in the fasta file
    # population dictionary with how many times that sequence appears
    for seq in set(concatSeq):
    	averageD[seq] = (concatSeq.count(seq)/numRecords)

    # determine the length of the alignment
    alignmentLen = 0
    alignmentLen = len(records[0].seq)

    # intiailize list to hold the RCV values per ith taxa 
    # that will later be summed
    indivRCVvalues = []

    # loop through records again and calculate RCV for 
    # each taxa and append to indivRCVvalues
    for record in records:
        # temp holds a temporary value of the numerator before appending
        # to numeratorRCVvalues and then is reassigned to 0 when it goes
        # through the loop again
        temp = 0
        # calculates the absolute value of the ith sequence letter minus the average
        for seqLetter in set(concatSeq):
            temp += abs(record.seq.count(seqLetter)-averageD[seqLetter])
        indivRCVvalues.append(temp/(numRecords*alignmentLen))

    # print the sum of all RCV values
    print(sum(indivRCVvalues))


def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    alignment = ''

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
            ## explanation
            print("\nCalculate relative nucleotide composition variability (RCV) as defined by Phillips and Penny")
            print("2003, specifically in the manuscript titled 'The root of the mammalian tree inferred from")
            print("whole mitochondrial genomes. Mol Phylogenet Evol. 28:171–185.'")
            print("https://www.ncbi.nlm.nih.gov/pubmed/12878457\n")
            print("The formula for RCV is as follows:")
            print("the sum of i=1^n of (|Ai - A*| + |Ti - T*| + |Ci - C*| + |Gi - G*|)/n•t")
            print("where, Ai, Ti, Ci, and Gi are the number of nucleotides for the ith taxon,")
            print("A*, T*, C*, and G* are the averages of each nucleotide across the n taxa")
            print("t is the number of sites and n is the number of taxa.")
            print("\nThe output will a number that is the RCV value for a given alignment.")
            ## options
            # alignment files list
            print("\n-i <alignment fasta file>")
            print("\tA multi-fasta sequence alignment file")
            sys.exit()

        elif opt == '-i':
            if os.path.isfile(arg):
                alignment = arg
            else:
                # error message
                print("\n\nThe specified alignment list (-i) file does not exist.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()

    # pass to read_config parameter
    RCV(
        alignment
        )

if __name__ == '__main__':
    main(sys.argv[1:])