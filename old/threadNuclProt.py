#!/usr/bin/env python

import sys
import getopt
import os.path
from Bio import SeqIO

def thread_dna(
    PROT, NUCL, STOP
    ):
    """
    threads DNA on top of protein alignment
    
    Parameters
    ----------
    argv: PROT
        protein fasta file
    argv: NUCL
        nucleotide fasta file
    argv: STOP
        specifies if the stop codon should included or not
    """

    # read in aligned protein fasta file
    format    = "fasta"
    handle    = open(PROT)
    PROTa     = list(SeqIO.parse(handle, format))

    # read in aligned protein fasta file
    format    = "fasta"
    handle    = open(NUCL)
    NUCLf     = list(SeqIO.parse(handle, format))

    # initialize dictionary for sequences
    pal2nal = {}

    # loop through genes and thread DNA over
    # protein alignment
    for entry in range(0, len(PROTa)):
        # save the gene ID to ID
        ID = ''
        ID = PROTa[entry].id
        # save protein sequence to Pseq
        Pseq = ''
        Pseq = PROTa[entry].seq
        # save nucleotide sequence to Nseq
        Nseq = ''
        Nseq = NUCLf[entry].seq

        pal2nal[ID] = ''
        GAPcnt = 0
        # loop through the sequence
        for AA in range(0, (int(len(Pseq))+1) - 1, 1):
            # if stops should be included
            if STOP == 'T':
                # if AA is a gap insert a codon of gaps
                if Pseq[AA] == "-":
                    pal2nal[ID] += ('---')
                    GAPcnt+=1
                # if AA is not a gap, insert the corresponding codon
                elif Pseq[AA] != "-":
                    NTwin=(AA-GAPcnt)*3
                    pal2nal[ID] += (Nseq[NTwin:NTwin+3])
            # if stops should not be included
            elif STOP == 'F':
                # if AA is a gap insert a codon of gaps
                if Pseq[AA] == "-":
                    pal2nal[ID] += ('---')
                    GAPcnt+=1
                # if AA is not a gap, insert the corresponding codon
                elif Pseq[AA] != "-":
                    # if AA is a stop or ambiguous insert a codon of gaps
                    if Pseq[AA] == 'X' or Pseq[AA] == '*':
                        pal2nal[ID] += ('---')
                    else:
                        NTwin=(AA-GAPcnt)*3
                        pal2nal[ID] += (Nseq[NTwin:NTwin+3])



                ## this commented code will check if the nucleotide window  
                ## translates to the corresponding codon
                # if Pseq[AA] != Nseq[NTwin:NTwin+3].translate():
                #     print("\nAmino acid position", AA, "(",Pseq[AA],")", "does not correspond to codon")
                #     print(Nseq[NTwin:NTwin+3], "in nucleotide window", NTwin,"-",NTwin+3)
                #     print(Nseq[NTwin:NTwin+3], "translates to", Nseq[NTwin:NTwin+3].translate())
                #     print("Nucleotides cannot be threaded ontop of the protein sequence.")
                #     print("Exiting now...\n")
                #     sys.exit()

    # print out threaded DNA alignment
    for entry in range(0, len(PROTa)):
        print(">{}\n{}".format(PROTa[entry].id, pal2nal[PROTa[entry].id]))

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # initialize argument variables
    PROT  = ''
    NUCL  = ''
    STOP  = ''

    try:
        opts, args = getopt.getopt(argv, "hp:n:s:")
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
            print("\nThis script threads a DNA codons onto a protein alignment much like")
            print("Pal2Nal or transAlign. The order of entries in the protein and nucleotide")
            print("fasta files must be the same. The output will be the DNA codon alignment")
            print("with the gene IDs that are used in the protein alignment file.")
            # protein alignment file explanation
            print("\n-p\tprotein alignment:")
            print("\tSingle or multiple protein fasta alignment. Genes should appear")
            print("\tin the same order as the genes in the nucleotide fasta file")
            print("\tspecified with the -n parameter")
            # nucleotide file explanation
            print("\n-n\tnucleotide file:")
            print("\tSingle or multiple nucleotide fasta. Genes should appear in the")
            print("\tsame order as the genes in the protein alignment fasta file specified")
            print("\twith the -p parameter")
            # include or exclude stop codon explanation
            print("\n-s\tInclude or exclude stop codon(s):")
            print("\tA boolean 'T' or 'F' for the inclusion of stop codons. If stop codons")
            print("\tshould be kept in the resulting alignment, use 'T'")
            print("\tIf stop codons should not be kept in the resulting alignment, use 'F'")
            print("\tThis option is included because numerous programs will not use alignments")
            print("\tthat contain stop codons\n")
            sys.exit()
        elif opt == '-p':
            if os.path.isfile(arg):
                PROT = arg
            else:
                # error message
                print("\n\nThe specified protein alignment fasta file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-n':
            if os.path.isfile(arg):
                NUCL = arg
            else:
                # error message
                print("\n\nThe specified nucleotide fasta file does not exist.\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()
        elif opt == '-s':
            STOP = arg
            if STOP in ('T', 'F'):
                1
            else:
                # error message
                print("\n\nThe -s argument has been incorrectly input. Input must be a boolean 'T' or 'F'\n")
                print("For detailed explanation use -h argument\n")
                sys.exit()

    # pass to thread_dna function
    thread_dna(
        PROT, NUCL, STOP
        )

if __name__ == '__main__':
    main(sys.argv[1:])
