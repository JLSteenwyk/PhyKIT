#!/usr/bin/env python

import sys
import getopt
import os.path
#import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from datetime import datetime

def create_concat(
    alignment_list, taxa_list, 
    prefix, sequence_type, startTime
    ):
    """
    Reads the alignment list and taxa list
    
    Parameters
    ----------
    argv: alignment_list
        single column file of busco output dirs
    argv: taxa_list
        taxon occupancy value between 0 and 1
    argv: prefix
        prefix of the output files
    argv: sequence_type
        sequence type of alignment fastas
    argv: startTime
        timing how long the script runs for
    """

    # read lists into python lists
    alignments  = [line.rstrip('\n') for line in open(alignment_list)]
    taxa        = [line.rstrip('\n') for line in open(taxa_list)]

    # dictionary to populate concat.fa with,
    # and assignment the output file names
    concatD        = {}
    file_partition = prefix + ".partition"
    fastaOUT       = prefix + ".fa"
    file_occupancy = prefix + ".occupancy"

    # print log message
    print("\n"+"-"*(len("- General features -")))
    print("| General features |")
    print("-"*(len("- General features -")))
    print("Total number of taxa:", len(taxa))
    print("Total number of alignments:", len(alignments))
    print("\n\n"+"-"*(len("- Output files -")))
    print("| Output files |")
    print("-"*(len("- Output files -")))
    print("partition file output:", file_partition)
    print("concatenated fasta output:", fastaOUT)
    print("occupancy report:", file_occupancy)

    # assigning placeholders for lengths
    firstLen    = 1
    secondLen   = 0

    # create a dictionary with keys as taxa and empty list values
    concatD = {key: [] for key in taxa}

    # if the sequence is protein or nucl then assign first column as AUTO or DNA
    # this is used when creating the partition file
    if sequence_type == 'prot':
        field1 = 'AUTO'
    elif sequence_type == 'nucl':
        field1 = 'DNA'

    # create an empty file_partition
    open(file_partition, "w")
    # create an empty occupancy per OG file
    open(file_occupancy, "w")
    # loop through USCO trimal files
    for fasta in alignments:
        # list for keeping track of taxa per USCO
        OGtaxa    = []
        # open alignment file
        records = list(SeqIO.parse(fasta, "fasta"))

        # create lists of USCOtaxa and missing_taxa
        for record in records:
            # append taxa per USCO to USCOtaxa
            OGtaxa.append(record.id)
            
        # initialize list and determine missing taxa
        missing_taxa = []
        missing_taxa = list(set(taxa)-set(OGtaxa))

        # create string for missing seq
        OGlen     = ''
        missing_seq = ''
        OGlen   = len(records[0].seq)
        missingSeq = OGlen*'?'

        # create missing taxa sequence based on either it is protein or nucleotide
        if sequence_type == 'prot':
            # create record for missing taxa
            for indiv in missing_taxa:
                indivRecord = SeqRecord(Seq(missingSeq,IUPAC.protein), 
                    id = indiv, name = indiv, description = indiv)
                concatD[indiv].append(indivRecord.seq)
        elif sequence_type == 'nucl':
            # create record for missing taxa
            for indiv in missing_taxa:
                indivRecord = SeqRecord(Seq(missingSeq,IUPAC.ambiguous_dna), 
                    id = indiv, name = indiv, description = indiv)
                concatD[indiv].append(indivRecord.seq)        

        # create record for present data
        for record in records:
            if record.id in OGtaxa:
                concatD[record.id].append(record.seq)

        # create partition file
        with open(file_partition, "a") as f:
            # second value in partition file
            secondLen += OGlen
            entry = field1+", "+str(fasta)+"="+str(firstLen)+"-"+str(secondLen)+"\n"
            f.write(str(entry))
            # add to first value for partition file
            firstLen += OGlen

        # create partition file
        with open(file_occupancy, "a") as f:
            # determine number of present taxa
            NumPresent  = len(OGtaxa)
            # determine number of missing taxa
            NumMissing  = len(missing_taxa)
            # determine percent occupied
            PercentOccu = NumPresent/(NumPresent+NumMissing)
            entry = str(fasta)+'\t'+str(NumPresent)+'\t'+str(NumMissing)+'\t'+str(PercentOccu)+"\n"
            f.write(str(entry))

    # create concatenated fasta file based on either it is protein or nucleotide
    if sequence_type == 'prot':
        # join seqs of genes in value (list) and write to concat_fa   
        with open(fastaOUT, "w") as final_fasta_file: 
            for x in concatD:
                concatenated = Seq("", IUPAC.protein)
                for s in concatD[x]:
                    concatenated += s
                concatD[x] = concatenated
                entry = '>'+x+"\n"+concatD[x]+"\n"
                final_fasta_file.write(str(entry))
    elif sequence_type == 'nucl':
        # join seqs of genes in value (list) and write to concat_fa   
        with open(fastaOUT, "w") as final_fasta_file: 
            for x in concatD:
                concatenated = Seq("", IUPAC.ambiguous_dna)
                for s in concatD[x]:
                    concatenated += s
                concatD[x] = concatenated
                entry = '>'+x+"\n"+concatD[x]+"\n"
                final_fasta_file.write(str(entry))                

    # print log message
    print("Complete!\n")
    print("Total time:", datetime.now() - startTime)

def main(
    argv
    ):
    """
    Reads arguments 
    """

    # start time of script
    startTime = datetime.now()

    # initialize argument variables
    alignment_list  = ''
    taxa_list       = ''
    prefix          = ''

    try:
        opts, args = getopt.getopt(argv, "ha:t:p:c:m:")
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
            print("\nCreate a concatenation and partition file for the concatenation")
            print("and gene-partitioned based phylogenetic/phylgenomic inference.")
            print("Three output files are made.")
            print("(1) A fasta file with '.fa' appended to the end of the prefix. This")
            print("is the concatenated data matrix.")
            print("(2) A partition file ready for input into RAxML or IQ-tree.")
            print("(3) An occupancy file that summarizes the taxon occupancy per sequence")
            print("file. The output is formatted with col1: the file name, col2: the number")
            print("of taxa present, col3: the number of taxa missing, and col4: the percent")
            print("of occupancy for that particular sequence file.")
            ## options
            # alignment files list
            print("\n-a <list of alignment files>")
            print("\tSingle column file of the alignment files that will be concatenated.")
            # specify if nucleotide or protein files list
            print("\n-c <alignments are either nucleotide or protein fastas>")
            print("\tArgument can be either 'prot' or 'nucl' to specify if the")
            print("\talignments contain protein or nucleotide sequences.")
            # taxa name list
            print("\n-t <list of taxa names>")
            print("\tSingle column file of taxa names across alignment files. Any taxa")
            print("\tnot included will be ignored by the script. Name of taxa specified")
            print("\tin this file must match the alignment file with no extra characters.")
            # output file name
            print("\n-p <prefix of output files>")
            print("\tPrefix of the output files. For example, if the output prefix is")
            print("\tspecified to be 'concat' a file named 'concat.fa' (the concatenated")
            print("\tmatrix) and 'concat.partition' (the associated partition file) will be made.")
            sys.exit()

        elif opt == '-a':
            if os.path.isfile(arg):
                alignment_list = arg
            else:
                # error message
                print("\n\nThe specified alignment list (-a) file does not exist.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()
        elif opt == '-t':
            if os.path.isfile(arg):
                taxa_list = arg 
            else:
                # error message
                print("\n\nThe specified taxon list (-t) does not exist.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()
        elif opt == '-p':
            if arg:
                prefix = arg
            else:
                # error message
                print("\n\nPlease specify an output file prefix (-p).\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()
        elif opt == '-c':
            if arg in ('prot', 'nucl'):
                sequence_type = arg 
            else:
                # error message
                print("\n\nThe specified input for (-c) should be 'prot' or 'nucl'.\n")
                print("For detailed explanation of configuration file use -h argument\n")
                sys.exit()

    # print log output
    print("\n"+"-"*(len("- Argument parameters -")))
    print("| Argument parameters |")
    print("-"*(len("- Argument parameters -")))
    print("-a", alignment_list)
    print("-c", sequence_type)
    print("-t", taxa_list)
    print("-p", prefix+"\n")

    # pass to read_config parameter
    create_concat(
        alignment_list, taxa_list, 
        prefix, sequence_type, startTime
        )

if __name__ == '__main__':
    main(sys.argv[1:])