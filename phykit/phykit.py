#!/usr/bin/env python

import getopt
import logging
import os.path
import sys
import textwrap
import time

from argparse import (
    ArgumentParser,
    RawTextHelpFormatter,
    SUPPRESS,
    RawDescriptionHelpFormatter,
)
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)


def execute(
    inFile: str,
    outFile: str,
    inFileFormat: FileFormat,
    outFileFormat: FileFormat,
    gaps: float,
    complement: bool,
    mode: TrimmingMode,
    use_log: bool,
):
    """
    Master execute Function                                      
    This function executes the main functions and calls other    
    subfunctions to trim the input file  
    """
    # logic for whether or not to create a log file
    if use_log:
        # write INFO level logging to file for user
        logger.setLevel(logging.DEBUG)
        fh = logging.FileHandler(f"{outFile}.log", mode="w")
        fh.setLevel(logging.DEBUG)
        logger.addHandler(fh)

    # create start time logger
    start_time = time.time()


def main(argv=None):
    """
    Function that parses and collects arguments              
    """

    # ASCII text art was obtained from the following website
    # http://patorjk.com/software/taag/#p=display&f=Big&t=PhyKIT
    parser = ArgumentParser(
        add_help=False,
        formatter_class=RawDescriptionHelpFormatter,
        usage=SUPPRESS,
        description=textwrap.dedent(
            """\
             _____  _           _  _______ _______ 
            |  __ \| |         | |/ /_   _|__   __|
            | |__) | |__  _   _| ' /  | |    | |   
            |  ___/| '_ \| | | |  <   | |    | |   
            | |    | | | | |_| | . \ _| |_   | |   
            |_|    |_| |_|\__, |_|\_\_____|  |_|   
                           __/ |                   
                          |___/   
            Citation: Steenwyk et al. Journal, journal info, link

            PhyKIT helps with all things phylogenetics and phylogenomics.

            Usage: phykit <input> [optional arguments]
            """
        ),
    )

    # if no arguments are given, print help and exit
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()

    ## required arguments
    required = parser.add_argument_group(
        "required arguments",
        description=textwrap.dedent(
            """\
        <function>                                  phykit function 
                                                    (must be the first argument)
        """
        ),
    )

    required.add_argument("input", type=str, help=SUPPRESS)

    ## optional arguments
    optional = parser.add_argument_group(
        "optional arguments",
        description=textwrap.dedent(
            """\
        -t, --tree <output_file_name>               output file name 
                                                    (default: input file named with '.clipkit' suffix)

        --------------------------------------------------------------
        | Detailed explanation of functions (alphabetically ordered) | 
        --------------------------------------------------------------
        List of functions (alphabetically ordered)
        
        dvmc (degree of violation of molecular clock)
            - calculates blah
            - See Liu et al (2017), PNAS; https://www.pnas.org/content/114/35/E7282

        rcv (relative composition variability)
            

        treeness
        """
        ),
    )

    # TODO: communicate to Thomas that optional arguments will largely depend on the 
    # function being called. E.g., dvmc and treeness take trees as input; in contrast, 
    # relative composition variability takes a fasta file. Maybe we should refactor
    # the argparsing based on the function being called
    optional.add_argument(
        "-t",
        "--tree",
        help=SUPPRESS
    )

    # parse and assign arguments
    args = parser.parse_args()

    # pass to master execute function
    execute(*process_args(args))


if __name__ == "__main__":
    main(sys.argv[1:])
