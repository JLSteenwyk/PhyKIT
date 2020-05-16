#!/usr/bin/env python

import logging
import sys
import textwrap
from argparse import (
    ArgumentParser,
    RawTextHelpFormatter,
    SUPPRESS,
    RawDescriptionHelpFormatter,
)

from .services.dna_threader import DNAThreader
from .services.tree.treeness import Treeness
from .services.tree.internode_labeler import InternodeLabeler

logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)


class Phykit(object):
    def __init__(self):
        parser = ArgumentParser(
            add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
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

                Usage: phykit <command> [optional command arguments]
                """
            ),
        )
        parser.add_argument("command")
        args = parser.parse_args(sys.argv[1:2])

        if not hasattr(self, args.command):
            print("Invalid command option. See help for more details.")
            parser.print_help()
            sys.exit(1)

        getattr(self, args.command)()

    def dvmc(self):
        raise NotImplementedError()

    def rcv(self):
        raise NotImplementedError()

    def treeness(self):
        parser = ArgumentParser()
        parser.add_argument("tree", type=str)
        args = parser.parse_args(sys.argv[2:])
        Treeness(args).run()

    def internode_labeler(self):
        parser = ArgumentParser()
        parser.add_argument("tree", type=str)
        # TODO: add output option?
        args = parser.parse_args(sys.argv[2:])
        InternodeLabeler(args).run()

    def thread_dna(self):
        parser = ArgumentParser()
        parser.add_argument(
            "-p",
            type=str,
            help="Single or multiple protein fasta alignment. Genes should appear in the same order as the genes in the nucleotide alignment fasta file specified with the -n parameter",
        )
        parser.add_argument(
            "-n",
            type=str,
            help="Single or multiple nucleotide fasta. Genes should appear in the same order as the genes in the protein alignment fasta file specified with the -p parameter",
        )
        parser.add_argument(
            "-s",
            type=bool,
            help="Should stop codons be kept in the resulting alignments",
        )
        args = parser.parse_args(sys.argv[2:])
        DNAThreader(args).run()


if __name__ == "__main__":
    Phykit()
