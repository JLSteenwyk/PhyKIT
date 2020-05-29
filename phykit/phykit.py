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

from .services.alignment.alignment_length import AlignmentLength
from .services.alignment.variable_sites import VariableSites

from .services.tree.treeness import Treeness
from .services.tree.total_tree_length import TotalTreeLength
from .services.tree.internode_labeler import InternodeLabeler
from .services.tree.lb_score import LBScore
from .services.tree.dvmc import DVMC

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

    ## Tree functions
    def dvmc(self):
        parser = ArgumentParser()
        parser.add_argument(
            "-r", "--root", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        args = parser.parse_args(sys.argv[2:])
        DVMC(args).run()

    def rcv(self):
        raise NotImplementedError()

    def treeness(self):
        parser = ArgumentParser()
        parser.add_argument("tree", type=str)
        args = parser.parse_args(sys.argv[2:])
        Treeness(args).run()

    def total_tree_length(self):
        parser = ArgumentParser()
        parser.add_argument("tree", type=str)
        args = parser.parse_args(sys.argv[2:])
        TotalTreeLength(args).run()
    
    def lb_score(self):
        parser = ArgumentParser()
        parser.add_argument("tree", type=str)
        args = parser.parse_args(sys.argv[2:])
        LBScore(args).run()

    def internode_labeler(self):
        parser = ArgumentParser()
        parser.add_argument("tree", type=str)
        # TODO: add output option?
        args = parser.parse_args(sys.argv[2:])
        InternodeLabeler(args).run()

    ### Alignment functions
    def alignment_length(self):
        parser = ArgumentParser()
        parser.add_argument("alignment", type=str)
        args = parser.parse_args(sys.argv[2:])
        AlignmentLength(args).run()

    def variable_sites(self):
        parser = ArgumentParser()
        parser.add_argument("alignment", type=str)
        args = parser.parse_args(sys.argv[2:])
        VariableSites(args).run()


if __name__ == "__main__":
    Phykit()
