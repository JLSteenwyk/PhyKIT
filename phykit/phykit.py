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

from .services.treeness import Treeness

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


if __name__ == "__main__":
    Phykit()
