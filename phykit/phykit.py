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
from .services.alignment.parsimony_informative_sites import ParsimonyInformative
from .services.alignment.variable_sites import VariableSites
from .services.alignment.alignment_length_no_gaps import AlignmentLengthNoGaps
from .services.alignment.rcv import RelativeCompositionVariability
from .services.alignment.dna_threader import DNAThreader

from .services.tree.bipartition_support_stats import BipartitionSupportStats
from .services.tree.treeness import Treeness
from .services.tree.total_tree_length import TotalTreeLength
from .services.tree.internode_labeler import InternodeLabeler
from .services.tree.lb_score import LBScore
from .services.tree.dvmc import DVMC
from .services.tree.internal_branch_stats import InternalBranchStats
from .services.tree.patristic_distances import PatristicDistances
from .services.tree.rf_distance import RobinsonFouldsDistance
from .services.tree.treeness_over_rcv import TreenessOverRCV
from .services.tree.spurious_sequence import SpuriousSequence

logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)


class Phykit(object):

    help_header = """
                 _____  _           _  _______ _______ 
                |  __ \| |         | |/ /_   _|__   __|
                | |__) | |__  _   _| ' /  | |    | |   
                |  ___/| '_ \| | | |  <   | |    | |   
                | |    | | | | |_| | . \ _| |_   | |   
                |_|    |_| |_|\__, |_|\_\_____|  |_|   
                               __/ |                   
                              |___/   
                            
                Citation: Steenwyk et al. Journal, journal info, link
    """

    def __init__(self):
        parser = ArgumentParser(
            add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                PhyKIT helps with all things phylogenetics and phylogenomics.

                Usage: phykit <command> [optional command arguments]

                Command specific help messages can be viewed by adding a 
                -h/--help argument after the command. For example, to see the
                to see the help message for the command 'treeness', execute
                phykit treeness -h or phykit treeness --help.

                Alignment-based commands
                ========================
                alignment_length
                    - calculates alignment length
                alignment_length_no_gaps
                    - calculates alignment length after removing sites with gaps
                parsimony_informative_sites
                    - calculates the number and percentage of parsimony
                      informative sites in an alignment
                rcv
                    - calculates relative composition variability in an alignment
                variable_sites
                    - calculates the number and percentage of variable sites
                      in an alignment

                Tree-based commands
                ===================
                bipartition_support_stats
                    - calculates summary statistics for bipartition support
                dvmc 
                    - reports the degree of violation of the molecular clock
                spurious_sequence
                    - identifies putatively spurious sequences by identifying
                      branch lengths that are atypically long
                treeness
                    - reports treeness or stemminess, a measure of signal-to-
                      noise ratio in a phylogeny
                total_tree_length
                    - calculates total tree length
                lb_score
                    - calculates lb (long branch) score for taxa in a phylogeny
                internal_branch_stats
                    - calculates summary statistics for internal branch lengths 
                internode_labeler
                    - create labels at internodes in a phylogeny

                Alignment- and tree-based commands
                ==================================
                treeness_over_rcv
                    - calculates treeness/rcv, treeness, and rcv

                Helper commands
                ===============
                thread_dna
                    - thread dna sequences over a protein alignment
                • create concatenation matrix
                
                
                

                                        
                """
            ),
        )
        parser.add_argument("command", help=SUPPRESS)
        args = parser.parse_args(sys.argv[1:2])

        if not hasattr(self, args.command):
            print("Invalid command option. See help for more details.")
            parser.print_help()
            sys.exit(1)

        getattr(self, args.command)()

        ### Alignment functions

    ## Alignment functions
    def alignment_length(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Longer alignments are associated with strong phylogenetic signal.

                Length of the input alignment is calculated using this function.
                Association between alignment length and phylogenetic signal
                was determined by Shen et al., Genome Biology and Evolution (2016),
                doi: 10.1093/gbe/evw179.

                Usage:
                phykit alignment_length <file>

                Options
                =====================================================
                <file>                      first argument after 
                                            function name should be
                                            an alignment file           
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        AlignmentLength(args).run()

    def alignment_length_no_gaps(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Longer alignments when excluding sites with gaps is
                associated with strong phylogenetic signal.

                PhyKIT reports three tab delimited values:
                col1: number of sites without gaps
                col2: total number of sites
                col3: percentage of sites without gaps

                Association between alignment length when excluding sites
                with gaps and phylogenetic signal was determined by Shen 
                et al., Genome Biology and Evolution (2016), 
                doi: 10.1093/gbe/evw179.

                Usage:
                phykit alignment_length_no_gaps <file>

                Options
                =====================================================
                <file>                      first argument after 
                                            function name should be
                                            an alignment file          
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        AlignmentLengthNoGaps(args).run()

    def parsimony_informative_sites(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                The number of parsimony informative sites in an alignment
                is associated with strong phylogenetic signal.

                PhyKIT reports three tab delimited values:
                col1: number of parsimony informative sites
                col2: total number of sites
                col3: percentage of parsimony informative sites

                Association between the number of parsimony informative
                sites and phylogenetic signal was determined by Shen 
                et al., Genome Biology and Evolution (2016), 
                doi: 10.1093/gbe/evw179.

                Usage:
                phykit parsimony_informative_sites <file>

                Options
                =====================================================
                <file>                      first argument after 
                                            function name should be
                                            an alignment file          
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        ParsimonyInformative(args).run()

    def rcv(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Lower RCV (relative composition variability) values are thought
                to be desirable because they represent a lower composition bias.

                RCV describes the average variability in composition between taxa. 

                Calculate RCV following Phillips and Penny, Molecular Phylogenetics
                and Evolution (2003), doi: 10.1016/S1055-7903(03)00057-5.

                Usage:
                phykit rcv <file>

                Options
                =====================================================
                <file>                      first argument after 
                                            function name should be
                                            an alignment file          
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        RelativeCompositionVariability(args).run()

    def variable_sites(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                The number of variable sites in an alignment is 
                associated with strong phylogenetic signal.

                PhyKIT reports three tab delimited values:
                col1: number of variable sites
                col2: total number of sites
                col3: percentage of variable sites

                Association between the number of variable sites and
                phylogenetic signal was determined by Shen et al.,
                Genome Biology and Evolution (2016), 
                doi: 10.1093/gbe/evw179.

                Usage:
                phykit variable_sites <file>

                Options
                =====================================================
                <file>                      first argument after 
                                            function name should be
                                            an alignment file          
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        VariableSites(args).run()


    ## Tree functions
    def bipartition_support_stats(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}
                
                High bipartition support values are thought to be desirable because
                they are indicative of greater certainty in the tree topology.

                Calculate summary statistics for bipartition support. Summary
                statistics include mean, median, 25th percentile, 75th percentile,
                minimum, maximum, standard deviation, and variance. 

                To obtain all bipartition support values, use the -v/--verbose option.

                Usage:
                phykit bipartition_support_stats <file> 

                Options
                =====================================================
                <file>                      first argument after 
                                            function name should be
                                            a tree file 
            
                -v, --verbose               optional argument to print
                                            all bipartition support
                                            values
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-v", "--verbose", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        BipartitionSupportStats(args).run()

    def dvmc(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Lower DVMC (degree of violation of the molecular clock) values are
                thought to be desirable because they are indicative of a lower degree of
                violation in the molecular clock assumption.

                Typically, outgroup taxa are not included in molecular clock analysis. Thus,
                prior to calculating DVMC from a single gene tree, outgroup taxa are pruned
                from the phylogeny. PhyKIT will prune taxa using tip names from a single column
                file, which is specified using the -r/--root file. If the tip name does not
                exist in the input tree, rather than raising an error/warning message, the tip  
                name is skipped. If the user wants to calculate DVMC for phylogenies with incomplete
                taxa representation, this will allow the user to use one <root file> for all trees. 

                Calculate degree of violation of the molecular clock (or DVMC) in a tree
                following Liu et al., PNAS (2017), doi: 10.1073/pnas.1616744114.

                Usage:
                phykit dvmc <file> 

                Options
                =====================================================
                -t, --tree <file>           input file tree name
            
                -r, --root <file>           single column file with
                                            tip names of root taxa
                """
            ),
        )
        parser.add_argument(
            "-r", "--root", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        args = parser.parse_args(sys.argv[2:])
        DVMC(args).run()

    def internal_branch_stats(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Internal branch lengths can be useful for tree diagnostics.

                Summary statistics of internal branch lengths include mean,
                median, 25th percentile, 75th percentile, minimum, maximum,
                standard deviation, and variance of per taxon LB scores is reported.
                To obtain all internal branch lengths, use the -v/--verbose option. 

                Usage:
                phykit internal_branch_stats <file> 

                Options
                =====================================================
                <file>                      first argument after 
                                            function name should be
                                            a tree file

                -v, --verbose               optional argument to print
                                            all LB score values           
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-v", "--verbose", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        InternalBranchStats(args).run()

    # TODO: fix documentation and finish writing function
    def internode_labeler(self):
        parser = ArgumentParser()
        parser.add_argument("tree", type=str)
        # TODO: add output option?
        args = parser.parse_args(sys.argv[2:])
        InternodeLabeler(args).run()

    def lb_score(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Lower LB (long branch) scores are thought to be desirable
                because they are indicative of taxa or trees that likely do
                not have issues with long branch attraction.

                LB score is calculated from patristic distances (or the sum 
                of branches between to two taxa). More specifically, it is
                mean pairwise patristic distance of taxon i compared to
                all other taxa over the average pairwise patristic distance.
                Summary statistics reported include mean, median, 25th
                percentile, 75th percentile, minimum, maximum, standard 
                deviation, and variance of per taxon LB scores is reported.
                To obtain LB scores for each taxa, use the -v/--verbose option. 

                LB scores are calculated following Struck, Evolutionary 
                Bioinformatics (2014), doi: 10.4137/EBO.S14239.

                Usage:
                phykit lb_score <file> 

                Options
                =====================================================
                <file>                      first argument after 
                                            function name should be
                                            a tree file

                -v, --verbose               optional argument to print
                                            all LB score values            
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-v", "--verbose", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        LBScore(args).run()

    def patristic_distances(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Patristic distances describes the distance from tip to tip.

                Summary statistics reported include mean, median, 25th
                percentile, 75th percentile, minimum, maximum, standard 
                deviation, and variance of patristic distances across taxa.
                To obtain all patristic distances, use the -v/--verbose option.
                With the -v option, the first column will have two taxon names
                separated by a '-' followed by the patristic distance. Features
                will be tab separated. 

                Usage:
                phykit patristic_distances <file> 

                Options
                =====================================================
                <file>                      first argument after 
                                            function name should be
                                            a tree file

                -v, --verbose               optional argument to print
                                            all patristic distances between
                                            taxa            
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-v", "--verbose", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        PatristicDistances(args).run() 

    def rf_distance(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Low (RF) Robinson-Foulds distances reflect greater similarity between
                two phylogenies. This function prints out two values, the plain
                RF value and the normalized RF value, which are separated by a tab.
                Normalized RF values are calculated by taking the plain RF value and
                dividing it by 2(n-3) where n is the number of tips in the phylogeny. 

                RF distances are calculated following Robinson & Foulds, Mathematical 
                Biosciences (1981), doi: 10.1016/0025-5564(81)90043-2.

                Usage:
                phykit rf_distance <tree_file_zero> <tree_file_one>

                Options
                =====================================================
                <tree_file_zero>            first argument after 
                                            function name should be
                                            an alignment file

                <tree_file_one>             first argument after 
                                            function name should be
                                            an alignment file           
                """
            ),
        )
        parser.add_argument("tree_zero", type=str, help=SUPPRESS)
        parser.add_argument("tree_one", type=str, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        RobinsonFouldsDistance(args).run()

    # TODO: create unit test for this function
    def spurious_sequence(self):
        parser = ArgumentParser(add_help=True,
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

                Identifies potentially spurious sequences and reports
                tips in the phylogeny that could possibly be removed
                from the underlying multiple sequence alignments. PhyKIT
                does so by identifying and reporting long terminal branches
                defined as branches that are equal to or 20 times the median
                length of internal and terminal branches.
                
                Using this method to identify potentially spurious sequences
                was, to my knowledge, first introduced by Shen et al., (2018)
                Cell doi: 10.1016/j.cell.2018.10.023.                

                Usage:
                phykit spurious_sequence <file> [-f 20]

                Options
                =====================================================
                <file>                      first argument after 
                                            function name should be
                                            an alignment file

                -f/--factor                 factor to multiply median
                                            branch length by to calculate
                                            the threshold of long branches.
                                            Default: 20
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument(
            "-f", "--factor",
            type=float, required=False,
            help=SUPPRESS
        )
        args = parser.parse_args(sys.argv[2:])
        SpuriousSequence(args).run()

    def total_tree_length(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Calculate total tree length, which is a sum of all branches. 

                Usage:
                phykit total_tree_length <file>

                Options
                =====================================================
                <file>                      first argument after 
                                            function name should be
                                            a tree file
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        TotalTreeLength(args).run()

    def treeness(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Higher treeness values are thought to be desirable because they
                represent a higher signal-to-noise ratio.

                Treeness describes the proportion of tree distance on internal
                branches. Treeness can be used as a measure of the signal-to-noise
                ratio in a phylogeny. 

                Calculate treeness (also referred to as stemminess) following
                Lanyon, The Auk (1988), doi: 10.1093/auk/105.3.565 and
                Phillips and Penny, Molecular Phylogenetics and Evolution
                (2003), doi: 10.1016/S1055-7903(03)00057-5.

                Usage:
                phykit treeness <file>

                Options
                =====================================================
                <file>                      first argument after 
                                            function name should be
                                            a tree file
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        Treeness(args).run()

    def treeness_over_rcv(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Higher treeness/RCV values are thought to be desirable because
                they harbor a high signal-to-noise ratio are least susceptible
                to composition bias.

                PhyKIT reports three tab delimited values:
                col1: treeness/RCV
                col2: treeness
                col3: RCV

                Treeness describes the proportion of tree distance on internal
                branches. RCV

                Calculate treeness/RCV following Phillips and Penny, Molecular 
                Phylogenetics and Evolution (2003), doi: 10.1016/S1055-7903(03)00057-5.

                Usage:
                phykit treeness_over_rcv <file>

                Options
                =====================================================
                -a/--alignment              an alignment file
                
                -t/--tree                   a tree file
                """
            ),
        )
        parser.add_argument(
            "-a", "--alignment", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        args = parser.parse_args(sys.argv[2:])
        TreenessOverRCV(args).run()

    ### Helper commands
    # TODO: thread DNA unit tests
    def thread_dna(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Thread DNA sequence onto a protein alignment to create a
                codon-based alignment. Note, sequences should occur in the
                same order in the protein and nucleotide alignment.

                Usage:
                phykit thread_dna -p <file> -n <file> [-s]

                Options
                =====================================================
                -p/--protein                protein alignment file

                -n/--nucleotide             nucleotide alignment file

                -s/--stop                   boolean for whether or not
                                            stop codons should be kept
                """
            ),
        )
        parser.add_argument(
            "-p", "--protein",
            type=str,
            help=SUPPRESS
        )
        parser.add_argument(
            "-n", "--nucleotide",
            type=str,
            help=SUPPRESS
        )
        parser.add_argument(
            "-s", "--stop",
            type=bool,
            help=SUPPRESS
        )
        args = parser.parse_args(sys.argv[2:])
        DNAThreader(args).run()

if __name__ == "__main__":
    Phykit()
