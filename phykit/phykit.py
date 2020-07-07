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

from .services.alignment import (
    AlignmentLength,
    AlignmentLengthNoGaps,
    CreateConcatenationMatrix,
    DNAThreader,
    GCContent,
    PairwiseIdentity,
    ParsimonyInformative,
    RelativeCompositionVariability,
    RenameFastaEntries,
    VariableSites
)

from .services.tree import (
    BipartitionSupportStats,
    BranchLengthMultiplier,
    CollapseBranches,
    CovaryingEvolutionaryRates,
    DVMC,
    InternalBranchStats,
    InternodeLabeler,
    LBScore,
    TotalTreeLength,
    PatristicDistances,
    PolytomyTest,
    PrintTree,
    PruneTree,
    RenameTreeTips,
    RobinsonFouldsDistance,
    Saturation,
    SpuriousSequence,
    TipLabels,
    Treeness,
    TreenessOverRCV
)

from .helpers.boolean_argument_parsing import str2bool


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

                PhyKIT helps with building and understanding phylogenetic trees.

                More specifically, PhyKIT is a collection of functions that perform
                various operations on multiple sequence alignments and newick tree
                files. Generally, all functions are designed to help understand
                the contents of alignments (e.g., gc content or the number of parsimony
                informative sites) and the shape of trees (e.g., treeness, degree of
                violation of a molecular clock).

                Usage: phykit <command> [optional command arguments]

                Command specific help messages can be viewed by adding a 
                -h/--help argument after the command. For example, to see the
                to see the help message for the command 'treeness', execute
                "phykit treeness -h" or "phykit treeness --help".

                Lastly, each function comes with aliases to save the user some
                key strokes. For example, to get the help message for the 'treeness'
                function, you can type "phykit tness -h".

                Alignment-based commands
                ========================
                alignment_length (alias: aln_len; al)
                    - calculates alignment length
                alignment_length_no_gaps (alias: aln_len_no_gaps; alng)
                    - calculates alignment length after removing sites with gaps
                create_concatenation_matrix (alias: create_concat; cc)
                    - create concatenation matrix from a set of alignments                    
                gc_content (alias: gc)
                    - calculate GC content of a fasta entries or entries thereof
                pairwise_identity (alias: pi)
                    - calculates average pairwise identify among sequences in
                      an alignment file. This is a proxy for evolutionary rate
                parsimony_informative_sites (alias: pis)
                    - calculates the number and percentage of parsimony
                      informative sites in an alignment
                relative_composition_variability (alias: rcv)
                    - calculates relative composition variability in an alignment
                rename_fasta_entries (alias: rename_fasta)
                    - rename entries in a fasta file
                thread_dna (alias: pal2nal; p2n)
                    - thread dna sequences over a protein alignment
                variable_sites (alias: vs)
                    - calculates the number and percentage of variable sites
                      in an alignment

                Tree-based commands
                ===================
                bipartition_support_stats (alias: bss)
                    - calculates summary statistics for bipartition support
                branch_length_multiplier (alias: blm)
                    - multiply all branch lengths by a specified factor
                collapse_branches (alias: collapse; cb)
                    - collapses branches according to bipartition support
                covarying_evolutionary_rates (alias: cover)
                    - calculates correlation in the evolutionary rate of two trees
                degree_of_violation_of_a_molecular_clock (alias: dvmc)
                    - reports the degree of violation of the molecular clock
                internal_branch_stats (alias: ibs)
                    - calculates summary statistics for internal branch lengths 
                internode_labeler (alias: il)
                    - create labels at internodes in a phylogeny
                long_branch_score (alias: lb_score; lbs)
                    - calculates lb (long branch) score for taxa in a phylogeny
                total_tree_length (alias: tree_len)
                    - calculates total tree length
                patristic_distances (alias: pd)
                    - calculate all pairwise distances between tips in a tree
                polytomy_test
                    - conducts a polytomy test using triplet and gene
                      support frequencies
                print_tree (alias: print)
                    - prints ascii tree
                prune_tree (alias: prune)
                    - prune taxa from a phylogeny
                rename_tree_tips (alias: rename_tree)
                    - renames tips in a phylogeny according to a file with
                      the desired new tip names
                robinson_foulds_distance (alias: rf_distance; rf)
                    - calculates Robinson-Foulds distance between two trees
                spurious_sequence (alias: spurious_seq; ss)
                    - identifies putatively spurious sequences by identifying
                      branch lengths that are atypically long
                tip_labels (alias: labels)
                    - print leaf names in a phylogeny
                treeness (alias: tness)
                    - reports treeness or stemminess, a measure of signal-to-
                      noise ratio in a phylogeny
             
                Alignment- and tree-based commands
                ==================================
                saturation (alias: sat)
                    - calculates saturation by examining the slope of
                      patristic distance and uncorrected distances
                treeness_over_rcv (alias: toverr)
                    - calculates treeness/rcv, treeness, and rcv
                """
            ),
        )
        parser.add_argument("command", help=SUPPRESS)
        args = parser.parse_args(sys.argv[1:2])

        # if command is part of the possible commands (i.e., the long form
        # commands, run). Otherwise, assume it is an alias and look to the
        # run_alias function
        if hasattr(self, args.command):
            getattr(self, args.command)()
        else:
            self.run_alias(args.command)

    ## Aliases
    def run_alias(self, command):
        # Alignment aliases
        if command in ['aln_len', 'al']:
            return self.alignment_length()
        elif command in ['aln_len_no_gaps', 'alng']:
            return self.alignment_length_no_gaps()
        elif command == 'gc':
            return self.gc_content()
        elif command == 'pi':
            return self.pairwise_identity()
        elif command == 'pis':
            return self.parsimony_informative_sites()
        elif command == 'relative_composition_variability':
            return self.rcv()
        elif command == 'rename_fasta':
            return self.rename_fasta_entries()
        elif command == 'vs':
            return self.variable_sites()
        # Tree aliases
        elif command == 'bss':
            return self.bipartition_support_stats()
        elif command == 'blm':
            return self.branch_length_multiplier()
        elif command in ['collapse', 'cb']:
            return self.collapse_branches()
        elif command == 'cover':
            return self.covarying_evolutionary_rates()
        elif command == 'degree_of_violation_of_a_molecular_clock':
            return self.dvmc()
        elif command == 'ibs':
            return self.internal_branch_stats()
        elif command == 'il':
            return self.internode_labeler()
        elif command in ['long_branch_score', 'lbs']:
            return self.lb_score()
        elif command == 'pd':
            return self.patristic_distances()
        elif command in ['polyt_test', 'ptt', 'polyt']:
            return self.polytomy_test()
        elif command in ['print', 'pt']:
            return self.print_tree()
        elif command == 'prune':
            return self.prune_tree()
        elif command in ['rename_tree', 'rename_tips']:
            return self.rename_tree_tips()
        elif command in ['robinson_foulds_distance', 'rf_dist', 'rf']:
            return self.rf_distance()
        elif command in ['spurious_seq', 'ss']:
            return self.spurious_sequence()
        elif command == 'labels':
            return self.tip_labels()
        elif command == 'tree_len':
            return self.total_tree_length
        elif command == 'tness':
            return self.treeness()
        # Alignment- and tree-based aliases
        elif command == 'sat':
            return self.saturation()
        elif command == 'toverr':
            return self.treeness_over_rcv()
        # Helper aliases
        elif command in ['create_concat', 'cc']:
            return self.create_concatenation_matrix()
        elif command in ['pal2nal', 'p2n']:
            return self.thread_dna()
        else:
            print("Invalid command option. See help for more details.")
            parser.print_help()
            sys.exit(1)

    # TODO: include alias in help messages
    ## Alignment functions
    def alignment_length(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Length of the input alignment is calculated using this function.

                Longer alignments are associated with strong phylogenetic signal.
                
                Association between alignment length and phylogenetic signal
                was determined by Shen et al., Genome Biology and Evolution (2016),
                doi: 10.1093/gbe/evw179.

                Aliases: aln_len, al

                Usage:
                phykit alignment_length <alignment>

                Options
                =====================================================
                <alignment>                 first argument after 
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

                Calculate alignment length excluding sites with gaps.

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

                Aliases: aln_len_no_gaps, alng

                Usage:
                phykit alignment_length_no_gaps <alignment>

                Options
                =====================================================
                <alignment>                 first argument after 
                                            function name should be
                                            an alignment file          
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        AlignmentLengthNoGaps(args).run()

    # TODO: write tests
    def gc_content(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}
                
                Calculate GC content of a fasta file.

                GC content is negatively correlated to phylogenetic signal.

                If there are multiple entries, use the -v/--verbose option
                to determine the GC content of each fasta entry separately.

                Association between GC content and phylogenetic signal was
                determined by Shen et al., Genome Biology and Evolution (2016), 
                doi: 10.1093/gbe/evw179.

                Alias: gc

                Usage:
                phykit gc_content <fasta> [-v/--verbose]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file 
            
                -v, --verbose               optional argument to print
                                            the GC content of each fasta
                                            entry
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-v", "--verbose", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        GCContent(args).run()

    # TODO: write unit tests
    # TODO: consider renaming to evolutionary_rate
    def pairwise_identity(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Calculate the average pairwise identity among sequences.
                
                Pairwise identities can be used as proxies for the
                evolutionary rate of sequences.

                Pairwise identity is defined as the number of identical
                columns between two aligned sequences divided by the
                number of columns in the alignment. Summary statistics
                are reported but with the verbose option, all pairwise
                identities will be reported.

                An example of pairwise identities being used as a proxy
                for evolutionary rate can be found here: Chen et al. 
                Genome Biology and Evolution (2017), doi: 10.1093/gbe/evx147.

                Alias: pi

                Usage:
                phykit pairwise_identity <alignment> [-v/--verbose]

                Options
                =====================================================
                <alignment>                 first argument after 
                                            function name should be
                                            an alignment file  

                -v, --verbose               optional argument to print
                                            identity per pair       
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        parser.add_argument("-v", "--verbose", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        PairwiseIdentity(args).run()

    def parsimony_informative_sites(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Calculate the number and percentage of parismony
                informative sites in an alignment.

                The number of parsimony informative sites in an alignment
                is associated with strong phylogenetic signal.

                PhyKIT reports three tab delimited values:
                col1: number of parsimony informative sites
                col2: total number of sites
                col3: percentage of parsimony informative sites

                Association between the number of parsimony informative
                sites and phylogenetic signal was determined by Shen 
                et al., Genome Biology and Evolution (2016), 
                doi: 10.1093/gbe/evw179

                Alias: pis

                Usage:
                phykit parsimony_informative_sites <alignment>

                Options
                =====================================================
                <alignment>                 first argument after 
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

                Calculate RCV (relative composition variability) for an alignment.

                Lower RCV values are thought to be desirable because they represent
                a lower composition bias in an alignment. Statistically, RCV describes
                the average variability in composition among taxa. 

                Calculate RCV following Phillips and Penny, Molecular Phylogenetics
                and Evolution (2003), doi: 10.1016/S1055-7903(03)00057-5.

                Alias: rcv

                Usage:
                phykit relative_composition_variability <alignment>

                Options
                =====================================================
                <alignment>                 first argument after 
                                            function name should be
                                            an alignment file          
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        RelativeCompositionVariability(args).run()

    # TODO: write test
    def rename_fasta_entries(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Renames fasta entries.

                Renaming fasta entries will follow the scheme of a tab-delimited
                file wherein the first column is the current fasta entry name and
                the second column is the new fasta entry name in the resulting 
                output alignment. 

                Alias: rename_fasta

                Usage:
                phykit rename_fasta_entries <fasta> -i/--idmap <idmap>
                    [-o/--output <output_file>]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file

                -i/--idmap                  identifier map of current tip
                                            names (col1) and desired tip
                                            names (col2)

                -o/--output                 optional argument to write
                                            the renamed fasta file to.
                                            Default output has the same 
                                            name as the input file with
                                            the suffix ".renamed.fa" added
                                            to it.
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-i","--idmap", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        RenameFastaEntries(args).run()

    def variable_sites(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Calculate the number of variable sites in an alignment.

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

                Alias: vs

                Usage:
                phykit variable_sites <alignment>

                Options
                =====================================================
                <alignment>                 first argument after 
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
                they are indicative of greater certainty in tree topology.

                Calculate summary statistics for bipartition support. Summary
                statistics include mean, median, 25th percentile, 75th percentile,
                minimum, maximum, standard deviation, and variance. 

                To obtain all bipartition support values, use the -v/--verbose option.

                Alias: bss

                Usage:
                phykit bipartition_support_stats <tree> [-v/--verbose]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file 
            
                -v, --verbose               optional argument to print
                                            all bipartition support
                                            values
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument(
            "-v",
            "--verbose",
            action="store_true",
            required=False,
            help=SUPPRESS
        )
        args = parser.parse_args(sys.argv[2:])
        BipartitionSupportStats(args).run()

    # TODO: create unit test for this function
    def branch_length_multiplier(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header} 

                Multiply branch lengths in a phylogeny by a given factor.
                
                This can help modify reference trees when conducting simulations
                or other analyses.              

                Alias: blm

                Usage:
                phykit branch_length_multiplier <tree> -f n [-o/--output <output_file>]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            an tree file

                -f/--factor                 factor to multiply branch 
                                            lengths by 

                -o/--output                 optional argument to name 
                                            the outputted tree file.
                                            Default output will have 
                                            the same name as the input
                                            file but with the suffix 
                                            ".factor_(n).tre"
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument(
            "-f", "--factor",
            type=float, required=True,
            help=SUPPRESS
        )
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        BranchLengthMultiplier(args).run()

    # TODO: create unit test for this function
    def collapse_branches(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header} 

                Collapse branches on a phylogeny according to bipartition support.
                Bipartitions will be collapsed if they are less than the user specified
                value.              

                Alias: collapse, cb

                Usage:
                phykit collapse_branches <tree> -s/--support n [-o/--output <output_file>]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            an tree file

                -s/--support                bipartitions with support less
                                            than this value will be collapsed

                -o/--output                 optional argument to name 
                                            the outputted tree file.
                                            Default output will have 
                                            the same name as the input
                                            file but with the suffix 
                                            ".collapsed_(support).tre""

                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument(
            "-s", "--support",
            type=float, required=True,
            help=SUPPRESS
        )
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        CollapseBranches(args).run()

    # write unit tests
    def covarying_evolutionary_rates(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Determine if two genes have a signature of coevolving with one another.

                Genes that have covarying evolutionary histories tend to have 
                similar functions and expression levels.

                Input two trees and calculate the correlation of branch lengths
                between the tree trees. The two input trees do not have to have
                the same taxa. This function will first prune both trees to have
                the same tips. Additionally, branch lengths must be corrected. 
                Branch length correction is typically done using the putative
                species tree's branch lengths. As recommended by the original
                method developers, outlier branches are removed. Outlier branches
                have a relative evolutionary rate greater than five.

                PhyKIT reports three tab delimited values:
                col1: correlation coefficient
                col2: p-value

                Method is empirically evaluated by Clark et al., Genome Research
                (2012), doi: 10.1101/gr.132647.111.

                Alias: cover

                Usage:
                phykit covarying_evolutionary_rates <tree_file_zero> <tree_file_one>
                    -r/--reference <reference_tree_file> [-v/--verbose]

                Options
                =====================================================
                <tree_file_zero>            first argument after 
                                            function name should be
                                            an alignment file

                <tree_file_one>             first argument after 
                                            function name should be
                                            an alignment file 

                -r/--reference              a tree to correct branch
                                            lengths by in the two input
                                            trees. Typically, this is a
                                            putative species tree.

                -v/--verbose                print out corrected branch
                                            lengths shared between
                                            tree 0 and tree 1
                """
            ),
        )
        parser.add_argument("tree_zero", type=str, help=SUPPRESS)
        parser.add_argument("tree_one", type=str, help=SUPPRESS)
        parser.add_argument(
            "-r", "--reference", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument("-v", "--verbose", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        CovaryingEvolutionaryRates(args).run()

    def dvmc(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Calculate DVMC (degree of violation of the molecular clock) of a tree.

                Lower DVMC values are thought to be desirable because they are indicative
                of a lower degree of violation in the molecular clock assumption.

                Typically, outgroup taxa are not included in molecular clock analysis. Thus,
                prior to calculating DVMC from a single gene tree, outgroup taxa are pruned
                from the phylogeny. PhyKIT will prune taxa using tip names from a single column
                file, which is specified using the -r/--root file. If the tip name does not
                exist in the input tree, rather than raising an error/warning message, the tip  
                name is skipped. If the user wants to calculate DVMC for phylogenies with incomplete
                taxa representation, this will allow the user to use one <root file> for all trees. 

                Calculate degree of violation of the molecular clock (or DVMC) in a tree
                following Liu et al., PNAS (2017), doi: 10.1073/pnas.1616744114.

                Alias: dvmc

                Usage:
                phykit degree_of_violation_of_a_molecular_clock -t/--tree <tree>
                    -r/--root <root_taxa>

                Options
                =====================================================
                -t, --tree                  input file tree name
            
                -r, --root                  single column file with
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

                Calculate summary statistics for internal branch lengths in a tree.

                Internal branch lengths can be useful for tree diagnostics.

                Summary statistics of internal branch lengths include mean,
                median, 25th percentile, 75th percentile, minimum, maximum,
                standard deviation, and variance of per taxon LB scores is reported.
                To obtain all internal branch lengths, use the -v/--verbose option. 

                Alias: ibs

                Usage:
                phykit internal_branch_stats <tree> [-v/--verbose]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -v, --verbose               optional argument to print
                                            all LB score values           
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument(
            "-v",
            "--verbose",
            action="store_true",
            required=False,
            help=SUPPRESS
        )
        args = parser.parse_args(sys.argv[2:])
        InternalBranchStats(args).run()

    # TODO: fix documentation and finish writing function
    def internode_labeler(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Appends numerical identifiers to bipartitions in place
                of support values. This is helpful for pointing to
                specific internodes in supplementary files or otherwise. 

                Alias: il

                Usage:
                phykit internode_labeler <file> [-o/--output <file>]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -o/--output                 optional argument to name 
                                            the outputted tree file
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        # TODO: write in functionality of using the output argument
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        InternodeLabeler(args).run()

    def lb_score(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Calculate LB (long branch) scores in a phylogeny.

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

                Alias: lb_score, lbs

                Usage:
                phykit long_branch_score <tree> [-v/--verbose]

                Options
                =====================================================
                <tree>                      first argument after 
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

                Calculate patristic distances in a phylogeny.

                Patristic distances describes the distance from tip to tip.

                Summary statistics reported include mean, median, 25th
                percentile, 75th percentile, minimum, maximum, standard 
                deviation, and variance of patristic distances across taxa.
                To obtain all patristic distances, use the -v/--verbose option.
                With the -v option, the first column will have two taxon names
                separated by a '-' followed by the patristic distance. Features
                will be tab separated. 

                Alias: pd

                Usage:
                phykit patristic_distances <tree> [-v/--verbose]

                Options
                =====================================================
                <tree>                      first argument after 
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

    # TODO: write unit tests
    def polytomy_test(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Conduct a polytomy test for three clades in a phylogeny.

                Polytomy tests can be used to identify putative radiations
                as well as identify well supported alternative topologies.

                The polytomy testing function takes as input a file with
                the three groups of taxa to test the relationships for and
                a single column file with the names of the desired tree files
                to use for polytomy testing. Next, the script to examine
                support for the grouping of the three taxa using triplets
                and gene support frequencies. 

                This function can account for uncertainty in gene trees - 
                that is, the input phylogenies can have collapsed bipartitions.

                Thereafter, a chi-squared test is conducted to determine if there
                is evidence to reject the null hypothesis wherein the null 
                hypothesis is that the three possible topologies among the three
                groups are equally supported. This test is done individually for
                triplets and gene support frequencies.

                Alias: polyt_test, polyt, ptt

                Usage:
                phykit polytomy_test -t/--trees <trees> -g/--groups <groups>

                Options
                =====================================================
                -t, --trees                 single column file with names
                                            of phylogenies to use for
                                            polytomy testing

                -g, --groups                a tab-delimited file with the
                                            grouping designations to test.
                                            Lines starting with comments 
                                            are not considered. Names
                                            of individual taxa should be
                                            separated by a semi-colon  ';'
                                            
                For example, the groups file
                #labels group0  group1  group2
                name_of_test    tip_name_A;tip_name_B   tip_name_C  tip_name_D;tip_name_E
                """
            ),
        )
        parser.add_argument("-t", "--trees", type=str, help=SUPPRESS)
        parser.add_argument("-g", "--groups", type=str, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        PolytomyTest(args).run() 

    # TODO: unit test
    def print_tree(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Print ascii tree of input phylogeny.

                Phylogeny can be printed with or without branch lengths.
                By default, the phylogeny will be printed with branch lengths
                but branch lengths can be removed using the -r/--remove argument.

                Alias: print, pt

                Usage:
                phykit print_tree <tree> [-r/--remove]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -r, --remove                optional argument to print
                                            the phylogeny without branch
                                            lengths       
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-r", "--remove", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        PrintTree(args).run()

    # TODO: unit test
    def prune_tree(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Prune tips from a phylogeny.

                Provide a single column file with the names of the tips
                in the input phylogeny you would like to prune from the
                tree.

                Alias: prune

                Usage:
                phykit prune_tree <tree> <list_of_taxa> [-o/--output <output_file>]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                <list_of_taxa>                single column file with the
                                            names of the tips to remove
                                            from the phylogeny 

                -o/--output                 name of output file for the
                                            pruned phylogeny. 
                                            Default output will have 
                                            the same name as the input
                                            file but with the suffix 
                                            ".pruned"    
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("list_of_taxa", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        PruneTree(args).run()

    # TODO: fix documentation and finish writing function
    def rename_tree_tips(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Renames tips in a phylogeny.

                Renaming tip files will follow the scheme of a tab-delimited
                file wherein the first column is the current tip name and the
                second column is the desired tip name in the resulting 
                phylogeny. 

                Alias: rename_tree, rename_tips

                Usage:
                phykit rename_tree_tips <tree> -i/--idmap <idmap.txt>
                    [-o/--output <output_file>] 

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -i/--idmap                  identifier map of current tip
                                            names (col1) and desired tip
                                            names (col2)

                -o/--output                 optional argument to write
                                            the renamed tree files to.
                                            Default output will have 
                                            the same name as the input
                                            file but with the suffix 
                                            ".renamed"          
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-i","--idmap", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        RenameTreeTips(args).run()

    def rf_distance(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Calculate (RF) Robinson-Foulds distance between two trees.

                Low RF distances reflect greater similarity between two phylogenies. 
                This function prints out two values, the plain RF value and the
                normalized RF value, which are separated by a tab. Normalized RF values
                are calculated by taking the plain RF value and dividing it by 2(n-3)
                where n is the number of tips in the phylogeny. 

                PhyKIT will print out both the plain RF distance (col 1) and the
                normalized RF distance (col 2).

                RF distances are calculated following Robinson & Foulds, Mathematical 
                Biosciences (1981), doi: 10.1016/0025-5564(81)90043-2.

                Alias: rf_distance, rf_dist, rf

                Usage:
                phykit robinson_foulds_distance <tree_file_zero> <tree_file_one>

                Options
                =====================================================
                <tree_file_zero>            first argument after 
                                            function name should be
                                            a tree file

                <tree_file_one>             second argument after 
                                            function name should be
                                            a tree file           
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
                f"""\
                {self.help_header} 

                Determines potentially spurious orthologs using branch lengths.

                Identifies potentially spurious sequences and reports
                tips in the phylogeny that could possibly be removed
                from the underlying multiple sequence alignments. PhyKIT
                does so by identifying and reporting long terminal branches
                defined as branches that are equal to or 20 times the median
                length of internal and terminal branches.
                
                Using this method to identify potentially spurious sequences
                was, to my knowledge, first introduced by Shen et al., (2018)
                Cell doi: 10.1016/j.cell.2018.10.023.                

                Alias: spurious_seq, ss

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
                                            (Default: 20)
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

    # TODO: unit test
    def tip_labels(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Prints the labels a phylogeny.

                Alias: labels

                Usage:
                phykit tip_labels <tree>

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        TipLabels(args).run()

    def total_tree_length(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Calculate total tree length, which is a sum of all branches. 

                Alias: tree_len

                Usage:
                phykit total_tree_length <tree>

                Options
                =====================================================
                <tree>                      first argument after 
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

                Calculate treeness statistic for a phylogeny.

                Higher treeness values are thought to be desirable because they
                represent a higher signal-to-noise ratio.

                Treeness describes the proportion of tree distance on internal
                branches. Treeness can be used as a measure of the signal-to-noise
                ratio in a phylogeny. 

                Calculate treeness (also referred to as stemminess) following
                Lanyon, The Auk (1988), doi: 10.1093/auk/105.3.565 and
                Phillips and Penny, Molecular Phylogenetics and Evolution
                (2003), doi: 10.1016/S1055-7903(03)00057-5.

                Alias: tness

                Usage:
                phykit treeness <tree>

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        Treeness(args).run()

    ## Alignment and tree functions
    # TODO: write unit tests
    def saturation(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Calculate saturation for a given tree and alignment.

                Saturation is defined as sequences in multiple sequence
                alignments that have undergone numerous substitutions such
                that the distances between taxa are underestimated.

                Data with no saturation will have a value of 1. Completely
                saturated data will have a value of 0.  

                Saturation is calculated following Philippe et al., PLoS Biology
                (2011), doi: 10.1371/journal.pbio.1000602.

                Alias: sat

                Usage:
                phykit saturation -a <alignment> -t <tree> [-v/--verbose]

                Options
                =====================================================
                -a/--alignment              an alignment file
                
                -t/--tree                   a tree file

                -v/--verbose                print out patristic distances
                                            and uncorrected distances used
                                            to determine saturation
                """
            ),
        )
        parser.add_argument(
            "-a", "--alignment", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument("-v", "--verbose", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(sys.argv[2:])
        Saturation(args).run()

    def treeness_over_rcv(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Calculate treeness/RCV for a given alignment and tree.

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

                Alias: toverr

                Usage:
                phykit treeness_over_rcv -a/--alignment <alignment> -t/--tree <tree>

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
    # TODO: write unit tests
    def create_concatenation_matrix(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Create a concatenated alignment file. This function is 
                used to help in the construction of multi-locus data
                matrices.

                PhyKIT will output three files:
                1) A fasta file with '.fa' appended to the prefix specified
                   with the -p/--prefix parameter.
                2) A partition file ready for input into RAxML or IQ-tree.
                3) An occupancy file that summarizes the taxon occupancy
                   per sequence.

                Alias: create_concat, cc

                Usage:
                phykit create_concatenation_matrix -a <file> -p <string>

                Options
                =====================================================
                -a/--alignment              alignment list file. File
                                            should contain a single
                                            column list of alignment
                                            sequence files to concatenate 
                                            into a single matrix. Provide
                                            path to files relative to
                                            working directory or provide
                                            absolute path.

                -p/--prefix                 prefix of output files
                """
            ),
        )
        parser.add_argument(
            "-a", "--alignment_list",
            type=str,
            help=SUPPRESS
        )
        parser.add_argument(
            "-p", "--prefix",
            type=str,
            help=SUPPRESS
        )
        args = parser.parse_args(sys.argv[2:])
        CreateConcatenationMatrix(args).run()

    # TODO: thread DNA unit tests
    def thread_dna(self):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                Thread DNA sequence onto a protein alignment to create a
                codon-based alignment. 
                
                This function requires input alignments are in fasta format.
                Codon alignments are then printed to stdout. Note, sequences
                should occur in the same order in the protein and nucleotide
                alignment.

                Alias: pal2nal, p2n

                Usage:
                phykit thread_dna -p <file> -n <file> [-s]

                Options
                =====================================================
                -p/--protein                protein alignment file

                -n/--nucleotide             nucleotide alignment file

                -s/--stop                   boolean for whether or not
                                            stop codons should be kept. 
                                            If used, stop codons will 
                                            be removed.
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
            type=str2bool, 
            nargs='?',
            default=True,
            help=SUPPRESS
        )
        args = parser.parse_args(sys.argv[2:])
        DNAThreader(args).run()

def main(argv=None):
    Phykit()
