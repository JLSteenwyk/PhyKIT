#!/usr/bin/env python

import logging
import sys
import textwrap

from .version import __version__

from argparse import (
    ArgumentParser,
    RawTextHelpFormatter,
    SUPPRESS,
    RawDescriptionHelpFormatter,
)

from .services.alignment import (
    AlignmentLength,
    AlignmentLengthNoGaps,
    ColumnScore,
    CreateConcatenationMatrix,
    DNAThreader,
    Faidx,
    GCContent,
    PairwiseIdentity,
    ParsimonyInformative,
    RelativeCompositionVariability,
    RenameFastaEntries,
    SumOfPairsScore,
    VariableSites
)

from .services.tree import (
    BipartitionSupportStats,
    BranchLengthMultiplier,
    CollapseBranches,
    CovaryingEvolutionaryRates,
    DVMC,
    EvolutionaryRate,
    HiddenParalogyCheck,
    InternalBranchStats,
    InternodeLabeler,
    LastCommonAncestorSubtree,
    LBScore,
    MonophylyCheck,
    NearestNeighborInterchange,
    PatristicDistances,
    PolytomyTest,
    PrintTree,
    PruneTree,
    RenameTreeTips,
    RobinsonFouldsDistance,
    RootTree,
    Saturation,
    SpuriousSequence,
    TipLabels,
    TipToTipDistance,
    TotalTreeLength,
    Treeness,
    TreenessOverRCV
)

from .helpers.boolean_argument_parsing import str2bool


logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)

help_header = f"""
                 _____  _           _  _______ _______ 
                |  __ \| |         | |/ /_   _|__   __|
                | |__) | |__  _   _| ' /  | |    | |   
                |  ___/| '_ \| | | |  <   | |    | |   
                | |    | | | | |_| | . \ _| |_   | |   
                |_|    |_| |_|\__, |_|\_\_____|  |_|   
                               __/ |                   
                              |___/   
                            
                Version: {__version__}
                Citation: Steenwyk et al. 2021, Bioinformatics. doi: 10.1093/bioinformatics/btab096
                Documentation link: https://jlsteenwyk.com/PhyKIT
                Publication link: https://academic.oup.com/bioinformatics/article-abstract/37/16/2325/6131675

"""

class Phykit(object):
    help_header = f"""
                 _____  _           _  _______ _______ 
                |  __ \| |         | |/ /_   _|__   __|
                | |__) | |__  _   _| ' /  | |    | |   
                |  ___/| '_ \| | | |  <   | |    | |   
                | |    | | | | |_| | . \ _| |_   | |   
                |_|    |_| |_|\__, |_|\_\_____|  |_|   
                               __/ |                   
                              |___/   
                            
                Version: {__version__}
                Citation: Steenwyk et al. 2021, Bioinformatics. doi: 10.1093/bioinformatics/btab096
                Documentation link: https://jlsteenwyk.com/PhyKIT
                Publication link: https://academic.oup.com/bioinformatics/article-abstract/37/16/2325/6131675

    """
    
    def __init__(self):
        parser = ArgumentParser(
            add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {self.help_header}

                PhyKIT helps process and analyze multiple sequence alignments and phylogenies.

                Generally, all functions are designed to help understand the contents of alignments
                (e.g., gc content or the number of parsimony informative sites) and the shape
                of trees (e.g., treeness, degree of violation of a molecular clock).

                Some help messages indicate that summary statistics are reported (e.g., 
                bipartition_support_stats). Summary statistics include mean, median, 25th percentile,
                75th percentile, minimum, maximum, standard deviation, and variance. These functions
                typically have a verbose option that allows users to get the underlying data
                used to calculate summary statistics. 

                Usage: phykit <command> [optional command arguments]

                Command specific help messages can be viewed by adding a 
                -h/--help argument after the command. For example, to see the
                to see the help message for the command 'treeness', execute
                "phykit treeness -h" or "phykit treeness --help".

                Lastly, each function comes with aliases to save the user some
                key strokes. For example, to get the help message for the 'treeness'
                function, you can type "phykit tness -h". All aliases are specified
                in parentheses after the long form of the function name. 

                Alignment-based commands
                ========================
                alignment_length (alias: aln_len; al)
                    - calculates alignment length
                alignment_length_no_gaps (alias: aln_len_no_gaps; alng)
                    - calculates alignment length after removing sites with gaps
                column_score (alias: cs)
                    - calculate column score between a reference and query alignment
                create_concatenation_matrix (alias: create_concat; cc)
                    - create concatenation matrix from a set of alignments
                faidx (alias: get_entry; ge)
                    - extract query fasta entry from multi-fasta file              
                gc_content (alias: gc)
                    - calculate GC content of a fasta entries or entries thereof
                pairwise_identity (alias: pairwise_id, pi)
                    - calculates average pairwise identify among sequences in
                      an alignment file. This is a proxy for evolutionary rate
                parsimony_informative_sites (alias: pis)
                    - calculates the number and percentage of parsimony
                      informative sites in an alignment
                relative_composition_variability (alias: rel_comp_var, rcv)
                    - calculates relative composition variability in an alignment
                rename_fasta_entries (alias: rename_fasta)
                    - rename entries in a fasta file
                sum_of_pairs_score (alias: sops; sop)
                    - calculate sum-of-pairs score between a reference and query alignment
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
                evolutionary_rate (alias: evo_rate)
                    - reports a tree-based estimation of evolutionary rate for a gene
                hidden_paralogy_check (alias: clan_check)
                    - check for monophyly of specific clades of taxa
                internal_branch_stats (alias: ibs)
                    - calculates summary statistics for internal branch lengths 
                internode_labeler (alias: il)
                    - create labels at internodes in a phylogeny
                last_common_ancestor_subtree (alias: lca_subtree)
                    - get last common ancestor of a set of taxa
                long_branch_score (alias: lb_score; lbs)
                    - calculates lb (long branch) score for taxa in a phylogeny
                monophyly_check (alias: is_monophyletic)
                    - determines if a set of tip names are monophyletic
                nearest_neighbor_interchange (alias: nni)
                    - make nearest neighbor interchange moves on a tree
                patristic_distances (alias: pd)
                    - calculate all pairwise distances between tips in a tree
                polytomy_test (alias: polyt_test; polyt; ptt)
                    - conducts a polytomy test using gene
                      support frequencies
                print_tree (alias: print; pt)
                    - prints ascii tree
                prune_tree (alias: prune)
                    - prune taxa from a phylogeny
                rename_tree_tips (alias: rename_tree; rename_tips)
                    - renames tips in a phylogeny according to a file with
                      the desired new tip names
                robinson_foulds_distance (alias: rf_distance; rf_dist; rf)
                    - calculates Robinson-Foulds distance between two trees
                root_tree (alias: root; rt)
                    - roots tree on user-specified taxa or taxon
                spurious_sequence (alias: spurious_seq; ss)
                    - identifies putatively spurious sequences by identifying
                      branch lengths that are atypically long
                tip_labels (alias: tree_labels; labels; tl)
                    - print leaf names in a phylogeny
                tip_to_tip_distance (alias: t2t_dist; t2t)
                    - calculate tip-to-tip distance in a phylogeny
                total_tree_length (alias: tree_len)
                    - calculates total tree length
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
        try:
            if hasattr(self, args.command):
                getattr(self, args.command)(sys.argv[2:])
            else:
                self.run_alias(args.command, sys.argv[2:])
        except NameError:
            sys.exit()

    ## Aliases
    def run_alias(self, command, argv):
        # version
        if command in ['version', 'v']:
            return self.version()
        # Alignment aliases
        if command in ['aln_len', 'al']:
            return self.alignment_length(argv)
        elif command in ['aln_len_no_gaps', 'alng']:
            return self.alignment_length_no_gaps(argv)
        elif command in 'cs':
            return self.column_score(argv)
        elif command in ['get_entry', 'ge']:
            return self.faidx(argv)
        elif command == 'gc':
            return self.gc_content(argv)
        elif command in ['pairwise_id', 'pi']:
            return self.pairwise_identity(argv)
        elif command == 'pis':
            return self.parsimony_informative_sites(argv)
        elif command in ['rel_comp_var', 'relative_composition_variability']:
            return self.rcv(argv)
        elif command == 'rename_fasta':
            return self.rename_fasta_entries(argv)
        elif command in ['sum_of_pairs_score', 'sops', 'sop']:
            return self.sum_of_pairs_score(argv)
        elif command == 'vs':
            return self.variable_sites(argv)
        # Tree aliases
        elif command == 'bss':
            return self.bipartition_support_stats(argv)
        elif command == 'blm':
            return self.branch_length_multiplier(argv)
        elif command in ['collapse', 'cb']:
            return self.collapse_branches(argv)
        elif command == 'cover':
            return self.covarying_evolutionary_rates(argv)
        elif command == 'degree_of_violation_of_a_molecular_clock':
            return self.dvmc(argv)
        elif command == 'evo_rate':
            return self.evolutionary_rate(argv)
        elif command == 'clan_check':
            return self.hidden_paralogy_check(argv)
        elif command == 'ibs':
            return self.internal_branch_stats(argv)
        elif command == 'il':
            return self.internode_labeler(argv)
        elif command in ['lca_subtree']:
            return self.last_common_ancestor_subtree(argv)
        elif command in ['long_branch_score', 'lbs']:
            return self.lb_score(argv)
        elif command == 'is_monophyletic':
            return self.monophyly_check(argv)
        elif command == 'nni':
            return self.nearest_neighbor_interchange(argv)
        elif command == 'pd':
            return self.patristic_distances(argv)
        elif command in ['polyt_test', 'ptt', 'polyt']:
            return self.polytomy_test(argv)
        elif command in ['print', 'pt']:
            return self.print_tree(argv)
        elif command == 'prune':
            return self.prune_tree(argv)
        elif command in ['rename_tree', 'rename_tips']:
            return self.rename_tree_tips(argv)
        elif command in ['robinson_foulds_distance', 'rf_dist', 'rf']:
            return self.rf_distance(argv)
        elif command in ['root', 'rt']:
            return self.root_tree(argv)
        elif command in ['spurious_seq', 'ss']:
            return self.spurious_sequence(argv)
        elif command in ['labels', 'tree_labels', 'tl']:
            return self.tip_labels(argv)
        elif command in ['t2t_dist', 't2t']:
            return self.tip_to_tip_distance(argv)
        elif command == 'tree_len':
            return self.total_tree_length(argv)
        elif command == 'tness':
            return self.treeness(argv)
        # Alignment- and tree-based aliases
        elif command == 'sat':
            return self.saturation(argv)
        elif command in ['toverr', 'tor']:
            return self.treeness_over_rcv(argv)
        # Helper aliases
        elif command in ['create_concat', 'cc']:
            return self.create_concatenation_matrix(argv)
        elif command in ['pal2nal', 'p2n']:
            return self.thread_dna(argv)
        else:
            print(textwrap.dedent(help_header))
            print("Invalid command option. See help for a complete list of commands and aliases.")
            sys.exit(1)

    ## print version
    def version(self):
        print(textwrap.dedent(
            f"""\
            {self.help_header}
            """
        ))

    ## Alignment functions
    @staticmethod
    def alignment_length(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Length of an input alignment is calculated using this function.

                Longer alignments are associated with strong phylogenetic signal.
                
                Association between alignment length and phylogenetic signal
                was determined by Shen et al., Genome Biology and Evolution (2016),
                doi: 10.1093/gbe/evw179.

                Aliases:
                  alignment_length, aln_len, al
                Command line interfaces:
                  pk_alignment_length, pk_aln_len, pk_al

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
        args = parser.parse_args(argv)
        AlignmentLength(args).run()

    @staticmethod
    def alignment_length_no_gaps(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

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

                Aliases:
                  alignment_length_no_gaps, aln_len_no_gaps, alng
                Command line interfaces: 
                  pk_alignment_length_no_gaps, pk_aln_len_no_gaps, pk_alng

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
        args = parser.parse_args(argv)
        AlignmentLengthNoGaps(args).run()

    @staticmethod
    def column_score(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculates column score.

                Column is an accuracy metric for a multiple alignment relative
                to a reference alignment. It is calculated by summing the correctly
                aligned columns over all columns in an alignment. Thus, values range
                from 0 to 1 and higher values indicate more accurate alignments.

                Aliases:
                  column_score, cs
                Command line interfaces: 
                  pk_column_score, pk_cs

                Usage:
                phykit column_score <fasta> -r/--reference <ref.aln>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be a
                                            query fasta alignment file
                                            to be scored for accuracy

                -r/--reference              reference fasta alignment to 
                                            compare query alignment to
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-r","--reference", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        ColumnScore(args).run()

    @staticmethod
    def faidx(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Extracts sequence entry from fasta file.

                This function works similarly to the faidx function 
                in samtools, but does not requiring an indexing function.

                Aliases:
                  faidx, get_entry; ge
                Command line interfaces: 
                  pk_faidx, pk_get_entry, pk_ge
                  

                Usage:
                phykit faidx <fasta> -e/--entry <fasta entry>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be a
                                            query fasta file

                -e/--entry                  entry name to be extracted
                                            from the inputted fasta file
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-e","--entry", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        Faidx(args).run()

    @staticmethod
    def gc_content(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                
                Calculate GC content of a fasta file.

                GC content is negatively correlated with phylogenetic signal.

                If there are multiple entries, use the -v/--verbose option
                to determine the GC content of each fasta entry separately.

                Association between GC content and phylogenetic signal was
                determined by Shen et al., Genome Biology and Evolution (2016), 
                doi: 10.1093/gbe/evw179.

                Aliases:
                  gc_content, gc
                Command line interfaces: 
                  pk_gc_content, pk_gc

                Usage:
                phykit gc_content <fasta> [-v/--verbose]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file 
            
                -v/--verbose                optional argument to print
                                            the GC content of each fasta
                                            entry
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-v", "--verbose", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        GCContent(args).run()

    @staticmethod
    def pairwise_identity(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate the average pairwise identity among sequences.
                
                Pairwise identities can be used as proxies for the
                evolutionary rate of sequences.

                Pairwise identity is defined as the number of identical
                columns (including gaps) between two aligned sequences divided
                by the number of columns in the alignment. Summary statistics
                are reported unless used with the verbose option in which
                all pairwise identities will be reported.

                An example of pairwise identities being used as a proxy
                for evolutionary rate can be found here: Chen et al. 
                Genome Biology and Evolution (2017), doi: 10.1093/gbe/evx147.

                Aliases:
                  pairwise_identity, pairwise_id, pi
                Command line interfaces:
                  pk_pairwise_identity, pk_pairwise_id, pk_pi

                Usage:
                phykit pairwise_identity <alignment> [-v/--verbose]

                Options
                =====================================================
                <alignment>                 first argument after 
                                            function name should be
                                            an alignment file  

                -v/--verbose                optional argument to print
                                            identity per pair       
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        parser.add_argument("-v", "--verbose", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        PairwiseIdentity(args).run()

    @staticmethod
    def parsimony_informative_sites(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

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
                doi: 10.1093/gbe/evw179 and Steenwyk et al., bioRxiv
                (2020), doi: 10.1101/2020.06.08.140384.

                Aliases: 
                  parsimony_informative_sites, pis
                Command line interfaces:
                  pk_parsimony_informative_sites, pk_pis

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
        args = parser.parse_args(argv)
        ParsimonyInformative(args).run()

    @staticmethod
    def rcv(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate RCV (relative composition variability) for an alignment.

                Lower RCV values are thought to be desirable because they represent
                a lower composition bias in an alignment. Statistically, RCV describes
                the average variability in sequence composition among taxa. 

                RCV is calculated following Phillips and Penny, Molecular Phylogenetics
                and Evolution (2003), doi: 10.1016/S1055-7903(03)00057-5.

                Aliases: 
                  relative_composition_variability, rel_comp_var, rcv
                Command line interfaces:
                  pk_relative_composition_variability, pk_rel_comp_var, pk_rcv

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
        args = parser.parse_args(argv)
        RelativeCompositionVariability(args).run()

    @staticmethod
    def rename_fasta_entries(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Renames fasta entries.

                Renaming fasta entries will follow the scheme of a tab-delimited
                file wherein the first column is the current fasta entry name and
                the second column is the new fasta entry name in the resulting 
                output alignment. 

                Aliases:
                  rename_fasta_entries, rename_fasta
                Command line interfaces: 
                  pk_rename_fasta_entries, pk_rename_fasta

                Usage:
                phykit rename_fasta_entries <fasta> -i/--idmap <idmap>
                    [-o/--output <output_file>]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file

                -i/--idmap                  identifier map of current FASTA
                                            names (col1) and desired FASTA
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
        args = parser.parse_args(argv)
        RenameFastaEntries(args).run()

    @staticmethod
    def sum_of_pairs_score(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculates sum-of-pairs score.

                Sum-of-pairs is an accuracy metric for a multiple alignment relative
                to a reference alignment. It is calculated by summing the correctly
                aligned residue pairs over all pairs of sequences. Thus, values range
                from 0 to 1 and higher values indicate more accurate alignments.

                Aliases:
                  sum_of_pairs_score, sops, sop 
                Command line interfaces: 
                  pk_sum_of_pairs_score, pk_sops, pk_sop

                Usage:
                phykit sum_of_pairs_score <fasta> -r/--reference <ref.aln>

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be a
                                            query fasta alignment file
                                            to be scored for accuracy

                -r/--reference              reference fasta alignment to 
                                            compare query alignment to
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-r","--reference", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        SumOfPairsScore(args).run()

    @staticmethod
    def variable_sites(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

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

                Aliases:
                  variable_sites, vs
                Command line interfaces: 
                  pk_variable_sites, pk_vs

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
        args = parser.parse_args(argv)
        VariableSites(args).run()


    ## Tree functions
    @staticmethod
    def bipartition_support_stats(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                Calculate summary statistics for bipartition support.

                High bipartition support values are thought to be desirable because
                they are indicative of greater certainty in tree topology.

                To obtain all bipartition support values, use the -v/--verbose option.

                Aliases:
                  bipartition_support_stats, bss
                Command line interfaces:
                  pk_bipartition_support_stats, pk_bss

                Usage:
                phykit bipartition_support_stats <tree> [-v/--verbose]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file 
            
                -v/--verbose                optional argument to print
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
        args = parser.parse_args(argv)
        BipartitionSupportStats(args).run()

    @staticmethod
    def branch_length_multiplier(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header} 

                Multiply branch lengths in a phylogeny by a given factor.
                
                This can help modify reference trees when conducting simulations
                or other analyses.              

                Alias:
                  branch_length_multiplier, blm
                Command line interfaces:
                  pk_branch_length_multiplier, pk_blm

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
        args = parser.parse_args(argv)
        BranchLengthMultiplier(args).run()

    @staticmethod
    def collapse_branches(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header} 

                Collapse branches on a phylogeny according to bipartition support.
                Bipartitions will be collapsed if they are less than the user specified
                value.              

                Aliases:
                  collapse_branches, collapse, cb
                Command line interfaces:
                  pk_collapse_branches, pk_collapse, pk_cb

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
                                            ".collapsed_(support).tre"

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
        args = parser.parse_args(argv)
        CollapseBranches(args).run()

    @staticmethod
    def covarying_evolutionary_rates(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Determine if two genes have a signature of covariation with one another.

                Genes that have covarying evolutionary histories tend to have 
                similar functions and expression levels.

                Input two phylogenies and calculate the correlation among relative 
                evolutionary rates between the two phylogenies. The two input trees 
                do not have to have the same taxa. This function will first prune both
                trees to have the same tips. To transform branch lengths into relative
                rates, PhyKIT uses the putative species tree's branch lengths, which is
                inputted by the user. As recommended by the original method developers,
                outlier branche lengths are removed. Outlier branches have a relative 
                evolutionary rate greater than five.

                PhyKIT reports two tab delimited values:
                col1: correlation coefficient
                col2: p-value

                Method is empirically evaluated by Clark et al., Genome Research
                (2012), doi: 10.1101/gr.132647.111. Normalization method using a 
                species tree follows Sato et al., Bioinformatics (2005), doi: 
                10.1093/bioinformatics/bti564. 


                Aliases:
                  covarying_evolutionary_rates, cover
                Command line interfaces:
                  pk_covarying_evolutionary_rates, pk_cover

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
        args = parser.parse_args(argv)
        CovaryingEvolutionaryRates(args).run()

    @staticmethod
    def dvmc(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate degree of violation of a molecular clock (or DVMC) in a phylogeny.

                Lower DVMC values are thought to be desirable because they are indicative
                of a lower degree of violation in the molecular clock assumption.

                Typically, outgroup taxa are not included in molecular clock analysis. Thus,
                prior to calculating DVMC from a single gene tree, outgroup taxa are pruned
                from the phylogeny. PhyKIT will prune taxa using tip names from a single column
                file, which is specified using the -r/--root file. If the tip name does not
                exist in the input tree, rather than raising an error/warning message, the tip  
                name is skipped. If the user wants to calculate DVMC for phylogenies with incomplete
                taxa representation, this will allow the user to use one <root file> for all trees.
                Lastly, an empty root file can be used if the user does not wish to prune outgroup taxa. 

                Calculate DVMC in a tree following Liu et al., PNAS (2017), doi: 10.1073/pnas.1616744114.

                Aliases:
                  degree_of_violation_of_a_molecular_clock, dvmc
                Command line interfaces:
                  pk_degree_of_violation_of_a_molecular_clock, pk_dvmc

                Usage:
                phykit degree_of_violation_of_a_molecular_clock -t/--tree <tree>
                    -r/--root <root_taxa>

                Options
                =====================================================
                -t/--tree                   input file tree name
            
                -r/--root                   single column file with
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
        args = parser.parse_args(argv)
        DVMC(args).run()

    @staticmethod
    def evolutionary_rate(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                Calculate a tree-based estimation of the evolutionary rate of a gene.

                Evolutionary rate is the total tree length divided by the number
                of terminals.

                Calculate evolutionary rate following Telford et al., Proceedings
                of the Royal Society B (2014).

                Aliases:
                  evolutionary_rate, evo_rate
                Command line interfaces:
                  pk_evolutionary_rate, pk_evo_rate

                Usage:
                phykit evolutionary_rate <tree>

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file 
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        EvolutionaryRate(args).run()

    @staticmethod
    def hidden_paralogy_check(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                Scan tree for evidence of hidden paralogy.

                This analysis can be used to identify hidden paralogy. 
                Specifically, this method will examine if a set of
                well known monophyletic taxa are, in fact, monophyletic.
                If they are not, the evolutionary history of the gene may
                be subject to hidden paralogy. This analysis is typically
                done with single-copy orthologous genes.

                Requires a clade file, which species which monophyletic
                lineages to check for. Multiple monophyletic
                lineages can be specified. Each lineage should
                be specified on a single line and each tip name 
                (or taxon name) should be separated by a space.
                For example, if it is anticipated that tips
                "A", "B", and "C" are monophyletic and "D",
                "E", and "F" are expected to be monophyletic, the
                clade file should be formatted as follows:
                "
                A B C
                D E F
                "
                Tip names not present in the tree will not be considered
                when assessing hidden paralogy.

                The output will have six columns and as many rows
                as clades were specified in the -c file. For example,
                if there were three rows of clades to examine the 
                monophyly of, there will be three rows in the output
                where the first row in the output corresponds to the 
                results of the first row in the clade file.
                col 1: if the clade was or wasn't monophyletic
                col 2: average bipartition support value in the clade of interest
                col 3: maximum bipartition support value in the clade of interest
                col 4: minimum bipartition support value in the clade of interest
                col 5: standard deviation of bipartition support values in the clade of interest
                col 6: tip names of taxa monophyletic with the lineage of interest
                       excluding those that are listed in the taxa_of_interest file

                The concept behind this analysis follows
                Siu-Ting et al., Molecular Biology and Evolution (2019).

                Aliases:
                  hidden_paralogy_check, clan_check
                Command line interfaces:
                  pk_hidden_paralogy_check, pk_clan_check

                Usage:
                phykit hidden_paralogy_check <tree> -c/--clade <clade_file>

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                <clade_file>                clade file that specifies
                                            what monophyletic clades
                                            to expect
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-c", "--clade", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        HiddenParalogyCheck(args).run()

    @staticmethod
    def internal_branch_stats(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate summary statistics for internal branch lengths in a phylogeny.

                Internal branch lengths can be useful for phylogeny diagnostics.

                To obtain all internal branch lengths, use the -v/--verbose option. 

                Aliases:
                  internal_branch_stats, ibs
                Command line interfaces:
                  pk_internal_branch_stats, pk_ibs

                Usage:
                phykit internal_branch_stats <tree> [-v/--verbose]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -v/--verbose                optional argument to print
                                            all internal branch lengths
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
        args = parser.parse_args(argv)
        InternalBranchStats(args).run()

    @staticmethod
    def internode_labeler(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Appends numerical identifiers to bipartitions in place
                of support values. This is helpful for pointing to
                specific internodes in supplementary files or otherwise. 

                Alias:
                  internode_labeler, il
                Command line interfaces: 
                  pk_internode_labeler, pk_il

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
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        InternodeLabeler(args).run()

    @staticmethod
    def last_common_ancestor_subtree(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Obtains subtree from a phylogeny by getting
                the last common ancestor from a list of taxa. 

                Alias:
                  last_common_ancestor_subtree, lca_subtree
                Command line interfaces: 
                  pk_last_common_ancestor_subtree, pk_lca_subtree

                Usage:
                phykit last_common_ancestor_subtree <file> <list_of_taxa> 
                [-o/--output <file>]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                <list_of_taxa>              list of taxa to get the last
                                            common ancestor subtree for

                -o/--output                 optional argument to name 
                                            the outputted tree file
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("list_of_taxa", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        LastCommonAncestorSubtree(args).run()

    @staticmethod
    def lb_score(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate long branch (LB) scores in a phylogeny.

                Lower LB scores are thought to be desirable because
                they are indicative of taxa or trees that likely do
                not have issues with long branch attraction.

                LB score is the mean pairwise patristic distance of
                taxon i compared to all other taxa over the average 
                pairwise patristic distance. 
        
                PhyKIT reports summary statistics. To obtain LB scores
                for each taxa, use the -v/--verbose option. 

                LB scores are calculated following Struck, Evolutionary 
                Bioinformatics (2014), doi: 10.4137/EBO.S14239.

                Aliases:
                  long_branch_score, lb_score, lbs
                Command line interfaces:
                  pk_long_branch_score, pk_lb_score, pk_lbs

                Usage:
                phykit long_branch_score <tree> [-v/--verbose]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -v/--verbose                optional argument to print
                                            all LB score values            
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-v", "--verbose", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        LBScore(args).run()

    @staticmethod
    def monophyly_check(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}
                Check for monophyly of a lineage.

                This analysis can be used to determine if a set of 
                taxa are monophyletic.

                Requires a taxa file, which species which tip names
                are expected to be monophyletic. File format is a
                single column file with tip names. Tip names not
                present in the tree will not be considered when
                examining monophyly.

                The output will have six columns.
                col 1: if the clade was or wasn't monophyletic
                col 2: average bipartition support value in the clade of interest
                col 3: maximum bipartition support value in the clade of interest
                col 4: minimum bipartition support value in the clade of interest
                col 5: standard deviation of bipartition support values in the clade of interest
                col 6: tip names of taxa monophyletic with the lineage of interest
                       excluding those that are listed in the taxa_of_interest file

                Aliases:
                  monophyly_check, is_monophyletic
                Command line interfaces:
                  pk_monophyly_check, pk_is_monophyletic

                Usage:
                phykit monophyly_check <tree> <list_of_taxa>

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                <list_of_taxa>              single column file with
                                            list of tip names to 
                                            examine the monophyly of
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("list_of_taxa", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        MonophylyCheck(args).run()

    @staticmethod
    def nearest_neighbor_interchange(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Generate all nearest neighbor interchange moves for a binary
                rooted tree.

                The output file will also include the original phylogeny.

                Aliases:
                  nearest_neighbor_interchange, nni
                Command line interfaces:
                  pk_nearest_neighbor_interchange, pk_nni

                Usage:
                phykit nearest_neighbor_interchange <tree> [-o/--output <output_file>]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -o/--output                 name of output file that will
                                            contain all trees with the
                                            nearest neighbor interchange
                                            moves. 
                                            Default output will have 
                                            the same name as the input
                                            file but with the suffix 
                                            ".NNIs"      
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        NearestNeighborInterchange(args).run()

    @staticmethod
    def patristic_distances(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate summary statistics among patristic distances in a phylogeny.

                Patristic distances are all tip-to-tip distances in a phylogeny.

                To obtain all patristic distances, use the -v/--verbose option.
                With the -v option, the first column will have two taxon names
                separated by a '-' followed by the patristic distance. Features
                will be tab separated. 

                Aliases:
                  patristic_distances, pd
                Command line interfaces: 
                  pk_patristic_distances, pk_pd

                Usage:
                phykit patristic_distances <tree> [-v/--verbose]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -v/--verbose                optional argument to print
                                            all patristic distances between
                                            taxa            
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-v", "--verbose", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        PatristicDistances(args).run() 

    @staticmethod
    def polytomy_test(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

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
                groups are equally supported. This test is done using gene support
                frequencies.

                Aliases:
                  polytomy_test, polyt_test, polyt, ptt
                Command line interfaces: 
                  pk_polytomy_test, pk_polyt_test, pk_polyt, pk_ptt

                Usage:
                phykit polytomy_test -t/--trees <trees> -g/--groups <groups>

                Options
                =====================================================
                -t/--trees                 single column file with names
                                            of phylogenies to use for
                                            polytomy testing

                -g/--groups                a tab-delimited file with the
                                            grouping designations to test.
                                            Lines starting with comments 
                                            are not considered. Names
                                            of individual taxa should be
                                            separated by a semi-colon  ';'
                                            
                For example, the groups file could look like the following:
                #labels group0  group1  group2
                name_of_test    tip_name_A;tip_name_B   tip_name_C  tip_name_D;tip_name_E
                """
            ),
        )
        parser.add_argument("-t", "--trees", type=str, help=SUPPRESS)
        parser.add_argument("-g", "--groups", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        PolytomyTest(args).run() 

    @staticmethod
    def print_tree(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Print ascii tree of input phylogeny.

                Phylogeny can be printed with or without branch lengths.
                By default, the phylogeny will be printed with branch lengths
                but branch lengths can be removed using the -r/--remove argument.

                Aliases:
                  print_tree, print, pt
                Command line interfaces:
                  pk_print_tree, pk_print, pk_pt

                Usage:
                phykit print_tree <tree> [-r/--remove]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -r/--remove                 optional argument to print
                                            the phylogeny without branch
                                            lengths       
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-r", "--remove", action="store_true", required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        PrintTree(args).run()

    @staticmethod
    def prune_tree(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Prune tips from a phylogeny.

                Provide a single column file with the names of the tips
                in the input phylogeny you would like to prune from the
                tree.

                Aliases: 
                  prune_tree, prune
                Command line interfaces: 
                  pk_prune_tree, pk_prune

                Usage:
                phykit prune_tree <tree> <list_of_taxa> [-o/--output <output_file>]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                <list_of_taxa>              single column file with the
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
        args = parser.parse_args(argv)
        PruneTree(args).run()

    @staticmethod
    def rename_tree_tips(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Renames tips in a phylogeny.

                Renaming tip files will follow the scheme of a tab-delimited
                file wherein the first column is the current tip name and the
                second column is the desired tip name in the resulting 
                phylogeny. 

                Aliases:
                  rename_tree_tips, rename_tree, rename_tips
                Command line interfaces: 
                  pk_rename_tree_tips, pk_rename_tree, pk_rename_tips

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
        args = parser.parse_args(argv)
        RenameTreeTips(args).run()

    @staticmethod
    def rf_distance(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate Robinson-Foulds (RF) distance between two trees.

                Low RF distances reflect greater similarity between two phylogenies. 
                This function prints out two values, the plain RF value and the
                normalized RF value, which are separated by a tab. Normalized RF values
                are calculated by taking the plain RF value and dividing it by 2(n-3)
                where n is the number of tips in the phylogeny. Prior to calculating
                an RF value, PhyKIT will first determine the number of shared tips
                between the two input phylogenies and prune them to a common set of
                tips. Thus, users can input trees with different topologies and 
                infer an RF value among subtrees with shared tips.

                PhyKIT will print out 
                col 1; the plain RF distance and 
                col 2: the normalized RF distance.

                RF distances are calculated following Robinson & Foulds, Mathematical 
                Biosciences (1981), doi: 10.1016/0025-5564(81)90043-2.

                Aliases:
                  robinson_foulds_distance, rf_distance, rf_dist, rf
                Command line interfaces: 
                  pk_robinson_foulds_distance, pk_rf_distance, pk_rf_dist, pk_rf

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
        args = parser.parse_args(argv)
        RobinsonFouldsDistance(args).run()

    @staticmethod
    def root_tree(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Roots phylogeny using user-specified taxa.

                A list of taxa to root the phylogeny on should be
                specified using the -r argument. The root_taxa file
                should be a single-column file with taxa names. The
                outputted file will have the same name as the inputted
                tree file but with the suffix ".rooted".

                Aliases:
                  root_tree, root, rt
                Command line interfaces: 
                  pk_root_tree, pk_root, pk_rt

                Usage:
                phykit root_tree <tree> -r/--root <root_taxa>
                    [-o/--output <output_file>] 

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -r/--root                   single column file with
                                            tip names of root taxa

                -o/--output                 optional argument to write
                                            the rooted tree file to.
                                            Default output will have 
                                            the same name as the input
                                            file but with the suffix 
                                            ".rooted"          
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-r", "--root", type=str, required=True, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        args = parser.parse_args(argv)
        RootTree(args).run()

    @staticmethod
    def spurious_sequence(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header} 

                Determines potentially spurious homologs using branch lengths.

                Identifies potentially spurious sequences and reports
                tips in the phylogeny that could possibly be removed
                from the associated multiple sequence alignment. PhyKIT
                does so by identifying and reporting long terminal branches
                defined as branches that are equal to or 20 times the median
                length of all branches.

                PhyKIT reports the following information
                col1: name of tip that is a putatively spurious sequence
                col2: length of branch leading to putatively spurious sequence
                col3: threshold used to identify putatively spurious sequences
                col4: median branch length in the phylogeny

                If there are no putatively spurious sequences, "None" is reported.
                
                Using this method to identify potentially spurious sequences
                was, to my knowledge, first introduced by Shen et al., (2018)
                Cell doi: 10.1016/j.cell.2018.10.023.                

                Aliases:
                  spurious_sequence, spurious_seq, ss
                Command line interfaces:
                  pk_spurious_sequence, pk_spurious_seq, pk_ss

                Usage:
                phykit spurious_sequence <file> [-f 20]

                Options
                =====================================================
                <file>                      first argument after 
                                            function name should be
                                            an tree file

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
        args = parser.parse_args(argv)
        SpuriousSequence(args).run()

    @staticmethod
    def tip_labels(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Prints the tip labels (or names) a phylogeny.

                Aliases:
                  tip_labels, tree_labels; labels; tl
                Command line interfaces: 
                  pk_tip_labels, pk_tree_labels; pk_labels; pk_tl

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
        args = parser.parse_args(argv)
        TipLabels(args).run()

    @staticmethod
    def tip_to_tip_distance(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate distance between two tips (or leaves) in a phylogeny.

                Distances are in substitutions per site.

                Aliases:
                  tip_to_tip_distance, t2t_dist, t2t
                Command line interfaces: 
                  pk_tip_to_tip_distance, pk_t2t_dist, pk_t2t

                Usage:
                phykit tip_to_tip_distance <tree_file> <tip_1> <tip_2>

                Options
                =====================================================
                <tree_file>                 first argument after 
                                            function name should be
                                            a tree file

                <tip_1>                     second argument after 
                                            function name should be
                                            one of the tip names

                <tip_2>                     third argument after 
                                            function name should be
                                            the second tip name      
                """
            ),
        )
        parser.add_argument("tree_zero", type=str, help=SUPPRESS)
        parser.add_argument("tip_1", type=str, help=SUPPRESS)
        parser.add_argument("tip_2", type=str, help=SUPPRESS)
        args = parser.parse_args(argv)
        TipToTipDistance(args).run()

    @staticmethod
    def total_tree_length(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate total tree length, which is a sum of all branches. 

                Aliases:
                  total_tree_length, tree_len
                Command line interfaces: 
                  pk_total_tree_length, pk_tree_len

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
        args = parser.parse_args(argv)
        TotalTreeLength(args).run()

    @staticmethod
    def treeness(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate treeness statistic for a phylogeny.

                Higher treeness values are thought to be desirable because they
                represent a higher signal-to-noise ratio.

                Treeness describes the proportion of the tree distance found on
                internal branches. Treeness can be used as a measure of the 
                signal-to-noise ratio in a phylogeny. 

                Calculate treeness (also referred to as stemminess) following
                Lanyon, The Auk (1988), doi: 10.1093/auk/105.3.565 and
                Phillips and Penny, Molecular Phylogenetics and Evolution
                (2003), doi: 10.1016/S1055-7903(03)00057-5.

                Aliases:
                  treeness, tness
                Command line interfaces:
                  pk_treeness, pk_tness

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
        args = parser.parse_args(argv)
        Treeness(args).run()

    ## Alignment and tree functions
    @staticmethod
    def saturation(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate saturation for a given tree and alignment.

                Saturation is defined as sequences in multiple sequence
                alignments that have undergone numerous substitutions such
                that the distances between taxa are underestimated.

                Data with no saturation will have a value of 1. Completely
                saturated data will have a value of 0.  

                Saturation is calculated following Philippe et al., PLoS 
                Biology (2011), doi: 10.1371/journal.pbio.1000602.

                Aliases: 
                  saturation, sat
                Command line interfaces:
                  pk_saturation, pk_sat

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
        args = parser.parse_args(argv)
        Saturation(args).run()

    @staticmethod
    def treeness_over_rcv(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate treeness/RCV for a given alignment and tree.

                Higher treeness/RCV values are thought to be desirable because
                they harbor a high signal-to-noise ratio are least susceptible
                to composition bias.

                PhyKIT reports three tab delimited values:
                col1: treeness/RCV
                col2: treeness
                col3: RCV

                Calculate treeness/RCV following Phillips and Penny, Molecular 
                Phylogenetics and Evolution (2003), doi: 10.1016/S1055-7903(03)00057-5.

                Aliases:
                  treeness_over_rcv, toverr, tor
                Command line interfaces:
                  pk_treeness_over_rcv, pk_toverr, pk_tor

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
        args = parser.parse_args(argv)
        TreenessOverRCV(args).run()

    ### Helper commands
    @staticmethod
    def create_concatenation_matrix(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Create a concatenated alignment file. This function is 
                used to help in the construction of multi-locus data
                matrices.

                PhyKIT will output three files:
                1) A fasta file with '.fa' appended to the prefix specified
                   with the -p/--prefix parameter.
                2) A partition file ready for input into RAxML or IQ-tree.
                3) An occupancy file that summarizes the taxon occupancy
                   per sequence.

                Aliases:
                  create_concatenation_matrix, create_concat, cc
                Command line interfaces: 
                  pk_create_concatenation_matrix, pk_create_concat, pk_cc

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
        args = parser.parse_args(argv)
        CreateConcatenationMatrix(args).run()

    @staticmethod
    def thread_dna(argv):
        parser = ArgumentParser(add_help=True,
            usage=SUPPRESS,
            formatter_class=RawDescriptionHelpFormatter,
            description=textwrap.dedent(
                f"""\
                {help_header}

                Thread DNA sequence onto a protein alignment to create a
                codon-based alignment. 
                
                This function requires input alignments are in fasta format.
                Codon alignments are then printed to stdout. Note, sequences
                are assumed to occur in the same order in the protein and 
                nucleotide alignment.

                Aliases:
                  thread_dna, pal2nal, p2n
                Command line interfaces:
                  pk_thread_dna, pk_pal2nal, pk_p2n

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
        args = parser.parse_args(argv)
        DNAThreader(args).run()

def main(argv=None):
    Phykit()

# Alignment-based functions
def alignment_length(argv=None):
    Phykit.alignment_length(sys.argv[1:])

def alignment_length_no_gaps(argv=None):
    Phykit.alignment_length_no_gaps(sys.argv[1:])

def column_score(argv=None):
    Phykit.column_score(sys.argv[1:])

def faidx(argv=None):
    Phykit.faidx(sys.argv[1:])

def gc_content(argv=None):
    Phykit.gc_content(sys.argv[1:])

def pairwise_identity(argv=None):
    Phykit.pairwise_identity(sys.argv[1:])

def parsimony_informative_sites(argv=None):
    Phykit.parsimony_informative_sites(sys.argv[1:])

def rcv(argv=None):
    Phykit.rcv(sys.argv[1:])

def rename_fasta_entries(argv=None):
    Phykit.rename_fasta_entries(sys.argv[1:])

def sum_of_pairs_score(argv=None):
    Phykit.sum_of_pairs_score(sys.argv[1:])

def variable_sites(argv=None):
    Phykit.variable_sites(sys.argv[1:])

# Tree-based functions
def bipartition_support_stats(argv=None):
    Phykit.bipartition_support_stats(sys.argv[1:])

def branch_length_multiplier(argv=None):
    Phykit.branch_length_multiplier(sys.argv[1:])

def collapse_branches(argv=None):
    Phykit.collapse_branches(sys.argv[1:])

def covarying_evolutionary_rates(argv=None):
    Phykit.covarying_evolutionary_rates(sys.argv[1:])

def dvmc(argv=None):
    Phykit.dvmc(sys.argv[1:])

def evolutionary_rate(argv=None):
    Phykit.evolutionary_rate(sys.argv[1:])

def hidden_paralogy_check(argv=None):
    Phykit.hidden_paralogy_check(sys.argv[1:])

def internal_branch_stats(argv=None):
    Phykit.internal_branch_stats(sys.argv[1:])

def internode_labeler(argv=None):
    Phykit.internode_labeler(sys.argv[1:])

def last_common_ancestor_subtree(argv=None):
    Phykit.last_common_ancestor_subtree(sys.argv[1:])

def lb_score(argv=None):
    Phykit.lb_score(sys.argv[1:])

def monophyly_check(argv=None):
    Phykit.monophyly_check(sys.argv[1:])

def nearest_neighbor_interchange(argv=None):
    Phykit.nearest_neighbor_interchange(sys.argv[1:])

def patristic_distances(argv=None):
    Phykit.patristic_distances(sys.argv[1:])

def polytomy_test(argv=None):
    Phykit.polytomy_test(sys.argv[1:])

def print_tree(argv=None):
    Phykit.print_tree(sys.argv[1:])

def prune_tree(argv=None):
    Phykit.prune_tree(sys.argv[1:])

def rename_tree_tips(argv=None):
    Phykit.rename_tree_tips(sys.argv[1:])

def rf_distance(argv=None):
    Phykit.rf_distance(sys.argv[1:])

def root_tree(argv=None):
    Phykit.root_tree(sys.argv[1:])

def spurious_sequence(argv=None):
    Phykit.spurious_sequence(sys.argv[1:])

def tip_labels(argv=None):
    Phykit.tip_labels(sys.argv[1:])

def tip_to_tip_distance(argv=None):
    Phykit.tip_to_tip_distance(sys.argv[1:])

def total_tree_length(argv=None):
    Phykit.total_tree_length(sys.argv[1:])

def treeness(argv=None):
    Phykit.treeness(sys.argv[1:])

# Alignment- and tree-based functions
def saturation(argv=None):
    Phykit.saturation(sys.argv[1:])

def treeness_over_rcv(argv=None):
    Phykit.treeness_over_rcv(sys.argv[1:])

# Helper functions
def create_concatenation_matrix(argv=None):
    Phykit.create_concatenation_matrix(sys.argv[1:])

def thread_dna(argv=None):
    Phykit.thread_dna(sys.argv[1:])

