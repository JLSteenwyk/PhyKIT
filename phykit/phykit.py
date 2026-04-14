#!/usr/bin/env python

import logging
import sys
import textwrap

from .version import __version__

from argparse import (
    ArgumentParser,
    SUPPRESS,
    RawDescriptionHelpFormatter,
)

from .helpers.boolean_argument_parsing import str2bool
from .helpers.plot_config import add_plot_arguments
from .cli_registry import ALIAS_TO_HANDLER
from .service_factories import SERVICE_FACTORIES
from .errors import PhykitUserError

# Expose legacy factory names used by static command handlers in this module.
globals().update(SERVICE_FACTORIES)


logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
logger.addHandler(ch)

BANNER = rf"""
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

# Backward-compatible module alias used by command-specific parser descriptions.
help_header = BANNER


def _new_parser(*, description: str) -> ArgumentParser:
    return ArgumentParser(
        add_help=True,
        usage=SUPPRESS,
        formatter_class=RawDescriptionHelpFormatter,
        description=description,
    )


def _add_json_argument(parser: ArgumentParser) -> None:
    parser.add_argument("--json", action="store_true", required=False, help=SUPPRESS)


def _run_service(parser: ArgumentParser, argv, service_factory) -> None:
    args = parser.parse_args(argv)
    try:
        service_factory(args).run()
    except PhykitUserError as err:
        try:
            for message in err.messages:
                print(message)
        except BrokenPipeError:
            pass
        raise SystemExit(err.code)


class Phykit:
    help_header = BANNER

    def __init__(self):
        parser = _new_parser(
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
                alignment_entropy (alias: aln_entropy; entropy)
                    - calculates site-wise alignment entropy
                alignment_recoding (alias: aln_recoding, recode)
                    - recode alignments using reduced character schemes
                alignment_subsample (alias: aln_subsample; subsample)
                    - randomly subsample genes, partitions, or sites
                dstatistic (alias: dstat; abba_baba)
                    - Patterson's D-statistic (ABBA-BABA test) for
                      detecting introgression/gene flow
                dfoil (alias: dfoil_test)
                    - DFOIL test (Pease & Hahn 2015) for detecting
                      and polarizing introgression in a 5-taxon
                      symmetric phylogeny
                alignment_outlier_taxa (alias: outlier_taxa; aot)
                    - identify potential outlier taxa and why they were flagged
                column_score (alias: cs)
                    - calculate column score between a reference and query alignment
                compositional_bias_per_site (alias: comp_bias_per_site; cbps)
                    - detects site-wise compositional biases in an alignment
                composition_per_taxon (alias: comp_taxon; comp_tax)
                    - calculates sequence composition per taxon
                create_concatenation_matrix (alias: create_concat; cc)
                    - create concatenation matrix from a set of alignments
                evolutionary_rate_per_site (alias: evo_rate_per_site; erps)
                    - estimate evolutionary per site in an alignment
                faidx (alias: get_entry; ge)
                    - extract query fasta entry from multi-fasta file
                gc_content (alias: gc)
                    - calculate GC content of a fasta entries or entries thereof
                mask_alignment (alias: mask_aln; mask)
                    - mask alignment sites based on thresholds
                plot_alignment_qc (alias: plot_qc; paqc)
                    - generate multi-panel alignment quality-control plot
                occupancy_per_taxon (alias: occupancy_taxon; occ_tax)
                    - calculates alignment occupancy per taxon
                pairwise_identity (alias: pairwise_id, pi)
                    - calculates average pairwise identify among sequences in
                      an alignment file. This is a proxy for evolutionary rate
                identity_matrix (alias: id_matrix, seqid)
                    - compute pairwise sequence identity matrix and plot as
                      a clustered heatmap
                parsimony_informative_sites (alias: pis)
                    - calculates the number and percentage of parsimony
                      informative sites in an alignment
                relative_composition_variability (alias: rel_comp_var, rcv)
                    - calculates relative composition variability in an alignment
                relative_composition_variability_taxon (alias: rel_comp_var_taxon, rcvt)
                    - calculates relative composition variability of each taxa in an alignment
                rename_fasta_entries (alias: rename_fasta)
                    - rename entries in a fasta file
                sum_of_pairs_score (alias: sops; sop)
                    - calculate sum-of-pairs score between a reference and query alignment
                thread_dna (alias: pal2nal; p2n)
                    - thread dna sequences over a protein alignment
                phylo_gwas (alias: pgwas)
                    - phylogenetic genome-wide association study
                      (Pease et al. 2016 approach)
                variable_sites (alias: vs)
                    - calculates the number and percentage of variable sites
                      in an alignment
                taxon_groups (alias: tgroups; shared_taxa)
                    - group tree or FASTA files by their taxon set
                occupancy_filter (alias: occ_filter; filter_occupancy)
                    - filter alignments/trees by cross-file taxon occupancy

                Tree-based commands
                ===================
                independent_contrasts (alias: pic; phylo_contrasts)
                    - Felsenstein's phylogenetically independent contrasts
                parsimony_score (alias: parsimony; pars)
                    - Fitch parsimony score of a tree given an alignment
                character_map (alias: charmap; synapomorphy_map)
                    - map discrete character changes on a tree (Fitch parsimony)
                ancestral_state_reconstruction (alias: asr; anc_recon)
                    - estimate ancestral states for continuous traits using
                      ML (fast or VCV-based) with optional contMap plot
                concordance_asr (alias: conc_asr; casr)
                    - concordance-aware ancestral state reconstruction
                      incorporating gene tree discordance
                chronogram (alias: chrono; time_tree)
                    - plot a time-calibrated tree with geological timescale
                dtt (alias: disparity_through_time)
                    - disparity through time analysis (Harmon et al. 2003)
                bipartition_support_stats (alias: bss)
                    - calculates summary statistics for bipartition support
                branch_length_multiplier (alias: blm)
                    - multiply all branch lengths by a specified factor
                collapse_branches (alias: collapse; cb)
                    - collapses branches according to bipartition support
                covarying_evolutionary_rates (alias: cover)
                    - calculates correlation in the evolutionary rate of two trees
                consensus_network (alias: consnet; splitnet; splits_network)
                    - extract bipartition splits from gene trees and visualize
                      conflicting phylogenetic signal as a splits network
                neighbor_net (alias: nnet)
                    - construct a NeighborNet phylogenetic network from
                      pairwise distances (alignment or CSV distance matrix)
                consensus_tree (alias: consensus; ctree)
                    - infer strict or majority-rule consensus from a tree collection
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
                faiths_pd (alias: faith_pd; fpd; phylo_diversity)
                    - calculate Faith's phylogenetic diversity for a
                      community of tips
                patristic_distances (alias: pd)
                    - calculate all pairwise distances between tips in a tree
                phylogenetic_signal (alias: phylo_signal; ps)
                    - calculate phylogenetic signal (Blomberg's K or Pagel's
                      lambda) for continuous trait data
                trait_correlation (alias: trait_corr; phylo_corr)
                    - compute phylogenetic correlations between all pairs
                      of traits and display as a heatmap
                phylo_impute (alias: impute; phylo_imp)
                    - impute missing trait values using phylogenetic
                      relationships and between-trait correlations
                phylogenetic_ordination (alias: phylo_ordination; ordination; ord;
                    phylo_pca; phyl_pca; ppca; phylo_dimreduce; dimreduce; pdr)
                    - phylogenetic ordination (PCA, t-SNE, or UMAP) on
                      continuous multi-trait data
                phylo_heatmap (alias: pheatmap; ph)
                    - phylogeny alongside a heatmap of numeric trait values
                phylomorphospace (alias: phylomorpho; phmo)
                    - plot raw traits with phylogeny overlaid via ancestral
                      reconstruction
                phylogenetic_regression (alias: phylo_regression; pgls)
                    - fit phylogenetic generalized least squares (PGLS) regression
                phylogenetic_glm (alias: phylo_glm; pglm)
                    - fit phylogenetic GLM for binary (logistic) or count (Poisson) data
                phylo_anova (alias: panova; phylo_manova; pmanova)
                    - phylogenetic ANOVA / MANOVA using RRPP (Adams & Collyer 2018)
                phylo_path (alias: ppath; phylopath)
                    - phylogenetic path analysis (von Hardenberg & Gonzalez-Voyer 2013)
                phylo_logistic (alias: phylo_logreg; plogreg)
                    - fit phylogenetic logistic regression (Ives & Garland 2010)
                stochastic_character_map (alias: simmap; scm)
                    - stochastic character mapping (SIMMAP) of discrete traits
                simmap_summary (alias: smsummary; describe_simmap)
                    - per-branch SIMMAP summary with node posteriors
                cont_map (alias: contmap; cmap)
                    - continuous trait map (contMap) visualization on a phylogeny
                trait_rate_map (alias: rate_map; branch_rates)
                    - per-branch evolutionary rate map for a continuous trait
                density_map (alias: densitymap; dmap)
                    - density map of posterior state probabilities on a phylogeny
                cophylo (alias: tanglegram; tangle)
                    - cophylogenetic (tanglegram) plot of two phylogenies
                phenogram (alias: traitgram; tg)
                    - phenogram (traitgram) visualization of trait evolution
                rate_heterogeneity (alias: brownie; rh)
                    - test for rate heterogeneity across tree regimes
                      using multi-rate Brownian motion (O'Meara et al. 2006)
                fit_continuous (alias: fitcontinuous; fc)
                    - compare continuous trait evolution models
                      (BM, OU, EB, Lambda, Delta, Kappa, White)
                fit_discrete (alias: fitdiscrete; fd)
                    - compare discrete trait evolution models
                      (ER, SYM, ARD)
                ouwie (alias: fit_ouwie; multi_regime_ou)
                    - fit multi-regime OU models (BM1, BMS, OU1,
                      OUM, OUMV, OUMA, OUMVA; Beaulieu et al. 2012)
                ou_shift_detection (alias: ou_shifts; l1ou; detect_shifts)
                    - automatic OU shift detection using LASSO
                      (Khabbazian et al. 2016)
                quartet_network (alias: quartet_net; qnet; nanuq)
                    - quartet-based network inference (NANUQ-style)
                      distinguishing ILS from hybridization
                quartet_pie (alias: qpie; quartet_pie_chart)
                    - phylogram with quartet concordance pie charts
                      at internal nodes
                ltt (alias: gamma_stat; gamma)
                    - lineage-through-time plot and Pybus & Harvey
                      gamma statistic for diversification rate testing
                network_signal (alias: netsig; net_signal)
                    - phylogenetic signal on a network
                threshold_model (alias: threshold; thresh; threshbayes)
                    - Felsenstein (2012) threshold model for
                      trait correlation via MCMC
                polytomy_test (alias: polyt_test; polyt; ptt)
                    - conducts a polytomy test using gene
                      support frequencies
                print_tree (alias: print; pt)
                    - prints ascii tree
                prune_tree (alias: prune)
                    - prune taxa from a phylogeny
                subtree_prune_regraft (alias: spr)
                    - generate all SPR rearrangements for a specified subtree
                transfer_annotations (alias: transfer_annot; annotate_tree)
                    - transfer node annotations between trees (e.g., wASTRAL to RAxML/IQ-TREE)
                relative_rate_test (alias: rrt; tajima_rrt)
                    - Tajima's relative rate test for equal evolutionary
                      rates between two ingroup lineages
                rename_tree_tips (alias: rename_tree; rename_tips)
                    - renames tips in a phylogeny according to a file with
                      the desired new tip names
                kuhner_felsenstein_distance (alias: kf_distance; kf_dist; kf)
                    - calculates Kuhner-Felsenstein (branch score) distance
                      between two trees
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
                tip_to_tip_node_distance (alias: t2t_node_dist; t2t_nd)
                    - calculate tip-to-tip node distance in a phylogeny
                total_tree_length (alias: tree_len)
                    - calculates total tree length
                treeness (alias: tness)
                    - reports treeness or stemminess, a measure of signal-to-
                      noise ratio in a phylogeny
                hybridization (alias: hybrid; reticulation)
                    - estimate minimum reticulation events and localize
                      hybridization on a species tree using gene tree
                      discordance asymmetry
                spectral_discordance (alias: spec_disc; sd)
                    - PCA + spectral clustering of gene tree space via
                      bipartition decomposition
                tree_space (alias: tspace; tree_landscape)
                    - visualize gene tree topology space via MDS, t-SNE,
                      or UMAP on pairwise tree distance matrices

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
        parser.add_argument("command", nargs="?", default=None, help=SUPPRESS)
        args = parser.parse_args(sys.argv[1:2])

        if args.command is None:
            parser.print_help()
            sys.exit(0)

        # if command is part of the possible commands (i.e., the long form
        # commands, run). Otherwise, assume it is an alias and look to the
        # run_alias function
        try:
            if hasattr(self, args.command):
                getattr(self, args.command)(sys.argv[2:])
            else:
                self.run_alias(args.command, sys.argv[2:])
        except SystemExit:
            # Re-raise SystemExit as-is to preserve exit code
            raise
        except NameError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(2)

    ## Aliases
    def run_alias(self, command, argv):
        handler_name = ALIAS_TO_HANDLER.get(command)
        if handler_name:
            handler = getattr(self, handler_name)
            if handler_name == "version":
                return handler()
            return handler(argv)

        print(textwrap.dedent(BANNER))
        print(
            "Invalid command option. See help for a complete list of commands and aliases."
        )
        sys.exit(1)

    ## print version
    def version(self):
        print(
            textwrap.dedent(
                f"""\
            {self.help_header}
            """
            )
        )

    ## Alignment functions
    @staticmethod
    def alignment_length(argv):
        parser = _new_parser(
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
                phykit alignment_length <alignment> [--json]

                Options
                =====================================================
                <alignment>                 first argument after 
                                            function name should be
                                            an alignment file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, AlignmentLength)

    @staticmethod
    def alignment_length_no_gaps(argv):
        parser = _new_parser(
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
                phykit alignment_length_no_gaps <alignment> [--json]

                Options
                =====================================================
                <alignment>                 first argument after 
                                            function name should be
                                            an alignment file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, AlignmentLengthNoGaps)

    @staticmethod
    def alignment_entropy(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate alignment entropy.

                Site-wise entropy is calculated using Shannon entropy.
                By default, this function prints the mean site entropy.
                With the -v/--verbose option, entropy is printed for each
                site in the alignment.

                Aliases:
                  alignment_entropy, aln_entropy, entropy
                Command line interfaces:
                  pk_alignment_entropy, pk_aln_entropy, pk_entropy

                Usage:
                phykit alignment_entropy <alignment> [-v/--verbose] [--plot] [--plot-output <path>]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                <alignment>                 first argument after
                                            function name should be
                                            an alignment file

                -v/--verbose                optional argument to print
                                            entropy for each site

                --plot                      optional argument to save a
                                            per-site alignment entropy plot

                --plot-output               output path for plot
                                            (default: alignment_entropy_plot.png)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument("--plot", action="store_true", required=False, help=SUPPRESS)
        parser.add_argument(
            "--plot-output",
            type=str,
            default="alignment_entropy_plot.png",
            required=False,
            help=SUPPRESS,
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, AlignmentEntropy)

    @staticmethod
    def alignment_recoding(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Recode alignments using reduced character states.

                Alignments can be recoded using established or
                custom recoding schemes. Recoding schemes are
                specified using the -c/--code argument. Custom
                recoding schemes can be used and should be formatted
                as a two column file wherein the first column is the
                recoded character and the second column is the character
                in the alignment.
                
                Aliases:
                  alignment_recoding, aln_recoding, recode
                Command line interfaces: 
                  pk_alignment_recoding, pk_aln_recoding, pk_recode

                Usage:
                phykit alignment_recoding <fasta> -c/--code <code> [--json]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file

                -c/--code                   recoding scheme to use

                --json                      optional argument to output
                                            results as JSON

                Codes for which recoding scheme to use
                =====================================================
                RY-nucleotide
                    R = purines (i.e., A and G) 
                    Y = pyrimidines (i.e., T and C)
                
                SandR-6
                    0 = A, P, S, and T
                    1 = D, E, N, and G
                    2 = Q, K, and R
                    3 = M, I, V, and L
                    4 = W and C
                    5 = F, Y, and H

                KGB-6
                    0 = A, G, P, and S
                    1 = D, E, N, Q, H, K, R, and T
                    2 = M, I, and L
                    3 = W
                    4 = F and Y
                    5 = C and V

                Dayhoff-6
                    0 = A, G, P, S, and T
                    1 = D, E, N, and Q
                    2 = H, K, and R
                    3 = I, L, M, and V
                    4 = F, W, and Y
                    5 = C

                Dayhoff-9
                    0 = D, E, H, N, and Q
                    1 = I, L, M, and V
                    2 = F and Y
                    3 = A, S, and T
                    4 = K and R
                    5 = G
                    6 = P
                    7 = C
                    8 = W

                Dayhoff-12
                    0 = D, E, and Q
                    1 = M, L, I, and V
                    2 = F and Y
                    3 = K, H, and R
                    4 = G
                    5 = A
                    6 = P
                    7 = S
                    8 = T
                    9 = N
                    A = W
                    B = C
                
                Dayhoff-15
                    0 = D, E, and Q
                    1 = M and L
                    2 = I and V
                    3 = F and Y
                    4 = G
                    5 = A
                    6 = P
                    7 = S
                    8 = T
                    9 = N
                    A = K
                    B = H
                    C = R
                    D = W
                    E = C
                
                Dayhoff-18
                    0 = F and Y
                    1 = M and L
                    2 = I
                    3 = V
                    4 = G
                    5 = A
                    6 = P
                    7 = S
                    8 = T
                    9 = D
                    A = E
                    B = Q
                    C = N
                    D = K
                    E = H
                    F = R
                    G = W
                    H = C
                """  # noqa
            ),
        )

        parser.add_argument("alignment", type=str, help=SUPPRESS)
        parser.add_argument("-c", "--code", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, AlignmentRecoding)

    @staticmethod
    def alignment_outlier_taxa(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Identify potential outlier taxa in an alignment.

                The following features are evaluated for each taxon:
                1) gap_rate: fraction of gap/ambiguous characters
                2) occupancy: fraction of valid characters
                3) composition_distance: Euclidean distance from median composition profile
                4) long_branch_proxy: mean pairwise sequence distance to other taxa
                5) rcvt: relative composition variability per taxon
                6) entropy_burden: mean site entropy over valid positions

                Taxa are flagged when one or more feature values exceed
                their feature-specific outlier threshold (or drop below the
                occupancy threshold).

                Aliases:
                  alignment_outlier_taxa, outlier_taxa, aot
                Command line interfaces:
                  pk_alignment_outlier_taxa, pk_outlier_taxa, pk_aot

                Usage:
                phykit alignment_outlier_taxa <alignment>
                  [--gap-z <float>] [--composition-z <float>] [--distance-z <float>]
                  [--rcvt-z <float>] [--occupancy-z <float>] [--entropy-z <float>]
                  [--json]

                Options
                =====================================================
                <alignment>                 first argument after
                                            function name should be
                                            an alignment file

                --gap-z                     z-threshold for gap_rate
                                            outlier detection

                --composition-z             z-threshold for
                                            composition_distance outlier
                                            detection

                --distance-z                z-threshold for
                                            long_branch_proxy outlier
                                            detection

                --rcvt-z                    z-threshold for
                                            rcvt outlier detection

                --occupancy-z               z-threshold for
                                            low-occupancy outlier detection

                --entropy-z                 z-threshold for
                                            entropy_burden outlier detection

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        parser.add_argument("--gap-z", type=float, default=3.0, required=False, help=SUPPRESS)
        parser.add_argument(
            "--composition-z", type=float, default=3.0, required=False, help=SUPPRESS
        )
        parser.add_argument(
            "--distance-z", type=float, default=3.0, required=False, help=SUPPRESS
        )
        parser.add_argument("--rcvt-z", type=float, default=3.0, required=False, help=SUPPRESS)
        parser.add_argument("--occupancy-z", type=float, default=3.0, required=False, help=SUPPRESS)
        parser.add_argument("--entropy-z", type=float, default=3.0, required=False, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, AlignmentOutlierTaxa)

    @staticmethod
    def column_score(argv):
        parser = _new_parser(
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
                phykit column_score <fasta> -r/--reference <ref.aln> [--json]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be a
                                            query fasta alignment file
                                            to be scored for accuracy

                -r/--reference              reference fasta alignment to 
                                            compare query alignment to

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-r", "--reference", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, ColumnScore)

    @staticmethod
    def compositional_bias_per_site(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculates compositional bias per site in an alignment.

                Site-wise chi-squared tests are conducted in an alignment to
                detect compositional biases. PhyKIT outputs four columns:
                col 1: index in alignment
                col 2: chi-squared statistic (higher values indicate greater bias)
                col 3: multi-test corrected p-value (Benjamini-Hochberg false discovery rate procedure)
                col 4: uncorrected p-value

                Aliases:
                  compositional_bias_per_site; comp_bias_per_site; cbps
                Command line interfaces:
                  pk_compositional_bias_per_site; pk_compositional_bias_per_site; pk_cbps

                Usage:
                phykit compositional_bias_per_site <alignment> [--plot] [--plot-output <path>]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                <alignment>                 first argument after the
                                            function name should be a
                                            fasta alignment file

                --plot                      optional argument to save a
                                            Manhattan-style plot of
                                            compositional bias per site

                --plot-output               output path for plot
                                            (default: compositional_bias_per_site_plot.png)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        parser.add_argument("--plot", action="store_true", required=False, help=SUPPRESS)
        parser.add_argument(
            "--plot-output",
            type=str,
            default="compositional_bias_per_site_plot.png",
            required=False,
            help=SUPPRESS,
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, CompositionalBiasPerSite)

    @staticmethod
    def composition_per_taxon(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate sequence composition per taxon in an alignment.

                Composition is reported as symbol:frequency values for each
                taxon, where frequencies are calculated from valid
                (non-gap/non-ambiguous) characters.

                Aliases:
                  composition_per_taxon, comp_taxon, comp_tax
                Command line interfaces:
                  pk_composition_per_taxon, pk_comp_taxon, pk_comp_tax

                Usage:
                phykit composition_per_taxon <alignment> [--json]

                Options
                =====================================================
                <alignment>                 first argument after
                                            function name should be
                                            an alignment file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, CompositionPerTaxon)

    @staticmethod
    def evolutionary_rate_per_site(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Estimate evolutionary rate per site.

                Evolutionary rate per site is one minus the sum of squared
                frequency of different characters at a given site. Values
                may range from 0 (slow evolving; no diversity at the given
                site) to 1 (fast evolving; all characters appear only once).

                PhyKIT prints out two columns of information.
                col 1: site in alignment
                col 2: estimated evolutionary rate

                Aliases:
                  evolutionary_rate_per_site; evo_rate_per_site; erps
                Command line interfaces:
                  pk_evolutionary_rate_per_site; pk_evo_rate_per_site; pk_erps
        

                Usage:
                phykit evo_rate_per_site <fasta> [--plot] [--plot-output <path>]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                <fasta>                     first argument after
                                            function name should be a
                                            query fasta file

                --plot                      optional argument to save a
                                            per-site evolutionary-rate plot

                --plot-output               output path for plot
                                            (default: evolutionary_rate_per_site_plot.png)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        parser.add_argument("--plot", action="store_true", required=False, help=SUPPRESS)
        parser.add_argument(
            "--plot-output",
            type=str,
            default="evolutionary_rate_per_site_plot.png",
            required=False,
            help=SUPPRESS,
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, EvolutionaryRatePerSite)

    @staticmethod
    def faidx(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Extracts sequence entry from fasta file.

                This function works similarly to the faidx function 
                in samtools, but does not requiring an indexing step.

                To obtain multiple entries, input multiple entries separated
                by a comma (,). For example, if you want entries 
                named "seq_0" and "seq_1", the string "seq_0,seq_1"
                should be associated with the -e argument.

                Aliases:
                  faidx, get_entry; ge
                Command line interfaces: 
                  pk_faidx, pk_get_entry, pk_ge
                  

                Usage:
                phykit faidx <fasta> -e/--entry <fasta entry> [--json]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be a
                                            query fasta file

                -e/--entry                  entry name to be extracted
                                            from the inputted fasta file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-e", "--entry", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, Faidx)

    @staticmethod
    def gc_content(argv):
        parser = _new_parser(
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
                phykit gc_content <fasta> [-v/--verbose] [--json]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be
                                            a fasta file 
            
                -v/--verbose                optional argument to print
                                            the GC content of each fasta
                                            entry

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        _add_json_argument(parser)
        _run_service(parser, argv, GCContent)

    @staticmethod
    def mask_alignment(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Mask alignment sites based on threshold criteria.

                Sites are retained when they pass all active thresholds:
                maximum gap fraction, minimum occupancy, and maximum
                site entropy.

                Aliases:
                  mask_alignment, mask_aln, mask
                Command line interfaces:
                  pk_mask_alignment, pk_mask_aln, pk_mask

                Usage:
                phykit mask_alignment <alignment> [-g/--max_gap <float>]
                  [-o/--min_occupancy <float>] [-e/--max_entropy <float>] [--json]

                Options
                =====================================================
                <alignment>                 first argument after
                                            function name should be
                                            an alignment file

                -g/--max_gap                maximum allowed fraction of
                                            missing/invalid characters
                                            at a site (default: 1.0)

                -o/--min_occupancy          minimum required occupancy
                                            at a site (default: 0.0)

                -e/--max_entropy            maximum allowed site entropy
                                            (default: no filter)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        parser.add_argument("-g", "--max_gap", type=float, default=1.0, help=SUPPRESS)
        parser.add_argument(
            "-o", "--min_occupancy", type=float, default=0.0, help=SUPPRESS
        )
        parser.add_argument(
            "-e", "--max_entropy", type=float, required=False, help=SUPPRESS
        )
        _add_json_argument(parser)
        _run_service(parser, argv, MaskAlignment)

    @staticmethod
    def plot_alignment_qc(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Generate a multi-panel alignment quality-control plot.

                The figure summarizes per-taxon occupancy and gap rates,
                composition-distance versus long-branch proxy, and counts
                of feature-based outlier flags.

                Features evaluated:
                gap_rate, occupancy, composition_distance, long_branch_proxy,
                rcvt, entropy_burden

                Aliases:
                  plot_alignment_qc, plot_qc, paqc
                Command line interfaces:
                  pk_plot_alignment_qc, pk_plot_qc, pk_paqc

                Usage:
                phykit plot_alignment_qc <alignment> [-o/--output <path>]
                  [--width <float>] [--height <float>] [--dpi <int>]
                  [--gap-z <float>] [--composition-z <float>] [--distance-z <float>]
                  [--rcvt-z <float>] [--occupancy-z <float>] [--entropy-z <float>]
                  [--json]

                Options
                =====================================================
                <alignment>                 first argument after
                                            function name should be
                                            an alignment file

                -o/--output                output image path
                                            (default: alignment_qc.png)

                --width                    figure width in inches
                                            (default: 14.0)

                --height                   figure height in inches
                                            (default: 10.0)

                --dpi                      image DPI (default: 300)

                --gap-z                    z-threshold for gap_rate outliers
                --composition-z            z-threshold for composition_distance outliers
                --distance-z               z-threshold for long_branch_proxy outliers
                --rcvt-z                   z-threshold for rcvt outliers
                --occupancy-z              z-threshold for low occupancy outliers
                --entropy-z                z-threshold for entropy_burden outliers

                --json                     optional argument to output
                                            metadata as JSON
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, default="alignment_qc.png", required=False, help=SUPPRESS)
        parser.add_argument("--width", type=float, default=14.0, required=False, help=SUPPRESS)
        parser.add_argument("--height", type=float, default=10.0, required=False, help=SUPPRESS)
        parser.add_argument("--dpi", type=int, default=300, required=False, help=SUPPRESS)
        parser.add_argument("--gap-z", type=float, default=3.0, required=False, help=SUPPRESS)
        parser.add_argument("--composition-z", type=float, default=3.0, required=False, help=SUPPRESS)
        parser.add_argument("--distance-z", type=float, default=3.0, required=False, help=SUPPRESS)
        parser.add_argument("--rcvt-z", type=float, default=3.0, required=False, help=SUPPRESS)
        parser.add_argument("--occupancy-z", type=float, default=3.0, required=False, help=SUPPRESS)
        parser.add_argument("--entropy-z", type=float, default=3.0, required=False, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, PlotAlignmentQC)

    @staticmethod
    def occupancy_per_taxon(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate occupancy per taxon in an alignment.

                Occupancy is the fraction of valid (non-gap/non-ambiguous)
                characters for each taxon.

                Aliases:
                  occupancy_per_taxon, occupancy_taxon, occ_tax
                Command line interfaces:
                  pk_occupancy_per_taxon, pk_occupancy_taxon, pk_occ_tax

                Usage:
                phykit occupancy_per_taxon <alignment> [--json]

                Options
                =====================================================
                <alignment>                 first argument after
                                            function name should be
                                            an alignment file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, OccupancyPerTaxon)

    @staticmethod
    def pairwise_identity(argv):
        parser = _new_parser(
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
                phykit pairwise_identity <alignment> [-v/--verbose] [-e/--exclude_gaps] [--plot] [--plot-output <file>]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                <alignment>                 first argument after 
                                            function name should be
                                            an alignment file  

                -v/--verbose                optional argument to print
                                            identity per pair

                -e/--exclude_gaps           if a site has a gap, ignore it

                --plot                      optional argument to save a clustered
                                            pairwise-identity heatmap

                --plot-output               output path for heatmap
                                            (default: pairwise_identity_heatmap.png)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument(
            "-e", "--exclude_gaps", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument("--plot", action="store_true", required=False, help=SUPPRESS)
        parser.add_argument(
            "--plot-output",
            type=str,
            required=False,
            default="pairwise_identity_heatmap.png",
            help=SUPPRESS,
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, PairwiseIdentity)

    @staticmethod
    def identity_matrix(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Compute a pairwise sequence identity matrix from an
                alignment and plot it as a clustered heatmap.

                For each pair of taxa, identity is defined as the
                fraction of non-gap, non-ambiguous columns that are
                identical. Gaps, '?', 'N', 'X', and '*' in either
                sequence cause a column to be skipped.

                The matrix can be displayed as identity (default) or
                p-distance (1 - identity) using --metric.

                Ordering can be by hierarchical clustering (default),
                tree tip order (--sort tree --tree <file>), or
                alphabetical (--sort alpha).

                Aliases:
                  identity_matrix, id_matrix, seqid
                Command line interfaces:
                  pk_identity_matrix, pk_id_matrix, pk_seqid

                Usage:
                phykit identity_matrix -a <alignment> -o <output>
                  [--metric identity|p-distance]
                  [--tree <file>] [--sort alpha|cluster|tree]
                  [--partition <file>] [--json]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]

                Options
                =====================================================
                -a/--alignment              alignment file (FASTA or
                                            other supported format)

                -o/--output                 output figure path
                                            (.png, .pdf, .svg)

                --metric                    'identity' (fraction matching)
                                            or 'p-distance' (1 - identity)
                                            (default: identity)

                --tree                      tree file for tree-guided
                                            ordering (Newick format)

                --sort                      ordering method: 'cluster'
                                            (hierarchical), 'tree'
                                            (requires --tree), or 'alpha'
                                            (alphabetical)
                                            (default: cluster)

                --partition                 RAxML-style partition file
                                            (reserved for future use)

                --json                      output structured JSON
                                            instead of plain text
                """
            ),
        )
        parser.add_argument(
            "-a", "--alignment", type=str, required=True, help=SUPPRESS
        )
        parser.add_argument(
            "-o", "--output", type=str, required=True, help=SUPPRESS
        )
        parser.add_argument(
            "--metric", type=str, default="identity",
            choices=["identity", "p-distance"],
            required=False, help=SUPPRESS,
        )
        parser.add_argument(
            "--tree", type=str, default=None,
            required=False, help=SUPPRESS,
        )
        parser.add_argument(
            "--sort", type=str, default="cluster",
            choices=["alpha", "cluster", "tree"],
            required=False, help=SUPPRESS,
        )
        parser.add_argument(
            "--partition", type=str, default=None,
            required=False, help=SUPPRESS,
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, IdentityMatrix)

    @staticmethod
    def parsimony_informative_sites(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate the number and percentage of parsimony
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
                phykit parsimony_informative_sites <alignment> [--json]

                Options
                =====================================================
                <alignment>                 first argument after 
                                            function name should be
                                            an alignment file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, ParsimonyInformative)

    @staticmethod
    def rcv(argv):
        parser = _new_parser(
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
                phykit relative_composition_variability <alignment> [--json]

                Options
                =====================================================
                <alignment>                 first argument after 
                                            function name should be
                                            an alignment file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, RelativeCompositionVariability)

    @staticmethod
    def rcvt(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate RCVT (relative composition variability, taxon) for an alignment.

                RCVT is the relative composition variability metric for individual taxa.
                This facilitates identifying specific taxa that may have compositional
                biases. Lower RCVT values are more desirable because they indicate
                a lower composition bias for a given taxon in an alignment.

                Aliases: 
                  relative_composition_variability_taxon, rel_comp_var_taxon, rcvt
                Command line interfaces:
                  pk_relative_composition_variability_taxon, pk_rel_comp_var_taxon, pk_rcvt

                Usage:
                phykit relative_composition_variability_taxon <alignment>
                  [--plot] [--plot-output <path>]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                <alignment>                 first argument after 
                                            function name should be
                                            an alignment file          

                --plot                      optional argument to output
                                            an RCVT barplot (PNG)

                --plot-output               output path for plot
                                            (default: rcvt_plot.png)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        parser.add_argument("--plot", action="store_true", required=False, help=SUPPRESS)
        parser.add_argument("--plot-output", type=str, default="rcvt_plot.png", required=False, help=SUPPRESS)
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, RelativeCompositionVariabilityTaxon)

    @staticmethod
    def rename_fasta_entries(argv):
        parser = _new_parser(
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
                    [-o/--output <output_file>] [--json]

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

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-i", "--idmap", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, RenameFastaEntries)

    @staticmethod
    def sum_of_pairs_score(argv):
        parser = _new_parser(
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
                phykit sum_of_pairs_score <fasta> -r/--reference <ref.aln> [--json]

                Options
                =====================================================
                <fasta>                     first argument after 
                                            function name should be a
                                            query fasta alignment file
                                            to be scored for accuracy

                -r/--reference              reference fasta alignment to 
                                            compare query alignment to

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("fasta", type=str, help=SUPPRESS)
        parser.add_argument("-r", "--reference", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, SumOfPairsScore)

    @staticmethod
    def variable_sites(argv):
        parser = _new_parser(
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
                phykit variable_sites <alignment> [--json]

                Options
                =====================================================
                <alignment>                 first argument after 
                                            function name should be
                                            an alignment file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("alignment", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, VariableSites)

    @staticmethod
    def phylo_gwas(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Phylogenetic genome-wide association study following the
                Pease et al. (2016) approach.

                Performs per-site association tests between alignment
                columns and a phenotype (categorical or continuous),
                applies Benjamini-Hochberg FDR correction, optionally
                classifies significant associations as monophyletic or
                polyphyletic using a phylogenetic tree, and produces a
                Manhattan plot.

                Categorical phenotypes use Fisher's exact test (2 groups)
                or chi-squared test (>2 groups). Continuous phenotypes use
                point-biserial correlation. Only biallelic sites are tested.

                Aliases:
                  phylo_gwas, pgwas
                Command line interfaces:
                  pk_phylo_gwas, pk_pgwas

                Usage:
                phykit phylo_gwas -a <alignment> -d <phenotype> -o <output>
                  [-t <tree>] [-p <partition>] [--alpha 0.05]
                  [--exclude-monophyletic] [--csv <file>] [--json]
                  [shared plot options]

                Options
                =====================================================
                -a/--alignment              FASTA alignment file

                -d/--phenotype              two-column TSV file:
                                            taxon<tab>phenotype

                -o/--output                 output Manhattan plot path

                -t/--tree                   optional Newick tree for
                                            monophyletic/polyphyletic
                                            classification

                -p/--partition              optional RAxML-style partition
                                            file for gene annotations

                --alpha                     FDR significance threshold
                                            (default: 0.05)

                --exclude-monophyletic      exclude monophyletic
                                            associations from results

                --dot-size                  scale factor for dot size
                                            in the Manhattan plot
                                            (default: 1.0; use 2.0
                                            for double, 0.5 for half)

                --csv                       output per-site results as CSV
                                            to the specified file

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("-a", "--alignment", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("-d", "--phenotype", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("-o", "--output", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("-t", "--tree", type=str, default=None, help=SUPPRESS, metavar="")
        parser.add_argument("-p", "--partition", type=str, default=None, help=SUPPRESS, metavar="")
        parser.add_argument("--alpha", type=float, default=0.05, help=SUPPRESS, metavar="")
        parser.add_argument("--exclude-monophyletic", action="store_true", default=False, help=SUPPRESS)
        parser.add_argument("--dot-size", type=float, default=1.0, help=SUPPRESS, metavar="")
        parser.add_argument("--csv", type=str, default=None, help=SUPPRESS, metavar="")
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, PhyloGwas)

    @staticmethod
    def phylo_anova(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Phylogenetic ANOVA / MANOVA using the Residual
                Randomization Permutation Procedure (RRPP) of
                Adams & Collyer (2018).

                Tests whether a continuous trait (ANOVA) or multiple
                traits (MANOVA) differ across discrete groups while
                accounting for phylogenetic non-independence.

                Auto-detects univariate vs multivariate based on the
                number of response trait columns. Override with
                --method anova or --method manova.

                Aliases:
                  phylo_anova, panova, phylo_manova, pmanova
                Command line interfaces:
                  pk_phylo_anova, pk_panova, pk_phylo_manova, pk_pmanova

                Usage:
                phykit phylo_anova -t <tree> --traits <traits_file>
                  [--group-column <name>] [--method auto|anova|manova]
                  [--permutations <int>] [--pairwise]
                  [--plot-output <file>] [--plot-type boxplot|phylomorphospace]
                  [--seed <int>] [--json]

                Options
                =====================================================
                -t/--tree                   species tree file (required)

                --traits                    TSV file with taxon, group
                                            column, and one or more
                                            response trait columns
                                            (required)

                --group-column              name of the categorical
                                            grouping column (default:
                                            first non-taxon column)

                --method                    analysis method: auto, anova,
                                            or manova (default: auto)

                --permutations              number of RRPP permutations
                                            (default: 1000)

                --pairwise                  include post-hoc pairwise
                                            group comparisons

                --plot-output               output figure path
                                            (.png, .pdf, .svg)

                --plot-type                 boxplot or phylomorphospace
                                            (default: auto — boxplot for
                                            ANOVA, phylomorphospace for
                                            MANOVA)

                --seed                      random seed for reproducible
                                            permutations

                --json                      output results as JSON
                """
            ),
        )
        parser.add_argument("-t", "--tree", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--traits", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--group-column", type=str, default=None, help=SUPPRESS, metavar="")
        parser.add_argument("--method", type=str, default="auto", choices=["auto", "anova", "manova"], help=SUPPRESS, metavar="")
        parser.add_argument("--permutations", type=int, default=1000, help=SUPPRESS, metavar="")
        parser.add_argument("--pairwise", action="store_true", help=SUPPRESS)
        parser.add_argument("--plot-output", type=str, default=None, help=SUPPRESS, metavar="")
        parser.add_argument("--plot-type", type=str, default="auto", choices=["auto", "boxplot", "phylomorphospace"], help=SUPPRESS, metavar="")
        parser.add_argument("--seed", type=int, default=None, help=SUPPRESS, metavar="")
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, PhyloAnova)

    @staticmethod
    def phylo_path(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Phylogenetic path analysis (von Hardenberg & Gonzalez-Voyer
                2013). Compare competing causal DAGs using d-separation
                tests via PGLS with Pagel's lambda, rank models by CICc,
                and estimate model-averaged path coefficients.

                Aliases:
                  phylo_path, ppath, phylopath
                Command line interfaces:
                  pk_phylo_path, pk_ppath, pk_phylopath

                Usage:
                phykit phylo_path -t <tree> --traits <traits_file>
                  --models <models_file> [--best-only]
                  [--plot-output <file>] [--csv <file>] [--json]

                Options
                =====================================================
                -t/--tree                   species tree (required)

                --traits                    TSV file with taxon and
                                            continuous trait columns
                                            (required)

                --models                    model definition file
                                            with candidate DAGs
                                            (required). Format:
                                            name: A->B, B->C, ...

                --best-only                 report only best model
                                            coefficients (default:
                                            model averaging)

                --plot-output               output DAG plot file

                --csv                       output CSV with model
                                            comparison and path
                                            coefficients

                --json                      output results as JSON
                """
            ),
        )
        parser.add_argument("-t", "--tree", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--traits", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--models", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--best-only", action="store_true", help=SUPPRESS)
        parser.add_argument("--plot-output", type=str, default=None, help=SUPPRESS, metavar="")
        parser.add_argument("--csv", type=str, default=None, help=SUPPRESS, metavar="")
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, PhyloPath)

    @staticmethod
    def alignment_subsample(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Randomly subsample genes, partitions, or sites from
                phylogenomic datasets.

                Three modes are available:
                  genes       — subsample alignment files from a list
                  partitions  — subsample partitions from a supermatrix
                  sites       — subsample columns from an alignment

                Aliases:
                  alignment_subsample, aln_subsample, subsample
                Command line interfaces:
                  pk_alignment_subsample, pk_aln_subsample, pk_subsample

                Usage:
                phykit alignment_subsample --mode <genes|partitions|sites>
                  [-a/--alignment <file>] [-l/--list <file>]
                  [-p/--partition <file>]
                  [--number N | --fraction F]
                  [--seed S] [--bootstrap]
                  [-o/--output <prefix>] [--json]

                Options
                =====================================================
                --mode                      subsampling mode: genes,
                                            partitions, or sites

                -a/--alignment              alignment file (FASTA).
                                            Required for partitions
                                            and sites modes.

                -l/--list                   file listing alignment
                                            paths (one per line).
                                            Required for genes mode.

                -p/--partition              RAxML-style partition file.
                                            Required for partitions
                                            mode.

                --number                    exact number of items to
                                            select

                --fraction                  fraction of items to select
                                            (0.0 to 1.0)

                --seed                      random seed for
                                            reproducibility

                --bootstrap                 sample with replacement

                -o/--output                 output file prefix
                                            (default: subsampled)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("--mode", type=str, required=True, choices=["genes", "partitions", "sites"])
        parser.add_argument("-a", "--alignment", type=str, default=None)
        parser.add_argument("-l", "--list", type=str, default=None)
        parser.add_argument("-p", "--partition", type=str, default=None)
        parser.add_argument("--number", type=int, default=None)
        parser.add_argument("--fraction", type=float, default=None)
        parser.add_argument("--seed", type=int, default=None)
        parser.add_argument("--bootstrap", action="store_true", default=False)
        parser.add_argument("-o", "--output", type=str, default="subsampled")
        _add_json_argument(parser)
        _run_service(parser, argv, AlignmentSubsample)

    @staticmethod
    def dstatistic(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Compute Patterson's D-statistic (ABBA-BABA test) for
                detecting introgression or gene flow.

                Two input modes:
                1) Site patterns from an alignment (-a)
                2) Quartet topologies from gene trees (-g)

                Species topology: (((P1, P2), P3), Outgroup).
                Under ILS alone, ABBA and BABA patterns (or
                discordant topologies) are equally frequent. A
                significant excess indicates introgression.

                D > 0: introgression between P2 and P3.
                D < 0: introgression between P1 and P3.
                D = 0: consistent with ILS alone.
                Note: D identifies which lineages exchanged genes
                but cannot determine direction of flow.

                Gene trees can have any number of taxa; only the
                quartet induced by the four specified taxa is
                evaluated from each tree.

                Aliases:
                  dstatistic, dstat, abba_baba
                Command line interfaces:
                  pk_dstatistic, pk_dstat, pk_abba_baba

                Usage:
                phykit dstatistic -a <alignment> --p1 <taxon>
                    --p2 <taxon> --p3 <taxon> --outgroup <taxon>
                    [--block-size 100] [--json]
                phykit dstatistic -g <gene_trees> --p1 <taxon>
                    --p2 <taxon> --p3 <taxon> --outgroup <taxon>
                    [--json]

                Options
                =====================================================
                -a/--alignment              FASTA alignment file
                                            (site-pattern mode)

                -g/--gene-trees             gene trees file, one
                                            Newick per line (gene-
                                            tree mode; trees can
                                            have any number of taxa)

                --p1                        taxon name for P1
                                            (sister to P2)

                --p2                        taxon name for P2
                                            (sister to P1; potential
                                            recipient of gene flow)

                --p3                        taxon name for P3
                                            (donor lineage)

                --outgroup                  outgroup taxon name

                --block-size                block size for jackknife
                                            estimation of standard
                                            error (default: 100;
                                            alignment mode only)

                --support                   minimum branch support
                                            threshold for gene trees;
                                            branches below this value
                                            are collapsed (treated as
                                            unresolved). Gene-tree
                                            mode only.

                --json                      output results as JSON
                """
            ),
        )
        parser.add_argument("-a", "--alignment", type=str, required=False, default=None, help=SUPPRESS, metavar="")
        parser.add_argument("-g", "--gene-trees", type=str, required=False, default=None, help=SUPPRESS, metavar="")
        parser.add_argument("--p1", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--p2", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--p3", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--outgroup", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--block-size", type=int, default=100, help=SUPPRESS, metavar="")
        parser.add_argument("--support", type=float, default=None, help=SUPPRESS, metavar="")
        _add_json_argument(parser)
        _run_service(parser, argv, Dstatistic)

    @staticmethod
    def dfoil(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Compute DFOIL statistics (Pease & Hahn 2015) for
                detecting and polarizing introgression in a 5-taxon
                symmetric phylogeny.

                Topology: ((P1, P2), (P3, P4), Outgroup)
                P1 and P2 are sister taxa; P3 and P4 are sister
                taxa; the two pairs are sister to each other with
                an outgroup rooting the tree.

                Four D-statistics are computed:
                  DFO (far-outer), DIL (inner-left),
                  DFI (far-inner), DOL (outer-left)

                The sign pattern of these four statistics maps to
                a specific introgression scenario via the lookup
                table from Pease & Hahn (2015).

                Aliases:
                  dfoil, dfoil_test
                Command line interfaces:
                  pk_dfoil, pk_dfoil_test

                Usage:
                phykit dfoil -a <alignment> --p1 <taxon>
                    --p2 <taxon> --p3 <taxon> --p4 <taxon>
                    --outgroup <taxon> [--json]

                Options
                =====================================================
                -a/--alignment              FASTA alignment file

                --p1                        taxon name for P1
                                            (sister to P2)

                --p2                        taxon name for P2
                                            (sister to P1)

                --p3                        taxon name for P3
                                            (sister to P4)

                --p4                        taxon name for P4
                                            (sister to P3)

                --outgroup                  outgroup taxon name

                --json                      output results as JSON
                """
            ),
        )
        parser.add_argument("-a", "--alignment", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--p1", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--p2", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--p3", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--p4", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--outgroup", type=str, required=True, help=SUPPRESS, metavar="")
        _add_json_argument(parser)
        _run_service(parser, argv, Dfoil)

    ## Tree functions
    @staticmethod
    def parsimony_score(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Compute the Fitch (1971) maximum parsimony score of a
                tree given an alignment.

                The parsimony score is the minimum number of character
                state changes required to explain the alignment on the
                given tree topology. Each site is scored independently
                using the Fitch downpass algorithm.

                Gap characters (-, N, X, ?) are treated as wildcards.
                Multifurcations are automatically resolved.

                Cross-validated against R's phangorn::parsimony().

                Aliases:
                  parsimony_score, parsimony, pars
                Command line interfaces:
                  pk_parsimony_score, pk_parsimony, pk_pars

                Usage:
                phykit parsimony_score -t <tree> -a <alignment>
                  [-v/--verbose] [--json]

                Options
                =====================================================
                -t/--tree                   tree file (required)

                -a/--alignment              alignment file in FASTA
                                            format (required)

                -v/--verbose                print per-site parsimony
                                            scores

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-a", "--alignment", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        _add_json_argument(parser)
        _run_service(parser, argv, ParsimonyScore)

    @staticmethod
    def character_map(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Map discrete character changes onto a phylogenetic tree
                using Fitch parsimony, classifying each change as a
                synapomorphy, convergence, or reversal.

                Two optimization strategies are supported:
                  ACCTRAN (accelerated transformation) — assigns changes
                  as close to the root as possible.
                  DELTRAN (delayed transformation) — assigns changes as
                  close to the tips as possible.

                Output is a cladogram (default) or phylogram (with
                --phylogram) annotated with colored circles indicating
                character state changes on each branch. Polytomies are
                automatically resolved by adding zero-length branches.

                Summary statistics include tree length (total parsimony
                steps), consistency index (CI), and retention index (RI).

                Aliases:
                  character_map, charmap, synapomorphy_map
                Command line interfaces:
                  pk_character_map, pk_charmap, pk_synapomorphy_map

                Usage:
                phykit character_map -t <tree> -d <data> -o <output>
                  [--optimization acctran|deltran] [--phylogram]
                  [--characters 0,1,3] [--verbose] [--json]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]

                Options
                =====================================================
                -t/--tree                   tree file (required)

                -d/--data                   character matrix file in TSV
                                            format: taxon<tab>char0<tab>
                                            char1<tab>... (required)

                -o/--output                 output figure path (required;
                                            supports .png, .pdf, .svg)

                --optimization              optimization strategy:
                                            acctran or deltran
                                            (default: acctran)

                --phylogram                 draw a phylogram instead of
                                            a cladogram

                --characters                comma-separated character
                                            indices to display on the
                                            plot (e.g., 0,1,3); all
                                            characters shown by default

                --verbose                   print per-character details
                                            including CI/RI and changes

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-o", "--output", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--optimization", type=str, default="acctran",
            required=False, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--phylogram", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument(
            "--characters", type=str, default=None,
            required=False, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, CharacterMap)

    @staticmethod
    def independent_contrasts(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Compute Felsenstein's (1985) phylogenetically independent
                contrasts (PIC) for a continuous trait on a phylogeny.

                Each internal node yields one standardized contrast,
                producing n-1 contrasts for n tips. Contrasts are
                computed by postorder traversal, dividing the trait
                difference between sister clades by the square root of
                the sum of their branch lengths.

                Multifurcations are automatically resolved by adding
                zero-length branches.

                Cross-validated against R's ape::pic().

                Aliases:
                  independent_contrasts, pic, phylo_contrasts
                Command line interfaces:
                  pk_independent_contrasts, pk_pic

                Usage:
                phykit independent_contrasts -t <tree> -d <trait_data>
                  [--json]

                Options
                =====================================================
                -t/--tree                   tree file (required)

                -d/--trait_data             trait data file, two columns:
                                            taxon<tab>value (required)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        _add_json_argument(parser)
        _run_service(parser, argv, IndependentContrasts)

    @staticmethod
    def ancestral_state_reconstruction(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Estimate ancestral states using maximum likelihood.

                Supports two trait types:
                - continuous (default): Brownian Motion model, analogous
                  to R's phytools::fastAnc() and ape::ace(type="ML").
                  Optionally produce a contMap plot.
                - discrete: Mk model with marginal posterior probabilities
                  at each internal node, analogous to ape::ace(type="discrete").
                  Optionally produce a pie-chart phylogeny plot.

                Continuous methods (--type continuous):
                - fast (default): Felsenstein's pruning/contrasts, O(n)
                - ml: full VCV-based ML with exact conditional CIs, O(n^3)

                Discrete models (--type discrete):
                - ER (default): equal rates
                - SYM: symmetric rates
                - ARD: all rates different

                Input trait data can be either:
                (1) A two-column file (taxon<tab>value) when -c is omitted
                (2) A multi-trait file with header row when -c specifies
                    which column to use

                Aliases:
                  ancestral_state_reconstruction, asr, anc_recon
                Command line interfaces:
                  pk_ancestral_state_reconstruction, pk_asr, pk_anc_recon

                Usage:
                phykit ancestral_state_reconstruction -t <tree> -d <trait_data> [-c <trait>] [--type <type>] [-m <method>] [--model <model>] [--ci] [--plot <output>]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait_data             trait data file (2-column or
                                            multi-trait with header)

                -c/--trait                  trait column name (required
                                            for multi-trait files)

                --type                      trait type: continuous or
                                            discrete (default: continuous)

                -m/--method                 method to use: fast or ml
                                            (continuous only; default: fast)

                --model                     Mk model: ER, SYM, or ARD
                                            (discrete only; default: ER)

                --ci                        include 95% confidence
                                            intervals (continuous only)

                --plot                      output path for plot

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --plot-ci                   draw confidence interval bars
                                            at internal nodes on the
                                            contMap plot (requires --ci
                                            and --plot)

                --ci-size                   scale factor for CI bar
                                            size (default: 1.0; use
                                            2.0 for larger, 0.5 for
                                            smaller)

                --json                      output results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-c", "--trait", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--type", type=str, required=False, default="continuous",
            choices=["continuous", "discrete"], help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-m", "--method", type=str, required=False, default="fast",
            choices=["fast", "ml"], help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--model", type=str, required=False, default="ER",
            choices=["ER", "SYM", "ARD"], help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--ci", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument(
            "--plot", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--plot-ci", action="store_true", required=False, default=False,
            help=SUPPRESS,
        )
        parser.add_argument(
            "--ci-size", type=float, required=False, default=1.0,
            help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, AncestralReconstruction)

    @staticmethod
    def concordance_asr(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Concordance-aware ancestral state reconstruction (ASR)
                that incorporates gene tree discordance into ancestral
                estimates. Standard ASR operates on a single species tree
                and ignores gene tree conflict. This command propagates
                topological uncertainty from gene tree discordance into
                ancestral state estimates using gene concordance factors.

                Two methods are available:
                - weighted (default): concordance-weighted ASR on the
                  species tree with NNI alternatives, using law of total
                  variance to separate topological vs parameter uncertainty
                - distribution: reconstruct on each gene tree independently,
                  map nodes by descendant-set identity, report concordance-
                  weighted means and percentile CIs

                Aliases:
                  concordance_asr, conc_asr, casr
                Command line interfaces:
                  pk_concordance_asr, pk_conc_asr, pk_casr

                Usage:
                phykit concordance_asr -t <species_tree> -g <gene_trees> -d <trait_data>
                    [-c <trait>] [-m weighted|distribution] [--ci]
                    [--plot <output>] [--missing-taxa error|shared]
                    [--fig-width <float>] [--fig-height <float>]
                    [--dpi <int>] [--no-title] [--title <str>]
                    [--legend-position <str>]
                    [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                    [--title-fontsize <float>] [--axis-fontsize <float>]
                    [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                    [--json]

                Options
                =====================================================
                -t/--tree                   species tree file

                -g/--gene-trees             file with gene trees (multi-
                                            Newick, one per line)

                -d/--trait_data             trait data file (2-column or
                                            multi-trait with header)

                -c/--trait                  trait column name (required
                                            for multi-trait files)

                -m/--method                 method: weighted or distribution
                                            (default: weighted)

                --ci                        include 95%% confidence
                                            intervals

                --plot                      output path for concordance
                                            ASR plot

                --plot-uncertainty          output path for uncertainty
                                            plot showing the distribution
                                            of ancestral estimates across
                                            gene trees (distribution
                                            method) or concordance
                                            sources (weighted method)
                                            as violin + boxplots

                --missing-taxa              how to handle taxa mismatches:
                                            shared (default) or error

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      output results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-g", "--gene-trees", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-c", "--trait", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-m", "--method", type=str, required=False, default="weighted",
            choices=["weighted", "distribution"], help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--ci", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument(
            "--plot", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--plot-uncertainty", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--missing-taxa", type=str, required=False, default="shared",
            choices=["error", "shared"], help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, ConcordanceAsr)

    @staticmethod
    def chronogram(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Plot a chronogram (time-calibrated phylogeny) with
                geological timescale bands. Requires an ultrametric
                tree and the root age in millions of years (Ma).

                Geological epoch/period/era bands are drawn behind
                the tree as colored stripes, with labels along the
                top. The time axis runs from past (left) to present
                (right).

                Aliases:
                  chronogram, chrono, time_tree
                Command line interfaces:
                  pk_chronogram, pk_chrono, pk_time_tree

                Usage:
                phykit chronogram -t <tree> --root-age <float>
                  --plot-output <file> [--timescale auto|epoch|period|era]
                  [--node-ages] [--circular] [--ladderize]
                  [--color-file <file>] [--json]

                Options
                =====================================================
                -t/--tree                   ultrametric tree file
                                            (required)

                --root-age                  age of the root in Ma
                                            (required)

                --plot-output               output figure path
                                            (required; .png, .pdf, .svg)

                --timescale                 timescale level: auto
                                            (default), epoch, period,
                                            or era. Auto selects based
                                            on root age.

                --node-ages                 label internal nodes with
                                            divergence times (Ma)

                --circular                  draw circular chronogram

                --ladderize                 ladderize the tree

                --color-file                color annotation file
                                            (iTOL-inspired TSV)

                --json                      output node ages as JSON
                """
            ),
        )
        parser.add_argument("-t", "--tree", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--root-age", type=float, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--plot-output", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--timescale", type=str, default="auto", choices=["auto", "epoch", "period", "era"], help=SUPPRESS, metavar="")
        parser.add_argument("--node-ages", action="store_true", help=SUPPRESS)
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, Chronogram)

    @staticmethod
    def dtt(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Disparity through time (DTT) analysis. Computes how
                morphological disparity partitions among subclades
                through time (Harmon et al. 2003).

                At each branching time, calculates the mean relative
                subclade disparity. Under Brownian motion, this declines
                linearly. The MDI (Morphological Disparity Index) is the
                area between the observed DTT and the BM null median.

                Positive MDI = late disparity accumulation
                Negative MDI = early disparity accumulation (radiation)

                Aliases:
                  dtt, disparity_through_time
                Command line interfaces:
                  pk_dtt, pk_disparity_through_time

                Usage:
                phykit dtt -t <tree> --traits <traits_file>
                  [--trait <column>] [--index avg_sq|avg_manhattan]
                  [--nsim <int>] [--seed <int>]
                  [--plot-output <file>] [--json]

                Options
                =====================================================
                -t/--tree                   ultrametric tree file
                                            (required)

                --traits                    TSV file with trait data
                                            (required)

                --trait                     specific trait column name
                                            (default: all traits)

                --index                     disparity index: avg_sq
                                            (average squared Euclidean
                                            distance, default) or
                                            avg_manhattan

                --nsim                      number of BM simulations
                                            for null DTT envelope and
                                            MDI p-value (default: 0,
                                            no simulations)

                --seed                      random seed for
                                            reproducibility

                --plot-output               output figure path

                --json                      output results as JSON
                """
            ),
        )
        parser.add_argument("-t", "--tree", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--traits", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--trait", type=str, default=None, help=SUPPRESS, metavar="")
        parser.add_argument("--index", type=str, default="avg_sq", choices=["avg_sq", "avg_manhattan"], help=SUPPRESS, metavar="")
        parser.add_argument("--nsim", type=int, default=0, help=SUPPRESS, metavar="")
        parser.add_argument("--seed", type=int, default=None, help=SUPPRESS, metavar="")
        parser.add_argument("--plot-output", type=str, default=None, help=SUPPRESS, metavar="")
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, Dtt)

    @staticmethod
    def bipartition_support_stats(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}
                Calculate summary statistics for bipartition support.

                High bipartition support values are thought to be desirable because
                they are indicative of greater certainty in tree topology.

                To obtain all bipartition support values, use the -v/--verbose option.
                In addition to support values for each node, the names of all terminal
                branches tips are also included. Each terminal branch name is separated
                with a semi-colon (;).

                Aliases:
                  bipartition_support_stats, bss
                Command line interfaces:
                  pk_bipartition_support_stats, pk_bss

                Usage:
                phykit bipartition_support_stats <tree>
                  [-v/--verbose]
                  [--thresholds <comma-separated-floats>]
                  [--json]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file 
            
                -v/--verbose                optional argument to print
                                            all bipartition support
                                            values

                --thresholds                optional comma-separated
                                            support cutoffs; reports
                                            count and fraction of
                                            bipartitions below each
                                            cutoff

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument("--thresholds", type=str, required=False, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, BipartitionSupportStats)

    @staticmethod
    def branch_length_multiplier(argv):
        parser = _new_parser(
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
                phykit branch_length_multiplier <tree> -f n [-o/--output <output_file>] [--json]

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

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-f", "--factor", type=float, required=True, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, BranchLengthMultiplier)

    @staticmethod
    def collapse_branches(argv):
        parser = _new_parser(
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
                phykit collapse_branches <tree> -s/--support n [-o/--output <output_file>] [--json]

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

                --json                      optional argument to output
                                            results as JSON

                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-s", "--support", type=float, required=True, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, CollapseBranches)

    @staticmethod
    def covarying_evolutionary_rates(argv):
        parser = _new_parser(
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
                    -r/--reference <reference_tree_file> [-v/--verbose] [--plot] [--plot-output <path>]
                    [--fig-width <float>] [--fig-height <float>]
                    [--dpi <int>] [--no-title] [--title <str>]
                    [--legend-position <str>]
                    [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                    [--title-fontsize <float>] [--axis-fontsize <float>]
                    [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                    [--json]

                Options
                =====================================================
                <tree_file_zero>            first argument after
                                            function name should be
                                            a tree file

                <tree_file_one>             second argument after
                                            function name should be
                                            a tree file

                -r/--reference              a tree to correct branch
                                            lengths by in the two input
                                            trees. Typically, this is a
                                            putative species tree.

                -v/--verbose                print out corrected branch
                                            lengths shared between
                                            tree 0 and tree 1

                --plot                      optional argument to save a
                                            covarying-rates scatter plot

                --plot-output               output path for plot
                                            (default: covarying_rates_plot.png)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree_zero", type=str, help=SUPPRESS)
        parser.add_argument("tree_one", type=str, help=SUPPRESS)
        parser.add_argument(
            "-r", "--reference", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument("--plot", action="store_true", required=False, help=SUPPRESS)
        parser.add_argument(
            "--plot-output",
            type=str,
            default="covarying_rates_plot.png",
            required=False,
            help=SUPPRESS,
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, CovaryingEvolutionaryRates)

    @staticmethod
    def dvmc(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate degree of violation of a molecular clock (or DVMC) in a phylogeny.

                Lower DVMC values are thought to be desirable because they are indicative
                of a lower degree of violation in the molecular clock assumption.

                Typically, outgroup taxa are not included in molecular clock analysis. Thus,
                prior to calculating DVMC from a single gene tree, users may want to prune
                outgroup taxa from the phylogeny. To prune tips from a phylogeny, see the 
                prune_tree function. 

                Calculate DVMC in a tree following Liu et al., PNAS (2017), doi: 10.1073/pnas.1616744114.

                Aliases:
                  degree_of_violation_of_a_molecular_clock, dvmc
                Command line interfaces:
                  pk_degree_of_violation_of_a_molecular_clock, pk_dvmc

                Usage:
                phykit degree_of_violation_of_a_molecular_clock <tree> [--json]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, DVMC)

    @staticmethod
    def evolutionary_rate(argv):
        parser = _new_parser(
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
                phykit evolutionary_rate <tree> [--json]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, EvolutionaryRate)

    @staticmethod
    def hidden_paralogy_check(argv):
        parser = _new_parser(
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

                The output will report if the specified taxa were monophyletic
                or not. The number of rows will reflect how many groups of taxa
                were checked for monophyly. For example, if there were three
                rows of clades in the -c file, there will be three rows in the
                output where the first row in the output corresponds to the 
                results of the first row in the clade file.

                The concept behind this analysis follows
                Siu-Ting et al., Molecular Biology and Evolution (2019).

                Aliases:
                  hidden_paralogy_check, clan_check
                Command line interfaces:
                  pk_hidden_paralogy_check, pk_clan_check

                Usage:
                phykit hidden_paralogy_check <tree> -c/--clade <clade_file> [--json]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                <clade_file>                clade file that specifies
                                            what monophyletic clades
                                            to expect

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-c", "--clade", type=str, required=False, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, HiddenParalogyCheck)

    @staticmethod
    def internal_branch_stats(argv):
        parser = _new_parser(
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
                phykit internal_branch_stats <tree> [-v/--verbose] [--json]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -v/--verbose                optional argument to print
                                            all internal branch lengths

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        _add_json_argument(parser)
        _run_service(parser, argv, InternalBranchStats)

    @staticmethod
    def internode_labeler(argv):
        parser = _new_parser(
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
                phykit internode_labeler <tree> [-o/--output <file>] [--json]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -o/--output                 optional argument to name 
                                            the outputted tree file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, InternodeLabeler)

    @staticmethod
    def last_common_ancestor_subtree(argv):
        parser = _new_parser(
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
                [-o/--output <file>] [--json]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                <list_of_taxa>              list of taxa to get the last
                                            common ancestor subtree for

                -o/--output                 optional argument to name 
                                            the outputted tree file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("list_of_taxa", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, LastCommonAncestorSubtree)

    @staticmethod
    def lb_score(argv):
        parser = _new_parser(
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
                phykit long_branch_score <tree> [-v/--verbose] [--json]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -v/--verbose                optional argument to print
                                            all LB score values

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        _add_json_argument(parser)
        _run_service(parser, argv, LBScore)

    @staticmethod
    def monophyly_check(argv):
        parser = _new_parser(
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
                phykit monophyly_check <tree> <list_of_taxa> [--json]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                <list_of_taxa>              single column file with
                                            list of tip names to 
                                            examine the monophyly of

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("list_of_taxa", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, MonophylyCheck)

    @staticmethod
    def nearest_neighbor_interchange(argv):
        parser = _new_parser(
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
                phykit nearest_neighbor_interchange <tree> [-o/--output <output_file>] [--json]

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

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, NearestNeighborInterchange)

    @staticmethod
    def faiths_pd(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate Faith's phylogenetic diversity (PD) for a
                community of tips on a phylogeny.

                Faith's PD is the sum of branch lengths in the minimum
                subtree that connects a set of taxa. By default, the
                path from the community's most recent common ancestor
                up to the tree root is included, matching Faith (1992)
                and picante::pd(..., include.root = TRUE). Use
                --exclude-root to sum only the branches of the induced
                subtree rooted at the MRCA, matching
                picante::pd(..., include.root = FALSE).

                Aliases:
                  faiths_pd, faith_pd, fpd, phylo_diversity
                Command line interfaces:
                  pk_faiths_pd, pk_faith_pd, pk_fpd, pk_phylo_diversity

                Usage:
                phykit faiths_pd <tree> -t/--taxa <taxa_file>
                  [--exclude-root] [--json]

                Options
                =====================================================
                <tree>                      first argument after
                                            function name should be
                                            a tree file

                -t/--taxa                   file with one tip label per
                                            line defining the community

                --exclude-root              sum only branches of the
                                            induced subtree rooted at
                                            the community MRCA; by
                                            default the path up to the
                                            tree root is included

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument(
            "-t", "--taxa", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--exclude-root",
            dest="exclude_root",
            action="store_true",
            help=SUPPRESS,
        )
        _add_json_argument(parser)
        _run_service(parser, argv, FaithsPD)

    @staticmethod
    def patristic_distances(argv):
        parser = _new_parser(
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
                phykit patristic_distances <tree> [-v/--verbose] [--json]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -v/--verbose                optional argument to print
                                            all patristic distances between
                                            taxa

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        _add_json_argument(parser)
        _run_service(parser, argv, PatristicDistances)

    @staticmethod
    def phylogenetic_signal(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate phylogenetic signal for continuous trait data.

                Supports two methods:
                  - blombergs_k: Blomberg's K statistic (Blomberg et al. 2003)
                    with permutation-based p-value.
                  - lambda: Pagel's lambda (Pagel 1999) with likelihood ratio
                    test p-value.

                Trait file should be tab-delimited with two columns:
                  taxon_name<tab>trait_value

                Lines starting with '#' are treated as comments.

                Output for blombergs_k: K_value<tab>p_value
                Output for lambda: lambda_value<tab>log_likelihood<tab>p_value

                Aliases:
                  phylogenetic_signal, phylo_signal, ps
                Command line interfaces:
                  pk_phylogenetic_signal, pk_phylo_signal, pk_ps

                Usage:
                phykit phylogenetic_signal -t <tree> -d <trait_data> [-m <method>] [-p <permutations>] [--json]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait_data             tab-delimited trait file
                                            (taxon_name<tab>trait_value)

                -m/--method                 method to use: blombergs_k
                                            or lambda (default: blombergs_k)

                -p/--permutations           number of permutations for
                                            blombergs_k (default: 1000)

                -g/--gene-trees             optional multi-Newick file of
                                            gene trees for discordance-aware
                                            VCV computation

                --multivariate              compute K_mult (Adams 2014)
                                            for multivariate traits using
                                            a multi-column TSV file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-m",
            "--method",
            type=str,
            required=False,
            default="blombergs_k",
            choices=["blombergs_k", "lambda"],
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "-p",
            "--permutations",
            type=int,
            required=False,
            default=1000,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "-g", "--gene-trees", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--multivariate", action="store_true", default=False,
            help=SUPPRESS,
        )
        _add_json_argument(parser)
        _run_service(parser, argv, PhylogeneticSignal)

    @staticmethod
    def trait_correlation(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Compute phylogenetic correlations between all pairs of
                traits and display them as a heatmap with significance
                indicators.

                Uses GLS-centering via the tree's variance-covariance
                matrix to account for phylogenetic non-independence.
                P-values are computed from the t-distribution.

                Significance stars in the heatmap:
                  ***  p < 0.001
                  **   p < 0.01
                  *    p < alpha (default 0.05)

                Input is a phylogenetic tree and a tab-delimited
                multi-trait file with a header row:
                taxon<tab>trait1<tab>trait2<tab>...

                Aliases:
                  trait_correlation, trait_corr, phylo_corr
                Command line interfaces:
                  pk_trait_correlation, pk_trait_corr, pk_phylo_corr

                Usage:
                phykit trait_correlation -t <tree> -d <trait_data> -o <output>
                  [--alpha 0.05] [--cluster] [-g <gene_trees>] [--json]

                Options
                =====================================================
                -t/--tree                   tree file (required)

                -d/--trait-data             multi-trait TSV with header
                                            row (required)

                -o/--output                 output figure path (required;
                                            supports .png, .pdf, .svg)

                --alpha                     significance threshold for
                                            marking p-values
                                            (default: 0.05)

                --cluster                   cluster traits by correlation
                                            similarity and display
                                            dendrograms

                -g/--gene-trees             optional multi-Newick file of
                                            gene trees for discordance-
                                            aware VCV computation

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait-data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-o", "--output", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--alpha", type=float, default=0.05, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--cluster", action="store_true", default=False, help=SUPPRESS
        )
        parser.add_argument(
            "-g", "--gene-trees", type=str, default=None, help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, TraitCorrelation)

    @staticmethod
    def phylo_impute(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Phylogenetic imputation of missing trait values using
                conditional multivariate normal distributions.

                Captures both phylogenetic relationships (via the
                tree's variance-covariance matrix) and between-trait
                correlations to predict missing values. Reports
                imputed values with standard errors and 95% CIs.

                Missing values in the input trait file may be marked
                as NA, na, ?, or left empty.

                Input is a phylogenetic tree and a tab-delimited
                multi-trait file with a header row:
                taxon<tab>trait1<tab>trait2<tab>...

                Aliases:
                  phylo_impute, impute, phylo_imp
                Command line interfaces:
                  pk_phylo_impute, pk_impute, pk_phylo_imp

                Usage:
                phykit phylo_impute -t <tree> -d <trait_data> -o <output>
                  [-g <gene_trees>] [--json]

                Options
                =====================================================
                -t/--tree                   tree file (required)

                -d/--trait-data             multi-trait TSV with header
                                            row; missing values marked
                                            as NA, ?, or empty
                                            (required)

                -o/--output                 output TSV file with
                                            imputed values (required)

                -g/--gene-trees             optional multi-Newick file
                                            of gene trees for
                                            discordance-aware VCV

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait-data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-o", "--output", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-g", "--gene-trees", type=str, default=None, help=SUPPRESS, metavar=""
        )
        _add_json_argument(parser)
        _run_service(parser, argv, PhyloImpute)

    @staticmethod
    def phylogenetic_ordination(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Perform phylogenetic ordination (PCA, t-SNE, or UMAP) on
                continuous multi-trait data while accounting for
                phylogenetic non-independence among species.

                Phylogenetic correction uses GLS-centering via the tree's
                variance-covariance matrix. For PCA, eigendecomposition
                of the evolutionary rate matrix is performed. For t-SNE
                and UMAP, nonlinear embedding is applied to the
                GLS-centered data.

                Input is a phylogenetic tree and a tab-delimited
                multi-trait file with a header row:
                taxon<tab>trait1<tab>trait2<tab>...

                Aliases:
                  phylogenetic_ordination, phylo_ordination, ordination, ord,
                  phylo_pca, phyl_pca, ppca, phylo_dimreduce, dimreduce, pdr
                Command line interfaces:
                  pk_phylogenetic_ordination, pk_phylo_ordination, pk_ordination,
                  pk_ord, pk_phylo_pca, pk_phyl_pca, pk_ppca,
                  pk_phylo_dimreduce, pk_dimreduce, pk_pdr

                Usage:
                phykit phylogenetic_ordination -t <tree> -d <trait_data> [--method <pca|tsne|umap>] [--correction <BM|lambda>] [--mode <cov|corr>] [--n-components <int>] [--perplexity <float>] [--n-neighbors <int>] [--min-dist <float>] [--seed <int>] [--plot] [--plot-tree] [--no-plot-tree] [--color-by <col_or_file>] [--tree-color-by <col_or_file>] [--plot-output <path>]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait_data             tab-delimited multi-trait file
                                            with header row

                --method                    ordination method: pca, tsne, or
                                            umap (default: pca)

                --correction                phylogenetic correction: BM or
                                            lambda (default: BM)

                --mode                      PCA mode: cov or corr
                                            (default: cov; PCA only)

                --n-components              number of embedding dimensions
                                            (default: 2; tsne/umap only)

                --perplexity                t-SNE perplexity (default: auto)

                --n-neighbors               UMAP n_neighbors (default: auto)

                --min-dist                  UMAP min_dist (default: 0.1)

                --seed                      random seed for reproducibility

                --plot                      optional argument to save a
                                            scatter plot

                --plot-tree                 overlay phylogeny edges via
                                            ancestral reconstruction
                                            (default for tsne/umap)

                --no-plot-tree              disable phylogeny overlay for
                                            tsne/umap plots

                --color-by                  color tip points by trait;
                                            specify a column name from the
                                            multi-trait file or a separate
                                            tab-delimited file (taxon<tab>value)

                --tree-color-by             color phylogeny edges by a trait;
                                            specify a column name or a file
                                            (default: distance from root)

                --plot-output               output path for plot
                                            (default: phylo_ordination_plot.png)

                -g/--gene-trees             optional multi-Newick file of
                                            gene trees for discordance-aware
                                            VCV computation

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--method",
            type=str,
            required=False,
            default="pca",
            choices=["pca", "tsne", "umap"],
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--correction",
            type=str,
            required=False,
            default="BM",
            choices=["BM", "lambda"],
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--mode",
            type=str,
            required=False,
            default="cov",
            choices=["cov", "corr"],
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--n-components",
            type=int,
            required=False,
            default=2,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--perplexity",
            type=float,
            required=False,
            default=None,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--n-neighbors",
            type=int,
            required=False,
            default=None,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--min-dist",
            type=float,
            required=False,
            default=0.1,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--seed",
            type=int,
            required=False,
            default=None,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument("--plot", action="store_true", required=False, help=SUPPRESS)
        parser.add_argument("--plot-tree", action="store_true", required=False, help=SUPPRESS)
        parser.add_argument("--no-plot-tree", action="store_true", required=False, help=SUPPRESS)
        parser.add_argument(
            "--color-by",
            type=str,
            default=None,
            required=False,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--tree-color-by",
            type=str,
            default=None,
            required=False,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--plot-output",
            type=str,
            default="phylo_ordination_plot.png",
            required=False,
            help=SUPPRESS,
        )
        parser.add_argument(
            "-g", "--gene-trees", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, PhylogeneticOrdination)

    @staticmethod
    def phylogenetic_pca(argv):
        Phykit.phylogenetic_ordination(argv)

    @staticmethod
    def phylogenetic_dimreduce(argv):
        Phykit.phylogenetic_ordination(argv)

    @staticmethod
    def phylo_heatmap(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Draw a phylogenetic heatmap: a phylogeny alongside a
                color-coded matrix of numeric trait values. Rows are
                aligned to tree tips.

                Analogous to R's phytools::phylo.heatmap().

                Aliases:
                  phylo_heatmap, pheatmap, ph
                Command line interfaces:
                  pk_phylo_heatmap, pk_pheatmap, pk_ph

                Usage:
                phykit phylo_heatmap -t <tree> -d <data> -o <output>
                  [--split 0.3] [--standardize] [--cmap viridis]
                  [--json]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]

                Options
                =====================================================
                -t/--tree                   tree file (required)

                -d/--data                   numeric data matrix in TSV
                                            format with header row
                                            (required)

                -o/--output                 output figure path (required;
                                            supports .png, .pdf, .svg)

                --split                     fraction of figure width for
                                            the tree panel (default: 0.3)

                --standardize               z-score each column before
                                            coloring

                --cmap                      matplotlib colormap name
                                            (default: viridis)

                --cluster-columns           cluster trait columns by
                                            similarity and display a
                                            dendrogram at the top

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --ylabel-fontsize           font size for taxon labels;
                                            0 to hide

                --xlabel-fontsize           font size for trait column
                                            labels; 0 to hide

                --json                      optional argument to output
                                            metadata as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-o", "--output", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--split", type=float, required=False, default=0.3,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--standardize", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument(
            "--cmap", type=str, required=False, default="viridis",
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--cluster-columns", action="store_true", required=False,
            default=False, help=SUPPRESS
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, PhyloHeatmap)

    @staticmethod
    def phylomorphospace(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Plot a phylomorphospace: two raw traits in trait space
                with the phylogeny overlaid via ML-reconstructed
                ancestral states at internal nodes.

                Tree edges are colored by distance from root (coolwarm
                colormap). Tip points are colored by --color-by if
                specified, otherwise default blue.

                Input is a phylogenetic tree and a tab-delimited
                multi-trait file with a header row:
                taxon<tab>trait1<tab>trait2<tab>...

                If the trait file has exactly 2 traits and --trait-x /
                --trait-y are omitted, the first two columns are used
                automatically.

                Aliases:
                  phylomorphospace, phylomorpho, phmo
                Command line interfaces:
                  pk_phylomorphospace, pk_phylomorpho, pk_phmo

                Usage:
                phykit phylomorphospace -t <tree> -d <trait_data> [--trait-x <name>] [--trait-y <name>] [--color-by <col_or_file>] [--plot-output <path>]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait_data             tab-delimited multi-trait file
                                            with header row

                --trait-x                   column name for x-axis trait

                --trait-y                   column name for y-axis trait

                --color-by                  color tip points by trait;
                                            specify a column name from the
                                            multi-trait file or a separate
                                            tab-delimited file

                --plot-output               output path for plot
                                            (default: phylomorphospace_plot.png)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--trait-x",
            type=str,
            default=None,
            required=False,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--trait-y",
            type=str,
            default=None,
            required=False,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--color-by",
            type=str,
            default=None,
            required=False,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--plot-output",
            type=str,
            default="phylomorphospace_plot.png",
            required=False,
            help=SUPPRESS,
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, Phylomorphospace)

    @staticmethod
    def phylogenetic_regression(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Fit a Phylogenetic Generalized Least Squares (PGLS)
                regression while accounting for phylogenetic non-independence
                among species, analogous to R's caper::pgls().

                Input is a phylogenetic tree and a tab-delimited
                multi-trait file with a header row:
                taxon<tab>trait1<tab>trait2<tab>...

                Two methods are available:
                - BM (default): Brownian motion (lambda fixed at 1)
                - lambda: jointly estimates Pagel's lambda via ML

                Output includes coefficient estimates, standard errors,
                t-values, p-values, R-squared, F-statistic, log-likelihood,
                and AIC.

                Aliases:
                  phylogenetic_regression, phylo_regression, pgls
                Command line interfaces:
                  pk_phylogenetic_regression, pk_phylo_regression, pk_pgls

                Usage:
                phykit phylogenetic_regression -t <tree> -d <trait_data> -y <response> -x <predictor1> [predictor2 ...] [-m <method>] [--json]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait_data             tab-delimited multi-trait file
                                            with header row

                -y/--response               response (dependent) variable
                                            column name

                -x/--predictors             one or more predictor column
                                            names

                -m/--method                 method to use: BM or lambda
                                            (default: BM)

                -g/--gene-trees             optional multi-Newick file of
                                            gene trees for discordance-aware
                                            VCV computation

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-y", "--response", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-x",
            "--predictors",
            type=str,
            nargs="+",
            required=True,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "-m",
            "--method",
            type=str,
            required=False,
            default="BM",
            choices=["BM", "lambda"],
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "-g", "--gene-trees", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        _add_json_argument(parser)
        _run_service(parser, argv, PhylogeneticRegression)

    @staticmethod
    def phylogenetic_glm(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Fit a Phylogenetic Generalized Linear Model (GLM) for binary
                or count response data while accounting for phylogenetic
                non-independence among species.

                Two families are supported:
                - binomial: logistic regression via Maximum Penalized Likelihood
                  Estimation (logistic_MPLE; Ives & Garland 2010)
                - poisson: Poisson regression via Generalized Estimating Equations
                  (poisson_GEE; Paradis & Claude 2002)

                Input is a phylogenetic tree and a tab-delimited multi-trait
                file with a header row:
                taxon<tab>trait1<tab>trait2<tab>...

                Output includes coefficient estimates, standard errors,
                z-values, p-values, log-likelihood, and AIC. For binomial
                models, the estimated phylogenetic signal parameter alpha is
                reported. For Poisson models, the overdispersion parameter
                phi is reported.

                Aliases:
                  phylogenetic_glm, phylo_glm, pglm
                Command line interfaces:
                  pk_phylogenetic_glm, pk_phylo_glm, pk_pglm

                Usage:
                phykit phylogenetic_glm -t <tree> -d <trait_data> -y <response> -x <predictor1> [predictor2 ...] --family <binomial|poisson> [--json]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait_data             tab-delimited multi-trait file
                                            with header row

                -y/--response               response (dependent) variable
                                            column name

                -x/--predictors             one or more predictor column
                                            names

                --family                    distribution family: binomial
                                            or poisson

                --method                    estimation method: logistic_MPLE
                                            or poisson_GEE (auto from family)

                --btol                      linear predictor bound for
                                            logistic model (default: 10)

                --log-alpha-bound           bound on log(alpha*Tmax) for
                                            logistic model (default: 4)

                -g/--gene-trees             optional multi-Newick file of
                                            gene trees for discordance-aware
                                            VCV computation

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-y", "--response", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-x",
            "--predictors",
            type=str,
            nargs="+",
            required=True,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--family",
            type=str,
            required=True,
            choices=["binomial", "poisson"],
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--method",
            type=str,
            required=False,
            default=None,
            choices=["logistic_MPLE", "poisson_GEE"],
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--btol",
            type=float,
            required=False,
            default=10,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--log-alpha-bound",
            type=float,
            required=False,
            default=4,
            dest="log_alpha_bound",
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "-g", "--gene-trees", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        _add_json_argument(parser)
        _run_service(parser, argv, PhylogeneticGLM)

    @staticmethod
    def phylo_logistic(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Fit a Phylogenetic Logistic Regression for binary (0/1)
                response data while accounting for phylogenetic
                non-independence among species (Ives & Garland 2010).

                Uses Maximum Penalized Likelihood Estimation (logistic_MPLE)
                with Firth's bias-correction penalty and jointly estimates
                the phylogenetic signal parameter alpha via the
                OU-transformed variance-covariance matrix.

                Input is a phylogenetic tree and a tab-delimited multi-trait
                file with a header row:
                taxon<tab>trait1<tab>trait2<tab>...

                Output includes coefficient estimates, standard errors,
                z-values, p-values, alpha, log-likelihood, penalized
                log-likelihood, and AIC.

                Aliases:
                  phylo_logistic, phylo_logreg, plogreg
                Command line interfaces:
                  pk_phylo_logistic, pk_phylo_logreg, pk_plogreg

                Usage:
                phykit phylo_logistic -t <tree> -d <trait_data> --response <column> --predictor <column> [--method logistic_MPLE|logistic_IG10] [--json]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait-data             tab-delimited multi-trait file
                                            with header row

                --response                  binary response column name
                                            (must contain only 0 and 1)

                --predictor                 predictor column name(s),
                                            comma-separated for multiple

                --method                    estimation method: logistic_MPLE
                                            or logistic_IG10
                                            (default: logistic_MPLE)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait-data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--response", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--predictor", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--method",
            type=str,
            required=False,
            default="logistic_MPLE",
            choices=["logistic_MPLE", "logistic_IG10"],
            help=SUPPRESS,
            metavar="",
        )
        _add_json_argument(parser)
        _run_service(parser, argv, PhyloLogistic)

    @staticmethod
    def stochastic_character_map(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Perform Stochastic Character Mapping (SIMMAP) of discrete
                traits onto a phylogeny (Huelsenbeck et al. 2003;
                Bollback 2006), analogous to R's phytools::make.simmap().

                Fits a continuous-time Markov chain (CTMC) rate matrix Q
                via maximum likelihood, then simulates character histories
                conditioned on tip states. Three models are available:
                ER (equal rates), SYM (symmetric), ARD (all rates differ).

                Input is a phylogenetic tree and a tab-delimited file
                with a header row: taxon<tab>trait_column<tab>...

                Output includes the fitted Q matrix, log-likelihood,
                mean dwelling times, and mean transition counts.

                Aliases:
                  stochastic_character_map, simmap, scm
                Command line interfaces:
                  pk_stochastic_character_map, pk_simmap, pk_scm

                Usage:
                phykit stochastic_character_map -t <tree> -d <trait_data> -c <trait_column> [-m <model>] [-n <nsim>] [--seed <seed>] [--plot <output.png>]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait_data             tab-delimited trait file
                                            with header row

                -c/--trait                  column name for discrete
                                            character trait

                -m/--model                  substitution model: ER,
                                            SYM, or ARD (default: ER)

                -n/--nsim                   number of stochastic
                                            mapping simulations
                                            (default: 100)

                --seed                      random seed for
                                            reproducibility

                --plot                      output plot file path

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-c", "--trait", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-m",
            "--model",
            type=str,
            required=False,
            default="ER",
            choices=["ER", "SYM", "ARD"],
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "-n", "--nsim", type=int, required=False, default=100,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--seed", type=int, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--plot", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, StochasticCharacterMap)

    @staticmethod
    def simmap_summary(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Run N stochastic character maps and summarize per-branch
                dwelling time proportions, expected transitions, and
                posterior state probabilities at each node.

                This extends stochastic_character_map by providing a
                detailed per-branch summary analogous to
                phytools::describe.simmap() in R.

                Aliases:
                  simmap_summary, smsummary, describe_simmap
                Command line interfaces:
                  pk_simmap_summary, pk_smsummary, pk_describe_simmap

                Usage:
                phykit simmap_summary -t <tree> -d <trait_data> -c <trait>
                  [-m/--model ER|SYM|ARD] [-n/--nsim <int>]
                  [--seed <int>] [--plot <file>] [--csv <file>] [--json]

                Options
                =====================================================
                -t/--tree                   phylogenetic tree file
                                            (required)

                -d/--trait_data             tab-delimited trait file
                                            with header row (required)

                -c/--trait                  column name for the
                                            discrete character trait
                                            (required)

                -m/--model                  substitution model: ER
                                            (equal rates), SYM
                                            (symmetric), or ARD (all
                                            rates different). Default: ER

                -n/--nsim                   number of stochastic maps
                                            to simulate (default: 100)

                --seed                      random seed for
                                            reproducibility

                --plot                      output plot file showing
                                            tree with posterior pie
                                            charts at nodes

                --csv                       output CSV file with
                                            per-branch dwelling
                                            proportions and node
                                            posteriors

                --json                      output results as JSON
                """
            ),
        )
        parser.add_argument("-t", "--tree", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("-c", "--trait", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("-m", "--model", type=str, required=False, default="ER", help=SUPPRESS, metavar="")
        parser.add_argument("-n", "--nsim", type=int, required=False, default=100, help=SUPPRESS, metavar="")
        parser.add_argument("--seed", type=int, required=False, default=None, help=SUPPRESS, metavar="")
        parser.add_argument("--plot", type=str, required=False, default=None, help=SUPPRESS, metavar="")
        parser.add_argument("--csv", type=str, required=False, default=None, help=SUPPRESS, metavar="")
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, SimmapSummary)

    @staticmethod
    def cont_map(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Continuous Trait Map (contMap) visualization.

                Runs ancestral state reconstruction internally (fast
                Felsenstein two-pass algorithm) and produces a contMap
                plot: a phylogram with branches colored by a continuous
                gradient (coolwarm colormap) representing inferred
                trait values.

                Input is a phylogenetic tree and a tab-delimited file
                with two columns: taxon_name<tab>trait_value (no header).

                Aliases:
                  cont_map, contmap, cmap
                Command line interfaces:
                  pk_cont_map, pk_contmap, pk_cmap

                Usage:
                phykit cont_map -t <tree> -d <trait_data> -o <output.png>
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait_data             tab-delimited trait file
                                            (taxon<tab>value, no header)

                -o/--output                 output plot file path
                                            (required)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to also
                                            output results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-o", "--output", type=str, required=True, help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, ContMap)

    @staticmethod
    def density_map(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Density Map visualization of posterior discrete state
                probabilities along each branch of a phylogeny.

                Runs stochastic character mapping internally (N
                simulations), then for each point along each branch
                computes the fraction of simulations in each state.
                Branches are colored by a gradient reflecting these
                probabilities. Analogous to R's phytools::densityMap().

                Input is a phylogenetic tree and a tab-delimited file
                with a header row: taxon<tab>trait_column<tab>...

                Aliases:
                  density_map, densitymap, dmap
                Command line interfaces:
                  pk_density_map, pk_densitymap, pk_dmap

                Usage:
                phykit density_map -t <tree> -d <trait_data> -c <trait> -o <output.png> [-n <nsim>] [--seed <seed>]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait_data             tab-delimited trait file
                                            with header row

                -c/--trait                  column name for discrete
                                            character trait

                -n/--nsim                   number of stochastic
                                            mapping simulations
                                            (default: 100)

                --seed                      random seed for
                                            reproducibility

                -o/--output                 output plot file path
                                            (required)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to also
                                            output results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-c", "--trait", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-n", "--nsim", type=int, required=False, default=100,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--seed", type=int, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-o", "--output", type=str, required=True, help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, DensityMap)

    @staticmethod
    def phenogram(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Plot a phenogram (traitgram) showing continuous trait
                evolution across a phylogeny. X-axis shows distance from
                root (time), Y-axis shows trait values. Tips are plotted
                at observed values, internal nodes at ML ancestral
                estimates. Analogous to R's phytools::phenogram().

                Aliases:
                  phenogram, traitgram, tg
                Command line interfaces:
                  pk_phenogram, pk_traitgram, pk_tg

                Usage:
                phykit phenogram -t <tree> -d <trait_data> -o <output.png>
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait_data             tab-delimited trait file
                                            (taxon<tab>value)

                -o/--output                 output plot file path
                                            (required)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-o", "--output", type=str, required=True, help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, Phenogram)

    @staticmethod
    def cophylo(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Cophylogenetic (tanglegram) plot of two phylogenies.
                Draws two trees facing each other with connecting lines
                between matching taxa, analogous to R's phytools::cophylo().

                By default, taxa are matched by identical tip names.
                A mapping file can be provided to match differently named
                taxa. Internal nodes of tree2 are rotated to minimize
                line crossings.

                Aliases:
                  cophylo, tanglegram, tangle
                Command line interfaces:
                  pk_cophylo, pk_tanglegram, pk_tangle

                Usage:
                phykit cophylo -t <tree1> -t2 <tree2> -o <output.png> [-m <mapping>]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                -t/--tree1                  first tree file

                -t2/--tree2                 second tree file

                -o/--output                 output plot file path
                                            (required)

                -m/--mapping                optional tab-delimited
                                            mapping file
                                            (taxon1<tab>taxon2)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree1", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-t2", "--tree2", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-o", "--output", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-m", "--mapping", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, Cophylo)

    @staticmethod
    def rate_heterogeneity(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Test for rate heterogeneity across phylogenetic regimes
                using multi-rate Brownian motion (O'Meara et al. 2006),
                analogous to R's phytools::brownie.lite().

                Fits single-rate vs. multi-rate BM models and performs a
                likelihood ratio test. Users specify a tree, continuous
                trait data, and a regime file mapping tips to regimes.

                Regime assignments to internal branches are inferred via
                Fitch parsimony. Per-regime VCV matrices are decomposed
                and per-regime sigma-squared values are estimated via ML.

                Aliases:
                  rate_heterogeneity, brownie, rh
                Command line interfaces:
                  pk_rate_heterogeneity, pk_brownie, pk_rh

                Usage:
                phykit rate_heterogeneity -t <tree> -d <trait_data> -r <regime_data> [-n <nsim>] [--seed <seed>] [--plot <output.png>]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait_data             tab-delimited trait file
                                            (taxon<tab>value)

                -r/--regime_data            tab-delimited regime file
                                            (taxon<tab>regime_label)

                -n/--nsim                   number of parametric
                                            bootstrap simulations
                                            (default: 0, no bootstrap)

                --seed                      random seed for
                                            reproducibility

                --plot                      output plot file path

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-r", "--regime_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-n", "--nsim", type=int, required=False, default=0,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--seed", type=int, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--plot", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, RateHeterogeneity)

    @staticmethod
    def fit_discrete(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Compare models of discrete trait evolution on a phylogeny.

                Fits ER (Equal Rates), SYM (Symmetric), and ARD (All Rates
                Different) Mk models of discrete character evolution via
                maximum likelihood. Compares models using AIC and BIC.

                Analogous to R's geiger::fitDiscrete().

                Aliases:
                  fit_discrete, fitdiscrete, fd
                Command line interfaces:
                  pk_fit_discrete, pk_fitdiscrete, pk_fd

                Usage:
                phykit fit_discrete -t <tree> -d <trait_data> -c <trait>
                  [--models ER,SYM,ARD] [--json]

                Options
                =====================================================
                -t/--tree                   tree file (required)

                -d/--trait_data             trait data file in TSV format
                                            (required)

                -c/--trait                  column name for the discrete
                                            trait in the data file
                                            (required)

                --models                    comma-separated list of models
                                            to fit (default: ER,SYM,ARD)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-c", "--trait", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--models", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        _add_json_argument(parser)
        _run_service(parser, argv, FitDiscrete)

    @staticmethod
    def fit_continuous(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Compare models of continuous trait evolution on a phylogeny.

                Fits up to 7 models (BM, OU, EB, Lambda, Delta, Kappa, White)
                and ranks them by AIC, BIC, and AIC weights — analogous to
                R's geiger::fitContinuous().

                Models:
                  BM      — Brownian motion (baseline)
                  OU      — Ornstein-Uhlenbeck (stabilizing selection)
                  EB      — Early Burst (Harmon et al. 2010)
                  Lambda  — Pagel's lambda (phylogenetic signal)
                  Delta   — Pagel's delta (tempo of evolution)
                  Kappa   — Pagel's kappa (punctuational vs gradual)
                  White   — White noise (no phylogenetic signal)

                Aliases:
                  fit_continuous, fitcontinuous, fc
                Command line interfaces:
                  pk_fit_continuous, pk_fitcontinuous, pk_fc

                Usage:
                phykit fit_continuous -t <tree> -d <trait_data> [--models BM,OU,Lambda] [--json]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait_data             tab-delimited trait file
                                            (taxon<tab>value)

                --models                    comma-separated list of models
                                            to fit (default: all 7)

                -g/--gene-trees             optional multi-Newick file of
                                            gene trees for discordance-aware
                                            VCV computation

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--models", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-g", "--gene-trees", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        _add_json_argument(parser)
        _run_service(parser, argv, FitContinuous)

    @staticmethod
    def ouwie(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Fit multi-regime Ornstein-Uhlenbeck models of continuous
                trait evolution (Beaulieu et al. 2012), analogous to
                R's OUwie package.

                Models:
                  BM1   — single-rate Brownian motion
                  BMS   — multi-rate Brownian motion (per-regime sigma2)
                  OU1   — single-regime Ornstein-Uhlenbeck
                  OUM   — multi-regime OU (per-regime optima)
                  OUMV  — OUM + per-regime sigma2
                  OUMA  — OUM + per-regime alpha
                  OUMVA — all parameters regime-specific

                Aliases:
                  ouwie, fit_ouwie, multi_regime_ou
                Command line interfaces:
                  pk_ouwie, pk_fit_ouwie, pk_multi_regime_ou

                Usage:
                phykit ouwie -t <tree> -d <trait_data> -r <regime_data> [--models BM1,OUM,OUMVA] [--json]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait_data             tab-delimited trait file
                                            (taxon<tab>value)

                -r/--regime_data            tab-delimited regime file
                                            (taxon<tab>regime_label)

                --models                    comma-separated list of models
                                            to fit (default: all 7)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-r", "--regime_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--models", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        _add_json_argument(parser)
        _run_service(parser, argv, OUwie)

    @staticmethod
    def ou_shift_detection(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Automatic OU shift detection using LASSO (l1ou approach).

                Discovers where on the phylogeny the adaptive optimum
                changed, using the LASSO-based approach from
                Khabbazian et al. (2016). No regime file is needed —
                only a tree and trait data.

                Aliases:
                  ou_shift_detection, ou_shifts, l1ou, detect_shifts
                Command line interfaces:
                  pk_ou_shift_detection, pk_ou_shifts, pk_l1ou, pk_detect_shifts

                Usage:
                phykit l1ou -t <tree> -d <trait_data> [--criterion pBIC] [--max-shifts N] [--json]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait_data             tab-delimited trait file
                                            (taxon<tab>value)

                --criterion                 model selection criterion:
                                            pBIC (default), BIC, or AICc

                --max-shifts                maximum number of shifts to
                                            consider (default: n/2)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--criterion", type=str, required=False, default="pBIC",
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--max-shifts", type=int, required=False, default=None,
            help=SUPPRESS, metavar="", dest="max_shifts"
        )
        _add_json_argument(parser)
        _run_service(parser, argv, OUShiftDetection)

    @staticmethod
    def polytomy_test(argv):
        parser = _new_parser(
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
                phykit polytomy_test -t/--trees <trees> -g/--groups <groups> [--json]

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

                --json                      optional argument to output
                                            results as JSON
                                            
                For example, the groups file could look like the following:
                #labels group0  group1  group2
                name_of_test    tip_name_A;tip_name_B   tip_name_C  tip_name_D;tip_name_E
                """
            ),
        )
        parser.add_argument("-t", "--trees", type=str, help=SUPPRESS)
        parser.add_argument("-g", "--groups", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, PolytomyTest)

    @staticmethod
    def print_tree(argv):
        parser = _new_parser(
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
                phykit print_tree <tree> [-r/--remove] [--json]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -r/--remove                 optional argument to print
                                            the phylogeny without branch
                                            lengths

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument(
            "-r", "--remove", action="store_true", required=False, help=SUPPRESS
        )
        _add_json_argument(parser)
        _run_service(parser, argv, PrintTree)

    @staticmethod
    def consensus_tree(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Infer a consensus tree from a collection of trees.

                Input can be either:
                1) a file with one Newick tree per line, or
                2) a file with one tree-file path per line.

                If input trees have different taxon sets, use
                --missing-taxa shared to prune each tree to the
                intersection of taxa before inferring consensus.

                Aliases:
                  consensus_tree, consensus, ctree
                Command line interfaces:
                  pk_consensus_tree, pk_consensus, pk_ctree

                Usage:
                phykit consensus_tree -t/--trees <trees>
                    [-m/--method strict|majority]
                    [--missing-taxa error|shared] [--json]

                Options
                =====================================================
                -t/--trees                 file containing trees or
                                           tree paths

                -m/--method                consensus method to infer
                                           (Default: majority)

                --missing-taxa             how to handle mismatched
                                           taxa across trees:
                                           error or shared
                                           (Default: error)

                --json                     optional argument to output
                                           results as JSON
                """
            ),
        )
        parser.add_argument("-t", "--trees", type=str, required=True, help=SUPPRESS)
        parser.add_argument(
            "-m",
            "--method",
            type=str,
            choices=["strict", "majority"],
            default="majority",
            required=False,
            help=SUPPRESS,
        )
        parser.add_argument(
            "--missing-taxa",
            type=str,
            choices=["error", "shared"],
            default="error",
            required=False,
            help=SUPPRESS,
        )
        _add_json_argument(parser)
        _run_service(parser, argv, ConsensusTree)

    @staticmethod
    def consensus_network(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Extract bipartition splits from a collection of gene trees
                and summarize conflicting phylogenetic signal.

                Counts how frequently each non-trivial bipartition appears
                across input trees and filters by a minimum frequency
                threshold.  Optionally draws a circular splits network
                diagram.

                Polytomies (collapsed branches) in input trees are
                handled conservatively: splits from polytomous nodes
                are excluded since they represent unresolved
                relationships. Trifurcating roots (standard unrooted
                Newick) are not affected.

                Input can be either:
                1) a file with one Newick tree per line, or
                2) a file with one tree-file path per line.

                Aliases:
                  consensus_network, consnet, splitnet, splits_network
                Command line interfaces:
                  pk_consensus_network, pk_consnet, pk_splitnet,
                  pk_splits_network

                Usage:
                phykit consensus_network -t/--trees <trees>
                    [--threshold 0.1]
                    [--missing-taxa error|shared]
                    [--plot-output <file>]
                    [--fig-width <float>] [--fig-height <float>]
                    [--dpi <int>] [--no-title] [--title <str>]
                    [--legend-position <str>]
                    [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                    [--title-fontsize <float>] [--axis-fontsize <float>]
                    [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                    [--json]

                Options
                =====================================================
                -t/--trees                 file containing trees or
                                           tree paths

                --threshold                minimum split frequency
                                           to include (0-1)
                                           (Default: 0.1)

                --missing-taxa             how to handle mismatched
                                           taxa across trees:
                                           error or shared
                                           (Default: error)

                --plot-output              output filename for the
                                           circular splits network
                                           plot (optional)

                --max-splits               maximum number of splits
                                           to include in the network
                                           graph (default: 30; higher
                                           splits are still reported
                                           in text/JSON output)

                --histogram                output filename for a split
                                           frequency histogram (optional)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                     optional argument to output
                                           results as JSON
                """
            ),
        )
        parser.add_argument("-t", "--trees", type=str, required=True, help=SUPPRESS)
        parser.add_argument(
            "--threshold",
            type=float,
            default=0.1,
            required=False,
            help=SUPPRESS,
        )
        parser.add_argument(
            "--missing-taxa",
            type=str,
            choices=["allow", "error", "shared"],
            default="allow",
            required=False,
            help=SUPPRESS,
        )
        parser.add_argument(
            "--plot-output",
            type=str,
            default=None,
            required=False,
            help=SUPPRESS,
        )
        parser.add_argument(
            "--max-splits",
            type=int,
            default=30,
            required=False,
            help=SUPPRESS,
        )
        parser.add_argument(
            "--histogram",
            type=str,
            default=None,
            required=False,
            help=SUPPRESS,
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, ConsensusNetwork)

    @staticmethod
    def neighbor_net(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Construct a NeighborNet phylogenetic network from pairwise
                distances and visualize it as a planar splits graph.

                Accepts either a FASTA alignment (from which pairwise
                distances are computed using the chosen metric) or a
                pre-computed CSV distance matrix.

                The implementation builds a Neighbor-Joining tree to obtain
                a circular ordering, enumerates all circular splits
                compatible with that ordering, and estimates split weights
                via non-negative least squares (NNLS).

                Aliases:
                  neighbor_net, nnet
                Command line interfaces:
                  pk_neighbor_net, pk_nnet

                Usage:
                phykit neighbor_net -a/--alignment <alignment>
                    -o/--output <output_figure>
                    [--distance-matrix <csv>]
                    [--metric identity|p-distance|jc]
                    [--max-splits 30]
                    [--json]
                    [shared plot options]

                Options
                =====================================================
                -a/--alignment             FASTA alignment file
                                           (computes distances
                                           internally)

                --distance-matrix          pre-computed distance
                                           matrix as CSV (taxon
                                           labels as row/column
                                           headers). One of -a or
                                           --distance-matrix is
                                           required.

                -o/--output                output figure path
                                           (required)

                --metric                   distance metric when
                                           using -a: identity,
                                           p-distance, or jc
                                           (Jukes-Cantor)
                                           (Default: p-distance)

                --max-splits               maximum number of splits
                                           for visualization
                                           (Default: 30)

                --json                     optional argument to
                                           output results as JSON
                """
            ),
        )
        parser.add_argument(
            "-a", "--alignment",
            type=str,
            default=None,
            required=False,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--distance-matrix",
            type=str,
            default=None,
            required=False,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "-o", "--output",
            type=str,
            required=True,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--metric",
            type=str,
            default="p-distance",
            choices=["identity", "p-distance", "jc"],
            required=False,
            help=SUPPRESS,
            metavar="",
        )
        parser.add_argument(
            "--max-splits",
            type=int,
            default=30,
            required=False,
            help=SUPPRESS,
            metavar="",
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, NeighborNet)

    @staticmethod
    def quartet_pie(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Draw a phylogram with pie charts at internal nodes showing
                quartet concordance proportions.

                In native mode (-g provided), computes gene concordance
                factors (gCF, gDF1, gDF2) from a species tree and gene
                trees via bipartition matching. In ASTRAL mode (no -g),
                parses q1/q2/q3 annotations from ASTRAL -t 2 output or
                wASTRAL --support 3 output.

                Pie slices show: concordant (blue), discordant alt 1
                (red), discordant alt 2 (gray).

                Aliases:
                  quartet_pie, qpie, quartet_pie_chart
                Command line interfaces:
                  pk_quartet_pie, pk_qpie

                Usage:
                phykit quartet_pie -t <tree> [-g <gene_trees>] -o <output>
                  [--annotate] [--csv <file>] [--json]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]

                Options
                =====================================================
                -t/--tree                   species tree file (required)

                -g/--gene-trees             gene trees file, one Newick
                                            tree per line (optional;
                                            if omitted, ASTRAL -t 2 or
                                            wASTRAL --support 3
                                            annotations are parsed)

                -o/--output                 output figure path (required;
                                            supports .png, .pdf, .svg)

                --annotate                  show gCF/gDF values as text
                                            near each pie chart

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for tip labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors for
                                            concordant, disc1, disc2
                                            (default: "#2b8cbe,#d62728,
                                            #969696")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --branch-labels             show concordant gene count
                                            above and LPP support below
                                            each internal branch
                                            (PhyTop-style)

                --csv                       output per-branch concordance
                                            values as a CSV file

                --pie-size                  scale factor for pie chart
                                            size (default: 1.0; use
                                            2.0 for double, 0.5 for
                                            half, etc.)

                --json                      optional argument to output
                                            per-node concordance as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-g", "--gene-trees", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-o", "--output", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--annotate", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument(
            "--branch-labels", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument(
            "--csv", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--pie-size", type=float, required=False, default=1.0,
            help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, QuartetPie)

    @staticmethod
    def quartet_network(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Quartet-based network inference (NANUQ-style).

                Computes quartet concordance factors from gene trees,
                classifies each quartet as tree-like, hybrid, or
                unresolved using two chi-squared tests (matching the
                NANUQ algorithm from MSCquartets), and optionally
                visualizes the result as a species tree with
                reticulation edges overlaid.

                Star test: Pearson chi-squared against uniform (1/3
                each). If p > beta, the quartet is unresolved.
                Tree test (T3): G-test under the MSC tree model.
                If p > alpha, the quartet is tree-like.
                Otherwise the quartet shows hybrid signal.

                Polytomies (collapsed branches) in gene trees are
                handled conservatively: bipartitions from polytomous
                nodes are excluded, so quartets spanning a polytomy
                are treated as unresolved rather than misclassified.
                Trifurcating roots (standard unrooted Newick) are
                not affected.

                Input can be either:
                1) a file with one Newick tree per line, or
                2) a file with one tree-file path per line.

                Aliases:
                  quartet_network, quartet_net, qnet, nanuq
                Command line interfaces:
                  pk_quartet_network, pk_quartet_net, pk_qnet,
                  pk_nanuq

                Usage:
                phykit quartet_network -t/--trees <trees>
                    [--alpha 0.05] [--beta 0.95]
                    [--missing-taxa error|shared]
                    [--plot-output <file>]
                    [--fig-width <float>] [--fig-height <float>]
                    [--dpi <int>] [--no-title] [--title <str>]
                    [--legend-position <str>]
                    [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                    [--title-fontsize <float>] [--axis-fontsize <float>]
                    [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                    [--json]

                Options
                =====================================================
                -t/--trees                 file containing trees or
                                           tree paths

                --alpha                    significance level for the
                                           T3 tree model test
                                           (Default: 0.05)

                --beta                     threshold for the star tree
                                           test; quartets with
                                           p_star > beta are called
                                           unresolved
                                           (Default: 0.95)

                --missing-taxa             how to handle mismatched
                                           taxa across trees:
                                           error or shared
                                           (Default: error)

                --plot-output              output filename for the
                                           quartet network plot
                                           (optional)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                     optional argument to output
                                           results as JSON
                """
            ),
        )
        parser.add_argument("-t", "--trees", type=str, required=True, help=SUPPRESS)
        parser.add_argument(
            "--alpha",
            type=float,
            default=0.05,
            required=False,
            help=SUPPRESS,
        )
        parser.add_argument(
            "--beta",
            type=float,
            default=0.95,
            required=False,
            help=SUPPRESS,
        )
        parser.add_argument(
            "--missing-taxa",
            type=str,
            choices=["error", "shared"],
            default="error",
            required=False,
            help=SUPPRESS,
        )
        parser.add_argument(
            "--plot-output",
            type=str,
            default=None,
            required=False,
            help=SUPPRESS,
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, QuartetNetwork)

    @staticmethod
    def network_signal(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Phylogenetic signal on a network.

                Measures phylogenetic signal (Bloomberg's K and/or
                Pagel's lambda) on a phylogenetic network by
                incorporating hybrid edges inferred from quartet
                concordance factors.

                Hybrid edges can be specified directly (--hybrid)
                or loaded from a quartet network JSON file
                (--quartet-json).

                Polytomies (collapsed branches) in the input tree
                are represented as star topologies in the network
                VCV, which correctly models unresolved relationships
                as equal covariance among all children.

                Aliases:
                  network_signal, netsig, net_signal
                Command line interfaces:
                  pk_network_signal, pk_netsig, pk_net_signal

                Usage:
                phykit network_signal -t <tree> -d <trait_data>
                    (--hybrid <P H1 H2 gamma> | --quartet-json <file>)
                    [--method both|blombergs_k|lambda]
                    [--permutations 1000] [-v/--verbose] [--json]

                Options
                =====================================================
                -t/--tree                   a phylogeny file
                                            (required)

                -d/--trait-data             tab-delimited trait data
                                            file (required)

                --hybrid                    hybrid specification:
                                            parent hybrid child1
                                            child2 gamma (nargs +)

                --quartet-json              path to quartet network
                                            JSON output file

                --method                    which signal measure to
                                            compute: both, blombergs_k,
                                            or lambda
                                            (Default: both)

                --permutations              number of permutations
                                            for significance testing
                                            (Default: 1000)

                -v/--verbose                print detailed output

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("-t", "--tree", type=str, required=True, help=SUPPRESS)
        parser.add_argument("-d", "--trait-data", type=str, required=True, help=SUPPRESS)
        hybrid_group = parser.add_mutually_exclusive_group(required=True)
        hybrid_group.add_argument(
            "--hybrid", nargs="+", type=str, help=SUPPRESS
        )
        hybrid_group.add_argument(
            "--quartet-json", type=str, help=SUPPRESS
        )
        parser.add_argument(
            "--method",
            type=str,
            choices=["both", "blombergs_k", "lambda"],
            default="both",
            required=False,
            help=SUPPRESS,
        )
        parser.add_argument(
            "--permutations",
            type=int,
            default=1000,
            required=False,
            help=SUPPRESS,
        )
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        _add_json_argument(parser)
        _run_service(parser, argv, NetworkSignal)

    @staticmethod
    def ltt(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Lineage-through-time plot and gamma statistic.

                Computes the Pybus & Harvey (2000) gamma statistic
                to test for temporal variation in diversification
                rates.  Under a constant-rate pure-birth process,
                gamma ~ N(0,1).  Negative values indicate early
                diversification (decelerating); positive values
                indicate late diversification (accelerating).

                Optionally generates a lineage-through-time plot.

                Aliases:
                  ltt, gamma_stat, gamma
                Command line interfaces:
                  pk_ltt, pk_gamma_stat, pk_gamma

                Usage:
                phykit ltt -t <tree> [-v/--verbose]
                    [--plot-output <file>]
                    [--fig-width <float>] [--fig-height <float>]
                    [--dpi <int>] [--no-title] [--title <str>]
                    [--legend-position <str>]
                    [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                    [--title-fontsize <float>] [--axis-fontsize <float>]
                    [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                    [--json]

                Options
                =====================================================
                -t/--tree                   a rooted phylogeny file
                                            with branch lengths
                                            (required)

                -v/--verbose                print branching times
                                            and LTT data points

                --plot-output               output filename for the
                                            lineage-through-time
                                            plot (optional)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("-t", "--tree", type=str, required=True, help=SUPPRESS)
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument(
            "--plot-output",
            type=str,
            default=None,
            required=False,
            help=SUPPRESS,
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, LTT)

    @staticmethod
    def prune_tree(argv):
        parser = _new_parser(
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
                phykit prune_tree <tree> <list_of_taxa> [-o/--output <output_file>
                -k/--keep] [--json]

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

                -k/--keep                   optional argument. If used
                                            instead of pruning taxa in
                                            <list_of_taxa>, keep them

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("list_of_taxa", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        parser.add_argument(
            "-k", "--keep", type=str2bool, nargs="?", default=False, help=SUPPRESS
        )
        _add_json_argument(parser)
        _run_service(parser, argv, PruneTree)

    @staticmethod
    def subtree_prune_regraft(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Generate all possible SPR (Subtree Pruning and Regrafting)
                rearrangements for a specified subtree on a tree.

                The subtree is identified by specifying one or more taxa
                whose MRCA defines the clade to prune. The pruned subtree
                is then regrafted onto every other branch in the remaining
                tree, producing one Newick tree per regraft position.

                Aliases:
                  subtree_prune_regraft, spr
                Command line interfaces:
                  pk_subtree_prune_regraft, pk_spr

                Usage:
                phykit subtree_prune_regraft -t <tree> --subtree <taxa>
                  [-o/--output <output_file>] [--json]

                Options
                =====================================================
                -t/--tree                   input tree file in Newick
                                            format (required)

                --subtree                   comma-separated list of
                                            taxa defining the subtree
                                            to prune (MRCA resolved),
                                            or a single-column file
                                            with one taxon per line
                                            (required)

                -o/--output                 output file for SPR trees
                                            (one Newick per line).
                                            If omitted, prints to
                                            stdout.

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("-t", "--tree", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--subtree", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("-o", "--output", type=str, default=None, help=SUPPRESS, metavar="")
        _add_json_argument(parser)
        _run_service(parser, argv, Spr)

    @staticmethod
    def transfer_annotations(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Transfer internal node annotations from one tree onto
                another. Matches nodes by bipartition (descendant taxa
                set) and copies the annotation labels.

                Typical use case: transfer wASTRAL support annotations
                (q1/q2/q3, pp1, f1, etc.) from an annotated ASTRAL
                tree onto a branch-length-optimized topology from
                RAxML-NG, IQ-TREE, or any other tool. The output tree
                has the target's branch lengths with the source's
                annotations.

                Aliases:
                  transfer_annotations, transfer_annot, annotate_tree
                Command line interfaces:
                  pk_transfer_annotations, pk_transfer_annot, pk_annotate_tree

                Usage:
                phykit transfer_annotations --source <annotated_tree>
                  --target <branch_length_tree> [-o/--output <file>]
                  [--json]

                Options
                =====================================================
                --source                    annotated tree file (e.g.,
                                            wASTRAL output with
                                            --support 3)

                --target                    target tree file with
                                            branch lengths to keep
                                            (e.g., RAxML-NG or
                                            IQ-TREE output)

                -o/--output                 output file for the
                                            annotated tree (default:
                                            target file + ".annotated")

                --json                      output results as JSON
                """
            ),
        )
        parser.add_argument("--source", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("--target", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("-o", "--output", type=str, default=None, help=SUPPRESS, metavar="")
        _add_json_argument(parser)
        _run_service(parser, argv, TransferAnnotations)

    @staticmethod
    def relative_rate_test(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Tajima's relative rate test.

                Tests whether two ingroup lineages evolve at
                equal rates relative to an outgroup.  The tree
                must be rooted with a single outgroup taxon.
                All pairwise ingroup comparisons are performed
                with Bonferroni and BH-FDR correction.

                Provide either a single alignment (-a) or a
                file listing multiple alignment paths (-l) for
                batch (multi-gene) analysis.

                Aliases:
                  relative_rate_test, rrt, tajima_rrt
                Command line interfaces:
                  pk_relative_rate_test, pk_rrt, pk_tajima_rrt

                Usage:
                phykit relative_rate_test (-a <alignment> | -l <alignment_list>) -t <tree> [-v/--verbose]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                -a/--alignment              a single alignment file

                -l/--alignment-list         a file listing alignment
                                            paths (one per line)

                -t/--tree                   a rooted tree file
                                            (required)

                -v/--verbose                print detailed output

                --plot-output               save pairwise p-value
                                            heatmap to this path

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        aln_group = parser.add_mutually_exclusive_group()
        aln_group.add_argument(
            "-a", "--alignment", type=str, required=False, help=SUPPRESS, metavar=""
        )
        aln_group.add_argument(
            "-l", "--alignment-list", dest="alignment_list", type=str, required=False, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument(
            "--plot-output", type=str, required=False, help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, RelativeRateTest)

    @staticmethod
    def threshold_model(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Felsenstein (2012) threshold model.

                Estimates the correlation between two traits (binary
                discrete and/or continuous) using a latent-liability
                Brownian motion model and MCMC sampling.  Binary
                discrete characters are modelled as arising from
                continuous liabilities that cross a threshold at 0.

                This is the equivalent of phytools::threshBayes in R.

                Aliases:
                  threshold_model, threshold, thresh, threshbayes,
                  thresh_bayes
                Command line interfaces:
                  pk_threshold_model, pk_threshold, pk_thresh,
                  pk_threshbayes, pk_thresh_bayes

                Usage:
                phykit threshold_model -t <tree> -d <trait_data>
                    --traits <trait1,trait2> --types <type1,type2>
                    [--ngen 100000] [--sample 100] [--burnin 0.2]
                    [--seed <int>] [--plot <file>]
                    [--fig-width <float>] [--fig-height <float>]
                    [--dpi <int>] [--no-title] [--title <str>]
                    [--legend-position <str>]
                    [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                    [--title-fontsize <float>] [--axis-fontsize <float>]
                    [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                    [--json]

                Options
                =====================================================
                -t/--tree                   a rooted phylogeny file
                                            with branch lengths
                                            (required)

                -d/--trait-data             tab-delimited trait file
                                            with header row (required)

                --traits                    comma-separated pair of
                                            trait column names

                --types                     comma-separated pair of
                                            trait types (discrete or
                                            continuous)

                --ngen                      number of MCMC generations
                                            (default: 100000)

                --sample                    sample frequency
                                            (default: 100)

                --burnin                    burn-in fraction
                                            (default: 0.2)

                --seed                      random seed for
                                            reproducibility

                --plot                      output filename for
                                            trace plots (optional)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("-t", "--tree", type=str, required=True, help=SUPPRESS)
        parser.add_argument(
            "-d", "--trait-data", dest="trait_data", type=str, required=True, help=SUPPRESS
        )
        parser.add_argument(
            "--traits", type=str, required=True, help=SUPPRESS
        )
        parser.add_argument(
            "--types", type=str, required=True, help=SUPPRESS
        )
        parser.add_argument(
            "--ngen", type=int, default=100000, required=False, help=SUPPRESS
        )
        parser.add_argument(
            "--sample", type=int, default=100, required=False, help=SUPPRESS
        )
        parser.add_argument(
            "--burnin", type=float, default=0.2, required=False, help=SUPPRESS
        )
        parser.add_argument(
            "--seed", type=int, default=None, required=False, help=SUPPRESS
        )
        parser.add_argument(
            "--plot", dest="plot_output", type=str, default=None, required=False, help=SUPPRESS
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, ThresholdModel)

    @staticmethod
    def rename_tree_tips(argv):
        parser = _new_parser(
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
                    [-o/--output <output_file>] [--json]

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

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-i", "--idmap", type=str, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, RenameTreeTips)

    @staticmethod
    def kf_distance(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate Kuhner-Felsenstein (KF) branch score distance
                between two trees.

                Unlike Robinson-Foulds distance which only considers topology,
                KF distance incorporates both topology and branch length
                differences. The KF distance is calculated as:
                KF = sqrt( sum_over_all_splits( (b1_i - b2_i)^2 ) )
                where b1_i and b2_i are branch lengths for split i in each
                tree. Splits absent from one tree use branch length 0.

                PhyKIT will print out
                col 1: the plain KF distance and
                col 2: the normalized KF distance.

                KF distances are calculated following Kuhner & Felsenstein,
                Journal of Computational Biology (1994),
                doi: 10.1089/cmb.1994.1.183.

                Aliases:
                  kuhner_felsenstein_distance, kf_distance, kf_dist, kf
                Command line interfaces:
                  pk_kuhner_felsenstein_distance, pk_kf_distance, pk_kf_dist,
                  pk_kf

                Usage:
                phykit kf_distance <tree_file_zero> <tree_file_one> [--json]

                Options
                =====================================================
                <tree_file_zero>            first argument after
                                            function name should be
                                            a tree file

                <tree_file_one>             second argument after
                                            function name should be
                                            a tree file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree_zero", type=str, help=SUPPRESS)
        parser.add_argument("tree_one", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, KuhnerFelsensteinDistance)

    @staticmethod
    def rf_distance(argv):
        parser = _new_parser(
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
                phykit robinson_foulds_distance <tree_file_zero> <tree_file_one> [--json]

                Options
                =====================================================
                <tree_file_zero>            first argument after 
                                            function name should be
                                            a tree file

                <tree_file_one>             second argument after 
                                            function name should be
                                            a tree file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree_zero", type=str, help=SUPPRESS)
        parser.add_argument("tree_one", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, RobinsonFouldsDistance)

    @staticmethod
    def root_tree(argv):
        parser = _new_parser(
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
                    [-o/--output <output_file>] [--json]

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

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-r", "--root", type=str, required=True, help=SUPPRESS)
        parser.add_argument("-o", "--output", type=str, required=False, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, RootTree)

    @staticmethod
    def spurious_sequence(argv):
        parser = _new_parser(
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
                phykit spurious_sequence <file> [-f 20] [--json]

                Options
                =====================================================
                <file>                      first argument after 
                                            function name should be
                                            an tree file

                -f/--factor                 factor to multiply median
                                            branch length by to calculate
                                            the threshold of long branches.
                                            (Default: 20)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument("-f", "--factor", type=float, required=False, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, SpuriousSequence)

    @staticmethod
    def terminal_branch_stats(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate summary statistics for terminal branch lengths in a phylogeny.

                Terminal branch lengths can be useful for phylogeny diagnostics.

                To obtain all terminal branch lengths, use the -v/--verbose option. 

                Aliases:
                  terminal_branch_stats, tbs
                Command line interfaces:
                  pk_terminal_branch_stats, pk_tbs

                Usage:
                phykit terminal_branch_stats <tree> [-v/--verbose] [--json]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                -v/--verbose                optional argument to print
                                            all terminal branch lengths

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        _add_json_argument(parser)
        _run_service(parser, argv, TerminalBranchStats)

    @staticmethod
    def tip_labels(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Prints the tip labels (or names) a phylogeny.

                Aliases:
                  tip_labels, tree_labels; labels; tl
                Command line interfaces: 
                  pk_tip_labels, pk_tree_labels; pk_labels; pk_tl

                Usage:
                phykit tip_labels <tree> [--json]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, TipLabels)

    @staticmethod
    def tip_to_tip_distance(argv):
        parser = _new_parser(
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
                phykit tip_to_tip_distance <tree_file> <tip_1> <tip_2> [--json]
                phykit tip_to_tip_distance <tree_file> --all-pairs [--plot] [--plot-output <path>]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

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

                --all-pairs                 optional argument to report
                                            all pairwise tip distances

                --plot                      optional argument to save a
                                            clustered distance heatmap
                                            (requires --all-pairs)

                --plot-output               output path for heatmap
                                            (default: tip_to_tip_distance_heatmap.png)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree_zero", type=str, help=SUPPRESS)
        parser.add_argument("tip_1", type=str, nargs="?", help=SUPPRESS)
        parser.add_argument("tip_2", type=str, nargs="?", help=SUPPRESS)
        parser.add_argument("--all-pairs", action="store_true", required=False, help=SUPPRESS)
        parser.add_argument("--plot", action="store_true", required=False, help=SUPPRESS)
        parser.add_argument(
            "--plot-output",
            type=str,
            default="tip_to_tip_distance_heatmap.png",
            required=False,
            help=SUPPRESS,
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, TipToTipDistance)

    @staticmethod
    def tip_to_tip_node_distance(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate distance between two tips (or leaves) in a phylogeny.

                Distance is measured by the number of nodes between one tip
                and another.

                Aliases:
                  tip_to_tip_node_distance, t2t_node_dist, t2t_nd
                Command line interfaces: 
                  pk_tip_to_tip_node_distance, pk_t2t_node_dist, pk_t2t_nd

                Usage:
                phykit tip_to_tip_node_distance <tree_file> <tip_1> <tip_2> [--json]

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

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree_zero", type=str, help=SUPPRESS)
        parser.add_argument("tip_1", type=str, help=SUPPRESS)
        parser.add_argument("tip_2", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, TipToTipNodeDistance)

    @staticmethod
    def total_tree_length(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate total tree length, which is a sum of all branches. 

                Aliases:
                  total_tree_length, tree_len
                Command line interfaces: 
                  pk_total_tree_length, pk_tree_len

                Usage:
                phykit total_tree_length <tree> [--json]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, TotalTreeLength)

    @staticmethod
    def treeness(argv):
        parser = _new_parser(
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
                phykit treeness <tree> [--json]

                Options
                =====================================================
                <tree>                      first argument after 
                                            function name should be
                                            a tree file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("tree", type=str, help=SUPPRESS)
        _add_json_argument(parser)
        _run_service(parser, argv, Treeness)

    ## Alignment and tree functions
    @staticmethod
    def saturation(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Calculate saturation for a given tree and alignment.

                Saturation is defined as sequences in multiple sequence
                alignments that have undergone numerous substitutions such
                that the distances between taxa are underestimated.

                Data with no saturation will have a value of 1. The closer
                the value is to 1, the less saturated the data.

                This function outputs two values (as of v1.19.9). The first
                value is the saturation value and the second column is the absolute
                value of saturation minus 1. Thus, lower values in the second column
                are indicative of values closer to one and, thus, less saturation.

                Saturation is calculated following Philippe et al., PLoS 
                Biology (2011), doi: 10.1371/journal.pbio.1000602.

                Aliases: 
                  saturation, sat
                Command line interfaces:
                  pk_saturation, pk_sat

                Usage:
                phykit saturation -a <alignment> -t <tree> [-v/--verbose] [-e/--exclude_gaps] [--plot] [--plot-output <path>]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                -a/--alignment              an alignment file
                
                -t/--tree                   a tree file

                -e/--exclude_gaps           if a site has a gap, ignore it

                -v/--verbose                print out patristic distances
                                            and uncorrected distances used
                                            to determine saturation

                --plot                      optional argument to save a
                                            saturation scatter plot

                --plot-output               output path for saturation plot
                                            (default: saturation_plot.png)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-a", "--alignment", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument(
            "-e", "--exclude_gaps", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument("--plot", action="store_true", required=False, help=SUPPRESS)
        parser.add_argument(
            "--plot-output",
            type=str,
            default="saturation_plot.png",
            required=False,
            help=SUPPRESS,
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, Saturation)

    @staticmethod
    def treeness_over_rcv(argv):
        parser = _new_parser(
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
                phykit treeness_over_rcv -a/--alignment <alignment> -t/--tree <tree> [--json]

                Options
                =====================================================
                -a/--alignment              an alignment file
                
                -t/--tree                   a tree file

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-a", "--alignment", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        _add_json_argument(parser)
        _run_service(parser, argv, TreenessOverRCV)

    @staticmethod
    def evo_tempo_map(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Detect rate-topology associations by comparing branch length
                distributions between concordant and discordant gene trees at
                each species tree branch.

                Under the multispecies coalescent, discordant gene trees should
                have shorter internal branches near the discordant node. Deviations
                suggest substitution rate heterogeneity correlated with topology
                (adaptive evolution, different selective pressures, or model
                misspecification).

                For each internal branch of the species tree, gene trees are
                classified as concordant or discordant via bipartition matching.
                The homologous branch length is extracted from each gene tree
                and the two groups are compared using Mann-Whitney U and
                permutation tests. P-values are corrected for multiple testing
                using Benjamini-Hochberg FDR.

                A global treeness (internal/total branch length ratio) comparison
                is also reported.

                Aliases:
                  evo_tempo_map, etm
                Command line interfaces:
                  pk_evo_tempo_map, pk_etm

                Usage:
                phykit evo_tempo_map -t/--tree <tree> -g/--gene-trees <gene_trees>
                    [--plot <output>] [-v/--verbose]
                    [--fig-width <float>] [--fig-height <float>]
                    [--dpi <int>] [--no-title] [--title <str>]
                    [--legend-position <str>]
                    [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                    [--title-fontsize <float>] [--axis-fontsize <float>]
                    [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                    [--json]

                Options
                =====================================================
                -t/--tree                   a species tree file

                -g/--gene-trees             multi-Newick file of gene trees
                                            with branch lengths

                --plot                      optional output path for
                                            box/strip plot (PNG)

                -v/--verbose                print per-gene-tree details

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-g", "--gene-trees", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--plot", dest="plot_output", type=str, required=False,
            default=None, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, EvoTempoMap)

    @staticmethod
    def discordance_asymmetry(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Estimate the minimum number of reticulation events and
                localize where hybridization likely occurred. Tests
                whether the two discordant NNI topologies at each branch
                are significantly asymmetric (indicating introgression
                rather than ILS). Note: identifies which lineages
                exchanged genes but not the direction of flow.

                Under incomplete lineage sorting (ILS) alone, the two minor
                NNI alternatives (gDF1 and gDF2) should be equally frequent.
                When they are significantly asymmetric, it suggests
                introgression or gene flow between specific lineages.

                For each internal branch, a two-sided binomial test (H0:
                P(alt1) = 0.5) is applied. P-values are corrected for
                multiple testing using Benjamini-Hochberg FDR.

                Aliases:
                  discordance_asymmetry, disc_asym, da
                Command line interfaces:
                  pk_discordance_asymmetry, pk_disc_asym, pk_da

                Usage:
                phykit discordance_asymmetry -t/--tree <tree> -g/--gene-trees <gene_trees>
                    [--plot <output>] [-v/--verbose]
                    [--fig-width <float>] [--fig-height <float>]
                    [--dpi <int>] [--no-title] [--title <str>]
                    [--legend-position <str>]
                    [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                    [--title-fontsize <float>] [--axis-fontsize <float>]
                    [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                    [--json]

                Options
                =====================================================
                -t/--tree                   a species tree file

                -g/--gene-trees             multi-Newick file of gene trees

                --plot                      optional output path for
                                            asymmetry phylogram (PNG)

                -v/--verbose                print per-branch details

                --annotate                  show gCF values on the
                                            plot near each branch

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-g", "--gene-trees", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--plot", dest="plot_output", type=str, required=False,
            default=None, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-v", "--verbose", action="store_true", required=False, help=SUPPRESS
        )
        parser.add_argument(
            "--annotate", action="store_true", required=False, default=False,
            help=SUPPRESS,
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, DiscordanceAsymmetry)

    @staticmethod
    def hybridization(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Estimate the minimum number of reticulation (hybridization)
                events and localize where hybridization likely occurred on
                a species tree.

                For each internal branch, the four-group decomposition is
                used to count concordant and two NNI-alternative topologies
                across gene trees. A two-sided binomial test detects
                asymmetric discordance (a hallmark of hybridization or
                introgression). P-values are corrected using Benjamini-
                Hochberg FDR. Branches with significant asymmetry are
                flagged as putative reticulation events.

                Aliases:
                  hybridization, hybrid, reticulation
                Command line interfaces:
                  pk_hybridization, pk_hybrid, pk_reticulation

                Usage:
                phykit hybridization -t/--tree <tree> -g/--gene-trees <gene_trees>
                    [--support <float>] [--alpha <float>]
                    [--plot <output>] [--json]
                    [--fig-width <float>] [--fig-height <float>]
                    [--dpi <int>] [--no-title] [--title <str>]
                    [--legend-position <str>]
                    [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                    [--title-fontsize <float>] [--axis-fontsize <float>]
                    [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]

                Options
                =====================================================
                -t/--tree                   a species tree file

                -g/--gene-trees             multi-Newick file of gene trees

                --support                   collapse gene tree branches
                                            below this support value
                                            before topology determination

                --alpha                     significance threshold for
                                            asymmetry tests after FDR
                                            correction (default: 0.05)

                --plot                      optional output path for
                                            hybridization score phylogram
                                            (PNG)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-g", "--gene-trees", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--support", type=float, default=None, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--alpha", type=float, default=0.05, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--plot", dest="plot_output", type=str, required=False,
            default=None, help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, Hybridization)

    @staticmethod
    def spectral_discordance(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Spectral discordance decomposition — decompose gene tree
                space via PCA on a bipartition presence/absence (or
                branch-length) matrix, with spectral clustering and
                automatic cluster detection via the eigengap heuristic.

                Each gene tree is encoded as a vector over the union of
                all bipartitions observed across gene trees. PCA reveals
                the axes of topological variation, with loading vectors
                identifying which bipartitions drive each PC. Spectral
                clustering groups genes sharing similar topologies.

                Two metrics are available:
                - nrf (default): binary presence/absence (normalized RF)
                - wrf: branch-length weighted

                Polytomies (collapsed branches) in gene trees are
                handled conservatively: splits from polytomous nodes
                are excluded from the bipartition matrix since they
                represent unresolved relationships. Trifurcating
                roots (standard unrooted Newick) are not affected.

                Aliases:
                  spectral_discordance, spec_disc, sd
                Command line interfaces:
                  pk_spectral_discordance, pk_spec_disc, pk_sd

                Usage:
                phykit spectral_discordance -g <gene_trees> [-t <tree>] [--metric nrf|wrf] [--clusters K] [--n-pcs N] [--top-loadings N] [--plot <prefix>]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]
                  [--json]

                Options
                =====================================================
                -g/--gene-trees             file of gene trees (one
                                            Newick per line, or file
                                            of filenames)

                -t/--tree                   species tree (optional; flags
                                            species-tree bipartitions in
                                            loading output)

                --metric                    distance metric: nrf or wrf
                                            (default: nrf)

                --clusters                  override auto-detected K

                --n-pcs                     number of PCs to report
                                            (default: min(10, G-1))

                --top-loadings              top bipartitions per PC
                                            (default: 5)

                --plot                      output prefix for plots
                                            (generates _scatter.png and
                                            _eigengap.png)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      output results as JSON
                """
            ),
        )
        parser.add_argument(
            "-g", "--gene-trees", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--metric", type=str, required=False, default="nrf",
            choices=["nrf", "wrf"], help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--clusters", type=int, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--n-pcs", type=int, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--top-loadings", type=int, required=False, default=5,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--plot", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, SpectralDiscordance)

    @staticmethod
    def trait_rate_map(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Trait Rate Map — estimate per-branch evolutionary
                rates for a continuous trait and display them as a
                branch-colored phylogram.

                Ancestral states are reconstructed via Felsenstein's
                weighted-average method (inverse-branch-length
                weighting, postorder traversal). Per-branch rate is
                the squared standardized contrast:
                  rate = (child_val - parent_val)^2 / branch_length

                Input is a phylogenetic tree and either:
                  (a) a two-column TSV (taxon<tab>value, no header), or
                  (b) a multi-column TSV with header (use --trait to
                      select a column)

                Aliases:
                  trait_rate_map, rate_map, branch_rates
                Command line interfaces:
                  pk_trait_rate_map, pk_rate_map, pk_branch_rates

                Usage:
                phykit trait_rate_map -t <tree> -d <trait_data> -o <output>
                  [--trait <column>] [--json]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]

                Options
                =====================================================
                -t/--tree                   a tree file

                -d/--trait_data             tab-delimited trait file
                                            (two-column: taxon<tab>value,
                                            no header; or multi-column
                                            with header when --trait is
                                            used)

                -o/--output                 output plot file path
                                            (required)

                --trait                     column name to use from a
                                            multi-column trait file
                                            (if omitted, two-column
                                            format is expected)

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to also
                                            output results as JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-d", "--trait_data", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-o", "--output", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--trait", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, TraitRateMap)

    @staticmethod
    def tree_space(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Tree space visualization — visualize how gene trees
                cluster in topology space using MDS, t-SNE, or UMAP
                on pairwise tree distance matrices.

                Gene trees are compared pairwise using Robinson-Foulds
                (RF) or Kuhner-Felsenstein (KF) distance. The resulting
                distance matrix is embedded into 2D using MDS, t-SNE,
                or UMAP. Spectral clustering with eigengap heuristic
                auto-detects topological groups.

                An optional species tree can be highlighted as a
                distinct marker in the plot.

                Aliases:
                  tree_space, tspace, tree_landscape
                Command line interfaces:
                  pk_tree_space, pk_tspace, pk_tree_landscape

                Usage:
                phykit tree_space -t <trees> -o <output>
                  [--metric rf|kf] [--method mds|tsne|umap]
                  [--species-tree <file>] [--k <int>] [--seed <int>]
                  [--heatmap] [--distance-matrix <file>]
                  [--json]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>]

                Options
                =====================================================
                -t/--trees              file with gene trees (one
                                        Newick per line, or one path
                                        per line)

                -o/--output             output figure path (.png,
                                        .pdf, .svg)

                --metric                distance metric: rf
                                        (Robinson-Foulds, default)
                                        or kf (Kuhner-Felsenstein)

                --method                dimensionality reduction:
                                        mds (default), tsne, or umap

                --species-tree          optional species tree to
                                        highlight in the plot

                --k                     number of clusters (auto-
                                        detected via eigengap if
                                        omitted)

                --seed                  random seed for
                                        reproducibility (t-SNE/UMAP)

                --heatmap               draw a clustered distance
                                        heatmap instead of a scatter
                                        plot

                --distance-matrix       output pairwise distance
                                        matrix as CSV file

                --json                  output structured JSON
                """
            ),
        )
        parser.add_argument(
            "-t", "--trees", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-o", "--output", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--metric", type=str, required=False, default="rf",
            choices=["rf", "kf"], help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--method", type=str, required=False, default="mds",
            choices=["mds", "tsne", "umap"], help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--species-tree", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--k", type=int, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--seed", type=int, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--heatmap", action="store_true", default=False,
            help=SUPPRESS,
        )
        parser.add_argument(
            "--distance-matrix", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, TreeSpace)

    @staticmethod
    def taxon_groups(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Determine which tree or FASTA files share the same set
                of taxa. Reads a file listing paths to gene trees or
                alignments and groups them by their taxon set (exact
                match). Reports groups sorted by size (largest first),
                with the taxa present in each group.

                Useful for identifying subsets of genes with identical
                taxon sampling for concatenation or comparative analysis.

                Aliases:
                  taxon_groups, tgroups, shared_taxa
                Command line interfaces:
                  pk_taxon_groups, pk_tgroups, pk_shared_taxa

                Usage:
                phykit taxon_groups -l <file> [-f trees|fasta] [--json]

                Options
                =====================================================
                -l/--list                   file listing paths to gene
                                            trees or FASTA files (one
                                            per line). Blank lines and
                                            lines starting with # are
                                            skipped.

                -f/--format                 input format: trees (Newick)
                                            or fasta (FASTA alignment).
                                            default: trees

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("-l", "--list", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("-f", "--format", type=str, default="trees", choices=["trees", "fasta"], help=SUPPRESS, metavar="")
        _add_json_argument(parser)
        _run_service(parser, argv, TaxonGroups)

    @staticmethod
    def occupancy_filter(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Filter alignments and/or trees by cross-file taxon
                occupancy. Counts how many files each taxon appears in
                and retains only taxa meeting a minimum threshold.
                Outputs filtered copies of each input file.

                For FASTA files, removes sequences of filtered taxa.
                For tree files, prunes tips of filtered taxa.

                Aliases:
                  occupancy_filter, occ_filter, filter_occupancy
                Command line interfaces:
                  pk_occupancy_filter, pk_occ_filter, pk_filter_occupancy

                Usage:
                phykit occupancy_filter -l <file_list>
                  [-f/--format fasta|trees] [-t/--threshold <int>]
                  [-o/--output-dir <dir>] [--suffix <str>] [--json]

                Options
                =====================================================
                -l/--list                   file listing paths to
                                            alignment or tree files,
                                            one per line (required)

                -f/--format                 input file format: fasta
                                            or trees (default: fasta)

                -t/--threshold              minimum occupancy to retain
                                            a taxon. Values between 0
                                            and 1 are treated as a
                                            fraction (e.g., 0.5 = 50%
                                            of files). Values >= 1 are
                                            treated as an absolute
                                            count. (default: 0.5)

                -o/--output-dir             directory for filtered
                                            output files (default:
                                            same directory as input)

                --suffix                    suffix added to output
                                            filenames before the
                                            extension (default:
                                            ".filtered")

                --json                      output results as JSON
                """
            ),
        )
        parser.add_argument("-l", "--list", type=str, required=True, help=SUPPRESS, metavar="")
        parser.add_argument("-f", "--format", type=str, default="fasta", choices=["fasta", "trees"], help=SUPPRESS, metavar="")
        parser.add_argument("-t", "--threshold", type=float, default=0.5, help=SUPPRESS, metavar="")
        parser.add_argument("-o", "--output-dir", type=str, default=None, help=SUPPRESS, metavar="")
        parser.add_argument("--suffix", type=str, default=".filtered", help=SUPPRESS, metavar="")
        _add_json_argument(parser)
        _run_service(parser, argv, OccupancyFilter)

    ### Helper commands
    @staticmethod
    def create_concatenation_matrix(argv):
        parser = _new_parser(
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
                        - column 1: alignment name
                        - column 2: # of taxa present
                        - column 3: # of taxa missing
                        - column 4: fraction of occupancy
                        - column 5: names of missing taxa (; separated)

                Aliases:
                  create_concatenation_matrix, create_concat, cc
                Command line interfaces: 
                  pk_create_concatenation_matrix, pk_create_concat, pk_cc

                Usage:
                phykit create_concatenation_matrix -a <file> -p <string>
                  [--threshold <float>] [--plot-occupancy]
                  [--plot-output <path>] [--json]
                  [--fig-width <float>] [--fig-height <float>]
                  [--dpi <int>] [--no-title] [--title <str>]
                  [--legend-position <str>]
                  [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
                  [--title-fontsize <float>] [--axis-fontsize <float>]
                  [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]

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

                --threshold                 minimum fraction of informative
                                            (non-gap, non-ambiguous) sites
                                            across the concatenated alignment
                                            for a taxon to be included.
                                            Set to 0 to disable filtering.
                                            default: 0

                --plot-occupancy            optional argument to generate
                                            occupancy map figure

                --plot-output               output path for occupancy
                                            figure (supports .png, .pdf,
                                            .svg, .jpg). default:
                                            <prefix>.occupancy.png

                --fig-width                 figure width in inches
                                            (auto-scaled if omitted)

                --fig-height                figure height in inches
                                            (auto-scaled if omitted)

                --dpi                       resolution in DPI
                                            (default: 300)

                --no-title                  hide the plot title

                --title                     custom title text

                --legend-position           legend location (e.g.,
                                            "upper right", "none")

                --ylabel-fontsize           font size for y-axis labels;
                                            0 to hide

                --xlabel-fontsize           font size for x-axis labels;
                                            0 to hide

                --title-fontsize            font size for the title

                --axis-fontsize             font size for axis labels

                --colors                    comma-separated colors
                                            (hex or named, e.g.,
                                            "#ff0000,blue,#00ff00")

                --ladderize                 ladderize (sort) the tree
                                            before plotting

                --cladogram                 draw cladogram (equal branch
                                            lengths, tips aligned)
                                            instead of phylogram

                --circular                  draw circular (radial/fan)
                                            phylogram instead of
                                            rectangular

                --color-file                color annotation file for
                                            tip labels, clade ranges,
                                            and branch colors (iTOL-
                                            inspired TSV format)

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("-a", "--alignment_list", type=str, help=SUPPRESS)
        parser.add_argument("-p", "--prefix", type=str, help=SUPPRESS)
        parser.add_argument("--threshold", type=float, required=False, default=0, help=SUPPRESS)
        parser.add_argument("--plot-occupancy", action="store_true", required=False, help=SUPPRESS)
        parser.add_argument("--plot-output", type=str, required=False, default=None, help=SUPPRESS)
        add_plot_arguments(parser)
        _add_json_argument(parser)
        _run_service(parser, argv, CreateConcatenationMatrix)

    @staticmethod
    def thread_dna(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Thread DNA sequence onto a protein alignment to create a
                codon-based alignment. 
                
                This function requires input alignments are in fasta format.
                Codon alignments are then printed to stdout. Note, paired
                sequences are assumed to have the same name between the 
                protein and nucleotide file. The order does not matter.

                To thread nucleotide sequences over a trimmed amino acid
                alignment, provide PhyKIT with a log file specifying which
                sites have been trimmed and which have been kept. The log
                file must be formatted the same as the log files outputted
                by the alignment trimming toolkit ClipKIT (see -l in ClipKIT
                documentation.) Details about ClipKIT can be seen here:
                https://github.com/JLSteenwyk/ClipKIT.

                If using a ClipKIT log file, the untrimmed protein alignment
                should be provided in the -p/--protein argument.

                Aliases:
                  thread_dna, pal2nal, p2n
                Command line interfaces:
                  pk_thread_dna, pk_pal2nal, pk_p2n

                Usage:
                phykit thread_dna -p <file> -n <file> [-c/--clipkit_log_file
                  <clipkit outputted log file> -s] [--json]

                Options
                =====================================================
                -p/--protein                protein alignment file

                -n/--nucleotide             nucleotide sequence file

                -c/--clipkit_log            clipkit outputted log file

                -s/--stop                   boolean for whether or not
                                            stop codons should be kept. 
                                            If used, stop codons will 
                                            be removed.

                --json                      optional argument to output
                                            results as JSON
                """
            ),
        )
        parser.add_argument("-p", "--protein", type=str, help=SUPPRESS)
        parser.add_argument("-n", "--nucleotide", type=str, help=SUPPRESS)
        parser.add_argument(
            "-c",
            "--clipkit_log_file",
            type=str,
            required=False,
            help=SUPPRESS,
        )
        parser.add_argument(
            "-s", "--stop", type=str2bool, nargs="?", default=True, help=SUPPRESS
        )
        _add_json_argument(parser)
        _run_service(parser, argv, DNAThreader)


def main(argv=None):
    Phykit()


# Alignment-based functions
def alignment_length(argv=None):
    Phykit.alignment_length(sys.argv[1:])


def alignment_length_no_gaps(argv=None):
    Phykit.alignment_length_no_gaps(sys.argv[1:])


def alignment_entropy(argv=None):
    Phykit.alignment_entropy(sys.argv[1:])


def alignment_outlier_taxa(argv=None):
    Phykit.alignment_outlier_taxa(sys.argv[1:])


def column_score(argv=None):
    Phykit.column_score(sys.argv[1:])


def compositional_bias_per_site(argv=None):
    Phykit.compositional_bias_per_site(sys.argv[1:])


def composition_per_taxon(argv=None):
    Phykit.composition_per_taxon(sys.argv[1:])


def evolutionary_rate_per_site(argv=None):
    Phykit.evolutionary_rate_per_site(sys.argv[1:])


def faidx(argv=None):
    Phykit.faidx(sys.argv[1:])


def gc_content(argv=None):
    Phykit.gc_content(sys.argv[1:])


def mask_alignment(argv=None):
    Phykit.mask_alignment(sys.argv[1:])


def plot_alignment_qc(argv=None):
    Phykit.plot_alignment_qc(sys.argv[1:])


def occupancy_per_taxon(argv=None):
    Phykit.occupancy_per_taxon(sys.argv[1:])


def pairwise_identity(argv=None):
    Phykit.pairwise_identity(sys.argv[1:])


def identity_matrix(argv=None):
    Phykit.identity_matrix(sys.argv[1:])


def parsimony_informative_sites(argv=None):
    Phykit.parsimony_informative_sites(sys.argv[1:])


def rcv(argv=None):
    Phykit.rcv(sys.argv[1:])


def rcvt(argv=None):
    Phykit.rcvt(sys.argv[1:])


def rename_fasta_entries(argv=None):
    Phykit.rename_fasta_entries(sys.argv[1:])


def sum_of_pairs_score(argv=None):
    Phykit.sum_of_pairs_score(sys.argv[1:])


def variable_sites(argv=None):
    Phykit.variable_sites(sys.argv[1:])


def alignment_subsample(argv=None):
    Phykit.alignment_subsample(sys.argv[1:])


def dstatistic(argv=None):
    Phykit.dstatistic(sys.argv[1:])


def phylo_gwas(argv=None):
    Phykit.phylo_gwas(sys.argv[1:])


def phylo_anova(argv=None):
    Phykit.phylo_anova(sys.argv[1:])


def phylo_path(argv=None):
    Phykit.phylo_path(sys.argv[1:])


def dfoil(argv=None):
    Phykit.dfoil(sys.argv[1:])


# Tree-based functions
def parsimony_score(argv=None):
    Phykit.parsimony_score(sys.argv[1:])


def character_map(argv=None):
    Phykit.character_map(sys.argv[1:])


def independent_contrasts(argv=None):
    Phykit.independent_contrasts(sys.argv[1:])


def ancestral_state_reconstruction(argv=None):
    Phykit.ancestral_state_reconstruction(sys.argv[1:])


def concordance_asr(argv=None):
    Phykit.concordance_asr(sys.argv[1:])


def chronogram(argv=None):
    Phykit.chronogram(sys.argv[1:])


def dtt(argv=None):
    Phykit.dtt(sys.argv[1:])


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


def faiths_pd(argv=None):
    Phykit.faiths_pd(sys.argv[1:])


def patristic_distances(argv=None):
    Phykit.patristic_distances(sys.argv[1:])


def phylogenetic_signal(argv=None):
    Phykit.phylogenetic_signal(sys.argv[1:])


def trait_correlation(argv=None):
    Phykit.trait_correlation(sys.argv[1:])


def phylogenetic_ordination(argv=None):
    Phykit.phylogenetic_ordination(sys.argv[1:])


def phylogenetic_pca(argv=None):
    Phykit.phylogenetic_ordination(sys.argv[1:])


def phylogenetic_dimreduce(argv=None):
    Phykit.phylogenetic_ordination(sys.argv[1:])


def phylo_heatmap(argv=None):
    Phykit.phylo_heatmap(sys.argv[1:])


def phylomorphospace(argv=None):
    Phykit.phylomorphospace(sys.argv[1:])


def phylogenetic_regression(argv=None):
    Phykit.phylogenetic_regression(sys.argv[1:])


def phylogenetic_glm(argv=None):
    Phykit.phylogenetic_glm(sys.argv[1:])


def phylo_logistic(argv=None):
    Phykit.phylo_logistic(sys.argv[1:])


def stochastic_character_map(argv=None):
    Phykit.stochastic_character_map(sys.argv[1:])


def simmap_summary(argv=None):
    Phykit.simmap_summary(sys.argv[1:])


def cont_map(argv=None):
    Phykit.cont_map(sys.argv[1:])


def density_map(argv=None):
    Phykit.density_map(sys.argv[1:])


def phenogram(argv=None):
    Phykit.phenogram(sys.argv[1:])


def cophylo(argv=None):
    Phykit.cophylo(sys.argv[1:])


def rate_heterogeneity(argv=None):
    Phykit.rate_heterogeneity(sys.argv[1:])


def fit_continuous(argv=None):
    Phykit.fit_continuous(sys.argv[1:])


def fit_discrete(argv=None):
    Phykit.fit_discrete(sys.argv[1:])


def ouwie(argv=None):
    Phykit.ouwie(sys.argv[1:])


def ou_shift_detection(argv=None):
    Phykit.ou_shift_detection(sys.argv[1:])


def polytomy_test(argv=None):
    Phykit.polytomy_test(sys.argv[1:])


def print_tree(argv=None):
    Phykit.print_tree(sys.argv[1:])


def consensus_network(argv=None):
    Phykit.consensus_network(sys.argv[1:])


def neighbor_net(argv=None):
    Phykit.neighbor_net(sys.argv[1:])


def quartet_network(argv=None):
    Phykit.quartet_network(sys.argv[1:])


def quartet_pie(argv=None):
    Phykit.quartet_pie(sys.argv[1:])


def ltt(argv=None):
    Phykit.ltt(sys.argv[1:])


def network_signal(argv=None):
    Phykit.network_signal(sys.argv[1:])


def consensus_tree(argv=None):
    Phykit.consensus_tree(sys.argv[1:])


def prune_tree(argv=None):
    Phykit.prune_tree(sys.argv[1:])


def spr(argv=None):
    Phykit.subtree_prune_regraft(sys.argv[1:])


def subtree_prune_regraft(argv=None):
    Phykit.subtree_prune_regraft(sys.argv[1:])


def transfer_annotations(argv=None):
    Phykit.transfer_annotations(sys.argv[1:])


def relative_rate_test(argv=None):
    Phykit.relative_rate_test(sys.argv[1:])


def threshold_model(argv=None):
    Phykit.threshold_model(sys.argv[1:])


def rename_tree_tips(argv=None):
    Phykit.rename_tree_tips(sys.argv[1:])


def kf_distance(argv=None):
    Phykit.kf_distance(sys.argv[1:])


def rf_distance(argv=None):
    Phykit.rf_distance(sys.argv[1:])


def root_tree(argv=None):
    Phykit.root_tree(sys.argv[1:])


def spurious_sequence(argv=None):
    Phykit.spurious_sequence(sys.argv[1:])


def terminal_branch_stats(argv=None):
    Phykit.terminal_branch_stats(sys.argv[1:])


def tip_labels(argv=None):
    Phykit.tip_labels(sys.argv[1:])


def tip_to_tip_distance(argv=None):
    Phykit.tip_to_tip_distance(sys.argv[1:])


def tip_to_tip_node_distance(argv=None):
    Phykit.tip_to_tip_node_distance(sys.argv[1:])


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


def evo_tempo_map(argv=None):
    Phykit.evo_tempo_map(sys.argv[1:])


def discordance_asymmetry(argv=None):
    Phykit.discordance_asymmetry(sys.argv[1:])


def hybridization(argv=None):
    Phykit.hybridization(sys.argv[1:])


def spectral_discordance(argv=None):
    Phykit.spectral_discordance(sys.argv[1:])


def tree_space(argv=None):
    Phykit.tree_space(sys.argv[1:])


def phylo_impute(argv=None):
    Phykit.phylo_impute(sys.argv[1:])


def trait_rate_map(argv=None):
    Phykit.trait_rate_map(sys.argv[1:])


def taxon_groups(argv=None):
    Phykit.taxon_groups(sys.argv[1:])


def occupancy_filter(argv=None):
    Phykit.occupancy_filter(sys.argv[1:])
