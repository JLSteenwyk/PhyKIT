.. _usage:

Usage
=====

PhyKIT provides 100+ functions for processing and analyzing multiple sequence
alignments and phylogenies. Functions span alignment quality assessment,
tree manipulation, phylogenetic comparative methods, trait evolution modeling,
introgression detection, and more.

Some help messages indicate that summary statistics are reported (e.g.,
bipartition_support_stats). Summary statistics include mean, median, 25th percentile,
75th percentile, minimum, maximum, standard deviation, and variance. These functions
typically have a verbose option that allows users to get the underlying data
used to calculate summary statistics.

|

Quick start
-----------

Here is a typical workflow showing a few common PhyKIT operations:

.. code-block:: shell

   # Check alignment quality
   phykit pis alignment.fa                  # count parsimony informative sites
   phykit aot alignment.fa --json           # flag outlier taxa

   # Summarize tree properties
   phykit treeness species.tre              # treeness (internal/total branch length)
   phykit dvmc species.tre                  # degree of violation of a molecular clock

   # Phylogenetic comparative methods
   phykit pgls -t species.tre -d traits.tsv \
       --response brain_size --predictor body_mass   # PGLS regression
   phykit panova -t species.tre \
       --traits traits.tsv --pairwise                # phylogenetic ANOVA

   # Visualize gene tree concordance
   phykit qpie -t species.tre -g gene_trees.nwk \
       -o concordance.png --branch-labels            # quartet pie chart

|

General usage
-------------

Calling functions
#################

.. code-block:: shell

   phykit <command> [optional command arguments]

Command specific help messages can be viewed by adding a 
-h/--help argument after the command. For example, to see the help message
for the command 'treeness', execute:

.. code-block:: shell

   phykit treeness -h 
   # or
   phykit treeness --help

|

Function aliases
################

Each function comes with aliases to save the user some
key strokes. For example, to get the help message for the 'treeness'
function, you can type:

.. code-block:: shell

   phykit tness -h

|

Command line interfaces
#######################

As of version 1.2.0, all functions (including aliases) can be executed using
a command line interface that starts with *pk_*. For example, instead of typing
the previous command to get the help message of the treeness function, you can type:

.. code-block:: shell

   pk_treeness -h
   # or
   pk_tness -h


All possible function names are specified at the top of each function section. 

|

Command reference
-----------------

The complete reference is organized by task on the :doc:`command reference page </reference/index>`.

Legacy command anchors
######################

These links preserve the command fragments used by earlier versions of the documentation.

.. raw:: html

   <span id="cmd-alignment_entropy"></span>

- :doc:`Alignment entropy </reference/commands/alignment_entropy>`

.. raw:: html

   <span id="cmd-alignment_length"></span>

- :doc:`Alignment length </reference/commands/alignment_length>`

.. raw:: html

   <span id="cmd-alignment_length_no_gaps"></span>

- :doc:`Alignment length no gaps </reference/commands/alignment_length_no_gaps>`

.. raw:: html

   <span id="cmd-alignment_outlier_regions"></span>

- :doc:`Alignment outlier regions </reference/commands/alignment_outlier_regions>`

.. raw:: html

   <span id="cmd-alignment_outlier_taxa"></span>

- :doc:`Alignment outlier taxa </reference/commands/alignment_outlier_taxa>`

.. raw:: html

   <span id="cmd-alignment_recoding"></span>

- :doc:`Alignment recoding </reference/commands/alignment_recoding>`

.. raw:: html

   <span id="cmd-alignment_subsample"></span>

- :doc:`Alignment subsampling </reference/commands/alignment_subsample>`

.. raw:: html

   <span id="cmd-ancestral_state_reconstruction"></span>

- :doc:`Ancestral state reconstruction </reference/commands/ancestral_state_reconstruction>`

.. raw:: html

   <span id="cmd-bipartition_support_stats"></span>

- :doc:`Bipartition support statistics </reference/commands/bipartition_support_stats>`

.. raw:: html

   <span id="cmd-branch_length_multiplier"></span>

- :doc:`Branch length multiplier </reference/commands/branch_length_multiplier>`

.. raw:: html

   <span id="cmd-character_map"></span>

- :doc:`Character map (synapomorphy/homoplasy mapping) </reference/commands/character_map>`

.. raw:: html

   <span id="cmd-chronogram"></span>

- :doc:`Chronogram </reference/commands/chronogram>`

.. raw:: html

   <span id="cmd-codon_dnds"></span>

- :doc:`Codon dN/dS </reference/commands/codon_dnds>`

.. raw:: html

   <span id="cmd-collapse_branches"></span>

- :doc:`Collapse bipartitions </reference/commands/collapse_branches>`

.. raw:: html

   <span id="cmd-column_score"></span>

- :doc:`Column score </reference/commands/column_score>`

.. raw:: html

   <span id="cmd-composition_per_taxon"></span>

- :doc:`Composition per taxon </reference/commands/composition_per_taxon>`

.. raw:: html

   <span id="cmd-compositional_bias_per_site"></span>

- :doc:`Compositional bias per site </reference/commands/compositional_bias_per_site>`

.. raw:: html

   <span id="cmd-concordance_asr"></span>

- :doc:`Concordance-aware ancestral state reconstruction </reference/commands/concordance_asr>`

.. raw:: html

   <span id="cmd-consensus_network"></span>

- :doc:`Consensus network </reference/commands/consensus_network>`

.. raw:: html

   <span id="cmd-consensus_tree"></span>

- :doc:`Consensus tree </reference/commands/consensus_tree>`

.. raw:: html

   <span id="cmd-cont_map"></span>

- :doc:`Continuous trait mapping (contMap) </reference/commands/cont_map>`

.. raw:: html

   <span id="cmd-cophylo"></span>

- :doc:`Cophylogenetic plot (tanglegram) </reference/commands/cophylo>`

.. raw:: html

   <span id="cmd-covarying_evolutionary_rates"></span>

- :doc:`Covarying evolutionary rates </reference/commands/covarying_evolutionary_rates>`

.. raw:: html

   <span id="cmd-create_concatenation_matrix"></span>

- :doc:`Create concatenation matrix </reference/commands/create_concatenation_matrix>`

.. raw:: html

   <span id="cmd-degree_of_violation_of_a_molecular_clock"></span>

- :doc:`Degree of violation of the molecular clock </reference/commands/degree_of_violation_of_a_molecular_clock>`

.. raw:: html

   <span id="cmd-density_map"></span>

- :doc:`Density map </reference/commands/density_map>`

.. raw:: html

   <span id="cmd-dfoil"></span>

- :doc:`DFOIL test (Pease & Hahn 2015) </reference/commands/dfoil>`

.. raw:: html

   <span id="cmd-discordance_asymmetry"></span>

- :doc:`Discordance asymmetry </reference/commands/discordance_asymmetry>`

.. raw:: html

   <span id="cmd-dstatistic"></span>

- :doc:`D-statistic (ABBA-BABA test) </reference/commands/dstatistic>`

.. raw:: html

   <span id="cmd-dtt"></span>

- :doc:`Disparity through time (DTT) </reference/commands/dtt>`

.. raw:: html

   <span id="cmd-evo_tempo_map"></span>

- :doc:`Evolutionary tempo mapping </reference/commands/evo_tempo_map>`

.. raw:: html

   <span id="cmd-evolutionary_rate"></span>

- :doc:`Evolutionary rate </reference/commands/evolutionary_rate>`

.. raw:: html

   <span id="cmd-evolutionary_rate_per_site"></span>

- :doc:`Evolutionary Rate per Site </reference/commands/evolutionary_rate_per_site>`

.. raw:: html

   <span id="cmd-faidx"></span>

- :doc:`Faidx </reference/commands/faidx>`

.. raw:: html

   <span id="cmd-faiths_pd"></span>

- :doc:`Faith's phylogenetic diversity </reference/commands/faiths_pd>`

.. raw:: html

   <span id="cmd-fit_continuous"></span>

- :doc:`Continuous trait evolution model comparison (fitContinuous) </reference/commands/fit_continuous>`

.. raw:: html

   <span id="cmd-fit_discrete"></span>

- :doc:`Discrete trait evolution model comparison (fitDiscrete) </reference/commands/fit_discrete>`

.. raw:: html

   <span id="cmd-gc_content"></span>

- :doc:`Guanine-cytosine (GC) content </reference/commands/gc_content>`

.. raw:: html

   <span id="cmd-hidden_paralogy_check"></span>

- :doc:`Hidden paralogy check </reference/commands/hidden_paralogy_check>`

.. raw:: html

   <span id="cmd-hybridization"></span>

- :doc:`Hybridization analysis </reference/commands/hybridization>`

.. raw:: html

   <span id="cmd-identity_matrix"></span>

- :doc:`Identity matrix </reference/commands/identity_matrix>`

.. raw:: html

   <span id="cmd-independent_contrasts"></span>

- :doc:`Independent contrasts (PIC) </reference/commands/independent_contrasts>`

.. raw:: html

   <span id="cmd-internal_branch_stats"></span>

- :doc:`Internal branch statistics </reference/commands/internal_branch_stats>`

.. raw:: html

   <span id="cmd-internode_labeler"></span>

- :doc:`Internode labeler </reference/commands/internode_labeler>`

.. raw:: html

   <span id="cmd-kf_distance"></span>

- :doc:`Kuhner-Felsenstein distance </reference/commands/kf_distance>`

.. raw:: html

   <span id="cmd-last_common_ancestor_subtree"></span>

- :doc:`Last common ancestor subtree </reference/commands/last_common_ancestor_subtree>`

.. raw:: html

   <span id="cmd-long_branch_score"></span>

- :doc:`Long branch score </reference/commands/long_branch_score>`

.. raw:: html

   <span id="cmd-ltt"></span>

- :doc:`Lineage-through-time plot and gamma statistic </reference/commands/ltt>`

.. raw:: html

   <span id="cmd-mask_alignment"></span>

- :doc:`Mask alignment </reference/commands/mask_alignment>`

.. raw:: html

   <span id="cmd-monophyly_check"></span>

- :doc:`Monophyly check </reference/commands/monophyly_check>`

.. raw:: html

   <span id="cmd-nearest_neighbor_interchange"></span>

- :doc:`Nearest neighbor interchange </reference/commands/nearest_neighbor_interchange>`

.. raw:: html

   <span id="cmd-neighbor_net"></span>

- :doc:`NeighborNet </reference/commands/neighbor_net>`

.. raw:: html

   <span id="cmd-network_signal"></span>

- :doc:`Network signal </reference/commands/network_signal>`

.. raw:: html

   <span id="cmd-occupancy_filter"></span>

- :doc:`Occupancy filter </reference/commands/occupancy_filter>`

.. raw:: html

   <span id="cmd-occupancy_per_taxon"></span>

- :doc:`Occupancy per taxon </reference/commands/occupancy_per_taxon>`

.. raw:: html

   <span id="cmd-ou_shift_detection"></span>

- :doc:`OU shift detection (l1ou) </reference/commands/ou_shift_detection>`

.. raw:: html

   <span id="cmd-ouwie"></span>

- :doc:`Multi-regime OU models (OUwie) </reference/commands/ouwie>`

.. raw:: html

   <span id="cmd-pairwise_identity"></span>

- :doc:`Pairwise identity </reference/commands/pairwise_identity>`

.. raw:: html

   <span id="cmd-parsimony_informative_sites"></span>

- :doc:`Parsimony informative sites </reference/commands/parsimony_informative_sites>`

.. raw:: html

   <span id="cmd-parsimony_score"></span>

- :doc:`Parsimony score </reference/commands/parsimony_score>`

.. raw:: html

   <span id="cmd-patristic_distances"></span>

- :doc:`Patristic distances </reference/commands/patristic_distances>`

.. raw:: html

   <span id="cmd-phenogram"></span>

- :doc:`Phenogram (traitgram) </reference/commands/phenogram>`

.. raw:: html

   <span id="cmd-phylo_anova"></span>

- :doc:`Phylogenetic ANOVA / MANOVA </reference/commands/phylo_anova>`

.. raw:: html

   <span id="cmd-phylo_gwas"></span>

- :doc:`Phylo GWAS </reference/commands/phylo_gwas>`

.. raw:: html

   <span id="cmd-phylo_heatmap"></span>

- :doc:`Phylogenetic heatmap </reference/commands/phylo_heatmap>`

.. raw:: html

   <span id="cmd-phylo_impute"></span>

- :doc:`Phylogenetic imputation </reference/commands/phylo_impute>`

.. raw:: html

   <span id="cmd-phylo_logistic"></span>

- :doc:`Phylogenetic Logistic Regression </reference/commands/phylo_logistic>`

.. raw:: html

   <span id="cmd-phylo_path"></span>

- :doc:`Phylogenetic path analysis </reference/commands/phylo_path>`

.. raw:: html

   <span id="cmd-phylogenetic_glm"></span>

- :doc:`Phylogenetic GLM </reference/commands/phylogenetic_glm>`

.. raw:: html

   <span id="cmd-phylogenetic_ordination"></span>

- :doc:`Phylogenetic Ordination </reference/commands/phylogenetic_ordination>`

.. raw:: html

   <span id="cmd-phylogenetic_regression"></span>

- :doc:`Phylogenetic regression (PGLS) </reference/commands/phylogenetic_regression>`

.. raw:: html

   <span id="cmd-phylogenetic_signal"></span>

- :doc:`Phylogenetic signal </reference/commands/phylogenetic_signal>`

.. raw:: html

   <span id="cmd-phylomorphospace"></span>

- :doc:`Phylomorphospace </reference/commands/phylomorphospace>`

.. raw:: html

   <span id="cmd-plot_alignment_qc"></span>

- :doc:`Plot alignment QC </reference/commands/plot_alignment_qc>`

.. raw:: html

   <span id="cmd-polytomy_test"></span>

- :doc:`Polytomy testing </reference/commands/polytomy_test>`

.. raw:: html

   <span id="cmd-print_tree"></span>

- :doc:`Print tree </reference/commands/print_tree>`

.. raw:: html

   <span id="cmd-prune_tree"></span>

- :doc:`Prune tree </reference/commands/prune_tree>`

.. raw:: html

   <span id="cmd-quartet_network"></span>

- :doc:`Quartet network </reference/commands/quartet_network>`

.. raw:: html

   <span id="cmd-quartet_pie"></span>

- :doc:`Quartet pie chart </reference/commands/quartet_pie>`

.. raw:: html

   <span id="cmd-rate_heterogeneity"></span>

- :doc:`Rate heterogeneity test (multi-rate Brownian motion) </reference/commands/rate_heterogeneity>`

.. raw:: html

   <span id="cmd-relative_composition_variability"></span>

- :doc:`Relative composition variability </reference/commands/relative_composition_variability>`

.. raw:: html

   <span id="cmd-relative_composition_variability_taxon"></span>

- :doc:`Relative composition variability, taxon </reference/commands/relative_composition_variability_taxon>`

.. raw:: html

   <span id="cmd-relative_rate_test"></span>

- :doc:`Relative rate test </reference/commands/relative_rate_test>`

.. raw:: html

   <span id="cmd-rename_fasta_entries"></span>

- :doc:`Rename FASTA entries </reference/commands/rename_fasta_entries>`

.. raw:: html

   <span id="cmd-rename_tree_tips"></span>

- :doc:`Rename tree tips </reference/commands/rename_tree_tips>`

.. raw:: html

   <span id="cmd-robinson_foulds_distance"></span>

- :doc:`Robinson-Foulds distance </reference/commands/robinson_foulds_distance>`

.. raw:: html

   <span id="cmd-root_tree"></span>

- :doc:`Root tree </reference/commands/root_tree>`

.. raw:: html

   <span id="cmd-saturation"></span>

- :doc:`Saturation </reference/commands/saturation>`

.. raw:: html

   <span id="cmd-simmap_summary"></span>

- :doc:`SIMMAP summary </reference/commands/simmap_summary>`

.. raw:: html

   <span id="cmd-spectral_discordance"></span>

- :doc:`Spectral discordance decomposition </reference/commands/spectral_discordance>`

.. raw:: html

   <span id="cmd-spurious_sequence"></span>

- :doc:`Spurious homolog identification </reference/commands/spurious_sequence>`

.. raw:: html

   <span id="cmd-stochastic_character_map"></span>

- :doc:`Stochastic character mapping (SIMMAP) </reference/commands/stochastic_character_map>`

.. raw:: html

   <span id="cmd-subtree_prune_regraft"></span>

- :doc:`Subtree pruning and regrafting </reference/commands/subtree_prune_regraft>`

.. raw:: html

   <span id="cmd-sum_of_pairs_score"></span>

- :doc:`Sum-of-pairs score </reference/commands/sum_of_pairs_score>`

.. raw:: html

   <span id="cmd-taxon_groups"></span>

- :doc:`Taxon groups </reference/commands/taxon_groups>`

.. raw:: html

   <span id="cmd-terminal_branch_stats"></span>

- :doc:`Terminal branch statistics </reference/commands/terminal_branch_stats>`

.. raw:: html

   <span id="cmd-thread_dna"></span>

- :doc:`Protein-to-nucleotide alignment </reference/commands/thread_dna>`

.. raw:: html

   <span id="cmd-threshold_model"></span>

- :doc:`Threshold model </reference/commands/threshold_model>`

.. raw:: html

   <span id="cmd-tip_labels"></span>

- :doc:`Tip labels </reference/commands/tip_labels>`

.. raw:: html

   <span id="cmd-tip_to_tip_distance"></span>

- :doc:`Tip-to-tip distance </reference/commands/tip_to_tip_distance>`

.. raw:: html

   <span id="cmd-tip_to_tip_node_distance"></span>

- :doc:`Tip-to-tip node distance </reference/commands/tip_to_tip_node_distance>`

.. raw:: html

   <span id="cmd-total_tree_length"></span>

- :doc:`Total tree length </reference/commands/total_tree_length>`

.. raw:: html

   <span id="cmd-trait_correlation"></span>

- :doc:`Trait correlation </reference/commands/trait_correlation>`

.. raw:: html

   <span id="cmd-trait_rate_map"></span>

- :doc:`Trait rate map </reference/commands/trait_rate_map>`

.. raw:: html

   <span id="cmd-transfer_annotations"></span>

- :doc:`Transfer annotations </reference/commands/transfer_annotations>`

.. raw:: html

   <span id="cmd-tree_space"></span>

- :doc:`Tree space visualization </reference/commands/tree_space>`

.. raw:: html

   <span id="cmd-treeness"></span>

- :doc:`Treeness </reference/commands/treeness>`

.. raw:: html

   <span id="cmd-treeness_over_rcv"></span>

- :doc:`Treeness over RCV </reference/commands/treeness_over_rcv>`

.. raw:: html

   <span id="cmd-variable_sites"></span>

- :doc:`Variable sites </reference/commands/variable_sites>`

.. raw:: html

   <span id="cmd-version"></span>

- :doc:`Version </reference/commands/version>`
