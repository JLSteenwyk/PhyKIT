PhyKIT is a collection of functions that help users explore alignments and
trees. For example, PhyKIT can calculate statistics and parameters of alignments
and trees that are associated with phylogenetic signal, assist in creating 
concatenated or codon-based alignments for tree inference, among other functions.
  


Functions that calculate measures of phylogenetic signal include: 
- alignment_length
- alignment_length_no_gaps
- bipartition_support_stats
- gc_content
- internal_branch_stats
- long_branch_score
- parsimony_informative_sites
- relative_composition_variability
- saturation 
- treeness
- treeness_over_rcv
- variable_sites

Functions that perform calculations on multiple trees include:
- covarying_evolutionary_rates
- polytomy_test
- robinson_foulds_distance

Functions that help edit/manipulative or view alignments or trees include:
- rename_fasta_entries
- branch_length_multiplier
- collapse_branches
- internode_labeler
- print_tree
- prune_tree
- rename_tree_tips

Functions that help build alignments for tree inference include:
- create_concatenation_matrix
- thread_dna

Other functions that help understand alignments or trees include:
- pairwise_identity
- patristic_distances
- spurious_sequence


=============================================
Consider writing a function that conducts objective desirability-based integration of multiple statistics