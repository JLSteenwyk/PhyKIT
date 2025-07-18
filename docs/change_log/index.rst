.. _change_log:


Change log
==========

^^^^^

Major changes to PhyKIT are summarized here.

**2.0.2**:
Fixed bug in dna threading associated with how gaps were introduced in codons.

**2.0.1**:
Added arguments to exclude sites with gaps in the pairwise identities and saturation functions.

**2.0.0**:
Codebase overhaul to make PhyKIT more mem efficient and faster. For example, using list comprehension when appropriate.

**1.21.0**:
The partition file outputted from the create_concat function has been updated to the following format:
- column 1: alignment name
- column 2: # of taxa present
- column 3: # of taxa missing
- column 4: fraction of occupancy
- column 5: names of missing taxa (; separated)

**1.20.0**:
Fixed bug for thread_dna function when using a ClipKIT log file. Input protein alignment must be the untrimmed alignment.

**1.19.9**:
Saturation function now also reports the absolute value of 1-saturation. Lower values are indicative of less saturation.

**1.19.4**:
Saturation function forces y-intercept to be zero when calculating slope

**1.19.3**:
Saturation function now uses uncorrected distances instead of pairwise identities

**1.19.2**:
Verbose pairwise identity reporting separates pairwise identities by tabs and not a dash

**1.19.0**:
Added function to test for site-wise compositional biases in an alignment. See function compositional_bias_per_site.

**1.18.0**:
Added function to estimate site-wise evolutionary rate in an alignment. See function evolutionary_rate_per_site.

**1.15.0**:
Added function to recode alignments based on 8 different recoding schemes (7 for amino acids;
1 for nucleotides). See function recode.

**1.14.0**:
Added an optional argument to the thread_dna function. Now, PhyKIT can thread nucleotide
sequences onto a trimmed amino acid alignment. To do so, point PhyKIT to the ClipKIT outputted log
file using the -c argument. The ClipKIT log file can be generated when trimming an alignment with 
ClipKIT by adding the -l argument (see here for more details: https://jlsteenwyk.com/ClipKIT/).

**1.12.6**: relative composition variability is now adapted for calculating compositional biases in
individual taxa. The new function in rcvt (relative composition variability, taxon).

**1.12.4**: calculations of pairwise identity in alignment now supports excluding pairwise 
combinations with gaps.

**1.12.3**: hidden paralogy check now simply looks for monophyly or lack thereof for a set of taxa. Hidden paralogy
check still reports insufficient taxon representation.

**1.12.2**: removed root.txt file from DVMC function. User's are now recommended to trim outgroup taxa beforehand

**1.11.3**: Added an optional argument to the prune_tree function wherein instead of pruning tips
specified in the input file, those tips will be kept.

**1.11.1**: Modified sum of pairs score to divide the correct number
of pairs by the number of pairs in the reference alignment rather
than the query alignment alignment

**1.11.0**: Added terminal_branch_stats (alias: tbs) function to examine terminal branch lengths

**1.10.1**: Modified column score and sum of pairs score to divide the correct number
of columns or pairs by the number of columns or pairs in the query alignment rather
than the reference alignment

**1.10.0**: Added tip_to_tip_node_distance (alias: t2t_node_dist; t2t_nd) function to calculate
the phylogenetic distance between two leaves in a phylogeny. Distance is measured in nodes between
two leaves

**1.9.0**: Added monophyly_check (alias: is_monophyletic) function to examine monophyly 
among a specified set of taxa

**1.8.0**: Added hidden_paralogy_check (alias: clan_check) function to examine phylogenetic
tree for issues of hidden paralogy

**1.7.0**: Added nearest_neighbor_interchange (alias: nni) function to generate all NNI moves
for a binary rooted phylogeny

**1.6.0**: Added tip_to_tip_distance (alias: t2t_dist; t2t) function to calculate phylogenetic distance
between two leaves in a phylogeny

**1.5.0**: Added root_tree (alias: root; rt) function to root a phylogenetic tree

**1.4.0**: PhyKIT is now Python version 3.9 and BioPython 1.79 compatible

**1.3.0**: Added function that estimates the evolutionary rate of a gene using tree-based
properties. Function name is 'evolutionary_rate' or 'evo_rate' 

**1.2.2**: added function to get the subtree of the last common ancestor among a set of taxa

**1.2.0**: added command line interfaces for all functions so that each command 
can easily be executed. For example, 'phykit aln_len -h' can now be
called using 'pk_aln_len -h'

**1.1.0**: added faidx (alias: get_entry; ge) function to extract fasta entries from a
multi-fasta file

**1.0.3**: added rooting procedure before calculating RF to handle comparing unrooted
and rooted trees

**1.0.2**: function that calculates Robinson Foulds distance (robinson_foulds_distance;
rf_distance; rf_dist; rf) now can take trees that differ in topology. PhyKIT
will first determine shared tips between the two trees and prune both trees
to a common set of tips. Next, PhyKIT will calculate the Robinson Foulds 
distance.

**0.1.3**: Added function (column_score; cs) to calculate the quality of
an alignment given an input query alignment and a reference
alignment to compare it to

**0.1.2**: Added function (sum_of_pairs_score; sops; sop) to calculate
the quality of an alignment given an input query alignment
and a reference alignment to compare it to

**0.0.9**: PhyKIT now handles error stemming from piping output
