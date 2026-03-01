.. _usage:

Usage
=====

PhyKIT helps process and analyze multiple sequence alignments and phylogenies.

Generally, all functions are designed to help understand the contents of alignments
(e.g., gc content or the number of parsimony informative sites) and the shape
of trees (e.g., treeness, degree of violation of a molecular clock).

Some help messages indicate that summary statistics are reported (e.g., 
bipartition_support_stats). Summary statistics include mean, median, 25th percentile,
75th percentile, minimum, maximum, standard deviation, and variance. These functions
typically have a verbose option that allows users to get the underlying data
used to calculate summary statistics. 

|

General usage
-------------

Calling functions
#################

.. code-block:: shell

   phykit <command> [optional command arguments]

Command specific help messages can be viewed by adding a 
-h/\\-\\-help argument after the command. For example, to see the help message
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

Alignment-based functions
-------------------------

Alignment length
################

.. image:: ../_static/docs_img/aln_len.png 
   :align: center
   :width: 75%

Function names: alignment_length; aln_len; al |br|
Command line interface: pk_alignment_length; pk_aln_len; pk_al

Length of an input alignment is calculated using this function.

Longer alignments are associated with strong phylogenetic signal.
   
Association between alignment length and phylogenetic signal
was determined by Shen et al., Genome Biology and Evolution (2016),
doi: 10.1093/gbe/evw179.

.. code-block:: shell

		phykit aln_len <alignment> [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--json*: optional argument to print results as JSON

|

Alignment length no gaps
########################

.. image:: ../_static/docs_img/aln_len_no_gaps.png 
   :align: center
   :width: 75%

Function names: alignment_length_no_gaps; aln_len_no_gaps; alng |br|
Command line interface: pk_alignment_length_no_gaps; pk_aln_len_no_gaps; pk_alng

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

.. code-block:: shell

		phykit aln_len_no_gaps <alignment> [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--json*: optional argument to print results as JSON

|

Alignment entropy
#################
Function names: alignment_entropy; aln_entropy; entropy |br|
Command line interface: pk_alignment_entropy; pk_aln_entropy; pk_entropy

Calculate alignment entropy.

Site-wise entropy is calculated using Shannon entropy. By default,
PhyKIT reports the mean entropy across all sites in the alignment.
With the `-v/--verbose` option, PhyKIT reports entropy for each site.

.. code-block:: shell

	phykit alignment_entropy <alignment> [-v/--verbose] [--plot] [--plot-output <path>] [--json]

Example output (default):

.. code-block:: shell

	0.657

Example output (`-v`):

.. code-block:: shell

	1	0.0
	2	1.0
	3	0.971

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*-v/\\-\\-verbose*: optional argument to print entropy for each site |br|
*--plot*: save a per-site alignment entropy plot |br|
*--plot-output*: output path for plot (default: ``alignment_entropy_plot.png``) |br|
*--json*: optional argument to print results as JSON

|

Alignment recoding 
##################
Function names: alignment_recoding; aln_recoding; recode |br|
Command line interface: pk_alignment_recoding; pk_aln_recoding; pk_recode

Recode alignments using reduced character states.

Alignments can be recoded using established or
custom recoding schemes. Recoding schemes are
specified using the -c/--code argument. Custom
recoding schemes can be used and should be formatted
as a two column file wherein the first column is the
recoded character and the second column is the character
in the alignment.

.. code-block:: shell

	phykit alignment_recoding <fasta> [-c/--code <code>] [--json]

Codes for which recoding scheme to use: |br|
**RY-nucleotide** |br|
R = purines (i.e., A and G) |br|
Y = pyrimidines (i.e., T and C) |br|

**SandR-6** |br|
0 = A, P, S, and T |br|
1 = D, E, N, and G |br|
2 = Q, K, and R |br|
3 = M, I, V, and L |br|
4 = W and C |br|
5 = F, Y, and H |br|

**KGB-6** |br|
0 = A, G, P, and S |br|
1 = D, E, N, Q, H, K, R, and T |br|
2 = M, I, and L |br|
3 = W |br|
4 = F and Y |br|
5 = C and V |br|

**Dayhoff-6** |br|
0 = A, G, P, S, and T |br|
1 = D, E, N, and Q |br|
2 = H, K, and R |br|
3 = I, L, M, and V |br|
4 = F, W, and Y |br|
5 = C |br|

**Dayhoff-9** |br|
0 = D, E, H, N, and Q |br|
1 = I, L, M, and V |br|
2 = F and Y |br|
3 = A, S, and T |br|
4 = K and R |br|
5 = G |br|
6 = P |br|
7 = C |br|
8 = W |br|

**Dayhoff-12** |br|
0 = D, E, and Q |br|
1 = M, L, I, and V |br|
2 = F and Y |br|
3 = K, H, and R |br|
4 = G |br|
5 = A |br|
6 = P |br|
7 = S |br|
8 = T |br|
9 = N |br|
A = W |br|
B = C |br|

**Dayhoff-15** |br|
0 = D, E, and Q |br|
1 = M and L |br|
2 = I and V |br|
3 = F and Y |br|
4 = G |br|
5 = A |br|
6 = P |br|
7 = S |br|
8 = T |br|
9 = N |br|
A = K |br|
B = H |br|
C = R |br|
D = W |br|
E = C |br|

**Dayhoff-18** |br|
0 = F and Y |br|
1 = M and L |br|
2 = I |br|
3 = V |br|
4 = G |br|
5 = A |br|
6 = P |br|
7 = S |br|
8 = T |br|
9 = D |br|
A = E |br|
B = Q |br|
C = N |br|
D = K |br|
E = H |br|
F = R |br|
G = W |br|
H = C |br|

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*-c/\-\-code*: argument to specify the recoding scheme to use |br|
*--json*: optional argument to print results as JSON

|

Column score
############

.. image:: ../_static/docs_img/column_score.png 
   :align: center
   :width: 75%

Function names: column_score; cs |br|
Command line interface: pk_column_score; pk_cs

Calculates column score.

Column is an accuracy metric for a multiple alignment relative
to a reference alignment. It is calculated by summing the correctly
aligned columns over all columns in an alignment. Thus, values range
from 0 to 1 and higher values indicate more accurate alignments.

Column score is calculated following Thompson et al., Nucleic
Acids Research (1999), doi: 10.1093/nar/27.13.2682.

.. code-block:: shell

	phykit column_score <alignment> --reference <reference_alignment> [--json]

Options: |br|
*<alignment>*: first argument after function name should be a query
fasta alignment file to be scored for accuracy |br|
*-r/\\-\\-reference*: reference alignment to compare the query alignment
to |br|
*--json*: optional argument to print results as JSON

|

Compositional bias per site
###########################

Function names: compositional_bias_per_site; comp_bias_per_site; cbps |br|
Command line interface: pk_compositional_bias_per_site; pk_comp_bias_per_site; pk_cbps

Calculates compositional bias per site in an alignment.

Site-wise chi-squared tests are conducted in an alignment to
detect compositional biases. PhyKIT outputs four columns: |br|
col 1: index in alignment |br|
col 2: chi-squared statistic (higher values indicate greater bias) |br|
col 3: multi-test corrected p-value (Benjamini-Hochberg false discovery rate procedure) |br|
col 4: uncorrected p-value

.. code-block:: shell

	phykit comp_bias_per_site <alignment> [--plot] [--plot-output <path>] [--json]

Options: |br|
*<alignment>*: first argument after function name should be a query
fasta alignment to calculate the site-wise compositional bias of |br|
*--plot*: save a Manhattan-style plot of site-wise compositional bias |br|
*--plot-output*: output path for plot (default: ``compositional_bias_per_site_plot.png``) |br|
*--json*: optional argument to print results as JSON

|

Composition per taxon
#####################
Function names: composition_per_taxon; comp_taxon; comp_tax |br|
Command line interface: pk_composition_per_taxon; pk_comp_taxon; pk_comp_tax

Calculate sequence composition per taxon in an alignment.

Composition is reported as semicolon-separated `symbol:frequency` pairs
for each taxon. Frequencies are calculated from valid (non-gap/non-ambiguous)
characters. Symbol order is alphabetical.

.. code-block:: shell

	phykit composition_per_taxon <alignment> [--json]

Example output:

.. code-block:: shell

	1	A:0.4;C:0.0;G:0.2;T:0.4
	2	A:0.5;C:0.0;G:0.25;T:0.25

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--json*: optional argument to print results as JSON

|

Create concatenation matrix
###########################

.. image:: ../_static/docs_img/create_concat_matrix.png 
   :align: center
   :width: 75%

Function names: create_concatenation_matrix, create_concat, cc |br|
Command line interface: pk_create_concatenation_matrix, pk_create_concat, pk_cc

Create a concatenated alignment file. This function is 
used to help in the construction of multi-locus data
matrices.

PhyKIT will output three files:
1) A fasta file with '.fa' appended to the prefix specified with the -p/\\-\\-prefix parameter.
2) A partition file ready for input into RAxML or IQ-tree.
3) An occupancy file that summarizes the taxon occupancy per sequence.

.. code-block:: shell

	phykit create_concat -a <file> -p <string> [--threshold <float>] [--plot-occupancy] [--plot-output <path>] [--json]

Options: |br|
*-a/\\-\\-alignment*: alignment list file. File should contain a single column list of alignment
sequence files to concatenate into a single matrix. Provide path to files relative to
working directory or provide absolute path. |br|
*-p/\\-\\-prefix*: prefix of output files |br|
*--threshold*: minimum fraction of informative (non-gap, non-ambiguous) sites across the
concatenated alignment for a taxon to be included. Taxa whose effective occupancy falls
below this value are excluded from the output. Set to 0 to disable filtering
(default: 0). |br|
*--plot-occupancy*: optional argument to output an occupancy map figure where
x-axis shows concatenated positions with gene boundaries and y-axis shows taxa
sorted by total occupancy. Colors denote represented characters, gap/ambiguous
characters in present genes, and fully absent gene blocks. |br|
*--plot-output*: optional custom output path for occupancy map figure
(default: ``<prefix>.occupancy.png``). |br|
*--json*: optional argument to print summary metadata as JSON

|

Evolutionary Rate per Site
##########################

Function names: evolutionary_rate_per_site; evo_rate_per_site; erps |br|
Command line interface: pk_evolutionary_rate_per_site; pk_evo_rate_per_site; pk_erps

Estimate evolutionary rate per site.

Evolutionary rate per site is one minus the sum of squared frequency of different
characters at a given site. Values may range from 0 (slow evolving; no diversity
at the given site) to 1 (fast evolving; all characters appear only once).

.. code-block:: shell

	phykit evo_rate_per_site <alignment> [--plot] [--plot-output <path>] [--json]

Options: |br|
*<alignment>*: first argument after function name should be a query
fasta alignment to calculate the site-wise evolutionary rate of |br|
*--plot*: save a per-site evolutionary-rate plot |br|
*--plot-output*: output path for plot (default: ``evolutionary_rate_per_site_plot.png``) |br|
*--json*: optional argument to print results as JSON

|

Faidx
#####

.. image:: ../_static/docs_img/faidx.png 
   :align: center
   :width: 75%

Function names: faidx; get_entry; ge |br|
Command line interface: pk_faidx; pk_get_entry; pk_ge

Extracts sequence entry from fasta file.

This function works similarly to the faidx function 
in samtools, but does not requiring an indexing step.

To obtain multiple entries, input multiple entries separated
by a comma (,). For example, if you want entries 
named "seq_0" and "seq_1", the string "seq_0,seq_1"
should be associated with the -e argument.

.. code-block:: shell

	phykit faidx <fasta> -e/--entry <fasta entry> [--json]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-e/\\-\\-entry*: entry name to be extracted from the inputted fasta file |br|
*--json*: optional argument to print results as JSON

|

Guanine-cytosine (GC) content
#############################

.. image:: ../_static/docs_img/gc_content.png 
   :align: center
   :width: 75%

Function names: gc_content; gc |br|
Command line interface: pk_gc_content; pk_gc

Calculate GC content of a fasta file.

GC content is negatively correlated with phylogenetic signal.

If there are multiple entries, use the -v/\\-\\-verbose option
to determine the GC content of each fasta entry separately.
Association between GC content and phylogenetic signal was
determined by Shen et al., Genome Biology and Evolution (2016), 
doi: 10.1093/gbe/evw179.

.. code-block:: shell

		phykit gc_content <fasta> [-v/--verbose] [--json]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-v/\\-\\-verbose*: optional argument to print the GC content of each fasta
entry |br|
*--json*: optional argument to print results as JSON

|

Mask alignment
##############
Function names: mask_alignment; mask_aln; mask |br|
Command line interface: pk_mask_alignment; pk_mask_aln; pk_mask

Mask alignment sites based on threshold criteria.

Sites are retained when they pass all active thresholds:
maximum gap fraction, minimum occupancy, and maximum site entropy.

.. code-block:: shell

	phykit mask_alignment <alignment> [-g/--max_gap <float>] [-o/--min_occupancy <float>] [-e/--max_entropy <float>] [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*-g/\\-\\-max_gap*: maximum allowed fraction of missing/invalid characters at a site (default: 1.0) |br|
*-o/\\-\\-min_occupancy*: minimum required occupancy at a site (default: 0.0) |br|
*-e/\\-\\-max_entropy*: maximum allowed site entropy (default: no filter) |br|
*--json*: optional argument to print results as JSON

|

Plot alignment QC
#################
Function names: plot_alignment_qc; plot_qc; paqc |br|
Command line interface: pk_plot_alignment_qc; pk_plot_qc; pk_paqc

Generate a multi-panel alignment quality-control plot.

The figure includes:
1) occupancy per taxon |br|
2) gap rate per taxon |br|
3) composition distance vs long-branch proxy scatter |br|
4) count of flagged outliers by feature

Outlier evaluation uses the same features as ``alignment_outlier_taxa``:
``gap_rate``, ``occupancy``, ``composition_distance``, ``long_branch_proxy``,
``rcvt``, and ``entropy_burden``.

.. code-block:: shell

	phykit plot_alignment_qc <alignment> [-o/--output <path>] [--width <float>] [--height <float>] [--dpi <int>] [--gap-z <float>] [--composition-z <float>] [--distance-z <float>] [--rcvt-z <float>] [--occupancy-z <float>] [--entropy-z <float>] [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*-o/\\-\\-output*: output image path (default: ``alignment_qc.png``) |br|
*--width*: figure width in inches (default: ``14.0``) |br|
*--height*: figure height in inches (default: ``10.0``) |br|
*--dpi*: output image DPI (default: ``300``) |br|
*--gap-z*: z-threshold for gap-rate outliers (default: ``3.0``) |br|
*--composition-z*: z-threshold for composition-distance outliers (default: ``3.0``) |br|
*--distance-z*: z-threshold for long-branch-proxy outliers (default: ``3.0``) |br|
*--rcvt-z*: z-threshold for RCVT outliers (default: ``3.0``) |br|
*--occupancy-z*: z-threshold for low-occupancy outliers (default: ``3.0``) |br|
*--entropy-z*: z-threshold for entropy-burden outliers (default: ``3.0``) |br|
*--json*: optional argument to print plot metadata and outlier summary as JSON

|

Occupancy per taxon
###################
Function names: occupancy_per_taxon; occupancy_taxon; occ_tax |br|
Command line interface: pk_occupancy_per_taxon; pk_occupancy_taxon; pk_occ_tax

Calculate occupancy per taxon in an alignment.

Occupancy is the fraction of valid (non-gap/non-ambiguous) characters
for each taxon.

.. code-block:: shell

	phykit occupancy_per_taxon <alignment> [--json]

Example output:

.. code-block:: shell

	1	0.8333
	2	0.6667

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--json*: optional argument to print results as JSON

|

Alignment outlier taxa
######################
Function names: alignment_outlier_taxa; outlier_taxa; aot |br|
Command line interface: pk_alignment_outlier_taxa; pk_outlier_taxa; pk_aot

Identify potential outlier taxa in an alignment and explicitly report why each taxon was flagged.

The following features are evaluated per taxon:
1) ``gap_rate``: fraction of gap/ambiguous symbols |br|
2) ``occupancy``: fraction of valid symbols |br|
3) ``composition_distance``: Euclidean distance from the median composition profile |br|
4) ``long_branch_proxy``: mean pairwise sequence distance to other taxa |br|
5) ``rcvt``: relative composition variability per taxon |br|
6) ``entropy_burden``: average site entropy over this taxon's valid positions

If a taxon exceeds one or more feature-specific thresholds, PhyKIT reports the exact
feature(s), observed value(s), threshold(s), and explanation(s) for the flag.

.. code-block:: shell

	phykit alignment_outlier_taxa <alignment> [--gap-z <float>] [--composition-z <float>] [--distance-z <float>] [--rcvt-z <float>] [--occupancy-z <float>] [--entropy-z <float>] [--json]

Example output:

.. code-block:: shell

	features_evaluated	gap_rate,occupancy,composition_distance,long_branch_proxy,rcvt,entropy_burden
	thresholds	gap_rate>0.0;occupancy<0.4;composition_distance>0.1;long_branch_proxy>0.4181;rcvt>0.0095;entropy_burden>0.35
	taxon_d	composition_distance=1.3454>0.1;long_branch_proxy=1.0>0.4181	Unusual sequence composition profile relative to other taxa. | High mean pairwise sequence distance to other taxa.
	taxon_e	gap_rate=0.6>0.0;occupancy=0.4<0.4	High fraction of gap/ambiguous symbols compared to other taxa. | Low fraction of valid symbols compared to other taxa.

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--gap-z*: z-threshold used for high-gap outlier detection (default: 3.0) |br|
*--composition-z*: z-threshold used for composition-distance outlier detection (default: 3.0) |br|
*--distance-z*: z-threshold used for long-branch-proxy outlier detection (default: 3.0) |br|
*--rcvt-z*: z-threshold used for RCVT outlier detection (default: 3.0) |br|
*--occupancy-z*: z-threshold used for low-occupancy outlier detection (default: 3.0) |br|
*--entropy-z*: z-threshold used for entropy-burden outlier detection (default: 3.0) |br|
*--json*: optional argument to print results as JSON

|

Pairwise identity
#################

.. image:: ../_static/docs_img/pairwise_identity.png 
   :align: center
   :width: 75%

Function names: pairwise_identity; pairwise_id, pi |br|
Command line interface: pk_pairwise_identity; pk_pairwise_id, pk_pi

Calculate the average pairwise identity among sequences.

Pairwise identities can be used as proxies for the evolutionary rate of sequences.

Pairwise identity is defined as the number of identical
columns (including gaps) between two aligned sequences divided
by the number of columns in the alignment. Summary statistics
are reported unless used with the verbose option in which
all pairwise identities will be reported.

An example of pairwise identities being used as a proxy
for evolutionary rate can be found here: Chen et al. 
Genome Biology and Evolution (2017), doi: 10.1093/gbe/evx147.

.. code-block:: shell

	phykit pairwise_identity <alignment> [-v/--verbose] [-e/--exclude_gaps] [--plot] [--plot-output <file>] [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*-v/\\-\\-verbose*: optional argument to print identity per pair|br|
*-e/\-\-exclude_gaps*: if a site has a gap, ignore it |br|
*--plot*: save a clustered pairwise-identity heatmap |br|
*--plot-output*: output path for heatmap (default: pairwise_identity_heatmap.png) |br|
*--json*: optional argument to print results as JSON

|

Parsimony informative sites
###########################
Function names: parsimony_informative_sites; pis |br|
Command line interface: pk_parsimony_informative_sites; pk_pis

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
doi: 10.1093/gbe/evw179 and Steenwyk et al., PLOS Biology
(2020), doi: 10.1371/journal.pbio.3001007.

.. code-block:: shell

		phykit parsimony_informative_sites <alignment> [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--json*: optional argument to print results as JSON

|

Protein-to-nucleotide alignment
###############################
Function names: thread_dna; pal2nal, p2n |br|
Command line interface: pk_thread_dna; pk_pal2nal, pk_p2n

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

.. code-block:: shell

   phykit thread_dna -p <file> -n <file> [-s] [--json]

Options: |br|
*-p/\\-\\-protein*: protein alignment file |br|
*-n/\\-\\-nucleotide*: nucleotide sequence file |br|
*-c/\\-\\-clipkit_log*: clipkit outputted log file |br|
*-s/\\-\\-stop*: boolean for whether or not stop codons should be kept. 
If used, stop codons will be removed. |br|
*--json*: optional argument to print results as JSON

|

Relative composition variability
################################
Function names: relative_composition_variability; rel_comp_var; rcv |br|
Command line interface: pk_relative_composition_variability; pk_rel_comp_var; pk_rcv

Calculate RCV (relative composition variability) for an alignment.

Lower RCV values are thought to be desirable because they represent
a lower composition bias in an alignment. Statistically, RCV describes
the average variability in sequence composition among taxa. 

RCV is calculated following Phillips and Penny, Molecular Phylogenetics
and Evolution (2003), doi: 10.1016/S1055-7903(03)00057-5.

RCV calculations are case-insensitive. Gap and ambiguous characters are
excluded from composition counts and correction terms, and each taxon is
normalized by its valid (non-excluded) sequence length.

.. code-block:: shell

		phykit relative_composition_variability <alignment> [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--json*: optional argument to print results as JSON

|

Relative composition variability, taxon
#######################################
Function names: relative_composition_variability_taxon; rel_comp_var_taxon; rcvt |br|
Command line interface: pk_relative_composition_variability_taxon; pk_rel_comp_var_taxon; pk_rcvt

Calculate RCVT (relative composition variability, taxon) for an alignment.

RCVT is the relative composition variability metric for individual taxa.
This facilitates identifying specific taxa that may have compositional
biases. Lower RCVT values are more desirable because they indicate
a lower composition bias for a given taxon in an alignment.

RCVT calculations are case-insensitive and exclude gap/ambiguous symbols
from composition counts and normalization.

.. code-block:: shell

	phykit relative_composition_variability_taxon <alignment> [--plot] [--plot-output <path>] [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--plot*: optional argument to generate an RCVT per-taxon barplot |br|
*--plot-output*: output path for the RCVT plot (default: ``rcvt_plot.png``) |br|
*--json*: optional argument to print results as JSON

|

Rename FASTA entries
####################
Function names: rename_fasta_entries; rename_fasta |br|
Command line interface: pk_rename_fasta_entries; pk_rename_fasta

Renames fasta entries.

Renaming fasta entries will follow the scheme of a tab-delimited
file wherein the first column is the current fasta entry name and
the second column is the new fasta entry name in the resulting 
output alignment. Note, the input fasta file does not need to be
an alignment file.

.. code-block:: shell

	phykit rename_fasta_entries <fasta> -i/--idmap <idmap> [-o/--output <output_file>] [--json]

Options: |br|
*<fasta>*: first argument after function name should be a FASTA file |br|
*-i/\\-\\-idmap*: identifier map of current FASTA names (col1) and desired FASTA names (col2) |br|
*--json*: optional argument to print results as JSON

|

Sum-of-pairs score
##################
Function names: sum_of_pairs_score; sops; sop |br|
Command line interface: pk_sum_of_pairs_score; pk_sops; pk_sop

Calculates sum-of-pairs score.

Sum-of-pairs is an accuracy metric for a multiple alignment relative
to a reference alignment. It is calculated by summing the correctly
aligned residue pairs over all pairs of sequences. Thus, values range
from 0 to 1 and higher values indicate more accurate alignments.

Column score is calculated following Thompson et al., Nucleic
Acids Research (1999), doi: 10.1093/nar/27.13.2682.

.. code-block:: shell

	phykit sum_of_pairs_score <alignment> --reference <reference_alignment> [--json]

Options: |br|
*<alignment>*: first argument after function name should be a query
fasta alignment file to be scored for accuracy |br|
*-r/\\-\\-reference*: reference alignment to compare the query alignment
to |br|
*--json*: optional argument to print results as JSON

|

Variable sites
##############
Function names: variable_sites; vs |br|
Command line interface: pk_variable_sites; pk_vs

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

.. code-block:: shell

   phykit variable_sites <alignment> [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--json*: optional argument to print results as JSON

|

Tree-based functions
--------------------

Ancestral state reconstruction
##############################
Function names: ancestral_state_reconstruction; asr; anc_recon |br|
Command line interface: pk_ancestral_state_reconstruction; pk_asr; pk_anc_recon

Estimate ancestral states for continuous traits using maximum likelihood,
analogous to R's ``phytools::fastAnc()`` and ``ape::ace(type="ML")``.
Optionally produce a contMap plot showing continuous trait values mapped
onto the phylogeny.

Two methods are available:

- **fast** (default): Felsenstein's pruning/contrasts shortcut, O(n) time
- **ml**: full VCV-based ML with exact conditional CIs, O(n^3)

Both methods produce identical point estimates; ``ml`` gives exact
conditional confidence intervals.

Input trait data can be either a two-column file (``taxon<tab>value``)
when ``-c`` is omitted, or a multi-trait file with header row when ``-c``
specifies which column to use.

.. code-block:: shell

   phykit ancestral_state_reconstruction -t <tree> -d <trait_data> [-c <trait>] [-m <method>] [--ci] [--plot <output>] [--json]

Options: |br|
*-t/\\-\\-tree*: a phylogenetic tree file |br|
*-d/\\-\\-trait_data*: trait data file (two-column or multi-trait with header) |br|
*-c/\\-\\-trait*: trait column name (required for multi-trait files) |br|
*-m/\\-\\-method*: method to use: ``fast`` or ``ml`` (default: ``fast``) |br|
*--ci*: include 95% confidence intervals |br|
*--plot*: output path for contMap plot (requires matplotlib) |br|
*--json*: output results as JSON

|

Bipartition support statistics
##############################
Function names: bipartition_support_stats; bss |br|
Command line interface: pk_bipartition_support_stats; pk_bss

Calculate summary statistics for bipartition support.

High bipartition support values are thought to be desirable because
they are indicative of greater certainty in tree topology.

To obtain all bipartition support values, use the -v/\\-\\-verbose option.
In addition to support values for each node, the names of all terminal
branches tips are also included. Each terminal branch name is separated
with a semi-colon (;).

.. code-block:: shell

   phykit bipartition_support_stats <tree> [-v/--verbose]
       [--thresholds <comma-separated-floats>] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-v/\\-\\-verbose*: optional argument to print all bipartition support values |br|
*--thresholds*: optional comma-separated support cutoffs; prints count and
fraction of bipartitions below each cutoff |br|
*--json*: optional argument to print results as JSON

Example JSON output (summary mode):

.. code-block:: shell

   phykit bipartition_support_stats test.tre --thresholds 70,90 --json
   {"summary": {"maximum": 100, "mean": 95.71428571428571, "median": 100, "minimum": 85, "seventy_fifth": 100.0, "standard_deviation": 7.319250547113999, "twenty_fifth": 92.5, "variance": 53.57142857142857}, "thresholds": [{"count_below": 0, "fraction_below": 0.0, "threshold": 70.0}, {"count_below": 2, "fraction_below": 0.2857142857142857, "threshold": 90.0}], "verbose": false}

Example JSON output (verbose mode):

.. code-block:: shell

   phykit bipartition_support_stats test.tre -v --json
   {"bipartitions": [{"support": 85, "terminals": ["taxon_a", "taxon_b"]}, {"support": 100, "terminals": ["taxon_c", "taxon_d"]}], "thresholds": [], "verbose": true}

|

Branch length multiplier
########################
Function names: branch_length_multiplier; blm |br|
Command line interface: pk_branch_length_multiplier; pk_blm

Multiply branch lengths in a phylogeny by a given factor.
                
This can help modify reference trees when conducting simulations
or other analyses.  

.. code-block:: shell

   phykit branch_length_multiplier <tree> -f n [-o/--output <output_file>] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-f/\\-\\-factor*: factor to multiply branch lengths by |br|
*-o/\\-\\-output*: optional argument to name the outputted tree file. Default 
output will have the same name as the input file but with the suffix ".factor_(n).tre" |br|
*--json*: optional argument to print results as JSON

|

Collapse bipartitions
#####################
Function names: collapse_branches, collapse, cb |br|
Command line interface: pk_collapse_branches, pk_collapse, pk_cb

Collapse branches on a phylogeny according to bipartition support.

Bipartitions will be collapsed if they are less than the user specified
value.    

.. code-block:: shell

   phykit collapse_branches <tree> -s/--support n [-o/--output <output_file>] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-s/\\-\\-support*: bipartitions with support less than this value will be 
collapsed |br|
*-o/\\-\\-output*: optional argument to name the outputted tree file. Default 
output will have the same name as the input file but with the suffix 
".collapsed_(support).tre" |br|
*--json*: optional argument to print results as JSON

|

Covarying evolutionary rates
############################
Function names: covarying_evolutionary_rates; cover |br|
Command line interface: pk_covarying_evolutionary_rates; pk_cover

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

.. code-block:: shell

   phykit covarying_evolutionary_rates <tree_file_zero> <tree_file_one> -r/--reference <reference_tree_file> [-v/--verbose] [--plot] [--plot-output <path>] [--json]

Options: |br|
*<tree_file_zero>*: first argument after function name should be an alignment file |br|
*<tree_file_one>*: first argument after function name should be an alignment file |br| 
*-r/\\-\\-reference*: a tree to correct branch lengths by in the two input trees. Typically, 
this is a putative species tree. |br|
*-v/\\-\\-verbose*: print out corrected branch lengths shared between tree 0 and tree 1 |br|
*--plot*: save a covarying-rates scatter plot |br|
*--plot-output*: output path for plot (default: ``covarying_rates_plot.png``) |br|
*--json*: optional argument to print results as JSON

|

Degree of violation of the molecular clock
##########################################
Function names: degree_of_violation_of_a_molecular_clock, dvmc |br|
Command line interface: pk_degree_of_violation_of_a_molecular_clock, pk_dvmc

Calculate degree of violation of a molecular clock (or DVMC) in a phylogeny.

Lower DVMC values are thought to be desirable because they are indicative
of a lower degree of violation in the molecular clock assumption.

Typically, outgroup taxa are not included in molecular clock analysis. Thus,
prior to calculating DVMC from a single gene tree, users may want to prune
outgroup taxa from the phylogeny. To prune tips from a phylogeny, see the 
prune_tree function. 

Calculate DVMC in a tree following Liu et al., PNAS (2017), doi: 10.1073/pnas.1616744114.

.. code-block:: shell

   phykit degree_of_violation_of_a_molecular_clock <tree> [--json]

Options: |br|
*<tree>*: input file tree name |br|
*--json*: optional argument to print results as JSON

|

Evolutionary rate
#################
Function names: evolutionary_rate, evo_rate |br|
Command line interface: pk_evolutionary_rate, pk_evo_rate

Calculate a tree-based estimation of the evolutionary rate of a gene.

Evolutionary rate is the total tree length divided by the number
of terminals.

Calculate evolutionary rate following Telford et al., Proceedings
of the Royal Society B (2014). 

.. code-block:: shell

   phykit evolutionary_rate <tree> [--json]

Options: |br|
*<tree>*: input file tree name |br|
*--json*: optional argument to print results as JSON

|

Hidden paralogy check
#####################
Function names: hidden_paralogy_check, clan_check |br|
Command line interface: pk_hidden_paralogy_check, pk_clan_check

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
clade file should be formatted as follows: |br|
" |br|
A B C |br|
D E F |br|
"

The output will report if the specified taxa were monophyletic
or not. The number of rows will reflect how many groups of taxa
were checked for monophyly. For example,
if there were three rows of clades in the -c file, there will be
three rows in the output
where the first row in the output corresponds to the 
results of the first row in the clade file. |br|

The concept behind this analysis follows
Siu-Ting et al., Molecular Biology and Evolution (2019),
doi: 10.1093/molbev/msz067.

.. code-block:: shell

   phykit hidden_paralogy_check <tree> -c/--clade <clade_file> [--json]

Options: |br|
*-t/\\-\\-tree*: input file tree name |br|
*-c/\\-\\-clade*: clade file detailing which monophyletic lineages should
be scanned for |br|
*--json*: optional argument to print results as JSON

|

Internal branch statistics
##########################
Function names: internal_branch_stats; ibs |br|
Command line interface: pk_internal_branch_stats; pk_ibs

Calculate summary statistics for internal branch lengths in a phylogeny.

Internal branch lengths can be useful for phylogeny diagnostics.

To obtain all internal branch lengths, use the -v/\\-\\-verbose option.   

.. code-block:: shell

   phykit internal_branch_stats <tree> [-v/--verbose] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-v/\\-\\-verbose*: optional argument to print all internal branch lengths |br|
*--json*: optional argument to print results as JSON

|

Internode labeler
#################
Function names: internode_labeler; il |br|
Command line interface: pk_internode_labeler; pk_il

Appends numerical identifiers to bipartitions in place of support values.
This is helpful for pointing to specific internodes in supplementary files
or otherwise.  

.. code-block:: shell

   phykit internode_labeler <tree> [-o/--output <file>] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-o/\\-\\-output*: optional argument to name the outputted tree file |br|
*--json*: optional argument to print results as JSON

|

Last common ancestor subtree
############################
Function names: last_common_ancestor_subtree; lca_subtree |br|
Command line interface: pk_last_common_ancestor_subtree; pk_lca_subtree

Obtains subtree from a phylogeny by getting the last common ancestor
from a list of taxa.

.. code-block:: shell

   phykit last_common_ancestor_subtree <file> <list_of_taxa> [-o/--output <file>] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*<list_of_taxa>*: second argument after function name should be a single column
file with the list of taxa to get the last common ancestor subtree for
*-o/\\-\\-output*: optional argument to name the outputted tree file |br|
*--json*: optional argument to print results as JSON

|

Long branch score
#################
Function names: long_branch_score; lb_score; lbs |br|
Command line interface: pk_long_branch_score; pk_lb_score; pk_lbs

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

.. code-block:: shell

   phykit long_branch_score <tree> [-v/--verbose] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-v/\\-\\-verbose*: optional argument to print all LB score values |br|
*--json*: optional argument to print results as JSON

|

Monophyly check
###############
Function names: monophyly_check; is_monophyletic |br|
Command line interface: pk_monophyly_check; pk_is_monophyletic

This analysis can be used to determine if a set of 
taxa are exclusively monophyletic. By exclusively monophyletic,
if other taxa are in the same clade, the lineage will not be
considered exclusively monophyletic.

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
col 6: tip names of taxa monophyletic with the lineage of interest excluding those that are listed in the taxa_of_interest file

.. code-block:: shell

   phykit monophyly_check <tree> <list_of_taxa> [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*<list_of_taxa>*: single column file with list of tip names to 
examine the monophyly of |br|
*--json*: optional argument to print results as JSON

|

Nearest neighbor interchange
############################
Function names: nearest_neighbor_interchange; nni |br|
Command line interface: pk_nearest_neighbor_interchange; pk_nni

Generate all nearest neighbor interchange moves for a binary
rooted tree.

By default, the output file will have the same name as the input
file but with the suffix ".nnis"

The output file will also include the original phylogeny.

.. code-block:: shell

   phykit nearest_neighbor_interchange <tree> [-o/--output <output_file>] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-o/\\-\\-output*: optional argument to specify output file name |br|
*--json*: optional argument to print summary metadata as JSON

|

Patristic distances
###################
Function names: patristic_distances; pd |br|
Command line interface: pk_patristic_distances; pk_pd

Calculate summary statistics among patristic distances in a phylogeny.

Patristic distances are all tip-to-tip distances in a phylogeny.

To obtain all patristic distances, use the -v/--verbose option.
With the -v option, the first column will have two taxon names
separated by a '-' followed by the patristic distance. Features
will be tab separated. 

.. code-block:: shell

   phykit patristic_distances <tree> [-v/--verbose] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-v/\\-\\-verbose*: optional argument to print all tip-to-tip distances |br|
*--json*: optional argument to print results as JSON

|

Phylogenetic signal
####################
Function names: phylogenetic_signal; phylo_signal; ps |br|
Command line interface: pk_phylogenetic_signal; pk_phylo_signal; pk_ps

Calculate phylogenetic signal for continuous trait data on a phylogeny.

Two methods are available:

- **Blomberg's K** (Blomberg et al. 2003): measures the degree of phylogenetic
  signal relative to expectation under Brownian motion. K = 1 indicates trait
  variation consistent with BM; K < 1 indicates less phylogenetic signal than
  expected; K > 1 indicates more. P-value is computed via permutation test.
- **Pagel's lambda** (Pagel 1999): a tree-scaling parameter estimated by
  maximum likelihood. Lambda = 0 indicates no phylogenetic signal; lambda = 1
  indicates trait evolution consistent with BM. P-value is computed via
  likelihood ratio test against lambda = 0.

The trait file should be tab-delimited with two columns (taxon_name<tab>trait_value).
Lines starting with '#' are treated as comments. If the tree and trait file have
different taxa, the intersection is used and warnings are printed to stderr.

Output for Blomberg's K: K_value<tab>p_value |br|
Output for Pagel's lambda: lambda_value<tab>log_likelihood<tab>p_value

Results have been validated against the R package phytools (``phylosig`` function)
across 95 simulated datasets spanning diverse tree sizes (5-50 tips), topologies
(pure-birth, coalescent), trait models (random, Brownian motion, known lambda),
and branch length scales. All metrics show Pearson r > 0.999 with phytools.

.. image:: ../_static/docs_img/phylogenetic_signal_validation.png
   :align: center

|

.. code-block:: shell

   phykit phylogenetic_signal -t <tree> -d <trait_data> [-m <method>] [-p <permutations>] [--json]

Options: |br|
*-t/\\-\\-tree*: a tree file in Newick format |br|
*-d/\\-\\-trait_data*: tab-delimited trait file (taxon_name<tab>trait_value) |br|
*-m/\\-\\-method*: method to use: ``blombergs_k`` or ``lambda`` (default: blombergs_k) |br|
*-p/\\-\\-permutations*: number of permutations for Blomberg's K (default: 1000) |br|
*--json*: optional argument to print results as JSON

|

Phylogenetic Ordination
#######################
Function names: phylogenetic_ordination; phylo_ordination; ordination; ord;
phylo_pca; phyl_pca; ppca; phylo_dimreduce; dimreduce; pdr |br|
Command line interface: pk_phylogenetic_ordination; pk_phylo_ordination;
pk_ordination; pk_ord; pk_phylo_pca; pk_phyl_pca; pk_ppca;
pk_phylo_dimreduce; pk_dimreduce; pk_pdr

Perform phylogenetic ordination (PCA, t-SNE, or UMAP) on continuous multi-trait
data while accounting for phylogenetic non-independence among species.

All methods use GLS-centering via the phylogenetic variance-covariance matrix.
For PCA, eigendecomposition of the evolutionary rate matrix is performed to
extract principal components. For t-SNE and UMAP, nonlinear embedding is
applied to the GLS-centered data.

Three ordination methods are available via ``--method``:

- **pca** (default): phylogenetic PCA (Revell 2009). Performs eigendecomposition
  of the evolutionary rate matrix.
- **tsne**: t-SNE embedding of the GLS-centered data. Perplexity is auto-set
  to ``min(30, (n-1)/3)``; requires at least 4 taxa.
- **umap**: UMAP embedding of the GLS-centered data. n_neighbors is auto-set
  to ``min(15, n-1)``; requires at least 3 taxa.

Two phylogenetic correction modes are available via ``--correction``:

- **BM** (default): assumes traits evolved under Brownian motion. The
  phylogenetic VCV matrix is used directly.
- **lambda**: jointly estimates a single Pagel's lambda parameter across all
  traits by maximum likelihood, then rescales the off-diagonal elements of the
  VCV matrix. This accounts for deviations from Brownian motion.

For PCA, two modes are available via ``--mode``:

- **cov** (default): PCA on the evolutionary rate (covariance) matrix.
- **corr**: PCA on the evolutionary rate correlation matrix, which
  standardizes traits to equal variance before extracting components.

The multi-trait input file should be tab-delimited with a header row:
``taxon<tab>trait1<tab>trait2<tab>...``
Lines starting with '#' are treated as comments. If the tree and trait file
have different taxa, the intersection is used and warnings are printed to
stderr.

PCA output includes eigenvalues with proportion of variance explained, trait
loadings, and taxon scores for each principal component. t-SNE/UMAP output
includes method parameters and embedding coordinates. When ``correction=lambda``
is used, the estimated lambda and log-likelihood are also reported.

PCA results have been benchmarked against the R package
`phytools <https://cran.r-project.org/package=phytools>`_
(``phyl.pca`` function; Revell 2012). Eigenvalues, loadings, and scores match
phytools across all method/mode combinations (BM+cov, BM+corr, lambda+cov,
lambda+corr) within numerical tolerance (1e-4).

.. image:: ../_static/docs_img/phylogenetic_pca_plot.png
   :align: center

|

.. image:: ../_static/docs_img/phylogenetic_tsne_plot.png
   :align: center

|

.. image:: ../_static/docs_img/phylogenetic_umap_plot.png
   :align: center

|

.. code-block:: shell

   phykit phylogenetic_ordination -t <tree> -d <trait_data> [--method <pca|tsne|umap>] [--correction <BM|lambda>] [--mode <cov|corr>] [--n-components <int>] [--perplexity <float>] [--n-neighbors <int>] [--min-dist <float>] [--seed <int>] [--plot] [--plot-tree] [--no-plot-tree] [--color-by <col_or_file>] [--tree-color-by <col_or_file>] [--plot-output <path>] [--json]

Options: |br|
*-t/\\-\\-tree*: a tree file in Newick format |br|
*-d/\\-\\-trait_data*: tab-delimited multi-trait file with header row |br|
*--method*: ordination method: ``pca``, ``tsne``, or ``umap`` (default: pca) |br|
*--correction*: phylogenetic correction: ``BM`` or ``lambda`` (default: BM) |br|
*--mode*: PCA mode: ``cov`` or ``corr`` (default: cov; PCA only) |br|
*--n-components*: number of embedding dimensions (default: 2; tsne/umap only) |br|
*--perplexity*: t-SNE perplexity (default: auto) |br|
*--n-neighbors*: UMAP n_neighbors (default: auto) |br|
*--min-dist*: UMAP min_dist (default: 0.1) |br|
*--seed*: random seed for reproducibility |br|
*--plot*: optional argument to save a scatter plot |br|
*--plot-tree*: overlay the phylogeny via ancestral reconstruction (default for tsne/umap; opt-in for pca) |br|
*--no-plot-tree*: disable the phylogeny overlay for tsne/umap plots |br|
*--color-by*: color tip points by a trait; specify a column name from the multi-trait file or a separate tab-delimited file (taxon<tab>value) for continuous or discrete coloring |br|
*--tree-color-by*: color phylogeny edges by a trait; specify a column name or a tab-delimited file (default: distance from root) |br|
*--plot-output*: output path for plot (default: phylo_ordination_plot.png) |br|
*--json*: optional argument to print results as JSON

|

Phylomorphospace
################
Function names: phylomorphospace; phylomorpho; phmo |br|
Command line interface: pk_phylomorphospace; pk_phylomorpho; pk_phmo

Plot a phylomorphospace: two raw traits in trait space with the phylogeny
overlaid via ML-reconstructed ancestral states at internal nodes. This differs
from the ``phylogenetic_ordination --plot-tree`` option, which plots in PC space;
``phylomorphospace`` plots raw trait values directly on the x and y axes.

Tree edges connect species through ML-reconstructed ancestral states and are
colored by distance from root (coolwarm colormap with colorbar). Tip points
default to blue, or can be colored by a continuous or discrete variable using
the ``--color-by`` option.

The multi-trait input file should be tab-delimited with a header row:
``taxon<tab>trait1<tab>trait2<tab>...``
Lines starting with '#' are treated as comments. If the tree and trait file
have different taxa, the intersection is used and warnings are printed to
stderr.

If the trait file has exactly 2 trait columns and ``--trait-x`` / ``--trait-y``
are omitted, the first two columns are selected automatically.

.. image:: ../_static/docs_img/phylomorphospace_plot.png
   :align: center

|

.. code-block:: shell

   phykit phylomorphospace -t <tree> -d <trait_data> [--trait-x <name>] [--trait-y <name>] [--color-by <col_or_file>] [--plot-output <path>] [--json]

Options: |br|
*-t/\\-\\-tree*: a tree file in Newick format |br|
*-d/\\-\\-trait_data*: tab-delimited multi-trait file with header row |br|
*--trait-x*: column name for x-axis trait |br|
*--trait-y*: column name for y-axis trait |br|
*--color-by*: color tip points by a trait; specify a column name from the multi-trait file or a separate tab-delimited file (taxon<tab>value) for continuous or discrete coloring |br|
*--plot-output*: output path for plot (default: phylomorphospace_plot.png) |br|
*--json*: optional argument to print results as JSON

|

Phylogenetic regression (PGLS)
##############################
Function names: phylogenetic_regression; phylo_regression; pgls |br|
Command line interface: pk_phylogenetic_regression; pk_phylo_regression; pk_pgls

Fit a Phylogenetic Generalized Least Squares (PGLS) regression while
accounting for phylogenetic non-independence among species, analogous to
R's ``caper::pgls()``.

The multi-trait input file should be tab-delimited with a header row:
``taxon<tab>trait1<tab>trait2<tab>...``
Lines starting with '#' are treated as comments. If the tree and trait file
have different taxa, the intersection is used and warnings are printed to
stderr.

Two methods are available:

- **BM** (default): Brownian motion with lambda fixed at 1
- **lambda**: jointly estimates Pagel's lambda via maximum likelihood

Output includes coefficient estimates, standard errors, t-values, p-values,
R-squared, adjusted R-squared, F-statistic, log-likelihood, and AIC.

The implementation uses the raw phylogenetic variance-covariance (VCV) matrix
for GLS estimation, matching the approach used by R's ``caper::pgls()``.
Note that this differs from ``nlme::gls()`` with ``corBrownian``, which
normalizes the VCV to a correlation matrix (ones on the diagonal). For
non-ultrametric trees the two approaches yield different coefficient estimates;
PhyKIT follows the ``caper`` convention.

All results have been validated against R 4.4.0 using manual GLS with
``ape::vcv()`` (raw VCV). Coefficients, standard errors, t-values, p-values,
R-squared, F-statistic, log-likelihood, and AIC match R to at least four
decimal places for both simple and multiple regression.

.. code-block:: shell

   phykit phylogenetic_regression -t <tree> -d <trait_data> -y <response> -x <predictor1> [predictor2 ...] [-m <method>] [--json]

Options: |br|
*-t/\\-\\-tree*: a tree file in Newick format |br|
*-d/\\-\\-trait_data*: tab-delimited multi-trait file with header row |br|
*-y/\\-\\-response*: response (dependent) variable column name |br|
*-x/\\-\\-predictors*: one or more predictor column names |br|
*-m/\\-\\-method*: method to use: BM or lambda (default: BM) |br|
*--json*: optional argument to print results as JSON

|

Phylogenetic GLM
################
Function names: phylogenetic_glm; phylo_glm; pglm |br|
Command line interface: pk_phylogenetic_glm; pk_phylo_glm; pk_pglm

Fit a Phylogenetic Generalized Linear Model (GLM) for binary or count
response data while accounting for phylogenetic non-independence among
species.

Two families are supported:

- **binomial**: logistic regression via Maximum Penalized Likelihood
  Estimation (logistic_MPLE; Ives & Garland 2010). Uses Firth's penalty
  to prevent bias from complete/quasi-complete separation, and jointly
  estimates the phylogenetic signal parameter alpha via a two-state
  continuous-time Markov chain on the tree.
- **poisson**: Poisson regression via Generalized Estimating Equations
  (poisson_GEE; Paradis & Claude 2002). Uses Fisher scoring with the
  phylogenetic correlation matrix and reports an overdispersion parameter.

The multi-trait input file should be tab-delimited with a header row:
``taxon<tab>trait1<tab>trait2<tab>...``

Output includes coefficient estimates, standard errors, z-values,
p-values, log-likelihood, and AIC.

.. code-block:: shell

   phykit phylogenetic_glm -t <tree> -d <trait_data> -y <response> -x <predictor1> [predictor2 ...] --family <binomial|poisson> [--json]

Options: |br|
*-t/\\-\\-tree*: a tree file in Newick format |br|
*-d/\\-\\-trait_data*: tab-delimited multi-trait file with header row |br|
*-y/\\-\\-response*: response (dependent) variable column name |br|
*-x/\\-\\-predictors*: one or more predictor column names |br|
*--family*: distribution family: binomial or poisson |br|
*--method*: estimation method: logistic_MPLE or poisson_GEE (auto from family) |br|
*--btol*: linear predictor bound for logistic model (default: 10) |br|
*--log-alpha-bound*: bound on log(alpha) for logistic model (default: 4) |br|
*--json*: optional argument to print results as JSON

|

Stochastic character mapping (SIMMAP)
#####################################
Function names: stochastic_character_map; simmap; scm |br|
Command line interface: pk_stochastic_character_map; pk_simmap; pk_scm

Perform Stochastic Character Mapping (SIMMAP) of discrete traits onto a
phylogeny (Huelsenbeck et al. 2003; Bollback 2006), analogous to R's
``phytools::make.simmap()``.

Fits a continuous-time Markov chain (CTMC) rate matrix Q via maximum
likelihood using Felsenstein's pruning algorithm. Three substitution models
are available:

- **ER** (equal rates): 1 free parameter
- **SYM** (symmetric rates): k(k-1)/2 free parameters
- **ARD** (all rates differ): k(k-1) free parameters

The input trait file should be tab-delimited with a header row:
``taxon<tab>trait_column<tab>...``
Lines starting with '#' are treated as comments. If the tree and trait file
have different taxa, the intersection is used and warnings are printed to
stderr.

Output includes the fitted Q matrix, log-likelihood, mean dwelling times,
mean transition counts, and posterior probabilities at internal nodes.

Optionally, a horizontal phylogram plot with branches colored by mapped
character state can be generated using the ``--plot`` argument.

All results have been validated against R 4.4.0 using ``phytools::fitMk``
and ``phytools::make.simmap``. ER log-likelihoods match R to within 0.002
(both optimizers converge to the same flat plateau). For the ARD model,
PhyKIT's multi-start optimizer finds a slightly better local optimum
than R (loglik -8.38 vs -8.43). Dwelling times sum to total tree length
exactly, and Q matrix structural properties (rows sum to zero, model
nesting) are verified.

.. code-block:: shell

   phykit stochastic_character_map -t <tree> -d <trait_data> -c <trait_column> [-m <model>] [-n <nsim>] [--seed <seed>] [--plot <output.png>] [--json]

Options: |br|
*-t/\\-\\-tree*: a tree file in Newick format |br|
*-d/\\-\\-trait_data*: tab-delimited trait file with header row |br|
*-c/\\-\\-trait*: column name for discrete character trait |br|
*-m/\\-\\-model*: substitution model: ER, SYM, or ARD (default: ER) |br|
*-n/\\-\\-nsim*: number of stochastic mapping simulations (default: 100) |br|
*--seed*: random seed for reproducibility |br|
*--plot*: output plot file path for phylogram with colored branches |br|
*--json*: optional argument to print results as JSON

|

Rate heterogeneity test (multi-rate Brownian motion)
####################################################
Function names: rate_heterogeneity; brownie; rh |br|
Command line interface: pk_rate_heterogeneity; pk_brownie; pk_rh

Test for rate heterogeneity across phylogenetic regimes using multi-rate
Brownian motion (O'Meara et al. 2006), analogous to R's
``phytools::brownie.lite()``.

Fits single-rate vs. multi-rate BM models and performs a likelihood ratio
test. Users specify a tree, continuous trait data, and a regime file
mapping tips to regimes (tab-delimited, ``taxon<tab>regime_label``).

Regime assignments to internal branches are inferred via Fitch parsimony.
Per-regime VCV matrices are decomposed and per-regime sigma-squared values
are estimated via maximum likelihood.

Optionally, a parametric bootstrap can be run to compute a simulated
p-value (``-n/--nsim``). A horizontal phylogram with branches colored by
regime can be generated using ``--plot``.

.. code-block:: shell

   phykit rate_heterogeneity -t <tree> -d <trait_data> -r <regime_data> [-n <nsim>] [--seed <seed>] [--plot <output.png>] [--json]

Options: |br|
*-t/\\-\\-tree*: a tree file in Newick format |br|
*-d/\\-\\-trait_data*: tab-delimited trait file (taxon<tab>value) |br|
*-r/\\-\\-regime_data*: tab-delimited regime file (taxon<tab>regime_label) |br|
*-n/\\-\\-nsim*: number of parametric bootstrap simulations (default: 0) |br|
*--seed*: random seed for reproducibility |br|
*--plot*: output plot file path for phylogram with colored branches |br|
*--json*: optional argument to print results as JSON

|

Continuous trait evolution model comparison (fitContinuous)
##########################################################
Function names: fit_continuous; fitcontinuous; fc |br|
Command line interface: pk_fit_continuous; pk_fitcontinuous; pk_fc

Compare models of continuous trait evolution on a phylogeny, analogous to
R's ``geiger::fitContinuous()``. Fits up to 7 models and ranks them by
AIC, BIC, and AIC weights.

Models:

- **BM** -- Brownian motion (baseline, 2 params)
- **OU** -- Ornstein-Uhlenbeck / stabilizing selection (3 params)
- **EB** -- Early Burst (Harmon et al. 2010) (3 params)
- **Lambda** -- Pagel's lambda / phylogenetic signal (3 params)
- **Delta** -- Pagel's delta / tempo of evolution (3 params)
- **Kappa** -- Pagel's kappa / punctuational vs gradual (3 params)
- **White** -- White noise / no phylogenetic signal (2 params)

.. code-block:: shell

   phykit fit_continuous -t <tree> -d <trait_data> [--models BM,OU,Lambda] [--json]

Options: |br|
*-t/\\-\\-tree*: a tree file in Newick format |br|
*-d/\\-\\-trait_data*: tab-delimited trait file (taxon<tab>value) |br|
*--models*: comma-separated list of models to fit (default: all 7) |br|
*--json*: optional argument to print results as JSON

Example output:

.. code-block:: text

   Model Comparison (fitContinuous)

   Number of tips: 8

   Model       Param     Value      Sigma2    z0        LL         AIC     dAIC    AICw    BIC     dBIC
   BM          -         -          0.0384    1.6447    -11.570    27.14   0.00    0.453   27.83   0.00
   OU          alpha     0.0012     0.0385    1.6420    -11.568    29.14   2.00    0.167   30.18   2.35
   ...

   Best model (AIC): BM
   Best model (BIC): BM

|

Multi-regime OU models (OUwie)
##############################
Function names: ouwie; fit_ouwie; multi_regime_ou |br|
Command line interface: pk_ouwie; pk_fit_ouwie; pk_multi_regime_ou

Fit multi-regime Ornstein-Uhlenbeck models of continuous trait evolution,
analogous to R's ``OUwie`` package (Beaulieu et al. 2012). Fits up to 7
models and ranks them by AICc, BIC, and AICc weights. Regime assignments
to internal branches are inferred via Fitch parsimony.

Models:

- **BM1** -- single-rate Brownian motion (2 params)
- **BMS** -- multi-rate Brownian motion with per-regime sigma2 (R+1 params)
- **OU1** -- single-regime Ornstein-Uhlenbeck (3 params)
- **OUM** -- multi-regime OU with per-regime trait optima (R+2 params)
- **OUMV** -- OUM + per-regime sigma2 (2R+1 params)
- **OUMA** -- OUM + per-regime alpha (2R+1 params)
- **OUMVA** -- all parameters regime-specific (3R params)

.. code-block:: shell

   phykit ouwie -t <tree> -d <trait_data> -r <regime_data> [--models BM1,OUM,OUMVA] [--json]

Options: |br|
*-t/\\-\\-tree*: a tree file in Newick format |br|
*-d/\\-\\-trait_data*: tab-delimited trait file (taxon<tab>value) |br|
*-r/\\-\\-regime_data*: tab-delimited regime file (taxon<tab>regime_label) |br|
*--models*: comma-separated list of models to fit (default: all 7) |br|
*--json*: optional argument to print results as JSON

The trait data file is a two-column tab-delimited file mapping taxon names
to continuous trait values:

.. code-block:: none

   dog	1.1
   bear	1.9
   raccoon	1.5
   seal	1.8
   sea_lion	1.8
   cat	0.5
   weasel	1.7
   monkey	0.3

The regime data file is a two-column tab-delimited file mapping taxon names
to discrete regime labels (e.g., habitat, diet category):

.. code-block:: none

   dog	terrestrial
   bear	terrestrial
   raccoon	terrestrial
   seal	aquatic
   sea_lion	aquatic
   cat	terrestrial
   weasel	terrestrial
   monkey	terrestrial

Example output:

.. code-block:: none

   OUwie Model Comparison
   ======================
   Regimes: aquatic, terrestrial

   Model      logLik     AICc      BIC     k  AICc_w  Params
   -----  ----------  -------  -------  ----  ------  ------
   OUMVA     -6.9859  27.9717  29.5459     6  0.0040  alpha={aquatic:0.38, terrestrial:0.38}, sigma2={aquatic:0.01, terrestrial:0.05}, theta={aquatic:1.80, terrestrial:1.24}
   OUMA      -6.9859  27.9717  29.0119     5  0.0040  alpha={aquatic:0.38, terrestrial:0.38}, sigma2=0.0384, theta={aquatic:1.80, terrestrial:1.24}
   OUMV      -6.9859  27.9717  29.0119     5  0.0040  alpha=0.3849, sigma2={aquatic:0.01, terrestrial:0.05}, theta={aquatic:1.80, terrestrial:1.24}
   OUM       -8.6297  25.2594  26.3276     4  0.0488  alpha=0.0706, sigma2=0.0329, theta={aquatic:1.80, terrestrial:1.33}
   OU1      -10.2890  27.2447  28.0459     3  0.0063  alpha=0.0398, sigma2=0.0363, theta=1.64
   BMS      -11.2046  29.0759  29.8771     3  0.0024  sigma2={aquatic:0.01, terrestrial:0.05}, z0=1.64
   BM1      -11.5697  27.1393  27.6735     2  0.0073  sigma2=0.0384, z0=1.64

   Best model (AICc): OUM
   Best model (BIC): OUM

|

Continuous trait mapping (contMap)
##################################
Function names: cont_map; contmap; cmap |br|
Command line interface: pk_cont_map; pk_contmap; pk_cmap

Plot a phylogram with branches colored by continuous trait values
(analogous to R's ``phytools::contMap()``). Ancestral states are
estimated via maximum-likelihood (two-pass Felsenstein algorithm)
and mapped onto branches using a color gradient (coolwarm colormap).

.. code-block:: shell

   phykit cont_map -t <tree> -d <trait_data> -o <output.png> [--json]

Options: |br|
*-t/\\-\\-tree*: a tree file in Newick format |br|
*-d/\\-\\-trait_data*: tab-delimited trait file (taxon<tab>value) |br|
*-o/\\-\\-output*: output plot file path (required) |br|
*--json*: optional argument to print results as JSON

.. image:: ../_static/img/contmap_example.png
   :align: center
   :width: 80%

|

Density map
###########
Function names: density_map; densitymap; dmap |br|
Command line interface: pk_density_map; pk_densitymap; pk_dmap

Plot a phylogram with branches colored by posterior probabilities of
discrete character states from stochastic character mapping (analogous
to R's ``phytools::densityMap()``). Runs N simulations of stochastic
character mapping internally and, for each point along each branch,
computes the fraction of simulations in each state.

.. code-block:: shell

   phykit density_map -t <tree> -d <trait_data> -c <trait_column> -o <output.png> [-n <nsim>] [--seed <seed>] [--json]

Options: |br|
*-t/\\-\\-tree*: a tree file in Newick format |br|
*-d/\\-\\-trait_data*: tab-delimited trait file (taxon<tab>state) |br|
*-c/\\-\\-trait*: column name of the trait to map |br|
*-o/\\-\\-output*: output plot file path (required) |br|
*-n/\\-\\-nsim*: number of stochastic mapping simulations (default: 100) |br|
*--seed*: random seed for reproducibility |br|
*--json*: optional argument to print results as JSON

.. image:: ../_static/img/densitymap_example.png
   :align: center
   :width: 80%

|

Phenogram (traitgram)
#####################
Function names: phenogram; traitgram; tg |br|
Command line interface: pk_phenogram; pk_traitgram; pk_tg

Plot a phenogram (traitgram) showing trait evolution along a phylogeny
(analogous to R's ``phytools::phenogram()``). The X-axis represents
distance from the root and the Y-axis represents trait values.
Ancestral states are reconstructed via maximum-likelihood.

.. code-block:: shell

   phykit phenogram -t <tree> -d <trait_data> -o <output.png> [--json]

Options: |br|
*-t/\\-\\-tree*: a tree file in Newick format |br|
*-d/\\-\\-trait_data*: tab-delimited trait file (taxon<tab>value) |br|
*-o/\\-\\-output*: output plot file path (required) |br|
*--json*: optional argument to print results as JSON

.. image:: ../_static/img/phenogram_example.png
   :align: center
   :width: 80%

|

Cophylogenetic plot (tanglegram)
################################
Function names: cophylo; tanglegram; tangle |br|
Command line interface: pk_cophylo; pk_tanglegram; pk_tangle

Plot a cophylogenetic tanglegram of two phylogenies (analogous to R's
``phytools::cophylo()``). Draws two trees facing each other with
connecting lines between matching taxa. By default, taxa are matched
by identical tip names. Internal nodes of tree2 are rotated to minimize
line crossings.

.. code-block:: shell

   phykit cophylo -t <tree1> -t2 <tree2> -o <output.png> [-m <mapping>] [--json]

Options: |br|
*-t/\\-\\-tree1*: first tree file in Newick format |br|
*-t2/\\-\\-tree2*: second tree file in Newick format |br|
*-o/\\-\\-output*: output plot file path (required) |br|
*-m/\\-\\-mapping*: optional tab-delimited mapping file (taxon1<tab>taxon2) |br|
*--json*: optional argument to print results as JSON

.. image:: ../_static/img/cophylo_example.png
   :align: center
   :width: 80%

|

Polytomy testing
################
Function names: polytomy_test; polyt_test; polyt; ptt |br|
Command line interface: pk_polytomy_test; pk_polyt_test; pk_polyt; pk_ptt

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

.. code-block:: shell

   phykit polytomy_test -t/--trees <trees> -g/--groups <groups> [--json]

Options: |br|
*-t/\\-\\-trees <trees>*: single column file with the names of 
phylogenies to use for polytomy testing |br|
*-g/\\-\\-groups*: a tab-delimited file with the grouping designations
to test. Lines starting with commetns are not considered. Names of
individual taxa should be separated by a semi-colon ';' |br|
*--json*: optional argument to print results as JSON

For example, the groups file could look like the following:

.. code-block:: shell

   #label group0  group1  group2
   name_of_test    tip_name_A;tip_name_B   tip_name_C  tip_name_D;tip_name_E

|

Print tree
##########
Function names: print_tree; print; pt |br|
Command line interface: pk_print_tree; pk_print; pk_pt

Print ascii tree of input phylogeny.

Phylogeny can be printed with or without branch lengths.
By default, the phylogeny will be printed with branch lengths
but branch lengths can be removed using the -r/--remove argument.

.. code-block:: shell

   phykit print_tree <tree> [-r/--remove] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-r/\\-\\-remove*: optional argument to print the phylogeny without branch
lengths |br|
*--json*: optional argument to print results as JSON

|

Consensus tree
##############
Function names: consensus_tree; consensus; ctree |br|
Command line interface: pk_consensus_tree; pk_consensus; pk_ctree

Infer a consensus tree from a collection of trees.

Input can be either:
1) a file with one Newick tree per line, or
2) a file with one tree-file path per line.

Consensus methods:
* ``majority``: majority-rule consensus (default) |br|
* ``strict``: strict consensus

Missing taxa handling:
* ``--missing-taxa error`` (default): exits if trees do not share identical tip sets |br|
* ``--missing-taxa shared``: prunes all trees to the intersection of taxa before inferring consensus

.. code-block:: shell

   phykit consensus_tree -t/--trees <trees> [-m/--method strict|majority] [--missing-taxa error|shared] [--json]

Options: |br|
*-t/\\-\\-trees*: file containing trees (one Newick per line) or tree-file paths (one per line) |br|
*-m/\\-\\-method*: consensus method (``strict`` or ``majority``; default: ``majority``) |br|
*--missing-taxa*: handling strategy for mismatched taxa (``error`` or ``shared``; default: ``error``) |br|
*--json*: optional argument to print results as JSON

|

Prune tree
##########
Function names: prune_tree; prune |br|
Command line interface: pk_prune_tree; pk_prune

Prune tips from a phylogeny.

Provide a single column file with the names of the tips
in the input phylogeny you would like to prune from the
tree.

.. code-block:: shell

   phykit prune_tree <tree> <list_of_taxa> [-o/--output <output_file>] [-k/--keep] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*<list_of_taxa>*: single column file with the names of the tips to remove
from the phylogeny |br|
*-o/\\-\\-output*: name of output file for the pruned phylogeny. 
Default output will have the same name as the input file but with the suffix 
".pruned" |br|
*-k/\-\-keep*: optional argument. If used instead of pruning taxa in <list_of_taxa>,
keep them |br|
*--json*: optional argument to print results as JSON
|

Rename tree tips
################
Function names: rename_tree_tips; rename_tree; rename_tips |br|
Command line interface: pk_rename_tree_tips; pk_rename_tree; pk_rename_tips

Renames tips in a phylogeny.

Renaming tip files will follow the scheme of a tab-delimited
file wherein the first column is the current tip name and the
second column is the desired tip name in the resulting 
phylogeny. 

.. code-block:: shell

   phykit rename_tree_tips <tree> -i/--idmap <idmap.txt> [-o/--output <output_file>] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-i/\\-\\-idmap*: identifier map of current tip names (col1) and desired
tip names (col2) |br|
*-o/\\-\\-output*: optional argument to write the renamed tree files to. Default
output will have the same name as the input file but with the suffix ".renamed" |br|
*--json*: optional argument to print results as JSON

|

Robinson-Foulds distance
########################
Function names: robinson_foulds_distance; rf_distance; rf_dist; rf |br|
Command line interface: pk_robinson_foulds_distance; pk_rf_distance; pk_rf_dist; pk_rf

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

.. code-block:: shell

   phykit robinson_foulds_distance <tree_file_zero> <tree_file_one> [--json]

Options: |br|
*<tree_file_zero>*: first argument after function name should be a tree file |br|
*<tree_file_one>*: second argument after function name should be a tree file |br|
*--json*: optional argument to print results as JSON

|

Root tree
#########
Function names: root_tree; root; rt |br|
Command line interface: pk_root_tree; pk_root; pk_rt

Roots phylogeny using user-specified taxa.

A list of taxa to root the phylogeny on should be specified using the -r
argument. The root_taxa file should be a single-column file with taxa names.
The outputted file will have the same name as the inputted tree file but with
the suffix ".rooted".

.. code-block:: shell

   phykit root_tree <tree> -r/--root <root_taxa> [-o/--output <output_file>] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file to root|br|
*-r/\\-\\-root*: single column file with taxa names to root the phylogeny on|br|
*-o/\\-\\-output*: optional argument to specify the name of the output file |br|
*--json*: optional argument to print results as JSON

|

Spurious homolog identification
###############################
Function names: spurious_sequence; spurious_seq; ss |br|
Command line interface: pk_spurious_sequence; pk_spurious_seq; pk_ss

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

.. code-block:: shell

   phykit spurious_seq <file> -f/\\-\\-factor [--json]

Options: |br|
*<file>*: first argument after function name should be a tree file |br|
*-f/\\-\\-factor*: factor to multiply median branch length by to calculate
the threshold of long branches. (Default: 20) |br|
*--json*: optional argument to print results as JSON

|

Terminal branch statistics
##########################
Function names: terminal_branch_stats; tbs |br|
Command line interface: pk_terminal_branch_stats; pk_tbs

Calculate summary statistics for terminal branch lengths in a phylogeny.

Terminal branch lengths can be useful for phylogeny diagnostics.

To obtain all terminal branch lengths, use the -v/\\-\\-verbose option.   

.. code-block:: shell

   phykit terminal_branch_stats <tree> [-v/--verbose] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-v/\\-\\-verbose*: optional argument to print all terminal branch lengths |br|
*--json*: optional argument to print results as JSON

|

Tip labels
##########
Function names: tip_labels; tree_labels; labels; tl |br|
Command line interface: pk_tip_labels; pk_tree_labels; pk_labels; pk_tl

Prints the tip labels (or names) a phylogeny.

.. code-block:: shell

   phykit tip_labels <tree> [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*--json*: optional argument to print results as JSON

|

Tip-to-tip distance
###################
Function names: tip_to_tip_distance; t2t_dist; t2t |br|
Command line interface: pk_tip_to_tip_distance; pk_t2t_dist; pk_t2t

Calculate distance between two tips (or leaves) in a phylogeny.

Distances are in substitutions per site.

.. code-block:: shell

   phykit tip_to_tip_distance <tree_file> <tip_1> <tip_2> [--json]
   phykit tip_to_tip_distance <tree_file> --all-pairs [--plot] [--plot-output <path>] [--json]

Options: |br|
*<tree_file>*: first argument after function name should be a tree file |br|
*<tip_1>*: second argument should be the name of the first tip of interest |br|
*<tip_2>*: third argument should be the name of the second tip of interest |br|
*--all-pairs*: optional argument to report all pairwise tip distances |br|
*--plot*: optional argument to save a clustered distance heatmap (requires ``--all-pairs``) |br|
*--plot-output*: output path for heatmap (default: ``tip_to_tip_distance_heatmap.png``) |br|
*--json*: optional argument to print results as JSON

|

Tip-to-tip node distance
########################
Function names: tip_to_tip_node_distance; t2t_node_dist; t2t_nd |br|
Command line interface: pk_tip_to_tip_node_distance; pk_t2t_node_dist; pk_t2t_nd

Calculate distance between two tips (or leaves) in a phylogeny.

Distance is measured by the number of nodes between one tip
and another.

.. code-block:: shell

   phykit tip_to_tip_node_distance <tree_file> <tip_1> <tip_2> [--json]

Options: |br|
*<tree_file>*: first argument after function name should be a tree file |br|
*<tip_1>*: second argument should be the name of the first tip of interest |br|
*<tip_2>*: third argument should be the name of the second tip of interest |br|
*--json*: optional argument to print results as JSON

|

Total tree length
#################
Function names: total_tree_length; tree_len |br|
Command line interface: pk_total_tree_length; pk_tree_len

Calculate total tree length, which is a sum of all branches.

.. code-block:: shell

   phykit total_tree_length <tree> [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*--json*: optional argument to print results as JSON

|

Treeness
########
Function names: treeness; tness |br|
Command line interface: pk_treeness; pk_tness

Calculate treeness statistic for a phylogeny.

Higher treeness values are thought to be desirable because they
represent a higher signal-to-noise ratio.

Treeness is the sum of internal branch lengths divided by the total
tree length. Therefore, values range from 0 to 1. Treeness can be
used as a measure of the signal-to-noise ratio in a phylogeny. 

Calculate treeness (also referred to as stemminess) following
Lanyon, The Auk (1988), doi: 10.1093/auk/105.3.565 and
Phillips and Penny, Molecular Phylogenetics and Evolution
(2003), doi: 10.1016/S1055-7903(03)00057-5.

.. code-block:: shell

   phykit treeness <tree> [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*--json*: optional argument to print results as JSON

|

Alignment- and tree-based functions
-----------------------------------

Saturation
##########
Function names: saturation; sat |br|
Command line interface: pk_saturation; pk_sat

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

.. code-block:: shell

   phykit saturation -a <alignment> -t <tree> [-v/--verbose] [-e/--exclude_gaps] [--plot] [--plot-output <path>] [--json]

Options: |br|
*-a/\\-\\-alignment*: an alignment file |br|
*-t/\\-\\-tree*: a tree file |br|
*-e/\-\-exclude_gaps*: if a site has a gap, ignore it |br|
*-v/\\-\\-verbose*: print out patristic distances and uncorrected |br|
distances used to determine saturation |br|
*--plot*: save a saturation scatter plot with fitted slope through origin |br|
*--plot-output*: output path for saturation plot (default: ``saturation_plot.png``) |br|
*--json*: optional argument to print results as JSON

Treeness over RCV
#################
Function names: treeness_over_rcv; toverr; tor |br|
Command line interface: pk_treeness_over_rcv; pk_toverr; pk_tor

Calculate treeness/RCV for a given alignment and tree.

Higher treeness/RCV values are thought to be desirable because
they harbor a high signal-to-noise ratio and are least susceptible
to composition bias.

PhyKIT reports three tab delimited values:
col1: treeness/RCV
col2: treeness
col3: RCV

Calculate treeness/RCV following Phillips and Penny, Molecular 
Phylogenetics and Evolution (2003), doi: 10.1016/S1055-7903(03)00057-5.

.. code-block:: shell

   phykit treeness_over_rcv -a/--alignment <alignment> -t/--tree <tree> [--json]

Options: |br|
*-a/\\-\\-alignment*: an alignment file |br|
*-t/\\-\\-tree*: a tree file |br|
*--json*: optional argument to print results as JSON

.. |br| raw:: html

  <br/>
