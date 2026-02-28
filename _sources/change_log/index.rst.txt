.. _change_log:


Change log
==========

Major changes to PhyKIT are summarized here.

**2.1.12**:
Added phylogenetic regression (PGLS):

* Added new ``phylogenetic_regression`` command (aliases: ``phylo_regression``, ``pgls``)
  for fitting Phylogenetic Generalized Least Squares regression while
  accounting for phylogenetic non-independence among species
* Supports Brownian motion (BM) and Pagel's lambda estimation methods
* Outputs coefficient estimates, standard errors, t-values, p-values,
  R-squared, F-statistic, log-likelihood, and AIC
* Supports multiple predictor variables
* JSON output support via ``--json``
* Results validated against R 4.4.0 (manual GLS with ``ape::vcv()``, matching
  ``caper::pgls()`` behavior); coefficients, standard errors, t-values,
  p-values, R-squared, F-statistic, log-likelihood, and AIC match to at least
  four decimal places
* Uses the raw phylogenetic VCV matrix (not the normalized correlation matrix
  used by ``nlme::gls`` with ``corBrownian``)
* Added new CLI entry points:
  ``pk_phylogenetic_regression``, ``pk_phylo_regression``, ``pk_pgls``

**2.1.11**:
Maintenance:

* Added ``matplotlib>=3.7.0`` as a required dependency

**2.1.10**:
Added phylomorphospace:

* Added new ``phylomorphospace`` command (aliases: ``phylomorpho``, ``phmo``)
  for plotting raw traits with the phylogeny overlaid via ML-reconstructed
  ancestral states at internal nodes
* Tree edges colored by distance from root (coolwarm colormap with colorbar)
* Optional ``--color-by`` for tip point coloring (continuous or discrete)
* Auto-selects first two trait columns when the file has exactly 2 traits
  and ``--trait-x`` / ``--trait-y`` are omitted
* JSON output support via ``--json``
* Added new CLI entry points:
  ``pk_phylomorphospace``, ``pk_phylomorpho``, ``pk_phmo``

**2.1.9**:
Added phylogenetic PCA:

* Added new ``phylogenetic_pca`` command (aliases: ``phylo_pca``, ``phyl_pca``, ``ppca``)
  implementing Revell (2009) phylogenetic PCA for multi-trait data
* Two methods: ``BM`` (Brownian motion) and ``lambda`` (joint Pagel's lambda
  estimation across all traits)
* Two modes: ``cov`` (covariance PCA) and ``corr`` (correlation PCA)
* Tab-delimited multi-trait file input with header row, comment/blank-line support,
  and automatic taxon-mismatch handling (intersection with stderr warnings)
* Optional ``--plot`` argument to generate PCA scatter plot (PC1 vs PC2) with
  taxon labels and variance-explained axes
* JSON output support via ``--json``
* Added new CLI entry points:
  ``pk_phylogenetic_pca``, ``pk_phylo_pca``, ``pk_phyl_pca``, ``pk_ppca``
* Benchmarked against R phytools ``phyl.pca()`` (Revell 2012): eigenvalues,
  loadings, and scores match across all method/mode combinations within 1e-4

**2.1.8**:
Added phylogenetic signal analysis:

* Added new ``phylogenetic_signal`` command (aliases: ``phylo_signal``, ``ps``)
  with support for Blomberg's K and Pagel's lambda
* Blomberg's K includes permutation-based p-value (configurable ``--permutations``)
* Pagel's lambda uses ML optimization with likelihood ratio test p-value
* Tab-delimited trait file input with comment/blank-line support and
  automatic taxon-mismatch handling (intersection with stderr warnings)
* JSON output support via ``--json``
* Added new CLI entry points:
  ``pk_phylogenetic_signal``, ``pk_phylo_signal``, ``pk_ps``
* Validated against R phytools ``phylosig()`` across 95 simulated datasets
  (5-50 tips; pure-birth, coalescent; random, BM, and known-lambda traits):
  Pearson r = 1.0 for K, r > 0.999 for lambda, log-likelihood, and LRT p-value

**2.1.7**:
Added a missing-taxa-aware consensus tree utility:

* Added new ``consensus_tree`` command (aliases: ``consensus``, ``ctree``)
  with strict and majority-rule modes
* Added ``--missing-taxa`` handling:
  ``error`` (default) or ``shared`` (prune all trees to the shared taxa set)
* Added JSON output support for consensus metadata and Newick output
* Added new CLI entry points:
  ``pk_consensus_tree``, ``pk_consensus``, ``pk_ctree``
* Added unit/integration test coverage for parsing, alias dispatch,
  and missing-taxa behavior
* Updated usage documentation and top-level command help text

**2.1.6**:
Expanded plotting support for alignment/tree QC workflows:

* Added optional plotting to:
  ``pairwise_identity``, ``saturation``, ``covarying_evolutionary_rates``,
  ``compositional_bias_per_site``, ``evolutionary_rate_per_site``, and
  ``alignment_entropy``
* Added ``tip_to_tip_distance --all-pairs`` mode with optional clustered
  heatmap output (``--plot``)
* Standardized plotting arguments and JSON metadata:
  ``--plot``, ``--plot-output``, and ``plot_output`` in JSON payloads
* Updated CLI help text and online usage documentation for all newly plotted commands
* Expanded integration test coverage for plotting and JSON+plot behavior and
  reran full regression and documentation build checks
* Removed ``cython`` from runtime dependencies; repository has no Cython build
  pipeline (no ``.pyx``/``cythonize`` usage), so this was unnecessary
* Sunset Python 3.9 and 3.10 support; CI, packaging classifiers, and
  ``python_requires`` now target Python 3.11+

**2.1.5**:
JSON output expansion and harmonization:

* Added ``--json`` support to the remaining CLI commands, including
  ``create_concatenation_matrix`` and ``nearest_neighbor_interchange``
* Standardized JSON metadata key naming for improved consistency across commands
* Added canonical ``rows`` payloads for list-style JSON outputs while preserving
  legacy keys for backward compatibility
* Expanded integration test coverage for JSON payloads and performed full
  regression verification

**2.1.4**:
New alignment utilities, masking support, and composition/RCV correctness updates:

* Added new alignment commands:
  - ``alignment_entropy`` (site entropy reporting)
  - ``occupancy_per_taxon`` (per-taxon valid-site occupancy)
  - ``composition_per_taxon`` (per-taxon symbol composition, excluding invalid symbols)
  - ``mask_alignment`` (column masking by gap fraction, occupancy, and optional entropy)
* Updated saturation slope calculation to use NumPy no-intercept least-squares
  (fit constrained through the origin), replacing sklearn in this code path
* Updated ``rcv`` and ``rcvt`` handling to be case-insensitive and to exclude
  gaps/ambiguous symbols from counts and normalization
* Clarified valid-length normalization in RCV calculations and related docs/help text
* Documentation maintenance updates to improve rendering consistency and remove
  duplicate changelog maintenance burden

**2.1.0**:
Major performance improvements, expanded Python support, and bug fixes:

* **Compatibility:**

  - Added support for Python 3.12 and 3.13
  - Maintains compatibility with Python 3.9, 3.10, and 3.11

* **Performance Optimizations:**

  - Added multiprocessing support for computationally intensive functions (up to 8x faster):
    patristic distances, pairwise identity, polytomy test, LB score, bipartition support stats,
    sum of pairs score, saturation analysis, hidden paralogy check, and covarying evolutionary rates
  - Implemented tree caching with LRU cache to avoid re-parsing files
  - Added NumPy vectorization for alignment operations (5-10x faster)
  - Optimized file I/O with streaming for large concatenation operations
  - Added pickle-based fast tree copying for NNI operations

* **Bug Fixes:**

  - Fixed tree caching side effects that caused tree modifications to persist
  - Fixed spurious sequence detection to correctly use only terminal branches
  - Fixed DNA threader array broadcasting issues
  - Standardized error exit codes to 2
  - Fixed test infrastructure issues
  - Updated saturation slope fitting to use NumPy no-intercept least-squares
    (fit constrained through the origin), replacing sklearn in this code path
  - Updated ``rcv`` and ``rcvt`` to be case-insensitive and to exclude
    gaps/ambiguous symbols from composition counts and normalization; RCV now
    normalizes each taxon by valid (non-excluded) sequence length

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
