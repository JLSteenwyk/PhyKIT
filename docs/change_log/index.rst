.. _change_log:


Change log
==========

Major changes to PhyKIT are summarized here.

**2.1.18**:
Added phylogenetic generalized linear models for binary and count data:

* Added new ``phylogenetic_glm`` command (aliases: ``phylo_glm``, ``pglm``)
  for fitting phylogenetic GLMs
* Binomial family: logistic regression via Maximum Penalized Likelihood
  Estimation (logistic_MPLE; Ives & Garland 2010) with Firth's penalty.
  Log-likelihood computed via pruning algorithm for a 2-state CTMC on the
  phylogeny. Fisher information computed via O(n) tree-based three-point
  algorithm.
* Poisson family: Poisson regression via Generalized Estimating Equations
  (poisson_GEE; Paradis & Claude 2002) with overdispersion estimation
* Jointly estimates phylogenetic signal parameter alpha (binomial) or
  overdispersion phi (Poisson)
* JSON output support via ``--json``
* Added new CLI entry points: ``pk_phylogenetic_glm``, ``pk_phylo_glm``,
  ``pk_pglm``
* Results validated against R 4.4.0 (``phylolm::phyloglm()``):

  Poisson GEE (``count_trait ~ body_mass``):

  .. list-table::
     :header-rows: 1
     :widths: 30 20 20 20

     * - Parameter
       - PhyKIT
       - R ``phyloglm``
       - Difference
     * - Intercept
       - 0.6741
       - 0.6741
       - < 1e-4
     * - body_mass
       - 0.5968
       - 0.5968
       - < 1e-4
     * - SE(Intercept)
       - 0.1678
       - 0.1678
       - < 1e-4
     * - SE(body_mass)
       - 0.0877
       - 0.0877
       - < 1e-4
     * - Overdispersion (phi)
       - 0.1730
       - 0.1730
       - < 1e-4

  Poisson GEE matches R to within numerical precision.

  Logistic MPLE (``binary_trait ~ body_mass``):

  .. list-table::
     :header-rows: 1
     :widths: 30 20 20 20

     * - Parameter
       - PhyKIT
       - R ``phyloglm``
       - Difference
     * - Intercept
       - -2.2374
       - -2.1210
       - 0.116
     * - body_mass
       - 2.2432
       - 2.2158
       - 0.027
     * - alpha
       - 0.0215
       - 0.0274
       - 0.006
     * - Log-likelihood
       - -1.833
       - -1.870
       - 0.037
     * - AIC
       - 9.665
       - 9.740
       - 0.075

  Logistic MPLE coefficients agree to within ~5%. Small differences arise
  from how R's ``ape::branching.times()`` computes node heights for
  non-ultrametric trees, which slightly affects the Ives & Garland branch
  length transformation for the Fisher information penalty. Both
  implementations use the same 2-state CTMC pruning log-likelihood
  (verified to match R exactly at -1.870 when evaluated at R's optimal
  parameters).

**2.1.17**:
Unified phylogenetic PCA and dimensionality reduction into a single
``phylogenetic_ordination`` command, and added continuous trait evolution
model comparison:

* Merged ``phylogenetic_pca`` and ``phylogenetic_dimreduce`` into a unified
  ``phylogenetic_ordination`` command (aliases: ``phylo_ordination``,
  ``ordination``, ``ord``) supporting PCA, t-SNE, and UMAP via ``--method``
* All previous aliases remain functional: ``phylo_pca``, ``phyl_pca``,
  ``ppca``, ``phylo_dimreduce``, ``dimreduce``, ``pdr``
* PCA's old ``-m/--method`` (BM/lambda) is now ``--correction``
* GLS-centering via phylogenetic VCV matrix for all methods
* Auto-adjusted parameters for small datasets (t-SNE/UMAP)
* Optional Pagel's lambda correction via ``--correction lambda``
* Optional scatter plot with phylogeny overlay (``--plot``, ``--plot-tree``)
* JSON output support via ``--json``
* Added new CLI entry points: ``pk_phylogenetic_ordination``,
  ``pk_phylo_ordination``, ``pk_ordination``, ``pk_ord``

Added continuous trait evolution model comparison:

* Added new ``fit_continuous`` command (aliases: ``fitcontinuous``, ``fc``)
  for comparing models of continuous trait evolution on a phylogeny,
  analogous to R's ``geiger::fitContinuous()``
* Fits 7 models: BM, OU, EB, Lambda, Delta, Kappa, White
* Ranks models by AIC, BIC, and AIC weights
* Optional subset of models via ``--models``
* JSON output support via ``--json``
* Added new CLI entry points: ``pk_fit_continuous``,
  ``pk_fitcontinuous``, ``pk_fc``

**2.1.16**:
Added visualization commands and rate heterogeneity test:

* Added new ``cont_map`` command (aliases: ``contmap``, ``cmap``)
  for plotting a phylogram with branches colored by continuous trait
  values via ML ancestral reconstruction (analogous to R's
  ``phytools::contMap()``)
* Added new ``density_map`` command (aliases: ``densitymap``, ``dmap``)
  for plotting posterior probabilities of discrete character states
  along phylogeny branches from stochastic character mapping (analogous
  to R's ``phytools::densityMap()``)
* Added new ``phenogram`` command (aliases: ``traitgram``, ``tg``)
  for plotting a phenogram (traitgram) showing trait evolution along a
  phylogeny with X-axis = distance from root, Y-axis = trait value
  (analogous to R's ``phytools::phenogram()``)
* Added new ``cophylo`` command (aliases: ``tanglegram``, ``tangle``)
  for plotting cophylogenetic tanglegrams of two phylogenies with
  connecting lines between matching taxa and node rotation to minimize
  crossings (analogous to R's ``phytools::cophylo()``)

Added rate heterogeneity test:

* Added new ``rate_heterogeneity`` command (aliases: ``brownie``, ``rh``)
  for testing whether continuous trait evolution rates differ across
  phylogenetic regimes using multi-rate Brownian motion (O'Meara et al.
  2006), analogous to R's ``phytools::brownie.lite()``
* Fits single-rate vs. multi-rate BM models and performs a likelihood
  ratio test (chi-squared)
* Optional parametric bootstrap via ``-n/--nsim``
* Regime assignments to internal branches inferred via Fitch parsimony
* Optional ``--plot`` argument for regime-colored phylogram
* JSON output support via ``--json``
* Added new CLI entry points: ``pk_rate_heterogeneity``,
  ``pk_brownie``, ``pk_rh``
* Results validated against R 4.4.0 (``phytools::brownie.lite`` with
  ``paintSubTree(stem=TRUE)``):

  PhyKIT uses Fitch parsimony for regime assignment, which matches
  R's ``paintSubTree(stem=TRUE)`` behavior.

  .. list-table::
     :header-rows: 1
     :widths: 30 20 20 20

     * - Parameter
       - PhyKIT
       - R ``brownie.lite``
       - Difference
     * - Single-rate sigma2
       - 0.03841
       - 0.03841
       - < 1e-11
     * - Single-rate LL
       - -11.56968
       - -11.56968
       - < 1e-14
     * - Single-rate anc. state
       - 1.64469
       - 1.64469
       - < 1e-15
     * - Multi-rate sigma2 (terrestrial)
       - 0.05002
       - 0.04998
       - 3.9e-05
     * - Multi-rate sigma2 (aquatic)
       - 0.00881
       - 0.00889
       - 8.1e-05
     * - Multi-rate LL
       - -11.20459
       - -11.20461
       - 1.6e-05
     * - Chi-squared p-value
       - 0.39283
       - 0.39284
       - 1.1e-05

  Single-rate model matches to machine precision. Multi-rate model
  matches to within optimizer convergence tolerance (both converge to
  the same flat likelihood plateau).

**2.1.15**:
Added ancestral state reconstruction:

* Added new ``ancestral_state_reconstruction`` command (aliases: ``asr``,
  ``anc_recon``) for estimating ancestral states of continuous traits using
  maximum likelihood, analogous to R's ``phytools::fastAnc()`` and
  ``ape::ace(type="ML")``
* Two methods: ``fast`` (two-pass Felsenstein's algorithm, O(n)) and ``ml``
  (full VCV-based ML with exact conditional CIs, O(n^3))
* Optional ``--ci`` flag to include 95% confidence intervals
* Optional ``--plot`` argument to generate a contMap plot showing continuous
  trait values mapped onto the phylogeny with a color gradient
* Supports both two-column single-trait files and multi-trait files with
  ``-c`` flag to select a trait column
* JSON output support via ``--json``
* Added new CLI entry points: ``pk_ancestral_state_reconstruction``,
  ``pk_asr``, ``pk_anc_recon``
* Results validated against R 4.4.0 (``phytools::fastAnc`` with
  ``vars=TRUE, CI=TRUE``):

  Point estimates — ``fast`` method vs R's ``phytools::fastAnc()``

  .. list-table::
     :header-rows: 1
     :widths: 15 30 15 15 15

     * - Node
       - Descendants
       - PhyKIT
       - R ``fastAnc``
       - Error
     * - N1 (root)
       - all 8 tips
       - 1.6446924
       - 1.6446924
       - 0.0000000
     * - N2
       - bear, raccoon
       - 1.7012405
       - 1.7012405
       - 0.0000000
     * - N3
       - 5 taxa
       - 1.4564597
       - 1.4564597
       - 0.0000000
     * - N4
       - sea_lion, seal
       - 1.8090745
       - 1.8090745
       - 0.0000000
     * - N5
       - cat, monkey, weasel
       - 1.2565917
       - 1.2565917
       - 0.0000000
     * - N6
       - cat, monkey
       - 0.9894725
       - 0.9894725
       - 0.0000000

  95% CIs — ``fast`` method vs R's ``fastAnc(CI=TRUE)``

  .. list-table::
     :header-rows: 1
     :widths: 12 22 22 12

     * - Node
       - PhyKIT CI
       - R CI
       - Error
     * - N1 (root)
       - [0.894, 2.396]
       - [0.894, 2.396]
       - 0.000
     * - N2
       - [0.970, 2.433]
       - [0.970, 2.433]
       - 0.000
     * - N3
       - [0.639, 2.274]
       - [0.639, 2.274]
       - 0.000
     * - N4
       - [0.976, 2.642]
       - [0.976, 2.642]
       - 0.000
     * - N5
       - [0.355, 2.158]
       - [0.355, 2.158]
       - 0.000
     * - N6
       - [-0.565, 2.544]
       - [-0.565, 2.544]
       - 0.000

  - ``fast`` and ``ml`` methods produce identical point estimates
    (within 1e-6)
  - Sigma-squared (BM rate) = 0.04389 (matches R's PIC-based estimate
    within 1e-3)

**2.1.13**:
Added stochastic character mapping (SIMMAP):

* Added new ``stochastic_character_map`` command (aliases: ``simmap``, ``scm``)
  for performing Stochastic Character Mapping of discrete traits onto a
  phylogeny (Huelsenbeck et al. 2003; Bollback 2006), analogous to R's
  ``phytools::make.simmap()``
* Fits a continuous-time Markov chain (CTMC) rate matrix Q via maximum
  likelihood using Felsenstein's pruning algorithm
* Three substitution models: ER (equal rates), SYM (symmetric), ARD (all
  rates differ)
* Simulates character histories conditioned on tip states via rejection
  sampling
* Reports mean dwelling times, mean transition counts, and posterior node
  probabilities across simulations
* Optional ``--plot`` argument to generate a horizontal phylogram with
  branches colored by mapped character state
* Reproducible simulations via ``--seed`` argument
* JSON output support via ``--json``
* Results validated against R 4.4.0 (``phytools::fitMk`` and
  ``phytools::make.simmap``):

  - ER log-likelihood: R = -8.7889, PhyKIT = -8.7874 (within 0.002);
    both are valid ML estimates on a flat likelihood surface
  - ARD log-likelihood: R = -8.4305, PhyKIT = -8.3845; PhyKIT's
    multi-start optimizer finds a slightly better local optimum
  - Total tree length conserved exactly (277.2772)
  - Dwelling times sum to total tree length across all simulations
  - Q matrix structural properties verified: rows sum to zero,
    off-diagonal elements positive, diagonal elements negative
  - Model nesting confirmed: ARD loglik >= SYM loglik >= ER loglik
* Added new CLI entry points:
  ``pk_stochastic_character_map``, ``pk_simmap``, ``pk_scm``

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
