.. _change_log:


Change log
==========

Major changes to PhyKIT are summarized here.

**2.1.84**:
Added ``occupancy_filter`` command:

* Added ``occupancy_filter`` (``occ_filter`` / ``filter_occupancy``):
  filter alignments and/or trees by cross-file taxon occupancy. Given
  a list of FASTA or tree files, counts how many files each taxon
  appears in and retains only taxa meeting a minimum threshold. The
  threshold can be a fraction (e.g., ``0.5`` for 50% of files, the
  default) or an absolute count (e.g., ``5``). Outputs filtered copies
  of each input file. For FASTA files, sequences of removed taxa are
  dropped; for tree files, tips are pruned. Supports ``--output-dir``,
  ``--suffix``, and ``--json`` output.

**2.1.83**:
Added ``phylo_path`` command for phylogenetic path analysis:

* Added ``phylo_path`` (``ppath`` / ``phylopath``): phylogenetic path
  analysis following von Hardenberg & Gonzalez-Voyer (2013). Compares
  competing causal DAGs using d-separation tests via PGLS with Pagel's
  lambda, ranks models by CICc, and estimates conditionally model-averaged
  path coefficients. Users define candidate models in a simple text file
  (``name: A->B, B->C``). Supports ``--best-only`` for single-model
  inference, ``--plot-output`` for DAG visualization with path coefficients,
  ``--csv``, and ``--json`` output. Validated against R ``phylopath``
  package (see ``tests/r_validation/validate_phylo_path.R``).

**2.1.82**:
Added ``simmap_summary`` command:

* Added ``simmap_summary`` (``smsummary`` / ``describe_simmap``): run N
  stochastic character maps and provide per-branch summaries of dwelling
  time proportions, expected transitions, and posterior state probabilities
  at each internal node — analogous to ``phytools::describe.simmap()`` in R.
  Supports ``--csv`` for per-branch output, ``--plot`` for posterior pie
  chart tree, and ``--json`` output. Q matrix and log-likelihood validated
  against ``phytools::fitMk()``
  (see ``tests/r_validation/validate_simmap_summary.R``).

**2.1.81**:
Added ``phylo_anova`` command for phylogenetic ANOVA / MANOVA:

* Added ``phylo_anova`` (``panova`` / ``phylo_manova`` / ``pmanova``):
  phylogenetic ANOVA (univariate) or MANOVA (multivariate) using the
  Residual Randomization Permutation Procedure (RRPP; Adams & Collyer
  2018). Auto-detects univariate vs multivariate from the data;
  override with ``--method anova`` or ``--method manova``. Supports
  ``--pairwise`` post-hoc comparisons, ``--plot-output`` for violin+boxplot
  (ANOVA) or phylomorphospace (MANOVA) visualizations, ``--seed`` for
  reproducibility, and ``--json`` output. Deterministic values (SS, MS, F,
  Pillai's trace) validated against R (see
  ``tests/r_validation/validate_phylo_anova.R``).

**2.1.80**:
Added ``--branch-labels`` option to ``quartet_pie``:

* Added ``--branch-labels`` flag to ``quartet_pie`` (``qpie``): displays
  the number of concordant gene trees (blue, above branch) and the
  local posterior probability (LPP; red, below branch) on each internal
  branch, in the style of PhyTop. In ASTRAL mode, values are parsed
  from ``f1`` and ``pp1`` annotations; in native mode, the concordant
  gene count is computed directly from gene trees.

**2.1.79**:
Added ``subtree_prune_regraft`` command:

* Added ``subtree_prune_regraft`` (``spr``): generate all possible SPR
  (Subtree Pruning and Regrafting) rearrangements for a user-specified
  subtree on a parent tree. The subtree is defined by one or more
  comma-separated taxa (resolved by MRCA). Each output Newick tree
  represents the result of pruning the subtree and regrafting it onto
  a different branch. Supports ``--json`` output and ``-o`` output file.

**2.1.78**:
Added ``taxon_groups`` utility command:

* Added ``taxon_groups`` (``tgroups`` / ``shared_taxa``): determine
  which tree or FASTA files share the same set of taxa. Groups files
  by identical taxon sets, reports groups sorted by size with the
  taxa present in each group. Useful for identifying subsets of genes
  with identical taxon sampling for concatenation or comparative
  analysis.

**2.1.77**:
Added ``hybridization`` command for reticulation analysis:

* Added ``hybridization`` (``hybrid`` / ``reticulation``): estimate
  the minimum number of reticulation events and localize where
  hybridization likely occurred on a species tree. For each branch,
  computes a binomial test for asymmetric discordance (FDR-corrected),
  a hybridization score, and identifies which pair of lineages
  likely exchanged genes (but not the direction of flow).
  Supports ``--support`` threshold for
  collapsing low-confidence gene tree branches before analysis.
  Visualizes as a branch-colored phylogram (gray = no signal,
  red = strong hybridization signal, stars at significant nodes).

**2.1.76**:
Added ``neighbor_net`` command (Bryant & Moulton 2004):

* Added ``neighbor_net`` (``nnet``): construct a NeighborNet
  phylogenetic network from pairwise distances. Infers a splits
  graph from a FASTA alignment or pre-computed distance matrix
  using NJ-based circular ordering + NNLS split weight estimation.
  Supports p-distance, identity, and Jukes-Cantor distance metrics.
  Visualizes as a planar Buneman splits graph.

**2.1.75**:
Improved ``consensus_network`` for large datasets and plot customization:

* Added ``--max-splits N`` (default 30) to cap the number of splits
  used in the Buneman graph visualization, avoiding exponential
  blowup with many splits. All splits are still reported in text/JSON.
* Added ``--histogram <file>`` to output a split frequency distribution
  plot showing how many splits occur at each frequency level.
* Added warning for datasets with >100 taxa recommending ``--histogram``
  over the network graph, which doesn't scale well to large trees.
* Non-boundary taxa (not involved in displayed splits) now draw short,
  faded pendant edges to reduce visual clutter in the network graph.
* ``--ylabel-fontsize`` now controls tip label size in all phylogram
  services; ``--ylabel-fontsize 0`` hides labels entirely. Fixed across
  all 11 services that previously hardcoded fontsize=9.
* ``discordance_asymmetry`` annotations are now opt-in with
  ``--annotate`` (values shown as numbers without "gCF=" prefix).
  ``--legend-position none`` now hides the colorbar too.

**2.1.73**:
Added ``--missing-taxa allow`` mode to ``consensus_network``:

* ``consensus_network`` now supports ``--missing-taxa allow`` (the new
  default), which handles gene trees with different taxon sets by
  extracting splits from each tree using its own taxon set. This is
  the standard approach for phylogenomic datasets with incomplete
  taxon sampling. Previously, ``--missing-taxa shared`` required all
  taxa to be present in all trees, which fails when no taxon is
  universal across all gene trees.

**2.1.72**:
Added ``phylo_impute`` command for missing data imputation:

* Added ``phylo_impute`` (``impute`` / ``phylo_imp``): impute missing
  continuous trait values using phylogenetic relationships and
  between-trait correlations under Brownian motion. For each missing
  value, computes the conditional expectation from (1) phylogenetic
  neighbors' observed values and (2) the taxon's own observed traits
  via trait covariance. Reports imputed values with standard errors
  and 95% CIs. Supports discordance-aware VCV via ``-g``.
  Cross-compared against R's Rphylopars::phylopars(); imputed values
  agree in direction and magnitude. R validation script provided in
  ``tests/r_validation/validate_phylo_impute.R``.

**2.1.71**:
Added ``phylo_gwas`` command (Pease et al. 2016):

* Added ``phylo_gwas`` (``pgwas``): phylogenetic genome-wide
  association study. Tests each alignment site for association with
  a categorical (Fisher's exact) or continuous (point-biserial)
  phenotype, with Benjamini-Hochberg FDR correction. Optionally
  classifies significant hits as monophyletic (potentially inherited)
  or polyphyletic (convergent, the interesting candidates) using a
  phylogenetic tree. Produces a Manhattan plot with color-coded
  significance, optional gene partition annotations, and CSV/JSON
  output of per-site results. Supports ``--exclude-monophyletic``
  to filter inherited associations. ``--dot-size`` controls Manhattan
  plot dot size (same scale-factor pattern as ``--pie-size``).

**2.1.70**:
Added multivariate K_mult to ``phylogenetic_signal`` (Adams 2014):

* ``phylogenetic_signal`` now supports ``--multivariate`` to compute
  K_mult, the generalized K statistic for multivariate trait data
  (Adams 2014). Uses distance-based formulation that works even when
  the number of traits exceeds the number of taxa. Significance
  assessed via permutation test (shuffling species, preserving trait
  correlations). Cross-validated against R's geomorph::physignal():
  K_mult matches exactly (0.563803); p-values differ only due to
  permutation stochasticity. R validation script provided in
  ``tests/r_validation/validate_kmult.R``.

**2.1.69**:
Added ``phylo_logistic`` command (Ives & Garland 2010):

* Added ``phylo_logistic`` (``phylo_logreg`` / ``plogreg``):
  phylogenetic logistic regression for binary (0/1) response
  variables using maximum penalized likelihood estimation (MPLE)
  with Firth's bias-correction penalty. Estimates regression
  coefficients, phylogenetic correlation parameter (alpha),
  standard errors, Z-scores, p-values, log-likelihood, and AIC.
  Jointly optimizes beta and alpha via L-BFGS-B. Cross-validated
  against R's phylolm::phyloglm() on a 50-taxon simulated dataset:
  predicted probabilities correlate at r = 0.9999, intercept
  within 2%, slope within 6%, log-likelihood within 0.1%.
  R validation script provided in
  ``tests/r_validation/validate_phylo_logistic_50taxa.R``.

**2.1.68**:
Added CI visualization to ``ancestral_state_reconstruction``:

* ``ancestral_state_reconstruction`` now supports ``--plot-ci`` to
  draw confidence interval bars at each internal node on the contMap
  phylogram. Shows a vertical bar with caps spanning the 95% CI and
  a dot for the point estimate. Requires ``--ci`` and ``--plot``.
* ``--ci-size`` scale factor controls bar size (default 1.0).
  Works in both rectangular and circular modes.

**2.1.67**:
Added ``dfoil`` command (Pease & Hahn 2015):

* Added ``dfoil`` (``dfoil_test``): compute the four DFOIL statistics
  (DFO, DIL, DFI, DOL) for detecting and polarizing introgression
  in a 5-taxon symmetric phylogeny ``((P1, P2), (P3, P4), Outgroup)``.
  Counts 16 binary site patterns, computes four D-statistics with
  chi-squared significance, and interprets the joint sign pattern
  to identify which lineages exchanged genes and the direction of
  gene flow. Based on Pease & Hahn (*Systematic Biology*, 2015).

**2.1.66**:
Added ``--pie-size`` to ``quartet_pie``:

* ``quartet_pie`` now supports ``--pie-size <float>`` to scale pie
  chart sizes relative to the default (1.0 = default, 2.0 = double,
  0.5 = half). Works in both rectangular and circular modes.

**2.1.65**:
Added gene-tree mode and support filtering to ``dstatistic`` (ABBA-BABA):

* ``dstatistic`` now supports ``-g/--gene-trees`` as an alternative
  to ``-a/--alignment``. Gene trees can have any number of taxa;
  the induced quartet for the four specified taxa is extracted from
  each tree. Significance via chi-squared test.
* Added ``--support`` threshold for gene-tree mode: branches with
  support values below the threshold are collapsed (treated as
  unresolved), accounting for uncertainty in gene tree topology.
  For example, ``--support 70`` ignores branches with bootstrap
  support below 70%.

**2.1.64**:
Added ``dstatistic`` (ABBA-BABA) command:

* Added ``dstatistic`` (``dstat`` / ``abba_baba``): compute Patterson's
  D-statistic for detecting introgression. Two modes:

  - **Site-pattern mode** (``-a``): counts ABBA/BABA from a FASTA
    alignment with block jackknife significance testing
  - **Gene-tree mode** (``-g``): counts discordant quartet topologies
    from gene trees (any number of taxa) with chi-squared significance.
    Extracts the induced quartet for the four specified taxa from each
    multi-taxon gene tree.

**2.1.63**:
Added ``trait_rate_map`` command:

* Added ``trait_rate_map`` (``rate_map`` / ``branch_rates``):
  estimate per-branch evolutionary rates for a continuous trait
  using squared standardized contrasts and display as a branch-
  colored phylogram. Branches with faster evolution appear in
  warmer colors. Supports rectangular and circular layouts,
  cladogram mode, and color file annotations.

**2.1.62**:
Added ``trait_correlation`` command (with dendrogram fix).

* Added ``trait_correlation`` (``trait_corr`` / ``phylo_corr``):
  compute phylogenetic correlations between all pairs of continuous
  traits using GLS-centered covariance, and display as a heatmap
  with significance stars (``*`` p<0.05, ``**`` p<0.01, ``***``
  p<0.001). Supports optional hierarchical clustering of traits
  (``--cluster``), custom significance threshold (``--alpha``), and
  discordance-aware VCV via gene trees (``-g``). Cross-validates
  against pairwise PGLS.

**2.1.60**:
Added ``--cluster-columns`` to ``phylo_heatmap`` and partition support
to ``identity_matrix``:

* ``phylo_heatmap`` now supports ``--cluster-columns`` to cluster trait
  columns by similarity and display a dendrogram at the top. Column
  labels move to the bottom when the dendrogram is shown.
* ``identity_matrix`` now supports ``--partition <file>`` with a
  RAxML-style partition file. A per-gene identity panel is displayed
  alongside the main heatmap.

**2.1.58**:
Added ``identity_matrix`` command:

* Added ``identity_matrix`` (``id_matrix`` / ``seqid``): compute
  pairwise sequence identity (or p-distance) matrix from an alignment
  and plot as a clustered heatmap with dendrograms. Supports three
  ordering modes: hierarchical clustering (default), tree-guided
  (``--sort tree --tree <file>``), or alphabetical (``--sort alpha``).
  Partition file support planned for a future release.

**2.1.57**:
Added ``--heatmap`` and ``--distance-matrix`` options to ``tree_space``:

* ``--heatmap``: draw a clustered distance heatmap (with dendrogram)
  instead of the default MDS/t-SNE/UMAP scatter plot
* ``--distance-matrix <file>``: export the raw pairwise distance
  matrix as a CSV file (works with both scatter and heatmap modes)

**2.1.56**:
Added ``tree_space`` command for gene tree topology visualization:

* Added ``tree_space`` (``tspace`` / ``tree_landscape``): visualize
  gene tree topology space via MDS, t-SNE, or UMAP on pairwise
  Robinson-Foulds or Kuhner-Felsenstein distance matrices. Includes
  spectral clustering with eigengap auto-detection and optional
  species tree highlighting as a distinct marker.
* ``phylogenetic_ordination`` already supports ``--color-by`` for
  coloring points by a continuous trait column or a categorical
  group file (taxon-to-group TSV). This complements the auto-
  clustering in ``tree_space``.

**2.1.55**:
Added ``--csv`` output to ``quartet_pie``:

* ``quartet_pie`` now supports ``--csv <file>`` to output per-branch
  concordance values (gCF, gDF1, gDF2, concordant/discordant counts)
  as a CSV table. Works in both native (gene tree) and ASTRAL modes.

**2.1.54**:
Added ``--color-file`` plot option and ``alignment_subsample`` command:

* Added ``alignment_subsample`` (``aln_subsample`` / ``subsample``):
  randomly subsample genes, partitions, or sites from phylogenomic
  datasets. Three modes: ``genes`` (subsample gene list), ``partitions``
  (subsample from supermatrix + partition file), ``sites`` (subsample
  alignment columns). Supports bootstrap resampling (``--bootstrap``),
  exact count (``--number``) or fraction (``--fraction``), and
  reproducibility via ``--seed``.
* Added ``--color-file`` plot option for clade and tip label coloring:

* ``--color-file``: iTOL-inspired TSV annotation file for coloring
  tip labels (``label`` type), highlighting clades with transparent
  background bands (``range`` type), and coloring clade branches
  (``clade`` type). Clades are defined by listing taxa whose MRCA
  determines the clade. Works with all phylogram-drawing commands
  in both rectangular and circular modes; ``clade`` branch coloring
  is silently skipped for trait-colored commands. Labeled ranges
  and clades appear in the figure legend.
* Added a plot customization tutorial to the documentation covering
  all shared plot options with example figures.

**2.1.52**:
Added ``--circular`` plot option:

* ``--circular``: draw circular (radial/fan) phylograms with the root
  at the center, branches radiating outward, and tips around the
  perimeter. Curved arcs connect sister clades (FigTree/iTOL style).
  Combinable with ``--cladogram`` and ``--ladderize``. Supported in
  all 11 phylogram-drawing commands.

**2.1.51**:
Added ``--cladogram`` plot option to all phylogram-drawing commands:

* When set, trees are drawn with equal branch lengths and tips aligned
  at the right edge (topological depth layout), instead of the default
  phylogram layout where branch lengths are proportional to evolutionary
  distance
* Supported in all 11 phylogram-drawing commands: ``phylo_heatmap``,
  ``cont_map``, ``density_map``, ``stochastic_character_map``,
  ``quartet_pie``, ``ancestral_state_reconstruction``,
  ``concordance_asr``, ``discordance_asymmetry``,
  ``rate_heterogeneity``, ``cophylo``, and ``character_map``
* Shared utility ``compute_node_x_cladogram`` ensures consistent
  cladogram layout across all commands

**2.1.50**:
Added character map command and polytomy handling for network commands.

**2.1.49**:
Added wASTRAL ``--support 3`` compatibility and ``--ladderize`` plot option:

* ``quartet_pie`` (``qpie``) now explicitly supports wASTRAL
  ``--mode 1 -R --support 3`` output trees in addition to standard
  ASTRAL ``-t 2`` output; the q1/q2/q3 annotations are parsed
  automatically from the extended wASTRAL node label format
  (CULength, f1–f3, localPP, pp1–pp3, q1, q2, q3)
* Added ``--ladderize`` flag to all plot-generating commands; when
  set, the tree is ladderized (sorted by number of descendant tips)
  before rendering, producing a cleaner visual layout
* Supported in all 10 phylogram-drawing commands: ``phylo_heatmap``,
  ``cont_map``, ``density_map``, ``stochastic_character_map``,
  ``quartet_pie``, ``ancestral_state_reconstruction``,
  ``concordance_asr``, ``discordance_asymmetry``,
  ``rate_heterogeneity``, and ``cophylo``
* Added polytomy (collapsed branch) handling to hybridization and
  network analysis commands; gene trees with collapsed low-support
  branches can now be used directly as input:

  - ``quartet_network``, ``consensus_network``, ``spectral_discordance``:
    bipartitions from polytomous nodes are excluded (treated as
    uninformative), so quartets spanning a polytomy are classified as
    unresolved rather than misclassified
  - ``network_signal``: polytomies are represented as star topologies
    in the network VCV, correctly modeling unresolved relationships
  - Trifurcating roots (standard unrooted Newick) are not affected
* Added character map command (``character_map`` / ``charmap`` /
  ``synapomorphy_map``): maps synapomorphies and homoplasies onto a
  phylogeny using Fitch parsimony with ACCTRAN or DELTRAN optimization

  - Color-coded circles on branches: blue (synapomorphy), red
    (convergence), gray (reversal)
  - Supports cladogram (default) and phylogram layouts
  - Reports consistency index (CI) and retention index (RI)
  - Optional ``--characters`` filter to display specific characters
  - Cross-validated against R's phangorn package (CI, RI, tree length
    match exactly)

**2.1.47**:
Added Fitch parsimony score command (``parsimony_score`` / ``pars``):

* Computes the minimum number of character state changes required to
  explain an alignment on a given tree topology (Fitch 1971)
* Uses the Fitch downpass algorithm, scoring each site independently
* Gap characters (-, N, X, ?) treated as wildcards
* Automatically resolves multifurcations
* Optional ``-v/--verbose`` flag prints per-site parsimony scores
* Supports ``--json`` output
* Cross-validated against R's phangorn::parsimony(method="fitch");
  exact match (score=4). R validation script provided in
  ``tests/r_validation/validate_parsimony.R``

**2.1.46**:
Added phylogenetically independent contrasts command
(``independent_contrasts`` / ``pic``):

* Computes Felsenstein's (1985) phylogenetically independent contrasts
  for continuous traits on a phylogeny
* Produces n-1 standardized contrasts for n tips via postorder traversal
* Automatically resolves multifurcations by adding zero-length branches
* Reports individual contrasts with associated tip groups, mean absolute
  contrast, and variance of contrasts
* Supports ``--json`` output with per-node contrast values
* Cross-validated against R's ape::pic(); sum of squared contrasts
  matches R exactly (0.307253). R validation script provided in
  ``tests/r_validation/validate_pic.R``

**2.1.45**:
Added phylogenetic heatmap command (``phylo_heatmap`` / ``pheatmap`` /
``ph``):

* Draws a phylogeny alongside a color-coded heatmap of numeric trait
  values, with rows aligned to tree tips
* Analogous to R's phytools::phylo.heatmap()
* Input: species tree + multi-column numeric TSV (header with trait
  names, one row per taxon)
* ``--split`` controls tree vs heatmap width ratio (default: 0.3)
* ``--standardize`` z-scores each column before coloring
* ``--cmap`` selects any matplotlib colormap (default: viridis)
* Supports all shared plot options and ``--json`` output
* Supports ``.png``, ``.pdf``, ``.svg`` output

**2.1.42**:
Added quartet pie chart visualization command (``quartet_pie`` / ``qpie``):

* Draws a phylogram with pie charts at internal nodes showing gene
  concordance (gCF) and discordance (gDF1, gDF2) proportions
* Native mode: computes quartet proportions from species tree +
  gene trees via bipartition matching (four-group decomposition)
* ASTRAL mode: parses q1/q2/q3 annotations from ASTRAL ``-t 2``
  output, supporting multiple annotation formats
* Optional ``--annotate`` flag adds numeric values near each pie
* Default colors: blue (concordant), red (discordant alt 1),
  gray (discordant alt 2); overridable via ``--colors``
* Supports all shared plot options (``--fig-width``, ``--dpi``,
  ``--no-title``, etc.) and ``--json`` output with per-node
  concordance counts
* gCF/gDF values validated against manual bipartition matching
  computation on sample data

**2.1.41**:
Added discrete trait model comparison command (``fit_discrete`` / ``fd``):

* Compares ER (Equal Rates), SYM (Symmetric), and ARD (All Rates
  Different) Mk models of discrete character evolution via maximum
  likelihood
* Reports log-likelihood, AIC, delta-AIC, Akaike weights, BIC, and
  number of parameters for each model
* Extracts shared Q-matrix fitting, Felsenstein pruning, and trait
  parsing code from ``stochastic_character_map`` and
  ``ancestral_reconstruction`` into ``phykit/helpers/discrete_models.py``,
  eliminating code duplication
* Supports ``--models`` flag to select a subset of models (e.g.,
  ``--models ER,ARD``)
* Supports ``--json`` output with full Q-matrix and rate parameters
* Cross-validated against R's geiger::fitDiscrete(); R validation
  script provided in ``tests/r_validation/validate_fit_discrete.R``

**2.1.40**:
Added Kuhner-Felsenstein (branch score) distance command
(``kf_distance`` / ``kf``):

* Computes the KF distance between two phylogenies, incorporating
  both topology and branch length differences (Kuhner & Felsenstein
  1994)
* Reports plain and normalized KF distance
* Includes both internal and terminal branch lengths in the
  computation, matching the standard definition
* Prunes to shared taxa when input trees have different tip sets
* Supports ``--json`` output
* Cross-validated against R's phangorn::KF.dist(); R validation
  script provided in ``tests/r_validation/validate_kf_distance.R``

**2.1.39**:
Added shared plot configuration system with user-customizable CLI
arguments for all 27 plotting commands:

* New ``PlotConfig`` system (``phykit/helpers/plot_config.py``) provides
  auto-scaling figure dimensions and font sizes based on dataset size
* All plotting commands now accept ``--fig-width``, ``--fig-height``,
  ``--dpi``, ``--no-title``, ``--title``, ``--legend-position``,
  ``--ylabel-fontsize``, ``--xlabel-fontsize``, ``--title-fontsize``,
  ``--axis-fontsize``, and ``--colors`` arguments
* Figure height and font sizes auto-scale for large datasets — labels
  shrink for 50-550 taxa and auto-hide beyond 800 taxa
* Output format determined by file extension: ``.png``, ``.pdf``,
  ``.svg``, and ``.jpg`` are all supported via ``--plot-output``
* Custom colors can partially override defaults using comma-separated
  values (e.g., ``--colors ",,#e41a1c"`` to change only the third color)
* Updated CLI help text and Sphinx documentation for all commands

**2.1.31**:
Added evolutionary tempo mapping command (``evo_tempo_map`` / ``etm``) for
detecting rate-topology associations in phylogenomic datasets:

* Compares branch length distributions between concordant and discordant
  gene trees at each species tree branch via bipartition matching
* Tests for differences using Mann-Whitney U and permutation tests with
  Benjamini-Hochberg FDR correction across branches
* Reports global treeness (internal/total branch length ratio) comparison
  between concordant and discordant gene trees
* Optional ``--plot`` flag generates a grouped box/strip plot showing
  concordant vs. discordant branch lengths per species tree branch
* Supports ``--json`` and ``-v/--verbose`` output modes
* Cross-validated against existing PhyKIT tools: bipartition matching
  matches ``concordance_asr`` gCF, treeness matches ``treeness`` command,
  FDR correction matches ``relative_rate_test``
* Updated tutorial 19 (gene tree discordance pipeline) with a new
  evolutionary tempo mapping step including expected output and figures
* Enhanced command reference documentation with example output and plot

**2.1.30**:
Added phylogenetic effect size (R² variance decomposition) to all
phylogenetic comparative method commands:

* ``phylogenetic_signal``: reports R²_phylo = 1 - (σ²_BM / σ²_WN),
  the fraction of trait variance explained by phylogenetic structure
* ``phylogenetic_regression`` (PGLS): reports three-way decomposition:
  R²_total (phylo + predictor), R²_pred (predictor given phylogeny),
  and R²_phylo (phylogeny's unique contribution)
* ``phylogenetic_glm``: reports McFadden's pseudo-R² for both Poisson
  GEE and Logistic MPLE, computed from full vs. intercept-only model
  log-likelihoods
* ``fit_continuous``: reports per-model R² = 1 - (σ²_model / σ²_White),
  measuring how much each evolutionary model reduces unexplained variance
  compared to white noise
* ``rate_heterogeneity``: reports R²_regime, the variance reduction from
  regime-specific rates vs. a single rate, weighted by tips per regime
* ``ouwie``: reports per-model R² = 1 - (σ²_model / σ²_BM1), measuring
  improvement over the simplest Brownian motion baseline
* All effect sizes appear in both text and JSON output
* R validation scripts provided in ``tests/r_validation/``

**2.1.29**:
Added discordance-aware VCV matrix support for phylogenetic comparative methods:

* When gene trees are provided via ``-g/--gene-trees``, all five phylogenetic
  comparative method commands now compute a genome-wide average VCV matrix from
  per-gene-tree VCVs instead of using the species tree alone
* This accounts for incomplete lineage sorting (ILS) and introgression, giving
  more accurate covariance estimates for downstream analyses
* Affected commands: ``phylogenetic_signal``, ``phylogenetic_regression`` (PGLS),
  ``fit_continuous``, ``phylogenetic_ordination``, and ``phylogenetic_glm``
* Algorithm: parse gene trees from a multi-Newick file, prune to shared taxa,
  build a VCV matrix from each gene tree, average them, and correct to nearest
  positive semi-definite matrix via eigenvalue clipping
* Auto-prunes gene trees to the intersection of taxa shared across the species
  tree and all gene trees; errors if fewer than 3 shared taxa
* JSON output includes ``vcv_metadata`` with number of gene trees used, number
  of shared taxa, and whether PSD correction was applied
* When no gene trees are provided, behavior is unchanged (full backward
  compatibility)
* Consolidated duplicated ``_build_vcv_matrix`` code from 5 service files into
  a shared ``vcv_utils`` module

**2.1.28**:
Added Felsenstein (2012) threshold model for trait correlation via MCMC:

* Added new ``threshold_model`` command (aliases: ``threshold``, ``thresh``,
  ``threshbayes``, ``thresh_bayes``) for estimating evolutionary correlations
  between binary discrete and/or continuous traits using a latent-liability
  Brownian motion model
* Implements the Gibbs/Metropolis-Hastings MCMC sampler following
  phytools::threshBayes (Revell 2014): binary characters are modelled as
  continuous liabilities crossing a threshold at 0
* Supports all trait combinations: discrete+continuous, discrete+discrete,
  and continuous+continuous
* Adaptive proposal tuning during burn-in targeting ~23% acceptance rate
* Output includes posterior summary (mean, median, 95% HPD) for the
  correlation (r), rate parameters (sigma2), and ancestral values
* Optional ``--plot`` generates a 3x2 figure: trace plots (left column)
  and posterior density histograms with 95% HPD shading (right column)
  for convergence diagnostics and posterior visualization
* JSON output (``--json``) with full posterior samples for custom analysis
* Reproducible results via ``--seed`` for random number generation
* Cross-validated against R's phytools::threshBayes: posterior means and
  95% intervals agree within expected MCMC sampling variation

**2.1.27**:
Added lineage-through-time (LTT) plot and Pybus & Harvey gamma statistic:

* Added new ``ltt`` command (aliases: ``gamma_stat``, ``gamma``)
  for testing temporal variation in diversification rates
* Implements the Pybus & Harvey (2000) gamma statistic: under
  constant-rate pure-birth, gamma ~ N(0,1); negative = decelerating
  diversification, positive = accelerating
* Optional ``--plot-output`` generates a step-function LTT plot with
  log-scaled y-axis showing lineage accumulation through time
* Verbose mode (``-v``) prints branching times and full LTT data
* JSON output support via ``--json``
* Validated against R's ``ape::gammaStat()`` (ape v5.8.1, R 4.4.0):
  gamma values match to 10 decimal places across 4 test topologies
  (balanced 8-tip: -1.4142135624, ladder 5-tip: -0.7142857143,
  recent burst 10-tip: 2.2824790785, early burst 7-tip: -3.5362021857)

**2.1.26**:
Added phylogenetic signal on networks:

* Added new ``network_signal`` command (aliases: ``netsig``, ``net_signal``)
  for computing Blomberg's K and Pagel's lambda on phylogenetic networks
  using the Bastide et al. (2018) variance-covariance algorithm
* Accounts for hybridization/introgression when estimating how strongly a
  continuous trait tracks evolutionary history
* Accepts explicit hybrid edge specifications (``--hybrid donor:recipient:gamma``)
  or auto-infers from ``quartet_network`` JSON output (``--quartet-json``)
* Blomberg's K on phylogenetic networks is a novel capability not available
  in any other tool; Pagel's lambda on networks was previously only
  available in Julia (PhyloNetworks.jl)
* K and lambda formulas are identical to ``phylogenetic_signal``; only the
  VCV matrix differs (network VCV vs tree VCV)

**2.1.25**:
Added Tajima's relative rate test:

* Added new ``relative_rate_test`` command (aliases: ``rrt``, ``tajima_rrt``)
  for testing whether two lineages evolved at equal rates
* Implements Tajima's (1993) chi-squared test on unique substitution counts
  (m1/m2) relative to an outgroup
* Single alignment mode (``-a``) and batch mode (``-l``) for multi-gene analysis
* Automatic outgroup inference from rooted tree
* Bonferroni and Benjamini-Hochberg FDR multiple testing correction
* JSON output support via ``--json``
* Validated against R's ``pegas::rr.test()``: chi-squared statistics and
  p-values match to machine precision (< 1e-8 difference)

**2.1.24**:
Added quartet-based network inference (NANUQ-style) for distinguishing ILS
from hybridization:

* Added new ``quartet_network`` command (aliases: ``quartet_net``, ``qnet``,
  ``nanuq``) for computing quartet concordance factors from gene trees and
  classifying each quartet as tree-like, hybrid, or unresolved
* Implements the NANUQ algorithm (Allman, Baños & Rhodes 2019): star test
  (Pearson chi-squared against uniform 1/3) followed by T3 tree model test
  (G-test / likelihood ratio with conservative chi-squared df=1 p-value)
* Separate ``--alpha`` (tree test threshold, default 0.05) and ``--beta``
  (star test threshold, default 0.95) parameters matching MSCquartets
* ``--plot-output`` option to generate a species tree with reticulation
  arcs overlaid for hybrid quartets
* ``--missing-taxa shared`` support for trees with different taxon sets
* JSON output support via ``--json``
* Added new CLI entry points: ``pk_quartet_network``, ``pk_quartet_net``,
  ``pk_qnet``, ``pk_nanuq``
* Validated against R's MSCquartets v3.2 ``NANUQ()`` function:

  - Star test (p_star) p-values match R exactly
  - T3 tree test (p_tree) p-values match R exactly for large samples;
    small-sample values are slightly conservative (e.g., 0.096 vs 0.214
    for counts 8,0,2) but yield identical classifications
  - All 15 quartets from the sample gene tree file classified identically
    to R's NANUQ (100% agreement)

  .. list-table::
     :header-rows: 1
     :widths: 20 15 15 15 15

     * - Counts
       - p_star (PK)
       - p_star (R)
       - p_tree (PK)
       - p_T3 (R)
     * - (70, 15, 15)
       - 0.0000
       - 0.0000
       - 1.0000
       - 1.000
     * - (45, 35, 20)
       - 0.0087
       - 0.0087
       - 0.0418
       - 0.042
     * - (10, 0, 0)
       - 0.0000
       - 0.0000
       - 1.0000
       - 1.000
     * - (8, 0, 2)
       - 0.0055
       - 0.0055
       - 0.0959
       - 0.214

  Classifications agree in all cases. The p_tree difference for small
  samples (8,0,2) is due to MSCquartets using a specialized T3 density
  integration for the p-value; the conservative chi-squared(df=1)
  approach is a well-established approximation that improves with sample
  size.

**2.1.22**:
Added consensus splits network for visualizing conflicting phylogenetic signal:

* Added new ``consensus_network`` command (aliases: ``consnet``, ``splitnet``,
  ``splits_network``) for extracting bipartition splits from gene trees and
  summarizing conflicting phylogenetic signal
* Counts frequency of each non-trivial bipartition across input trees
* ``--threshold`` option to filter splits by minimum frequency (default: 0.1)
* ``--plot-output`` option to generate a circular splits network diagram
* ``--missing-taxa shared`` support for trees with different taxon sets
* JSON output support via ``--json``
* Added new CLI entry points: ``pk_consensus_network``, ``pk_consnet``,
  ``pk_splitnet``, ``pk_splits_network``

**2.1.21**:
Added automatic OU shift detection (l1ou):

* Added new ``ou_shift_detection`` command (aliases: ``ou_shifts``, ``l1ou``,
  ``detect_shifts``) for automatic detection of adaptive optimum shifts on
  a phylogeny using the LASSO-based approach of Khabbazian et al. (2016)
* No regime file needed — only a tree and trait data
* Model selection via pBIC (default), BIC, or AICc
* ``--max-shifts`` option to limit number of candidate shifts
* JSON output support via ``--json``
* Added new CLI entry points: ``pk_ou_shift_detection``, ``pk_ou_shifts``,
  ``pk_l1ou``, ``pk_detect_shifts``
* Validated against R's l1ou package: same shift count, alpha, and pBIC
  on a 100-tip Anolis dataset

**2.1.20**:
Added taxon occupancy threshold to ``create_concatenation_matrix``:

* New ``--threshold`` option (default 0) excludes taxa whose effective
  representation (fraction of informative, non-gap/non-ambiguous characters)
  falls below the specified value
* Filtering is disabled by default; set ``--threshold 0.5`` (for example)
  to exclude poorly represented taxa
* Excluded taxa are reported to stderr with their effective occupancy
* JSON output includes ``threshold`` and ``excluded_taxa`` fields

**2.1.19**:
Added multi-regime Ornstein-Uhlenbeck models (OUwie):

* Added new ``ouwie`` command (aliases: ``fit_ouwie``, ``multi_regime_ou``)
  for fitting multi-regime OU models of continuous trait evolution
* Seven models: BM1, BMS, OU1, OUM, OUMV, OUMA, OUMVA
  (Beaulieu et al. 2012)
* Regime assignments to internal branches via Fitch parsimony
* Model comparison via AIC, AICc, BIC, and AICc weights
* JSON output support via ``--json``
* Added new CLI entry points: ``pk_ouwie``, ``pk_fit_ouwie``,
  ``pk_multi_regime_ou``
* Results validated against R 4.4.0 (``OUwie`` v2.10 with
  ``root.station=FALSE``):

  .. list-table::
     :header-rows: 1
     :widths: 12 18 18 12 30

     * - Model
       - PhyKIT LL
       - R OUwie LL
       - Diff
       - Notes
     * - BM1
       - -11.5697
       - -11.5697
       - < 1e-4
       - Exact match
     * - BMS
       - -11.2046
       - -11.1357
       - 0.069
       - Rooting artifact (R adds 1e-6 branch)
     * - OU1
       - -10.2890
       - -11.5697
       - 1.281
       - R stuck at alpha=0 (BM boundary)
     * - OUM
       - -8.6297
       - -10.9823
       - 2.353
       - R stuck at alpha=0 (BM boundary)
     * - OUMV
       - -6.9859
       - -10.2705
       - 3.285
       - R stuck at alpha=0 (BM boundary)
     * - OUMA
       - -6.9859
       - -6.9892
       - 0.003
       - Excellent match
     * - OUMVA
       - -6.9859
       - -7.0063
       - 0.020
       - Very close

  BM1 matches R to machine precision. OUMA and OUMVA agree within
  0.003-0.02 log-likelihood units. For OU1, OUM, and OUMV, R's OUwie
  optimizer converges to alpha=0 (the Brownian motion boundary),
  while PhyKIT's multi-interval search finds genuinely better OU optima
  with positive alpha. BMS shows a small difference (0.07 LL units)
  attributable to R's ``resolve.root=TRUE`` adding a 1e-6 length branch
  and optimizer convergence differences.

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
