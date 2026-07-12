.. _cmd-phylogenetic_ordination:
.. _command-phylogenetic_ordination:

Phylogenetic Ordination
=======================

Ordination incorporating phylogenetic structure (supports discordance-aware VCV with ``-g``)

Command identity
----------------

:Canonical command: ``phylogenetic_ordination``
:Handler: ``phylogenetic_ordination``
:Aliases: dimreduce, ord, ordination, pdr, phyl_pca, phylo_dimreduce, phylo_ordination, phylo_pca, ppca
:Standalone executables: pk_phylogenetic_ordination, pk_dimreduce, pk_ord, pk_ordination, pk_pdr, pk_phyl_pca, pk_phylo_dimreduce, pk_phylo_ordination, pk_phylo_pca, pk_ppca
:Categories: Phylogenetic comparative methods

Runtime interface
-----------------

.. include:: /_generated/commands/phylogenetic_ordination.inc

Guidance, interpretation, and examples
--------------------------------------

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

.. image:: /_static/docs_img/phylogenetic_pca_plot.png
   :align: center

|

.. image:: /_static/docs_img/phylogenetic_tsne_plot.png
   :align: center

|

.. image:: /_static/docs_img/phylogenetic_umap_plot.png
   :align: center

|

.. code-block:: shell

   phykit phylogenetic_ordination -t <tree> -d <trait_data> [--method <pca|tsne|umap>] [--correction <BM|lambda>] [--mode <cov|corr>] [--n-components <int>] [--perplexity <float>] [--n-neighbors <int>] [--min-dist <float>] [--seed <int>] [--plot] [--plot-tree] [--no-plot-tree] [--color-by <col_or_file>] [--tree-color-by <col_or_file>] [--plot-output <path>] [-g <gene_trees>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-t/--tree*: a tree file in Newick format |br|
*-d/--trait_data*: tab-delimited multi-trait file with header row |br|
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
*-g/--gene-trees*: optional multi-Newick file of gene trees; when provided, uses a discordance-aware VCV (genome-wide average) instead of the species-tree VCV |br|
*--fig-width*: figure width in inches (auto-scaled if omitted) |br|
*--fig-height*: figure height in inches (auto-scaled if omitted) |br|
*--dpi*: resolution in DPI (default: 300) |br|
*--no-title*: hide the plot title |br|
*--title*: custom title text |br|
*--legend-position*: legend location (e.g., "upper right", "none" to hide) |br|
*--ylabel-fontsize*: font size for y-axis labels; 0 to hide |br|
*--xlabel-fontsize*: font size for x-axis labels; 0 to hide |br|
*--title-fontsize*: font size for the title |br|
*--axis-fontsize*: font size for axis labels |br|
*--colors*: comma-separated colors (hex or named) |br|
*--ladderize*: ladderize (sort) the tree before plotting |br|
*--cladogram*: draw cladogram (equal branch lengths, tips aligned) instead of phylogram |br|
*--circular*: draw circular (radial/fan) phylogram instead of rectangular |br|
*--color-file*: color annotation file for tip labels, clade ranges, and branch colors (iTOL-inspired TSV format) |br|
*--json*: optional argument to print results as JSON
