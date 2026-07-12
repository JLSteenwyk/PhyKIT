.. _cmd-phylogenetic_regression:
.. _command-phylogenetic_regression:

Phylogenetic regression (PGLS)
==============================

Phylogenetic generalized least squares regression (supports discordance-aware VCV with ``-g``)

Command identity
----------------

:Canonical command: ``phylogenetic_regression``
:Handler: ``phylogenetic_regression``
:Aliases: pgls, phylo_regression
:Standalone executables: pk_phylogenetic_regression, pk_pgls, pk_phylo_regression
:Categories: Phylogenetic comparative methods

Runtime interface
-----------------

.. include:: /_generated/commands/phylogenetic_regression.inc

Guidance, interpretation, and examples
--------------------------------------

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

A three-way variance decomposition is also reported: R²_total (variance
explained by phylogeny + predictors combined), R²_pred (predictor contribution
given phylogeny, = standard R²), and R²_phylo (phylogeny's unique contribution).
R²_phylo + R²_pred = R²_total.

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

   phykit phylogenetic_regression -t <tree> -d <trait_data> -y <response> -x <predictor1> [predictor2 ...] [-m <method>] [-g <gene_trees>] [--json]

Options: |br|
*-t/--tree*: a tree file in Newick format |br|
*-d/--trait_data*: tab-delimited multi-trait file with header row |br|
*-y/--response*: response (dependent) variable column name |br|
*-x/--predictors*: one or more predictor column names |br|
*-m/--method*: method to use: BM or lambda (default: BM) |br|
*-g/--gene-trees*: optional multi-Newick file of gene trees; when provided, uses a discordance-aware VCV (genome-wide average) instead of the species-tree VCV |br|
*--json*: optional argument to print results as JSON

**R validation:** Validated against ``caper``, ``nlme`` in R
(see ``tests/r_validation/validate_pgls_r2.R``).
