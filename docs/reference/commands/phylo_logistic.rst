.. _cmd-phylo_logistic:
.. _command-phylo_logistic:

Phylogenetic Logistic Regression
================================

Phylogenetic logistic regression for binary traits (Ives & Garland 2010)

Command identity
----------------

:Canonical command: ``phylo_logistic``
:Handler: ``phylo_logistic``
:Aliases: phylo_logreg, plogreg
:Standalone executables: pk_phylo_logistic, pk_phylo_logreg, pk_plogreg
:Categories: Phylogenetic comparative methods

Runtime interface
-----------------

.. include:: /_generated/commands/phylo_logistic.inc

Guidance, interpretation, and examples
--------------------------------------

Fit a Phylogenetic Logistic Regression for binary (0/1) response data while
accounting for phylogenetic non-independence among species (Ives & Garland 2010).

Uses Maximum Penalized Likelihood Estimation (logistic_MPLE) with Firth's
bias-correction penalty and jointly estimates the phylogenetic signal parameter
alpha via the OU-transformed variance-covariance matrix.

The multi-trait input file should be tab-delimited with a header row:
``taxon<tab>trait1<tab>trait2<tab>...``

Output includes coefficient estimates, standard errors, z-values,
p-values, alpha, log-likelihood, penalized log-likelihood, and AIC.

.. code-block:: shell

   phykit phylo_logistic -t <tree> -d <trait_data> --response <column> --predictor <column> [--method logistic_MPLE|logistic_IG10] [--json]

Options: |br|
*-t/--tree*: a tree file in Newick format |br|
*-d/--trait-data*: tab-delimited multi-trait file with header row |br|
*--response*: binary response column name (must contain only 0 and 1) |br|
*--predictor*: predictor column name(s), comma-separated for multiple |br|
*--method*: estimation method: logistic_MPLE or logistic_IG10 (default: logistic_MPLE) |br|
*--json*: optional argument to print results as JSON

**R validation:** Validated against ``phylolm`` in R
(see ``tests/r_validation/validate_phylo_logistic.R``).
