.. _cmd-phylogenetic_glm:
.. _command-phylogenetic_glm:

Phylogenetic GLM
================

Phylogenetic generalized linear model (supports discordance-aware VCV with ``-g``)

Command identity
----------------

:Canonical command: ``phylogenetic_glm``
:Handler: ``phylogenetic_glm``
:Aliases: pglm, phylo_glm
:Standalone executables: pk_phylogenetic_glm, pk_pglm, pk_phylo_glm
:Categories: Phylogenetic comparative methods

Runtime interface
-----------------

.. include:: /_generated/commands/phylogenetic_glm.inc

Guidance, interpretation, and examples
--------------------------------------

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
p-values, log-likelihood, AIC, and McFadden's pseudo-R² (computed from
full vs. intercept-only model log-likelihoods).

.. code-block:: shell

   phykit phylogenetic_glm -t <tree> -d <trait_data> -y <response> -x <predictor1> [predictor2 ...] --family <binomial|poisson> [-g <gene_trees>] [--json]

Options: |br|
*-t/--tree*: a tree file in Newick format |br|
*-d/--trait_data*: tab-delimited multi-trait file with header row |br|
*-y/--response*: response (dependent) variable column name |br|
*-x/--predictors*: one or more predictor column names |br|
*--family*: distribution family: binomial or poisson |br|
*--method*: estimation method: logistic_MPLE or poisson_GEE (auto from family) |br|
*--btol*: linear predictor bound for logistic model (default: 10) |br|
*--log-alpha-bound*: bound on log(alpha) for logistic model (default: 4) |br|
*-g/--gene-trees*: optional multi-Newick file of gene trees; when provided, uses a discordance-aware VCV (genome-wide average) instead of the species-tree VCV |br|
*--json*: optional argument to print results as JSON

**R validation:** Validated against ``phylolm`` in R
(see ``tests/r_validation/validate_glm_pseudo_r2.R``).
