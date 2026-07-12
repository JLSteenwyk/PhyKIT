.. _cmd-fit_continuous:
.. _command-fit_continuous:

Continuous trait evolution model comparison (fitContinuous)
===========================================================

Compare continuous trait evolution models (supports discordance-aware VCV with ``-g``)

Command identity
----------------

:Canonical command: ``fit_continuous``
:Handler: ``fit_continuous``
:Aliases: fc, fitcontinuous
:Standalone executables: pk_fit_continuous, pk_fc, pk_fitcontinuous
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/fit_continuous.inc

Guidance, interpretation, and examples
--------------------------------------

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

Each model reports R² = 1 - (σ²_model / σ²_White), measuring how much
variance is explained relative to the white noise baseline. White serves
as the reference (R² = 0).

.. code-block:: shell

   phykit fit_continuous -t <tree> -d <trait_data> [--models BM,OU,Lambda] [-g <gene_trees>] [--json]

Options: |br|
*-t/--tree*: a tree file in Newick format |br|
*-d/--trait_data*: tab-delimited trait file (taxon<tab>value) |br|
*--models*: comma-separated list of models to fit (default: all 7) |br|
*-g/--gene-trees*: optional multi-Newick file of gene trees; when provided, uses a discordance-aware VCV (genome-wide average) instead of the species-tree VCV |br|
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

**R validation:** Validated against ``geiger`` in R
(see ``tests/r_validation/validate_fit_continuous_r2.R``).
