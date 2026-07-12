.. _cmd-ouwie:
.. _command-ouwie:

Multi-regime OU models (OUwie)
==============================

Multi-regime Ornstein-Uhlenbeck models

Command identity
----------------

:Canonical command: ``ouwie``
:Handler: ``ouwie``
:Aliases: fit_ouwie, multi_regime_ou
:Standalone executables: pk_ouwie, pk_fit_ouwie, pk_multi_regime_ou
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/ouwie.inc

Guidance, interpretation, and examples
--------------------------------------

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

Each model reports R² = 1 - (σ²_model / σ²_BM1), measuring improvement
over the simplest Brownian motion baseline. For multi-regime models with
per-regime σ² values, the average is used.

.. code-block:: shell

   phykit ouwie -t <tree> -d <trait_data> -r <regime_data> [--models BM1,OUM,OUMVA] [--json]

Options: |br|
*-t/--tree*: a tree file in Newick format |br|
*-d/--trait_data*: tab-delimited trait file (taxon<tab>value) |br|
*-r/--regime_data*: tab-delimited regime file (taxon<tab>regime_label) |br|
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

**R validation:** Validated against ``OUwie`` in R
(see ``tests/r_validation/validate_ouwie_r2.R``).
