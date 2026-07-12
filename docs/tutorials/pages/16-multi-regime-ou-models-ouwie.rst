.. _tutorial-16:
.. _tutorial-multi-regime-ou-models-ouwie:

Tutorial 16: Multi-regime OU models (OUwie)
===========================================

Objectives
----------

- Complete the multi-regime ou models (ouwie) workflow.
- Interpret the reported values and generated artifacts in their scientific context.
- Identify the canonical command references for each analysis step.

Prerequisites and working directory
-----------------------------------

Install the current PhyKIT release and create a dedicated working directory.
Download the data linked in this tutorial into that directory before running
the commands. All paths below are relative to this directory.

.. code-block:: shell

   mkdir phykit-tutorial-16
   cd phykit-tutorial-16

Related command references
--------------------------

- :doc:`Ouwie </reference/commands/ouwie>`

Workflow
--------

PhyKIT's ``ouwie`` command (aliases: ``fit_ouwie``, ``multi_regime_ou``)
fits multi-regime Ornstein-Uhlenbeck models, analogous to R's OUwie
package (Beaulieu et al. 2012). These models allow different clades to
evolve toward different trait optima with potentially different rates.

.. centered::
   Download test data:
   :download:`Mammal phylogeny </data/tree_simple.tre>`;
   :download:`Continuous trait data </data/tree_simple_traits.tsv>`;
   :download:`Regime assignments </data/tree_simple_regimes.tsv>`

|

Step 1: Prepare input files
*****************************

You need three files:

1. **Newick tree file** -- a phylogenetic tree in Newick format. The tree
   should include all taxa present in the trait and regime files.

2. **Trait data file** -- a tab-delimited file with two columns: taxon name
   and continuous trait value (no header row). Lines starting with ``#``
   are treated as comments.

   .. code-block:: none

      dog	1.1
      bear	1.9
      raccoon	1.5
      seal	1.8
      sea_lion	1.8
      cat	0.5
      weasel	1.7
      monkey	0.3

3. **Regime data file** -- a tab-delimited file with two columns: taxon
   name and discrete regime label (no header row). Regime labels can be
   any string (e.g., ``aquatic`` / ``terrestrial``, ``herbivore`` /
   ``carnivore``). Internal branch regimes are inferred via Fitch
   parsimony.

   .. code-block:: none

      dog	terrestrial
      bear	terrestrial
      raccoon	terrestrial
      seal	aquatic
      sea_lion	aquatic
      cat	terrestrial
      weasel	terrestrial
      monkey	terrestrial

|

Step 2: Run OUwie with all models
************************************

.. code-block:: shell

   phykit ouwie -t tree_simple.tre -d tree_simple_traits.tsv -r tree_simple_regimes.tsv

Expected output:

.. code-block:: text

   OUwie Model Comparison (Multi-Regime OU)

   Number of tips: 8
   Regimes: 2 (aquatic, terrestrial)

   Model   k    LL          AIC       AICc      dAICc    AICcW    BIC       dBIC     R2
   BM1     2    -11.570     27.14     29.54     0.00     0.759    27.30     2.93     0.000
   OU1     3    -10.289     26.58     32.58     3.04     0.166    26.82     2.45     -30.335
   BMS     3    -11.205     28.41     34.41     4.87     0.067    28.65     4.28     -0.034
   OUM     4    -8.630      25.26     38.59     9.05     0.008    25.58     1.21     -19.695
   OUMA    5    -6.986      23.97     53.97     24.43    0.000    24.37     0.00     -143.555
   OUMV    5    -6.986      23.97     53.97     24.43    0.000    24.37     0.00     -59.882
   OUMVA   6    -6.986      25.97     109.97    80.43    0.000    26.45     2.08     -861.136

   Best model (AICc): BM1
   Best model (BIC):  OUMA

**Interpretation.** By AICc (preferred for small samples), BM1 is the best model
with an AICc weight of 0.76, meaning a single Brownian motion rate adequately
explains the data. The more complex OU models are penalized by the AICc correction
for only 8 taxa. By BIC, OUMA wins — but this should be treated cautiously given
the small sample size. This illustrates why model selection criteria matter: AICc
is more conservative with few taxa.

|

Step 3: Run with a subset of models
**************************************

To compare only specific models, use ``--models``:

.. code-block:: shell

   phykit ouwie -t tree.nwk -d traits.tsv -r regimes.tsv --models BM1,OUM,OUMVA

|

Step 4: Interpret the output
******************************

The output table includes:

- **logLik** -- log-likelihood (higher = better fit)
- **AICc** -- small-sample corrected AIC (lower = better)
- **BIC** -- Bayesian Information Criterion (lower = better)
- **k** -- number of free parameters
- **AICc_w** -- AICc weight (relative support, sums to 1 across models)
- **Params** -- estimated parameter values

Key parameters:

- **sigma2** -- diffusion rate (BM rate of evolution)
- **alpha** -- strength of selection toward the optimum (OU pull-back
  parameter; larger values = stronger constraint)
- **theta** -- trait optimum for each regime (the value the trait is
  "pulled" toward)
- **z0** -- root ancestral state (BM models only)

|

Step 5: Model selection guidance
**********************************

- If **BM1** is best, a single Brownian motion rate explains the data.
- If **BMS** is best, different regimes evolve at different rates but
  without directional selection.
- If **OUM** is best, different regimes have different trait optima --
  the most common finding in comparative studies.
- If **OUMV** or **OUMA** is best, regimes differ in both optima and
  either rates (OUMV) or selection strength (OUMA).
- If **OUMVA** is best, all parameters are regime-specific. Be cautious
  with this model on small datasets as it has the most parameters (3R).

Use AICc for small datasets (n < 40 * k) and BIC when you prefer a
more conservative penalty on model complexity.

|

Step 6: JSON output for downstream analysis
**********************************************

.. code-block:: shell

   phykit ouwie -t tree.nwk -d traits.tsv -r regimes.tsv --json

The JSON output includes all model results, the best model by AICc and
BIC, regime labels, and parameter estimates. This is useful for
programmatic downstream analysis or integration into pipelines.

|

Comparison with R's OUwie
***************************

PhyKIT's OUwie implementation has been validated against R's OUwie
package (v2.10). BM1 matches R exactly. OUMA and OUMVA agree within
0.003-0.02 log-likelihood units. For OU1, OUM, and OUMV, PhyKIT's
multi-interval optimizer finds better OU optima than R's default
optimizer, which can get stuck at the BM boundary (alpha=0) on some
datasets.

|

Expected artifacts
------------------

Each step identifies its expected terminal output or generated files. Confirm
that those artifacts exist before continuing to the next step; filenames are
relative to the tutorial working directory unless an absolute path is shown.

Troubleshooting
---------------

- Run ``phykit <command> --help`` to compare an invocation with the live interface.
- Confirm that downloaded files are in the current working directory and retain
  the filenames shown in the tutorial.
- For parsing errors, compare taxon names exactly across alignments, trees, and
  trait tables, including capitalization and underscores.
- See :doc:`Troubleshooting </troubleshooting/index>` for installation, format,
  and error-reporting guidance.
