.. _tutorial-15:
.. _tutorial-comparing-continuous-trait-evolution-models:

Tutorial 15: Comparing continuous trait evolution models
========================================================

Objectives
----------

- Complete the comparing continuous trait evolution models workflow.
- Interpret the reported values and generated artifacts in their scientific context.
- Identify the canonical command references for each analysis step.

Prerequisites and working directory
-----------------------------------

Install the current PhyKIT release and create a dedicated working directory.
Download the data linked in this tutorial into that directory before running
the commands. All paths below are relative to this directory.

.. code-block:: shell

   mkdir phykit-tutorial-15
   cd phykit-tutorial-15

Related command references
--------------------------

- :doc:`Fit Continuous </reference/commands/fit_continuous>`

Workflow
--------

A common analysis in comparative methods is determining which model of
continuous trait evolution best explains observed trait variation on a
phylogeny. PhyKIT's ``fit_continuous`` command (aliases: ``fitcontinuous``,
``fc``) fits up to 7 models and ranks them by AIC, BIC, and AIC weights,
analogous to R's ``geiger::fitContinuous()``.

.. centered::
   Download test data:
   :download:`Mammal phylogeny </data/tree_simple.tre>`;
   :download:`Continuous trait data </data/tree_simple_traits.tsv>`

|

Step 0: Prepare data
*********************

You need a Newick tree file and a tab-delimited trait file
(``taxon<tab>value``). For example:

.. code-block:: text

   raccoon	1.04
   bear	2.39
   sea_lion	2.30
   seal	1.88

|

Step 1: Run fit_continuous with all models
*******************************************

.. code-block:: shell

   phykit fit_continuous -t tree_simple.tre -d tree_simple_traits.tsv

Expected output:

.. code-block:: text

   Model Comparison (fitContinuous)

   Number of tips: 8

   Model       Param     Value      Sigma2    z0        LL         AIC      dAIC     AICw     BIC      dBIC     R2
   White       -         -          0.7667    1.2062    -10.289    24.58    0.00     0.304    24.74    0.00     0.000
   EB          a         -0.0785    0.0854    1.4827    -9.595     25.19    0.61     0.224    25.43    0.69     0.889
   Kappa       kappa     0.0100     0.3428    1.3230    -9.722     25.44    0.87     0.197    25.68    0.94     0.553
   OU          alpha     0.7848     1.2035    1.2063    -10.289    26.58    2.00     0.112    26.82    2.08     -0.570
   BM          -         -          0.0384    1.6447    -11.570    27.14    2.56     0.084    27.30    2.56     0.950
   Delta       delta     0.5188     0.1968    1.4939    -11.128    28.26    3.68     0.048    28.49    3.76     0.743
   Lambda      lambda    1.0000     0.0384    1.6447    -11.570    29.14    4.56     0.031    29.38    4.64     0.950

   Best model (AIC): White
   Best model (BIC): White

|

Step 2: Interpret the AIC/BIC table
*************************************

The output table shows each model's parameter estimate, sigma-squared,
ancestral state (z0), log-likelihood, AIC, delta-AIC, AIC weight, BIC,
and delta-BIC. Lower AIC/BIC values and higher AIC weights indicate
better-fitting models.

In this example, the White model (no phylogenetic structure) has the
lowest AIC and BIC, suggesting that with only 8 taxa the data are too
sparse to distinguish phylogenetic from non-phylogenetic models. However,
the R² column shows that BM (0.95) and Lambda (0.95) explain most of the
trait variance relative to the White model — with more taxa, these models
would likely be preferred. The EB model (R² = 0.89) also fits well,
consistent with early rapid evolution of body mass.

|

Step 3: Run with a subset of models
**************************************

.. code-block:: shell

   phykit fc -t tree.nwk -d traits.tsv --models BM,OU,Lambda

|

Step 4: JSON output for downstream analysis
*********************************************

.. code-block:: shell

   phykit fc -t tree.nwk -d traits.tsv --json

The JSON output includes all model results, the best model by AIC and
BIC, and the number of tips.

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
