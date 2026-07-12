.. _tutorial-17:
.. _tutorial-automatic-detection-of-adaptive-shifts-on-a-phylogeny:

Tutorial 17: Automatic detection of adaptive shifts on a phylogeny
==================================================================

Objectives
----------

- Complete the automatic detection of adaptive shifts on a phylogeny workflow.
- Interpret the reported values and generated artifacts in their scientific context.
- Identify the canonical command references for each analysis step.

Prerequisites and working directory
-----------------------------------

Install the current PhyKIT release and create a dedicated working directory.
Download the data linked in this tutorial into that directory before running
the commands. All paths below are relative to this directory.

.. code-block:: shell

   mkdir phykit-tutorial-17
   cd phykit-tutorial-17

Related command references
--------------------------

- :doc:`Ou Shift Detection </reference/commands/ou_shift_detection>`

Workflow
--------

OUwie (tutorial 15, above) requires specifying regimes *a priori*. But
what if you don't know where on the tree the trait optimum changed?
PhyKIT's ``ou_shift_detection`` command (aliases: ``l1ou``, ``ou_shifts``,
``detect_shifts``) answers this question automatically using the
LASSO-based approach of
`Khabbazian et al. (2016) <https://doi.org/10.1093/sysbio/syw062>`_.

**When to use this command:** You have a continuous trait and a phylogeny
and want to discover which lineages experienced shifts in their adaptive
optimum — for example, identifying which lizard clades adapted to
different body sizes, or which mammal lineages evolved distinct metabolic
rates. Unlike OUwie, no regime file is needed.

.. centered::
   Download test data:
   :download:`Mammal phylogeny </data/tree_simple.tre>`;
   :download:`Continuous trait data </data/tree_simple_traits.tsv>`

|

Step 1: Prepare input files
*****************************

You need two files:

1. **Newick tree file** -- a phylogenetic tree with branch lengths.

2. **Trait data file** -- a tab-delimited file with two columns: taxon
   name and continuous trait value (no header row). Lines starting with
   ``#`` are treated as comments.

   .. code-block:: none

      raccoon	1.04
      bear	2.39
      sea_lion	2.30
      seal	1.88
      cat	0.5
      weasel	1.7
      monkey	0.3

|

Step 2: Run OU shift detection
*********************************

.. code-block:: shell

   phykit l1ou -t tree_simple.tre -d tree_simple_traits.tsv

Expected output:

.. code-block:: text

   ============================================================
   OU Shift Detection (l1ou)
   ============================================================
   Number of tips:       8
   Number of shifts:     1
   Selection criterion:  pBIC
   Alpha (OU strength):  0.784768
   Sigma² (BM rate):     0.407049
   Root optimum (θ₀):    1.758001
   Log-likelihood:       -5.9531
   pBIC:                 26.7533
   BIC:                  22.3034
   AICc:                 51.9062

   Detected shifts:
   ------------------------------------------------------------
     Shift 1: stem of (cat, monkey, weasel)
              New optimum: 0.286667
   ============================================================

**Interpretation.** The algorithm detected one adaptive shift on the stem
branch leading to the (cat, monkey, weasel) clade, with the trait optimum
shifting from 1.76 (root) to 0.29. This means these three species are
evolving toward a much lower body mass optimum than the rest of the tree.
Alpha = 0.78 indicates moderate pull-back strength toward the optima.
The pBIC of 26.75 is the most conservative criterion; BIC and AICc give
lower values, suggesting the shift is well-supported.

|

Step 3: Interpret the output
******************************

The output reports:

- **Number of shifts** -- how many adaptive regime changes were detected
- **Alpha** -- OU selection strength (higher = trait is pulled more
  strongly toward the optimum)
- **Sigma²** -- Brownian rate of stochastic variation
- **Root optimum (θ₀)** -- ancestral/background trait optimum
- **Shifts** -- each shift includes a description of where it occurred
  and the new optimum value for that regime

A shift described as "stem of (taxon_A, taxon_B, ... +N more)" means the
adaptive optimum changed on the branch leading to the clade containing
those taxa. A "terminal branch to taxon_X" means only that single species
shifted to a different optimum.

|

Step 4: Choose a model selection criterion
********************************************

Three criteria are available:

.. code-block:: shell

   phykit l1ou -t tree.nwk -d traits.tsv --criterion pBIC   # default, most conservative
   phykit l1ou -t tree.nwk -d traits.tsv --criterion BIC    # standard BIC
   phykit l1ou -t tree.nwk -d traits.tsv --criterion AICc   # corrected AIC

- **pBIC** (phylogenetic BIC) accounts for the combinatorial cost of
  choosing which edges shifted and is recommended for most analyses. It
  tends to be conservative, avoiding false positives.
- **BIC** and **AICc** may detect more shifts, but carry a higher risk of
  overfitting, especially on large trees.

|

Step 5: Limit the search
**************************

For large trees, you can cap the maximum number of shifts to speed up
the analysis:

.. code-block:: shell

   phykit l1ou -t tree.nwk -d traits.tsv --max-shifts 10

|

Step 6: JSON output for downstream analysis
**********************************************

.. code-block:: shell

   phykit l1ou -t tree.nwk -d traits.tsv --json

The JSON output includes all model parameters, information criterion
scores, and a list of detected shift edges with their descriptions and
optima.

|

Biological example: lizard body size evolution
*************************************************

Consider a phylogeny of 100 Anolis lizard species with log-transformed
body size measurements. Running:

.. code-block:: shell

   phykit l1ou -t lizards.tre -d lizard_body_size.tsv

might reveal 8 adaptive shifts, indicating that 8 lineages independently
evolved toward different body size optima. The output shows which clades
shifted and their new optimal values, enabling you to map adaptive
radiation without specifying regimes in advance.

|

Comparison with R's l1ou
**************************

PhyKIT's OU shift detection implementation has been validated against
R's l1ou package (Khabbazian et al. 2016). On a 100-tip Anolis dataset,
PhyKIT detects the same 8 shifts with matching alpha (0.607) and
comparable pBIC scores (PhyKIT: 17.6, R: 16.8). The small numerical
difference arises from optimizer precision, not algorithmic differences.

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
