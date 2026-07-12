.. _cmd-ou_shift_detection:
.. _command-ou_shift_detection:

OU shift detection (l1ou)
=========================

Detect OU regime shifts on a phylogeny

Command identity
----------------

:Canonical command: ``ou_shift_detection``
:Handler: ``ou_shift_detection``
:Aliases: detect_shifts, l1ou, ou_shifts
:Standalone executables: pk_ou_shift_detection, pk_detect_shifts, pk_l1ou, pk_ou_shifts
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/ou_shift_detection.inc

Guidance, interpretation, and examples
--------------------------------------

Automatic OU shift detection using the LASSO-based approach from
Khabbazian et al. (2016). Discovers where on the phylogeny the adaptive
optimum changed without requiring an a priori regime assignment. Only a
tree and continuous trait data are needed.

The algorithm:

1. Fits a single-regime OU model to estimate alpha (selection strength)
2. Builds a design matrix with one column per candidate shift edge
3. Uses Cholesky transformation to remove phylogenetic correlation
4. Runs a LASSO path to identify candidate shift configurations
5. Selects the best model using pBIC, BIC, or AICc

.. code-block:: shell

   phykit l1ou -t <tree> -d <trait_data> [--criterion pBIC] [--max-shifts N] [--json]

Options: |br|
*-t/--tree*: a tree file in Newick format |br|
*-d/--trait_data*: tab-delimited trait file (taxon<tab>value) |br|
*--criterion*: model selection criterion: pBIC (default), BIC, or AICc |br|
*--max-shifts*: maximum number of shifts to consider (default: n/2) |br|
*--json*: optional argument to print results as JSON

Example output (no shifts detected):

.. code-block:: none

   ============================================================
   OU Shift Detection (l1ou)
   ============================================================
   Number of tips:       8
   Number of shifts:     0
   Selection criterion:  pBIC
   Alpha (OU strength):  0.784803
   Sigma² (BM rate):     1.203455
   Root optimum (θ₀):    1.206251
   Log-likelihood:       -10.2890
   pBIC:                 26.8163
   BIC:                  26.8163
   AICc:                 32.5780

   No shifts detected — single-regime OU is best.
   ============================================================

Example output (shifts detected):

.. code-block:: none

   ============================================================
   OU Shift Detection (l1ou)
   ============================================================
   Number of tips:       100
   Number of shifts:     8
   Selection criterion:  pBIC
   Alpha (OU strength):  0.606894
   Sigma² (BM rate):     0.062519
   Root optimum (θ₀):    0.248810
   Log-likelihood:       48.6896
   pBIC:                 17.6266
   BIC:                  -9.8811
   AICc:                 -49.8793

   Detected shifts:
   ------------------------------------------------------------
     Shift 1: terminal branch to valencienni
              New optimum: -0.564678
     Shift 2: terminal branch to insolitus
              New optimum: -0.876398
     Shift 3: stem of (barbatus, porcus, ... +2 more)
              New optimum: -0.635087
     Shift 4: stem of (altitudinalis, oporinus, ... +13 more)
              New optimum: -0.462944
   ============================================================

Results have been validated against R's l1ou package
(`Khabbazian et al. 2016 <https://doi.org/10.1093/sysbio/syw062>`_).
On a 100-tip lizard dataset, PhyKIT recovers the same 8 adaptive shifts
with matching alpha (0.607) and pBIC (17.6 vs R's 16.8).
