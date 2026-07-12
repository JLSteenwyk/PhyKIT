.. _cmd-simmap_summary:
.. _command-simmap_summary:

SIMMAP summary
==============

Per-branch SIMMAP summary with node posteriors (describe.simmap)

Command identity
----------------

:Canonical command: ``simmap_summary``
:Handler: ``simmap_summary``
:Aliases: describe_simmap, smsummary
:Standalone executables: pk_simmap_summary, pk_describe_simmap, pk_smsummary
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/simmap_summary.inc

Guidance, interpretation, and examples
--------------------------------------

Run N stochastic character maps (SIMMAP) and summarize the results
per branch, analogous to ``phytools::describe.simmap()`` in R.

For each branch, reports:

- **Dwelling time proportions**: the average fraction of time spent in
  each character state across all simulations
- **Expected transitions**: the mean number of state changes on that branch

For each internal node, reports:

- **Posterior probabilities**: the proportion of simulations in which
  each state was sampled at that node

**Interpreting the output:**

- A branch with dwelling proportions ``carnivore=0.72, herbivore=0.17,
  omnivore=0.11`` spent 72% of its evolutionary time in the carnivore
  state (averaged across N maps). This indicates high confidence the
  lineage was carnivorous along that branch.
- A node with posteriors ``carnivore=0.56, herbivore=0.25, omnivore=0.19``
  was reconstructed as carnivore in 56% of simulations. Low values
  across all states indicate genuine ancestral uncertainty.
- ``E[trans]`` is the expected number of state changes on each branch.
  High values on long branches are expected. Values > 1 on short branches
  suggest rapid trait evolution.

.. image:: /_static/simmap_summary_posteriors.png
   :alt: PhyKIT simmap summary posteriors figure
   :align: center
   :width: 90%

*Tree with pie charts at internal nodes showing posterior state
probabilities from 100 stochastic maps. Branches colored by the
dominant state.*

.. code-block:: shell

   phykit simmap_summary -t <tree> -d <trait_data> -c <trait>
       [-m/--model ER|SYM|ARD] [-n/--nsim <int>]
       [--seed <int>] [--plot <file>] [--csv <file>] [--json]

Options: |br|
*-t/--tree*: phylogenetic tree file (required) |br|
*-d/--trait_data*: tab-delimited trait file with header row (required) |br|
*-c/--trait*: column name for the discrete character trait (required) |br|
*-m/--model*: substitution model — ``ER`` (equal rates, default), ``SYM`` (symmetric), or ``ARD`` (all rates different) |br|
*-n/--nsim*: number of stochastic maps to simulate (default: 100) |br|
*--seed*: random seed for reproducibility |br|
*--plot*: output plot file showing tree with posterior pie charts at nodes (.png, .pdf, .svg) |br|
*--csv*: output CSV file with per-branch dwelling proportions and node posteriors |br|
*--json*: output results as JSON

**R validation:** Q matrix and log-likelihood validated against
``phytools::fitMk()`` and ``phytools::make.simmap()``
(see ``tests/r_validation/validate_simmap_summary.R``).
