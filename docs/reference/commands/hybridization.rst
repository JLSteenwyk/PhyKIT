.. _cmd-hybridization:
.. _command-hybridization:

Hybridization analysis
======================

Estimate reticulation events and localize hybridization on a species tree

Command identity
----------------

:Canonical command: ``hybridization``
:Handler: ``hybridization``
:Aliases: hybrid, reticulation
:Standalone executables: pk_hybridization, pk_hybrid, pk_reticulation
:Categories: Introgression & gene flow

Runtime interface
-----------------

.. include:: /_generated/commands/hybridization.inc

Guidance, interpretation, and examples
--------------------------------------

Estimate the minimum number of reticulation events and localize where
hybridization likely occurred on a species tree. For each internal branch,
tests whether the two discordant NNI topologies are significantly
asymmetric (binomial test with FDR correction). Asymmetric discordance
suggests introgression/gene flow rather than ILS.

Each branch receives:

- **p-value**: binomial test of gDF1 vs gDF2 equality (H0: ILS alone)
- **FDR-corrected p-value**: Benjamini-Hochberg correction across all branches
- **Hybridization score**: asymmetry ratio (0.5 = no signal, 1.0 = strong signal)
  for significant branches; 0 for non-significant
- **Favored alternative**: which NNI topology is in excess, identifying
  which pair of lineages likely exchanged genes (but not the direction
  of gene flow)

The estimated reticulation count is the number of branches with significant
asymmetric discordance. Gene tree uncertainty can be accounted for with
``--support`` to collapse low-support branches before topology determination.

.. code-block:: shell

   phykit hybridization -t <species_tree> -g <gene_trees>
       [--support <float>] [--alpha 0.05] [--plot <output>] [--json]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>]
       [--ladderize] [--cladogram] [--circular] [--ylabel-fontsize <float>]

Options: |br|
*-t/--tree*: species tree file (required) |br|
*-g/--gene-trees*: gene trees file, one Newick per line (required) |br|
*--support*: collapse gene tree branches below this support value before analysis |br|
*--alpha*: FDR significance threshold (default: 0.05) |br|
*--plot*: output hybridization score phylogram (branches colored by score, stars at significant nodes) |br|
*--json*: optional argument to print results as JSON
