.. _cmd-stochastic_character_map:
.. _command-stochastic_character_map:

Stochastic character mapping (SIMMAP)
=====================================

Stochastic character mapping on a phylogeny

Command identity
----------------

:Canonical command: ``stochastic_character_map``
:Handler: ``stochastic_character_map``
:Aliases: scm, simmap
:Standalone executables: pk_stochastic_character_map, pk_scm, pk_simmap
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/stochastic_character_map.inc

Guidance, interpretation, and examples
--------------------------------------

Perform Stochastic Character Mapping (SIMMAP) of discrete traits onto a
phylogeny (Huelsenbeck et al. 2003; Bollback 2006), analogous to R's
``phytools::make.simmap()``.

Fits a continuous-time Markov chain (CTMC) rate matrix Q via maximum
likelihood using Felsenstein's pruning algorithm. Three substitution models
are available:

- **ER** (equal rates): 1 free parameter
- **SYM** (symmetric rates): k(k-1)/2 free parameters
- **ARD** (all rates differ): k(k-1) free parameters

The input trait file should be tab-delimited with a header row:
``taxon<tab>trait_column<tab>...``
Lines starting with '#' are treated as comments. If the tree and trait file
have different taxa, the intersection is used and warnings are printed to
stderr.

Output includes the fitted Q matrix, log-likelihood, mean dwelling times,
mean transition counts, and posterior probabilities at internal nodes.

Optionally, a horizontal phylogram plot with branches colored by mapped
character state can be generated using the ``--plot`` argument.

All results have been validated against R 4.4.0 using ``phytools::fitMk``
and ``phytools::make.simmap``. ER log-likelihoods match R to within 0.002
(both optimizers converge to the same flat plateau). For the ARD model,
PhyKIT's multi-start optimizer finds a slightly better local optimum
than R (loglik -8.38 vs -8.43). Dwelling times sum to total tree length
exactly, and Q matrix structural properties (rows sum to zero, model
nesting) are verified.

.. code-block:: shell

   phykit stochastic_character_map -t <tree> -d <trait_data> -c <trait_column> [-m <model>] [-n <nsim>] [--seed <seed>] [--plot <output.png>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-t/--tree*: a tree file in Newick format |br|
*-d/--trait_data*: tab-delimited trait file with header row |br|
*-c/--trait*: column name for discrete character trait |br|
*-m/--model*: substitution model: ER, SYM, or ARD (default: ER) |br|
*-n/--nsim*: number of stochastic mapping simulations (default: 100) |br|
*--seed*: random seed for reproducibility |br|
*--plot*: output plot file path for phylogram with colored branches |br|
*--fig-width*: figure width in inches (auto-scaled if omitted) |br|
*--fig-height*: figure height in inches (auto-scaled if omitted) |br|
*--dpi*: resolution in DPI (default: 300) |br|
*--no-title*: hide the plot title |br|
*--title*: custom title text |br|
*--legend-position*: legend location (e.g., "upper right", "none" to hide) |br|
*--ylabel-fontsize*: font size for y-axis labels; 0 to hide |br|
*--xlabel-fontsize*: font size for x-axis labels; 0 to hide |br|
*--title-fontsize*: font size for the title |br|
*--axis-fontsize*: font size for axis labels |br|
*--colors*: comma-separated colors (hex or named) |br|
*--ladderize*: ladderize (sort) the tree before plotting |br|
*--cladogram*: draw cladogram (equal branch lengths, tips aligned) instead of phylogram |br|
*--circular*: draw circular (radial/fan) phylogram instead of rectangular |br|
*--color-file*: color annotation file for tip labels, clade ranges, and branch colors (iTOL-inspired TSV format) |br|
*--json*: optional argument to print results as JSON
