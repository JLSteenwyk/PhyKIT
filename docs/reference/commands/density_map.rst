.. _cmd-density_map:
.. _command-density_map:

Density map
===========

Posterior density of stochastic character maps

Command identity
----------------

:Canonical command: ``density_map``
:Handler: ``density_map``
:Aliases: densitymap, dmap
:Standalone executables: pk_density_map, pk_densitymap, pk_dmap
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/density_map.inc

Guidance, interpretation, and examples
--------------------------------------

Plot a phylogram with branches colored by posterior probabilities of
discrete character states from stochastic character mapping (analogous
to R's ``phytools::densityMap()``). Runs N simulations of stochastic
character mapping internally and, for each point along each branch,
computes the fraction of simulations in each state.

.. code-block:: shell

   phykit density_map -t <tree> -d <trait_data> -c <trait_column> -o <output.png> [-n <nsim>] [--seed <seed>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-t/--tree*: a tree file in Newick format |br|
*-d/--trait_data*: tab-delimited trait file (taxon<tab>state) |br|
*-c/--trait*: column name of the trait to map |br|
*-o/--output*: output plot file path (required) |br|
*-n/--nsim*: number of stochastic mapping simulations (default: 100) |br|
*--seed*: random seed for reproducibility |br|
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

.. image:: /_static/img/densitymap_example.png
   :align: center
   :width: 80%
