.. _cmd-tip_to_tip_distance:
.. _command-tip_to_tip_distance:

Tip-to-tip distance
===================

Distance between two tips in a tree

Command identity
----------------

:Canonical command: ``tip_to_tip_distance``
:Handler: ``tip_to_tip_distance``
:Aliases: t2t, t2t_dist
:Standalone executables: pk_tip_to_tip_distance, pk_t2t, pk_t2t_dist
:Categories: Tree summary statistics

Runtime interface
-----------------

.. include:: /_generated/commands/tip_to_tip_distance.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate distance between two tips (or leaves) in a phylogeny.

Distances are in substitutions per site.

.. code-block:: shell

   phykit tip_to_tip_distance <tree_file> <tip_1> <tip_2> [--json]
   phykit tip_to_tip_distance <tree_file> --all-pairs [--plot] [--plot-output <path>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*<tree_file>*: first argument after function name should be a tree file |br|
*<tip_1>*: second argument should be the name of the first tip of interest |br|
*<tip_2>*: third argument should be the name of the second tip of interest |br|
*--all-pairs*: optional argument to report all pairwise tip distances |br|
*--plot*: optional argument to save a clustered distance heatmap (requires ``--all-pairs``) |br|
*--plot-output*: output path for heatmap (default: ``tip_to_tip_distance_heatmap.png``) |br|
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
