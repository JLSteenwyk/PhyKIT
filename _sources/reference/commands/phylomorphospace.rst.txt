.. _cmd-phylomorphospace:
.. _command-phylomorphospace:

Phylomorphospace
================

Phylomorphospace visualization

Command identity
----------------

:Canonical command: ``phylomorphospace``
:Handler: ``phylomorphospace``
:Aliases: phmo, phylomorpho
:Standalone executables: pk_phylomorphospace, pk_phmo, pk_phylomorpho
:Categories: Phylogenetic comparative methods

Runtime interface
-----------------

.. include:: /_generated/commands/phylomorphospace.inc

Guidance, interpretation, and examples
--------------------------------------

Plot a phylomorphospace: two raw traits in trait space with the phylogeny
overlaid via ML-reconstructed ancestral states at internal nodes. This differs
from the ``phylogenetic_ordination --plot-tree`` option, which plots in PC space;
``phylomorphospace`` plots raw trait values directly on the x and y axes.

Tree edges connect species through ML-reconstructed ancestral states and are
colored by distance from root (coolwarm colormap with colorbar). Tip points
default to blue, or can be colored by a continuous or discrete variable using
the ``--color-by`` option.

The multi-trait input file should be tab-delimited with a header row:
``taxon<tab>trait1<tab>trait2<tab>...``
Lines starting with '#' are treated as comments. If the tree and trait file
have different taxa, the intersection is used and warnings are printed to
stderr.

If the trait file has exactly 2 trait columns and ``--trait-x`` / ``--trait-y``
are omitted, the first two columns are selected automatically.

.. image:: /_static/docs_img/phylomorphospace_plot.png
   :alt: PhyKIT phylomorphospace plot figure
   :align: center

|

.. code-block:: shell

   phykit phylomorphospace -t <tree> -d <trait_data> [--trait-x <name>] [--trait-y <name>] [--color-by <col_or_file>] [--plot-output <path>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-t/--tree*: a tree file in Newick format |br|
*-d/--trait_data*: tab-delimited multi-trait file with header row |br|
*--trait-x*: column name for x-axis trait |br|
*--trait-y*: column name for y-axis trait |br|
*--color-by*: color tip points by a trait; specify a column name from the multi-trait file or a separate tab-delimited file (taxon<tab>value) for continuous or discrete coloring |br|
*--plot-output*: output path for plot (default: phylomorphospace_plot.png) |br|
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
