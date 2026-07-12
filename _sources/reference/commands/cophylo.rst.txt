.. _cmd-cophylo:
.. _command-cophylo:

Cophylogenetic plot (tanglegram)
================================

Tanglegram for comparing two trees

Command identity
----------------

:Canonical command: ``cophylo``
:Handler: ``cophylo``
:Aliases: tangle, tanglegram
:Standalone executables: pk_cophylo, pk_tangle, pk_tanglegram
:Categories: Tree comparison & consensus

Runtime interface
-----------------

.. include:: /_generated/commands/cophylo.inc

Guidance, interpretation, and examples
--------------------------------------

Plot a cophylogenetic tanglegram of two phylogenies (analogous to R's
``phytools::cophylo()``). Draws two trees facing each other with
connecting lines between matching taxa. By default, taxa are matched
by identical tip names. Internal nodes of tree2 are rotated to minimize
line crossings.

.. code-block:: shell

   phykit cophylo -t <tree1> -t2 <tree2> -o <output.png> [-m <mapping>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-t/--tree1*: first tree file in Newick format |br|
*-t2/--tree2*: second tree file in Newick format |br|
*-o/--output*: output plot file path (required) |br|
*-m/--mapping*: optional tab-delimited mapping file (taxon1<tab>taxon2) |br|
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

.. image:: /_static/img/cophylo_example.png
   :alt: PhyKIT cophylo example figure
   :align: center
   :width: 80%
