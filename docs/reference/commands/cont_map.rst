.. _cmd-cont_map:
.. _command-cont_map:

Continuous trait mapping (contMap)
==================================

Map continuous traits onto a phylogeny

Command identity
----------------

:Canonical command: ``cont_map``
:Handler: ``cont_map``
:Aliases: cmap, contmap
:Standalone executables: pk_cont_map, pk_cmap, pk_contmap
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/cont_map.inc

Guidance, interpretation, and examples
--------------------------------------

Plot a phylogram with branches colored by continuous trait values
(analogous to R's ``phytools::contMap()``). Ancestral states are
estimated via maximum-likelihood (two-pass Felsenstein algorithm)
and mapped onto branches using a color gradient (coolwarm colormap).

.. code-block:: shell

   phykit cont_map -t <tree> -d <trait_data> -o <output.png>
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-t/--tree*: a tree file in Newick format |br|
*-d/--trait_data*: tab-delimited trait file (taxon<tab>value) |br|
*-o/--output*: output plot file path (required) |br|
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

.. image:: /_static/img/contmap_example.png
   :alt: PhyKIT contmap example figure
   :align: center
   :width: 80%
