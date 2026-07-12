.. _cmd-phenogram:
.. _command-phenogram:

Phenogram (traitgram)
=====================

Phenogram visualizing trait evolution

Command identity
----------------

:Canonical command: ``phenogram``
:Handler: ``phenogram``
:Aliases: tg, traitgram
:Standalone executables: pk_phenogram, pk_tg, pk_traitgram
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/phenogram.inc

Guidance, interpretation, and examples
--------------------------------------

Plot a phenogram (traitgram) showing trait evolution along a phylogeny
(analogous to R's ``phytools::phenogram()``). The X-axis represents
distance from the root and the Y-axis represents trait values.
Ancestral states are reconstructed via maximum-likelihood.

.. code-block:: shell

   phykit phenogram -t <tree> -d <trait_data> -o <output.png>
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

.. image:: /_static/img/phenogram_example.png
   :align: center
   :width: 80%
