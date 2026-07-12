.. _cmd-pairwise_identity:
.. _command-pairwise_identity:

Pairwise identity
=================

Pairwise sequence identity in an alignment

Command identity
----------------

:Canonical command: ``pairwise_identity``
:Handler: ``pairwise_identity``
:Aliases: pairwise_id, pi
:Standalone executables: pk_pairwise_identity, pk_pairwise_id, pk_pi
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/pairwise_identity.inc

Guidance, interpretation, and examples
--------------------------------------

.. image:: /_static/docs_img/pairwise_identity.png 
   :alt: PhyKIT pairwise identity figure
   :align: center
   :width: 75%


Calculate the average pairwise identity among sequences.

Pairwise identities can be used as proxies for the evolutionary rate of sequences.

Pairwise identity is defined as the number of identical
columns (including gaps) between two aligned sequences divided
by the number of columns in the alignment. Summary statistics
are reported unless used with the verbose option in which
all pairwise identities will be reported.

An example of pairwise identities being used as a proxy
for evolutionary rate can be found here: Chen et al. 
Genome Biology and Evolution (2017), doi: 10.1093/gbe/evx147.

.. code-block:: shell

	phykit pairwise_identity <alignment> [-v/--verbose] [-e/--exclude_gaps] [--plot] [--plot-output <file>]
		[--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
		[--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
		[--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*-v/--verbose*: optional argument to print identity per pair |br|
*-e/--exclude_gaps*: if a site has a gap, ignore it |br|
*--plot*: save a clustered pairwise-identity heatmap |br|
*--plot-output*: output path for heatmap (default: pairwise_identity_heatmap.png) |br|
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
