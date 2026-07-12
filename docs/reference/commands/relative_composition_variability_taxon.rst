.. _cmd-relative_composition_variability_taxon:
.. _command-relative_composition_variability_taxon:

Relative composition variability, taxon
=======================================

Per-taxon relative composition variability

Command identity
----------------

:Canonical command: ``relative_composition_variability_taxon``
:Handler: ``rcvt``
:Aliases: rcvt, rel_comp_var_taxon
:Standalone executables: pk_relative_composition_variability_taxon, pk_rcvt, pk_rel_comp_var_taxon
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/relative_composition_variability_taxon.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate RCVT (relative composition variability, taxon) for an alignment.

RCVT is the relative composition variability metric for individual taxa.
This facilitates identifying specific taxa that may have compositional
biases. Lower RCVT values are more desirable because they indicate
a lower composition bias for a given taxon in an alignment.

RCVT calculations are case-insensitive and exclude gap/ambiguous symbols
from composition counts and normalization.

.. code-block:: shell

	phykit relative_composition_variability_taxon <alignment> [--plot] [--plot-output <path>]
		[--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
		[--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
		[--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--plot*: optional argument to generate an RCVT per-taxon barplot |br|
*--plot-output*: output path for the RCVT plot (default: ``rcvt_plot.png``) |br|
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
