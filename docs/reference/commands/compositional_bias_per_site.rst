.. _cmd-compositional_bias_per_site:
.. _command-compositional_bias_per_site:

Compositional bias per site
===========================

Detect compositional bias across sites

Command identity
----------------

:Canonical command: ``compositional_bias_per_site``
:Handler: ``compositional_bias_per_site``
:Aliases: cbps, comp_bias_per_site
:Standalone executables: pk_compositional_bias_per_site, pk_cbps, pk_comp_bias_per_site
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/compositional_bias_per_site.inc

Guidance, interpretation, and examples
--------------------------------------

Calculates compositional bias per site in an alignment.

Site-wise chi-squared tests are conducted in an alignment to
detect compositional biases. PhyKIT outputs four columns: |br|
col 1: index in alignment |br|
col 2: chi-squared statistic (higher values indicate greater bias) |br|
col 3: multi-test corrected p-value (Benjamini-Hochberg false discovery rate procedure) |br|
col 4: uncorrected p-value

.. code-block:: shell

	phykit comp_bias_per_site <alignment> [--plot] [--plot-output <path>]
		[--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
		[--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
		[--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*<alignment>*: first argument after function name should be a query
fasta alignment to calculate the site-wise compositional bias of |br|
*--plot*: save a Manhattan-style plot of site-wise compositional bias |br|
*--plot-output*: output path for plot (default: ``compositional_bias_per_site_plot.png``) |br|
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
