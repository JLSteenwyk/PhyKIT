.. _cmd-evolutionary_rate_per_site:
.. _command-evolutionary_rate_per_site:

Evolutionary Rate per Site
==========================

Site-specific evolutionary rate estimation

Command identity
----------------

:Canonical command: ``evolutionary_rate_per_site``
:Handler: ``evolutionary_rate_per_site``
:Aliases: erps, evo_rate_per_site
:Standalone executables: pk_evolutionary_rate_per_site, pk_erps, pk_evo_rate_per_site
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/evolutionary_rate_per_site.inc

Guidance, interpretation, and examples
--------------------------------------

Estimate evolutionary rate per site.

Evolutionary rate per site is one minus the sum of squared frequency of different
characters at a given site. Values may range from 0 (slow evolving; no diversity
at the given site) to 1 (fast evolving; all characters appear only once).

.. code-block:: shell

	phykit evo_rate_per_site <alignment> [--plot] [--plot-output <path>]
		[--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
		[--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
		[--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*<alignment>*: first argument after function name should be a query
fasta alignment to calculate the site-wise evolutionary rate of |br|
*--plot*: save a per-site evolutionary-rate plot |br|
*--plot-output*: output path for plot (default: ``evolutionary_rate_per_site_plot.png``) |br|
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
