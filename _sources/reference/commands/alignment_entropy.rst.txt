.. _cmd-alignment_entropy:
.. _command-alignment_entropy:

Alignment entropy
=================

Shannon entropy across alignment sites

Command identity
----------------

:Canonical command: ``alignment_entropy``
:Handler: ``alignment_entropy``
:Aliases: aln_entropy, entropy
:Standalone executables: pk_alignment_entropy, pk_aln_entropy, pk_entropy
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/alignment_entropy.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate alignment entropy.

Site-wise entropy is calculated using Shannon entropy. By default,
PhyKIT reports the mean entropy across all sites in the alignment.
With the `-v/--verbose` option, PhyKIT reports entropy for each site.

.. code-block:: shell

	phykit alignment_entropy <alignment> [-v/--verbose] [--plot] [--plot-output <path>]
		[--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
		[--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
		[--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Example output (default):

.. code-block:: shell

	0.657

Example output (`-v`):

.. code-block:: shell

	1	0.0
	2	1.0
	3	0.971

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*-v/--verbose*: optional argument to print entropy for each site |br|
*--plot*: save a per-site alignment entropy plot |br|
*--plot-output*: output path for plot (default: ``alignment_entropy_plot.png``) |br|
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
