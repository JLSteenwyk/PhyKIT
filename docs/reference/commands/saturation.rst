.. _cmd-saturation:
.. _command-saturation:

Saturation
==========

Test for substitution saturation

Command identity
----------------

:Canonical command: ``saturation``
:Handler: ``saturation``
:Aliases: sat
:Standalone executables: pk_saturation, pk_sat
:Categories: Saturation & model adequacy

Runtime interface
-----------------

.. include:: /_generated/commands/saturation.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate saturation for a given tree and alignment.

Saturation is defined as sequences in multiple sequence
alignments that have undergone numerous substitutions such
that the distances between taxa are underestimated.

Data with no saturation will have a value of 1. The closer
the value is to 1, the less saturated the data.

This function outputs two values (as of v1.19.9). The first
value is the saturation value and the second column is the absolute
value of saturation minus 1. Thus, lower values in the second column
are indicative of values closer to one and, thus, less saturation.

Saturation is calculated following Philippe et al., PLoS 
Biology (2011), doi: 10.1371/journal.pbio.1000602.

.. code-block:: shell

   phykit saturation -a <alignment> -t <tree> [-v/--verbose] [-e/--exclude_gaps] [--plot] [--plot-output <path>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-a/--alignment*: an alignment file |br|
*-t/--tree*: a tree file |br|
*-e/--exclude_gaps*: if a site has a gap, ignore it |br|
*-v/--verbose*: print out patristic distances and uncorrected |br|
distances used to determine saturation |br|
*--plot*: save a saturation scatter plot with fitted slope through origin |br|
*--plot-output*: output path for saturation plot (default: ``saturation_plot.png``) |br|
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
