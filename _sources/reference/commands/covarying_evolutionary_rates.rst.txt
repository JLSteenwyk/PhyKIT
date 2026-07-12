.. _cmd-covarying_evolutionary_rates:
.. _command-covarying_evolutionary_rates:

Covarying evolutionary rates
============================

Detect covariation in evolutionary rates

Command identity
----------------

:Canonical command: ``covarying_evolutionary_rates``
:Handler: ``covarying_evolutionary_rates``
:Aliases: cover
:Standalone executables: pk_covarying_evolutionary_rates, pk_cover
:Categories: Evolutionary rate analysis

Runtime interface
-----------------

.. include:: /_generated/commands/covarying_evolutionary_rates.inc

Guidance, interpretation, and examples
--------------------------------------

Determine if two genes have a signature of covariation with one another.
Genes that have covarying evolutionary histories tend to have 
similar functions and expression levels.

Input two phylogenies and calculate the correlation among relative 
evolutionary rates between the two phylogenies. The two input trees 
do not have to have the same taxa. This function will first prune both
trees to have the same tips. To transform branch lengths into relative
rates, PhyKIT uses the putative species tree's branch lengths, which is
input by the user. As recommended by the original method developers,
outlier branch lengths are removed. Outlier branches have a relative 
evolutionary rate greater than five.

PhyKIT reports two tab delimited values:
col1: correlation coefficient
col2: p-value

Method is empirically evaluated by Clark et al., Genome Research
(2012), doi: 10.1101/gr.132647.111. Normalization method using a 
species tree follows Sato et al., Bioinformatics (2005), doi: 
10.1093/bioinformatics/bti564.  

.. code-block:: shell

   phykit covarying_evolutionary_rates <tree_file_zero> <tree_file_one> -r/--reference <reference_tree_file> [-v/--verbose] [--plot] [--plot-output <path>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*<tree_file_zero>*: first argument after function name should be a tree file |br|
*<tree_file_one>*: second argument after function name should be a tree file |br| 
*-r/--reference*: a tree to correct branch lengths by in the two input trees. Typically, 
this is a putative species tree. |br|
*-v/--verbose*: print out corrected branch lengths shared between tree 0 and tree 1 |br|
*--plot*: save a covarying-rates scatter plot |br|
*--plot-output*: output path for plot (default: ``covarying_rates_plot.png``) |br|
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
