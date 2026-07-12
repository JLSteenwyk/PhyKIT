.. _cmd-create_concatenation_matrix:
.. _command-create_concatenation_matrix:

Create concatenation matrix
===========================

Concatenate multiple alignments into a supermatrix

Command identity
----------------

:Canonical command: ``create_concatenation_matrix``
:Handler: ``create_concatenation_matrix``
:Aliases: cc, create_concat
:Standalone executables: pk_create_concatenation_matrix, pk_cc, pk_create_concat
:Categories: Alignment & dataset utilities

Runtime interface
-----------------

.. include:: /_generated/commands/create_concatenation_matrix.inc

Guidance, interpretation, and examples
--------------------------------------

.. image:: /_static/docs_img/create_concat_matrix.png 
   :alt: PhyKIT create concat matrix figure
   :align: center
   :width: 75%


Create a concatenated alignment file. This function is 
used to help in the construction of multi-locus data
matrices.

PhyKIT will output three files:
1) A fasta file with '.fa' appended to the prefix specified with the -p/--prefix parameter.
2) A partition file ready for input into RAxML or IQ-tree.
3) An occupancy file that summarizes the taxon occupancy per sequence.

.. code-block:: shell

	phykit create_concat -a <file> -p <string> [--threshold <float>] [--plot-occupancy] [--plot-output <path>]
		[--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
		[--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
		[--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-a/--alignment*: alignment list file. ``--alignment-list`` and the legacy
``--alignment_list`` spelling are also accepted. The file should contain a single column list of alignment
sequence files to concatenate into a single matrix. Provide path to files relative to
working directory or provide absolute path. |br|
*-p/--prefix*: prefix of output files |br|
*--threshold*: minimum fraction of informative (non-gap, non-ambiguous) sites across the
concatenated alignment for a taxon to be included. Taxa whose effective occupancy falls
below this value are excluded from the output. Set to 0 to disable filtering
(default: 0). |br|
*--plot-occupancy*: optional argument to output an occupancy map figure where
x-axis shows concatenated positions with gene boundaries and y-axis shows taxa
sorted by total occupancy. Colors denote represented characters, gap/ambiguous
characters in present genes, and fully absent gene blocks. |br|
*--plot-output*: optional custom output path for occupancy map figure
(default: ``<prefix>.occupancy.png``). |br|
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
*--json*: optional argument to print summary metadata as JSON
