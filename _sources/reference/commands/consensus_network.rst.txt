.. _cmd-consensus_network:
.. _command-consensus_network:

Consensus network
=================

Consensus network from multiple trees

Command identity
----------------

:Canonical command: ``consensus_network``
:Handler: ``consensus_network``
:Aliases: consnet, splitnet, splits_network
:Standalone executables: pk_consensus_network, pk_consnet, pk_splitnet, pk_splits_network
:Categories: Tree comparison & consensus

Runtime interface
-----------------

.. include:: /_generated/commands/consensus_network.inc

Guidance, interpretation, and examples
--------------------------------------

Extract bipartition splits from a collection of gene trees and summarize
conflicting phylogenetic signal. Counts how frequently each non-trivial
bipartition appears across input trees and filters by a minimum frequency
threshold. Optionally draws a circular splits network diagram.

Polytomies (collapsed branches) in input trees are handled conservatively:
splits from polytomous nodes are excluded since they represent unresolved
relationships. Trifurcating roots (standard unrooted Newick) are not
affected. This allows gene trees with collapsed low-support branches to be
used directly as input.

Input can be either:
1) a file with one Newick tree per line, or
2) a file with one tree-file path per line.

.. code-block:: shell

   phykit consensus_network -t/--trees <trees> [--threshold 0.1] [--missing-taxa allow|error|shared] [--plot-output <file>]
       [--max-splits <int>] [--histogram <file>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-t/--trees*: file containing trees (one Newick per line) or tree-file paths (one per line) |br|
*--threshold*: minimum split frequency to include, between 0 and 1 (default: ``0.1``) |br|
*--missing-taxa*: handling strategy for mismatched taxa: ``allow`` keeps each
tree's available taxa, ``shared`` prunes to the intersection, and ``error``
rejects mismatched sets (default: ``allow``) |br|
*--plot-output*: output filename for the circular splits network plot (optional) |br|
*--max-splits*: maximum number of splits drawn in the network graph; all
qualifying splits remain in text and JSON output (default: ``30``) |br|
*--histogram*: optional output filename for a split-frequency histogram |br|
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

When ``--plot-output`` is specified, a circular splits network diagram is produced.
Taxa are arranged at equal angles around a circle. Each split is drawn as a chord
connecting the boundary points between the two sides of the bipartition. Chord
thickness and opacity scale with split frequency — thicker, darker lines indicate
splits supported by more gene trees.

.. image:: /_static/img/consensus_network_example.png
   :align: center
   :width: 80%
