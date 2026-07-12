.. _cmd-quartet_network:
.. _command-quartet_network:

Quartet network
===============

Quartet-based network visualization

Command identity
----------------

:Canonical command: ``quartet_network``
:Handler: ``quartet_network``
:Aliases: nanuq, qnet, quartet_net
:Standalone executables: pk_quartet_network, pk_nanuq, pk_qnet, pk_quartet_net
:Categories: Tree comparison & consensus

Runtime interface
-----------------

.. include:: /_generated/commands/quartet_network.inc

Guidance, interpretation, and examples
--------------------------------------

Quartet-based network inference (NANUQ-style) for distinguishing incomplete
lineage sorting (ILS) from hybridization/gene flow using quartet concordance
factors from gene trees.

For each 4-taxon subset, counts how many gene trees display each of the 3
possible unrooted topologies and applies two hypothesis tests:

1. **Star test** (Pearson chi-squared): tests whether the three topology
   counts are consistent with a star tree (equal probabilities 1/3 each).
   If p_star > beta (default 0.95), the quartet is classified as
   *unresolved*.

2. **T3 tree model test** (G-test / likelihood ratio): tests whether the
   counts are consistent with any resolved quartet tree under the
   multispecies coalescent.  If p_tree > alpha (default 0.05), the quartet
   is classified as *tree-like* (conflict is due to ILS).

3. If both tests reject their null hypotheses, the quartet is classified as
   *hybrid* (asymmetric discordance indicating gene flow or hybridization).

This algorithm matches the NANUQ method of Allman, Baños & Rhodes (2019),
implemented in R's MSCquartets package.

Polytomies (collapsed branches) in gene trees are handled conservatively:
bipartitions from polytomous nodes are excluded, so quartets spanning a
polytomy are treated as unresolved rather than misclassified. Trifurcating
roots (standard unrooted Newick) are not affected. This allows gene trees
with collapsed low-support branches to be used directly as input.

Input can be either:
1) a file with one Newick tree per line, or
2) a file with one tree-file path per line.

.. code-block:: shell

   phykit quartet_network -t/--trees <trees> [--alpha 0.05] [--beta 0.95] [--missing-taxa error|shared] [--plot-output <file>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-t/--trees*: file containing trees (one Newick per line) or tree-file paths (one per line) |br|
*--alpha*: significance level for the T3 tree model test (default: ``0.05``) |br|
*--beta*: threshold for the star tree test; quartets with p_star > beta are called unresolved (default: ``0.95``) |br|
*--missing-taxa*: handling strategy for mismatched taxa (``error`` or ``shared``; default: ``error``) |br|
*--plot-output*: output filename for the quartet network plot (optional) |br|
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

When ``--plot-output`` is specified, a NANUQ-style splits graph is drawn from
the quartet distance matrix using Neighbor-Joining and circular split
decomposition.  Tree-like relationships appear as simple branching, while
hybridization / reticulation produces characteristic box (parallelogram)
structures — the same style of output produced by R's MSCquartets +
NeighborNet pipeline.

**No hybridization signal** — all quartets are tree-like, so the splits graph
is a clean unrooted tree:

.. image:: /_static/img/quartet_network_tree.png
   :alt: PhyKIT quartet network tree figure
   :align: center
   :width: 80%

|

**Hybridization signal present** — hybrid quartets introduce conflicting
splits that appear as boxes in the network, indicating reticulation among
C, D, E, and F:

.. image:: /_static/img/quartet_network_hybrid.png
   :alt: PhyKIT quartet network hybrid figure
   :align: center
   :width: 80%
