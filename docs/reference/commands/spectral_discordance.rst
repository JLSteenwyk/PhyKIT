.. _cmd-spectral_discordance:
.. _command-spectral_discordance:

Spectral discordance decomposition
==================================

PCA ordination and clustering of gene tree topologies

Command identity
----------------

:Canonical command: ``spectral_discordance``
:Handler: ``spectral_discordance``
:Aliases: sd, spec_disc
:Standalone executables: pk_spectral_discordance, pk_sd, pk_spec_disc
:Categories: Tree comparison & consensus

Runtime interface
-----------------

.. include:: /_generated/commands/spectral_discordance.inc

Guidance, interpretation, and examples
--------------------------------------

Decompose gene tree space via PCA on a bipartition presence/absence (or
branch-length) matrix, with spectral clustering and automatic cluster
detection via the eigengap heuristic.

Each gene tree is encoded as a vector over the union of all observed
bipartitions. PCA reveals the axes of topological variation, with loading
vectors directly identifying which bipartitions drive each PC. Spectral
clustering groups genes sharing similar topologies; the number of clusters
is auto-detected from the eigengap of the graph Laplacian.

Two metrics are available:

- **nrf** (default): binary presence/absence, consistent with normalized Robinson-Foulds distance
- **wrf**: branch-length weighted

Polytomies (collapsed branches) in gene trees are handled conservatively:
splits from polytomous nodes are excluded from the bipartition matrix since
they represent unresolved relationships. Trifurcating roots (standard
unrooted Newick) are not affected. This allows gene trees with collapsed
low-support branches to be used directly as input.

.. code-block:: shell

   phykit spectral_discordance -g <gene_trees> [-t <tree>] [--metric nrf|wrf] [--clusters K] [--n-pcs N] [--top-loadings N] [--plot <prefix>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-g/--gene-trees*: file of gene trees (one Newick per line, or file of filenames) |br|
*-t/--tree*: species tree (optional; flags species-tree bipartitions in loadings) |br|
*--metric*: distance metric: ``nrf`` or ``wrf`` (default: ``nrf``) |br|
*--clusters*: override auto-detected number of clusters |br|
*--n-pcs*: number of PCs to report (default: min(10, G-1)) |br|
*--top-loadings*: top bipartitions per PC to display (default: 5) |br|
*--plot*: output prefix for plots (generates ``_scatter.png`` and ``_eigengap.png``) |br|
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
*--json*: output results as JSON

.. image:: /_static/img/spectral_discordance_example_scatter.png
   :align: center
   :width: 80%

|

.. image:: /_static/img/spectral_discordance_example_eigengap.png
   :align: center
   :width: 80%
