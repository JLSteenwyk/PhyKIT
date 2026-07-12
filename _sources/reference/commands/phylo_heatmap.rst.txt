.. _cmd-phylo_heatmap:
.. _command-phylo_heatmap:

Phylogenetic heatmap
====================

Phylogeny alongside a heatmap of numeric trait values

Command identity
----------------

:Canonical command: ``phylo_heatmap``
:Handler: ``phylo_heatmap``
:Aliases: ph, pheatmap
:Standalone executables: pk_phylo_heatmap, pk_ph, pk_pheatmap
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/phylo_heatmap.inc

Guidance, interpretation, and examples
--------------------------------------

Draw a phylogenetic heatmap: a phylogeny alongside a color-coded matrix
of numeric trait values, with rows aligned to tree tips. Analogous to
R's phytools::phylo.heatmap().

.. image:: /_static/phylo_heatmap_example.png
   :alt: PhyKIT phylo heatmap example figure
   :align: center
   :width: 90%

.. code-block:: shell

	phykit phylo_heatmap -t <tree> -d <data> -o <output>
		[--split 0.3] [--standardize] [--cmap viridis] [--cluster-columns] [--json]

Options: |br|
*-t/--tree*: tree file (required) |br|
*-d/--data*: numeric data matrix in TSV format with header row (required) |br|
*-o/--output*: output figure path (required; supports .png, .pdf, .svg) |br|
*--split*: fraction of figure width for the tree panel (default: 0.3) |br|
*--standardize*: z-score each column before coloring |br|
*--cmap*: matplotlib colormap name (default: ``viridis``) |br|
*--cluster-columns*: cluster trait columns by similarity and display a dendrogram at the top |br|
*--json*: optional argument to output metadata as JSON
