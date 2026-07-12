.. _cmd-tree_space:
.. _command-tree_space:

Tree space visualization
========================

Visualize gene tree topology space via MDS, t-SNE, UMAP, or clustered distance heatmap

Command identity
----------------

:Canonical command: ``tree_space``
:Handler: ``tree_space``
:Aliases: tree_landscape, tspace
:Standalone executables: pk_tree_space, pk_tree_landscape, pk_tspace
:Categories: Tree comparison & consensus

Runtime interface
-----------------

.. include:: /_generated/commands/tree_space.inc

Guidance, interpretation, and examples
--------------------------------------

Visualize how gene trees cluster in topology space by computing pairwise
distances (Robinson-Foulds or Kuhner-Felsenstein) and projecting into
2D via MDS, t-SNE, or UMAP. Points are colored by spectral clustering
with automatic cluster detection via the eigengap heuristic.

Alternatively, use ``--heatmap`` to draw a clustered distance heatmap
with a dendrogram showing hierarchical relationships among gene trees.
Use ``--distance-matrix`` to export the raw pairwise distance matrix
as a CSV file.

Optionally highlight a species tree as a distinct marker to see where
it sits relative to the gene tree cloud.

.. code-block:: shell

   phykit tree_space -t <trees> -o <output>
       [--metric rf|kf] [--method mds|tsne|umap]
       [--species-tree <file>] [--k <int>] [--seed <int>]
       [--heatmap] [--distance-matrix <file>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-t/--trees*: file with gene trees (one Newick per line, or one file path per line) (required) |br|
*-o/--output*: output figure path (.png, .pdf, .svg) (required) |br|
*--metric*: distance metric: rf (normalized Robinson-Foulds, default) or kf (Kuhner-Felsenstein) |br|
*--method*: dimensionality reduction: mds (default), tsne, or umap |br|
*--species-tree*: optional species tree to highlight in the plot as a star marker |br|
*--k*: number of clusters (auto-detected via eigengap if omitted) |br|
*--seed*: random seed for reproducibility (t-SNE/UMAP) |br|
*--heatmap*: draw a clustered distance heatmap instead of a scatter plot |br|
*--distance-matrix*: output the pairwise distance matrix as a CSV file |br|
*--json*: optional argument to print results as JSON
