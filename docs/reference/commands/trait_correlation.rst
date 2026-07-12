.. _cmd-trait_correlation:
.. _command-trait_correlation:

Trait correlation
=================

Phylogenetic correlations between all pairs of traits with heatmap visualization

Command identity
----------------

:Canonical command: ``trait_correlation``
:Handler: ``trait_correlation``
:Aliases: phylo_corr, trait_corr
:Standalone executables: pk_trait_correlation, pk_phylo_corr, pk_trait_corr
:Categories: Phylogenetic signal

Runtime interface
-----------------

.. include:: /_generated/commands/trait_correlation.inc

Guidance, interpretation, and examples
--------------------------------------

Compute phylogenetic correlations between all pairs of traits and display
them as a heatmap with significance indicators.

Uses GLS-centering via the tree's variance-covariance matrix to account for
phylogenetic non-independence among species. P-values are computed from the
t-distribution and displayed as significance stars on the heatmap cells:

- ``***`` p < 0.001
- ``**`` p < 0.01
- ``*`` p < alpha (default 0.05)

The multi-trait input file should be tab-delimited with a header row:
``taxon<tab>trait1<tab>trait2<tab>...``
At least 2 traits are required. Lines starting with '#' are treated as
comments. If the tree and trait file have different taxa, the intersection
is used and warnings are printed to stderr.

The ``--cluster`` option clusters traits by correlation similarity
(using ``1 - |r|`` as distance) and draws dendrograms alongside the
heatmap.

.. code-block:: shell

   phykit trait_correlation -t <tree> -d <trait_data> -o <output>
       [--alpha 0.05] [--cluster] [-g <gene_trees>] [--json]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>]
       [--no-title] [--title <str>]

Options: |br|
*-t/--tree*: a tree file in Newick format |br|
*-d/--trait-data*: tab-delimited multi-trait file with header row (at least 2 traits) |br|
*-o/--output*: output figure path (supports .png, .pdf, .svg) |br|
*--alpha*: significance threshold for marking p-values (default: 0.05) |br|
*--cluster*: cluster traits by correlation similarity and display dendrograms |br|
*-g/--gene-trees*: optional multi-Newick file of gene trees for discordance-aware VCV |br|
*--json*: optional argument to print results as JSON
