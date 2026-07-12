.. _cmd-faiths_pd:
.. _command-faiths_pd:

Faith's phylogenetic diversity
==============================

Sum of branch lengths spanning a community of tips

Command identity
----------------

:Canonical command: ``faiths_pd``
:Handler: ``faiths_pd``
:Aliases: faith_pd, fpd, phylo_diversity
:Standalone executables: pk_faiths_pd, pk_faith_pd, pk_fpd, pk_phylo_diversity
:Categories: Tree summary statistics

Runtime interface
-----------------

.. include:: /_generated/commands/faiths_pd.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate Faith's phylogenetic diversity (PD) for a community of tips
on a phylogeny (Faith, *Biological Conservation*, 1992;
doi: 10.1016/0006-3207(92)91201-3).

Faith's PD is the sum of branch lengths in the minimum subtree that
connects a set of taxa. By default, the path from the community's most
recent common ancestor up to the tree root is included, matching the
original definition and ``picante::pd(..., include.root = TRUE)``. Pass
``--exclude-root`` to sum only the branches of the induced subtree
rooted at the MRCA, matching ``picante::pd(..., include.root = FALSE)``.

The community is supplied as a single-column text file of tip labels,
one per line.

.. code-block:: shell

   phykit faiths_pd <tree> -t/--taxa <taxa_file> [--exclude-root] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-t/--taxa*: file with one tip label per line defining the community |br|
*--exclude-root*: sum only branches of the induced subtree rooted at the community MRCA |br|
*--json*: optional argument to print results as JSON
