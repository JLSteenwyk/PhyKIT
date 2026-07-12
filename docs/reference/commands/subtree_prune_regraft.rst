.. _cmd-subtree_prune_regraft:
.. _command-subtree_prune_regraft:

Subtree pruning and regrafting
==============================

Generate all SPR rearrangements for a specified subtree

Command identity
----------------

:Canonical command: ``subtree_prune_regraft``
:Handler: ``subtree_prune_regraft``
:Aliases: spr
:Standalone executables: pk_subtree_prune_regraft, pk_spr
:Categories: Tree manipulation & utilities

Runtime interface
-----------------

.. include:: /_generated/commands/subtree_prune_regraft.inc

Guidance, interpretation, and examples
--------------------------------------

Generate all possible SPR (Subtree Pruning and Regrafting) rearrangements
for a specified subtree on a tree.

The subtree is identified by specifying one or more comma-separated taxa
whose MRCA (most recent common ancestor) defines the clade to prune. The
pruned subtree is then regrafted onto every other branch in the remaining
tree, producing one Newick tree per regraft position.

For a subtree defined by a single taxon (a single leaf), the leaf is pruned
and reattached to each branch in the remaining tree. For a clade (two or
more taxa), the entire clade is pruned and regrafted. Each output Newick
contains all original taxa.

.. code-block:: shell

   phykit subtree_prune_regraft -t <tree> --subtree <taxa> [-o/--output <output_file>] [--json]

Options: |br|
*-t/--tree*: input tree file in Newick format (required) |br|
*--subtree*: comma-separated list of taxa defining the subtree to prune (MRCA resolved), or a single-column file with one taxon per line (required) |br|
*-o/--output*: optional output file for SPR trees (one Newick per line); if omitted, prints to stdout |br|
*--json*: optional argument to print summary metadata and trees as JSON
