.. _cmd-transfer_annotations:
.. _command-transfer_annotations:

Transfer annotations
====================

Transfer node annotations between trees (e.g., wASTRAL to RAxML/IQ-TREE)

Command identity
----------------

:Canonical command: ``transfer_annotations``
:Handler: ``transfer_annotations``
:Aliases: annotate_tree, transfer_annot
:Standalone executables: pk_transfer_annotations, pk_annotate_tree, pk_transfer_annot
:Categories: Tree manipulation & utilities

Runtime interface
-----------------

.. include:: /_generated/commands/transfer_annotations.inc

Guidance, interpretation, and examples
--------------------------------------

Transfer internal node annotations from one tree onto another. Matches
nodes by bipartition (the set of descendant taxa) and copies the
annotation labels from the source to the target tree.

**Typical use case:** Transfer wASTRAL support annotations (q1/q2/q3,
pp1, f1, localPP, etc.) from a wASTRAL ``--support 3`` output tree onto
a branch-length-optimized topology from RAxML-NG, IQ-TREE, or any other
tool. The output tree has the target's branch lengths with the source's
annotations, ready for use with ``quartet_pie`` or other visualization
tools.

This function works with any Newick node annotations — not just wASTRAL.
Any labels or comments attached to internal nodes in the source tree will
be transferred to the matching nodes in the target tree (e.g., ASTRAL
``-t 2`` annotations, bootstrap support values, Bayesian posterior
probabilities, or custom labels).

**Workflow:**

1. Run wASTRAL to get an annotated species tree (mode 1 or 4, ``--support 3``)
2. Extract the unannotated s0 topology
3. Optimize branch lengths on s0 with RAxML-NG or IQ-TREE
4. Transfer annotations from step 1 onto the optimized tree from step 3:

.. code-block:: shell

   phykit transfer_annotations \
       --source wastral_annotated.tre \
       --target raxml_optimized.tre \
       -o combined.tre

   # Now use the combined tree with quartet_pie
   phykit qpie -t combined.tre -o concordance.png --branch-labels

.. code-block:: shell

   phykit transfer_annotations --source <annotated_tree> --target <branch_length_tree>
       [-o/--output <file>] [--json]

Options: |br|
*--source*: annotated tree file, e.g., wASTRAL output with ``--support 3`` or any tree with node labels/comments (required) |br|
*--target*: target tree file with branch lengths to keep, e.g., RAxML-NG or IQ-TREE output (required) |br|
*-o/--output*: output file for the annotated tree (default: target file + ``.annotated``) |br|
*--json*: output results as JSON
