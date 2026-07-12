.. _cmd-consensus_tree:
.. _command-consensus_tree:

Consensus tree
==============

Consensus tree from multiple trees

Command identity
----------------

:Canonical command: ``consensus_tree``
:Handler: ``consensus_tree``
:Aliases: consensus, ctree
:Standalone executables: pk_consensus_tree, pk_consensus, pk_ctree
:Categories: Tree comparison & consensus

Runtime interface
-----------------

.. include:: /_generated/commands/consensus_tree.inc

Guidance, interpretation, and examples
--------------------------------------

Infer a consensus tree from a collection of trees.

Input can be either:
1) a file with one Newick tree per line, or
2) a file with one tree-file path per line.

Consensus methods:
* ``majority``: majority-rule consensus (default) |br|
* ``strict``: strict consensus

Missing taxa handling:
* ``--missing-taxa error`` (default): exits if trees do not share identical tip sets |br|
* ``--missing-taxa shared``: prunes all trees to the intersection of taxa before inferring consensus

.. code-block:: shell

   phykit consensus_tree -t/--trees <trees> [-m/--method strict|majority] [--missing-taxa error|shared] [--json]

Options: |br|
*-t/--trees*: file containing trees (one Newick per line) or tree-file paths (one per line) |br|
*-m/--method*: consensus method (``strict`` or ``majority``; default: ``majority``) |br|
*--missing-taxa*: handling strategy for mismatched taxa (``error`` or ``shared``; default: ``error``) |br|
*--json*: optional argument to print results as JSON
