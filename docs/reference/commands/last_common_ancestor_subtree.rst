.. _cmd-last_common_ancestor_subtree:
.. _command-last_common_ancestor_subtree:

Last common ancestor subtree
============================

Extract subtree from LCA of specified taxa

Command identity
----------------

:Canonical command: ``last_common_ancestor_subtree``
:Handler: ``last_common_ancestor_subtree``
:Aliases: lca_subtree
:Standalone executables: pk_last_common_ancestor_subtree, pk_lca_subtree
:Categories: Tree manipulation & utilities

Runtime interface
-----------------

.. include:: /_generated/commands/last_common_ancestor_subtree.inc

Guidance, interpretation, and examples
--------------------------------------

Obtains subtree from a phylogeny by getting the last common ancestor
from a list of taxa.

.. code-block:: shell

   phykit last_common_ancestor_subtree <tree> <list_of_taxa> [-o/--output <file>] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*<list_of_taxa>*: second argument after function name should be a single column
file with the list of taxa to get the last common ancestor subtree for |br|
*-o/--output*: optional argument to name the outputted tree file |br|
*--json*: optional argument to print results as JSON
