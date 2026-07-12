.. _cmd-tip_labels:
.. _command-tip_labels:

Tip labels
==========

Print tip labels of a tree

Command identity
----------------

:Canonical command: ``tip_labels``
:Handler: ``tip_labels``
:Aliases: labels, tl, tree_labels
:Standalone executables: pk_tip_labels, pk_labels, pk_tl, pk_tree_labels
:Categories: Tree manipulation & utilities

Runtime interface
-----------------

.. include:: /_generated/commands/tip_labels.inc

Guidance, interpretation, and examples
--------------------------------------

Prints the tip labels (or names) of a phylogeny.

.. code-block:: shell

   phykit tip_labels <tree> [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*--json*: optional argument to print results as JSON
