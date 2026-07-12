.. _cmd-tip_to_tip_node_distance:
.. _command-tip_to_tip_node_distance:

Tip-to-tip node distance
========================

Node distance between two tips

Command identity
----------------

:Canonical command: ``tip_to_tip_node_distance``
:Handler: ``tip_to_tip_node_distance``
:Aliases: t2t_nd, t2t_node_dist
:Standalone executables: pk_tip_to_tip_node_distance, pk_t2t_nd, pk_t2t_node_dist
:Categories: Tree summary statistics

Runtime interface
-----------------

.. include:: /_generated/commands/tip_to_tip_node_distance.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate distance between two tips (or leaves) in a phylogeny.

Distance is measured by the number of nodes between one tip
and another.

.. code-block:: shell

   phykit tip_to_tip_node_distance <tree_file> <tip_1> <tip_2> [--json]

Options: |br|
*<tree_file>*: first argument after function name should be a tree file |br|
*<tip_1>*: second argument should be the name of the first tip of interest |br|
*<tip_2>*: third argument should be the name of the second tip of interest |br|
*--json*: optional argument to print results as JSON
