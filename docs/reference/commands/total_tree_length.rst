.. _cmd-total_tree_length:
.. _command-total_tree_length:

Total tree length
=================

Sum of all branch lengths

Command identity
----------------

:Canonical command: ``total_tree_length``
:Handler: ``total_tree_length``
:Aliases: tree_len
:Standalone executables: pk_total_tree_length, pk_tree_len
:Categories: Tree summary statistics

Runtime interface
-----------------

.. include:: /_generated/commands/total_tree_length.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate total tree length, which is a sum of all branches.

.. code-block:: shell

   phykit total_tree_length <tree> [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*--json*: optional argument to print results as JSON
