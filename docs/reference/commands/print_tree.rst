.. _cmd-print_tree:
.. _command-print_tree:

Print tree
==========

Print ASCII representation of a tree

Command identity
----------------

:Canonical command: ``print_tree``
:Handler: ``print_tree``
:Aliases: print, pt
:Standalone executables: pk_print_tree, pk_print, pk_pt
:Categories: Tree manipulation & utilities

Runtime interface
-----------------

.. include:: /_generated/commands/print_tree.inc

Guidance, interpretation, and examples
--------------------------------------

Print ascii tree of input phylogeny.

Phylogeny can be printed with or without branch lengths.
By default, the phylogeny will be printed with branch lengths
but branch lengths can be removed using the -r/--remove argument.

.. code-block:: shell

   phykit print_tree <tree> [-r/--remove] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-r/--remove*: optional argument to print the phylogeny without branch
lengths |br|
*--json*: optional argument to print results as JSON
