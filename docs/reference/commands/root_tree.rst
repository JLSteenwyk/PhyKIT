.. _cmd-root_tree:
.. _command-root_tree:

Root tree
=========

Root or reroot a tree

Command identity
----------------

:Canonical command: ``root_tree``
:Handler: ``root_tree``
:Aliases: root, rt
:Standalone executables: pk_root_tree, pk_root, pk_rt
:Categories: Tree manipulation & utilities

Runtime interface
-----------------

.. include:: /_generated/commands/root_tree.inc

Guidance, interpretation, and examples
--------------------------------------

Roots phylogeny using user-specified taxa.

A list of taxa to root the phylogeny on should be specified using the -r
argument. The root_taxa file should be a single-column file with taxa names.
The output file will have the same name as the input tree file but with
the suffix ".rooted".

.. code-block:: shell

   phykit root_tree <tree> -r/--root <root_taxa> [-o/--output <output_file>] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file to root |br|
*-r/--root*: single column file with taxa names to root the phylogeny on |br|
*-o/--output*: optional argument to specify the name of the output file |br|
*--json*: optional argument to print results as JSON
