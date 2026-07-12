.. _cmd-rename_tree_tips:
.. _command-rename_tree_tips:

Rename tree tips
================

Rename tip labels in a tree

Command identity
----------------

:Canonical command: ``rename_tree_tips``
:Handler: ``rename_tree_tips``
:Aliases: rename_tips, rename_tree
:Standalone executables: pk_rename_tree_tips, pk_rename_tips, pk_rename_tree
:Categories: Tree manipulation & utilities

Runtime interface
-----------------

.. include:: /_generated/commands/rename_tree_tips.inc

Guidance, interpretation, and examples
--------------------------------------

Renames tips in a phylogeny.

Renaming tips will follow the scheme of a tab-delimited
file wherein the first column is the current tip name and the
second column is the desired tip name in the resulting 
phylogeny. 

.. code-block:: shell

   phykit rename_tree_tips <tree> -i/--idmap <idmap.txt> [-o/--output <output_file>] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-i/--idmap*: identifier map of current tip names (col1) and desired
tip names (col2) |br|
*-o/--output*: optional argument to write the renamed tree files to. Default
output will have the same name as the input file but with the suffix ".renamed" |br|
*--json*: optional argument to print results as JSON
