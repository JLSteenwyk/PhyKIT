.. _cmd-monophyly_check:
.. _command-monophyly_check:

Monophyly check
===============

Test monophyly of a group of taxa

Command identity
----------------

:Canonical command: ``monophyly_check``
:Handler: ``monophyly_check``
:Aliases: is_monophyletic
:Standalone executables: pk_monophyly_check, pk_is_monophyletic
:Categories: Tree manipulation & utilities

Runtime interface
-----------------

.. include:: /_generated/commands/monophyly_check.inc

Guidance, interpretation, and examples
--------------------------------------

This analysis can be used to determine if a set of 
taxa are exclusively monophyletic. By exclusively monophyletic,
if other taxa are in the same clade, the lineage will not be
considered exclusively monophyletic.

Requires a taxa file, which specifies which tip names
are expected to be monophyletic. File format is a
single column file with tip names. Tip names not
present in the tree will not be considered when
examining monophyly.

The output will have six columns.
col 1: if the clade was or wasn't monophyletic
col 2: average bipartition support value in the clade of interest
col 3: maximum bipartition support value in the clade of interest
col 4: minimum bipartition support value in the clade of interest
col 5: standard deviation of bipartition support values in the clade of interest
col 6: tip names of taxa monophyletic with the lineage of interest excluding those that are listed in the taxa_of_interest file

.. code-block:: shell

   phykit monophyly_check <tree> <list_of_taxa> [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*<list_of_taxa>*: single column file with list of tip names to 
examine the monophyly of |br|
*--json*: optional argument to print results as JSON
