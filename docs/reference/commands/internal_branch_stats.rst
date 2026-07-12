.. _cmd-internal_branch_stats:
.. _command-internal_branch_stats:

Internal branch statistics
==========================

Summary statistics of internal branch lengths

Command identity
----------------

:Canonical command: ``internal_branch_stats``
:Handler: ``internal_branch_stats``
:Aliases: ibs
:Standalone executables: pk_internal_branch_stats, pk_ibs
:Categories: Tree summary statistics

Runtime interface
-----------------

.. include:: /_generated/commands/internal_branch_stats.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate summary statistics for internal branch lengths in a phylogeny.

Internal branch lengths can be useful for phylogeny diagnostics.

To obtain all internal branch lengths, use the -v/--verbose option.   

.. code-block:: shell

   phykit internal_branch_stats <tree> [-v/--verbose] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-v/--verbose*: optional argument to print all internal branch lengths |br|
*--json*: optional argument to print results as JSON
