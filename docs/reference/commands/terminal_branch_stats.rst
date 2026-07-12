.. _cmd-terminal_branch_stats:
.. _command-terminal_branch_stats:

Terminal branch statistics
==========================

Summary statistics of terminal branch lengths

Command identity
----------------

:Canonical command: ``terminal_branch_stats``
:Handler: ``terminal_branch_stats``
:Aliases: tbs
:Standalone executables: pk_terminal_branch_stats, pk_tbs
:Categories: Tree summary statistics

Runtime interface
-----------------

.. include:: /_generated/commands/terminal_branch_stats.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate summary statistics for terminal branch lengths in a phylogeny.

Terminal branch lengths can be useful for phylogeny diagnostics.

To obtain all terminal branch lengths, use the -v/--verbose option.   

.. code-block:: shell

   phykit terminal_branch_stats <tree> [-v/--verbose] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-v/--verbose*: optional argument to print all terminal branch lengths |br|
*--json*: optional argument to print results as JSON
