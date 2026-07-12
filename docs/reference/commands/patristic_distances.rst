.. _cmd-patristic_distances:
.. _command-patristic_distances:

Patristic distances
===================

Pairwise patristic distances between taxa

Command identity
----------------

:Canonical command: ``patristic_distances``
:Handler: ``patristic_distances``
:Aliases: pd
:Standalone executables: pk_patristic_distances, pk_pd
:Categories: Tree summary statistics

Runtime interface
-----------------

.. include:: /_generated/commands/patristic_distances.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate summary statistics among patristic distances in a phylogeny.

Patristic distances are all tip-to-tip distances in a phylogeny.

To obtain all patristic distances, use the -v/--verbose option.
With the -v option, the first column will have two taxon names
separated by a '-' followed by the patristic distance. Features
will be tab separated. 

.. code-block:: shell

   phykit patristic_distances <tree> [-v/--verbose] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-v/--verbose*: optional argument to print all tip-to-tip distances |br|
*--json*: optional argument to print results as JSON
