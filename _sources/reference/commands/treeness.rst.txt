.. _cmd-treeness:
.. _command-treeness:

Treeness
========

Ratio of internal to total branch lengths

Command identity
----------------

:Canonical command: ``treeness``
:Handler: ``treeness``
:Aliases: tness
:Standalone executables: pk_treeness, pk_tness
:Categories: Tree summary statistics

Runtime interface
-----------------

.. include:: /_generated/commands/treeness.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate treeness statistic for a phylogeny.

Higher treeness values are thought to be desirable because they
represent a higher signal-to-noise ratio.

Treeness is the sum of internal branch lengths divided by the total
tree length. Therefore, values range from 0 to 1. Treeness can be
used as a measure of the signal-to-noise ratio in a phylogeny. 

Calculate treeness (also referred to as stemminess) following
Lanyon, The Auk (1988), doi: 10.1093/auk/105.3.565 and
Phillips and Penny, Molecular Phylogenetics and Evolution
(2003), doi: 10.1016/S1055-7903(03)00057-5.

.. code-block:: shell

   phykit treeness <tree> [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*--json*: optional argument to print results as JSON
