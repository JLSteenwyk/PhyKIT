.. _cmd-treeness_over_rcv:
.. _command-treeness_over_rcv:

Treeness over RCV
=================

Treeness over relative composition variability

Command identity
----------------

:Canonical command: ``treeness_over_rcv``
:Handler: ``treeness_over_rcv``
:Aliases: tor, toverr
:Standalone executables: pk_treeness_over_rcv, pk_tor, pk_toverr
:Categories: Saturation & model adequacy

Runtime interface
-----------------

.. include:: /_generated/commands/treeness_over_rcv.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate treeness/RCV for a given alignment and tree.

Higher treeness/RCV values are thought to be desirable because
they indicate a high signal-to-noise ratio and suggest the data are least susceptible
to composition bias.

PhyKIT reports three tab delimited values:
col1: treeness/RCV
col2: treeness
col3: RCV

Calculate treeness/RCV following Phillips and Penny, Molecular 
Phylogenetics and Evolution (2003), doi: 10.1016/S1055-7903(03)00057-5.

.. code-block:: shell

   phykit treeness_over_rcv -a/--alignment <alignment> -t/--tree <tree> [--json]

Options: |br|
*-a/--alignment*: an alignment file |br|
*-t/--tree*: a tree file |br|
*--json*: optional argument to print results as JSON
