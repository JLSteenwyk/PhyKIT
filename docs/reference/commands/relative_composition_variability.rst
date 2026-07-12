.. _cmd-relative_composition_variability:
.. _command-relative_composition_variability:

Relative composition variability
================================

Composition variability across taxa

Command identity
----------------

:Canonical command: ``relative_composition_variability``
:Handler: ``rcv``
:Aliases: rcv, rel_comp_var
:Standalone executables: pk_relative_composition_variability, pk_rcv, pk_rel_comp_var
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/relative_composition_variability.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate RCV (relative composition variability) for an alignment.

Lower RCV values are thought to be desirable because they represent
a lower composition bias in an alignment. Statistically, RCV describes
the average variability in sequence composition among taxa. 

RCV is calculated following Phillips and Penny, Molecular Phylogenetics
and Evolution (2003), doi: 10.1016/S1055-7903(03)00057-5.

RCV calculations are case-insensitive. Gap and ambiguous characters are
excluded from composition counts and correction terms, and each taxon is
normalized by its valid (non-excluded) sequence length.

.. code-block:: shell

		phykit relative_composition_variability <alignment> [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--json*: optional argument to print results as JSON
