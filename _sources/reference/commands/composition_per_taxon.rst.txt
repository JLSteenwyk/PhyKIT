.. _cmd-composition_per_taxon:
.. _command-composition_per_taxon:

Composition per taxon
=====================

Nucleotide or amino acid composition per taxon

Command identity
----------------

:Canonical command: ``composition_per_taxon``
:Handler: ``composition_per_taxon``
:Aliases: comp_tax, comp_taxon
:Standalone executables: pk_composition_per_taxon, pk_comp_tax, pk_comp_taxon
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/composition_per_taxon.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate sequence composition per taxon in an alignment.

Composition is reported as semicolon-separated `symbol:frequency` pairs
for each taxon. Frequencies are calculated from valid (non-gap/non-ambiguous)
characters. Symbol order is alphabetical.

.. code-block:: shell

	phykit composition_per_taxon <alignment> [--json]

Example output:

.. code-block:: shell

	1	A:0.4;C:0.0;G:0.2;T:0.4
	2	A:0.5;C:0.0;G:0.25;T:0.25

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--json*: optional argument to print results as JSON
