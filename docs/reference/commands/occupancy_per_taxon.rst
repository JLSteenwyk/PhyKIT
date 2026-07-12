.. _cmd-occupancy_per_taxon:
.. _command-occupancy_per_taxon:

Occupancy per taxon
===================

Taxon occupancy in alignment columns

Command identity
----------------

:Canonical command: ``occupancy_per_taxon``
:Handler: ``occupancy_per_taxon``
:Aliases: occ_tax, occupancy_taxon
:Standalone executables: pk_occupancy_per_taxon, pk_occ_tax, pk_occupancy_taxon
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/occupancy_per_taxon.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate occupancy per taxon in an alignment.

Occupancy is the fraction of valid (non-gap/non-ambiguous) characters
for each taxon.

.. code-block:: shell

	phykit occupancy_per_taxon <alignment> [--json]

Example output:

.. code-block:: shell

	1	0.8333
	2	0.6667

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--json*: optional argument to print results as JSON
