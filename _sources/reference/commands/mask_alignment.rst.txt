.. _cmd-mask_alignment:
.. _command-mask_alignment:

Mask alignment
==============

Mask sites in an alignment

Command identity
----------------

:Canonical command: ``mask_alignment``
:Handler: ``mask_alignment``
:Aliases: mask, mask_aln
:Standalone executables: pk_mask_alignment, pk_mask, pk_mask_aln
:Categories: Alignment & dataset utilities

Runtime interface
-----------------

.. include:: /_generated/commands/mask_alignment.inc

Guidance, interpretation, and examples
--------------------------------------

Mask alignment sites based on threshold criteria.

Sites are retained when they pass all active thresholds:
maximum gap fraction, minimum occupancy, and maximum site entropy.

.. code-block:: shell

	phykit mask_alignment <alignment> [-g/--max_gap <float>] [-o/--min_occupancy <float>] [-e/--max_entropy <float>] [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*-g/--max_gap*: maximum allowed fraction of missing/invalid characters at a site (default: 1.0) |br|
*-o/--min_occupancy*: minimum required occupancy at a site (default: 0.0) |br|
*-e/--max_entropy*: maximum allowed site entropy (default: no filter) |br|
*--json*: optional argument to print results as JSON
