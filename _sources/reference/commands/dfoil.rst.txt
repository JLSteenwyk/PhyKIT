.. _cmd-dfoil:
.. _command-dfoil:

DFOIL test (Pease & Hahn 2015)
==============================

DFOIL test for detecting and polarizing introgression in a 5-taxon symmetric phylogeny

Command identity
----------------

:Canonical command: ``dfoil``
:Handler: ``dfoil``
:Aliases: dfoil_test
:Standalone executables: pk_dfoil, pk_dfoil_test
:Categories: Introgression & gene flow

Runtime interface
-----------------

.. include:: /_generated/commands/dfoil.inc

Guidance, interpretation, and examples
--------------------------------------

Compute DFOIL statistics (Pease & Hahn 2015) for detecting and polarizing
introgression in a 5-taxon symmetric phylogeny with topology
``((P1, P2), (P3, P4), Outgroup)``.

P1 and P2 are sister taxa; P3 and P4 are sister taxa; the two pairs are
sister to each other, with an outgroup rooting the tree.

Four D-statistics are computed: DFO (far-outer), DIL (inner-left),
DFI (far-inner), and DOL (outer-left). Each is tested for significance
using a chi-squared test (1 df). The sign pattern of the four statistics
maps to a specific introgression scenario via the lookup table from
Pease & Hahn (2015).

.. code-block:: shell

   phykit dfoil -a <alignment> --p1 <taxon> --p2 <taxon> --p3 <taxon> \
       --p4 <taxon> --outgroup <taxon> [--json]

Options: |br|
*-a/--alignment*: FASTA alignment file (required) |br|
*--p1*: taxon name for P1, sister to P2 (required) |br|
*--p2*: taxon name for P2, sister to P1 (required) |br|
*--p3*: taxon name for P3, sister to P4 (required) |br|
*--p4*: taxon name for P4, sister to P3 (required) |br|
*--outgroup*: outgroup taxon name (required) |br|
*--json*: optional argument to print results as JSON
