.. _cmd-gc_content:
.. _command-gc_content:

Guanine-cytosine (GC) content
=============================

GC content of an alignment

Command identity
----------------

:Canonical command: ``gc_content``
:Handler: ``gc_content``
:Aliases: gc
:Standalone executables: pk_gc_content, pk_gc
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/gc_content.inc

Guidance, interpretation, and examples
--------------------------------------

.. image:: /_static/docs_img/gc_content.png 
   :alt: PhyKIT gc content figure
   :align: center
   :width: 75%


Calculate GC content of a fasta file.

GC content is negatively correlated with phylogenetic signal.

If there are multiple entries, use the -v/--verbose option
to determine the GC content of each fasta entry separately.
The association between GC content and phylogenetic signal was
determined by Shen et al., Genome Biology and Evolution (2016), 
doi: 10.1093/gbe/evw179.

.. code-block:: shell

		phykit gc_content <fasta> [-v/--verbose] [--json]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-v/--verbose*: optional argument to print the GC content of each fasta
entry |br|
*--json*: optional argument to print results as JSON
