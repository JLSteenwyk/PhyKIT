.. _cmd-alignment_length:
.. _command-alignment_length:

Alignment length
================

Length of an input alignment

Command identity
----------------

:Canonical command: ``alignment_length``
:Handler: ``alignment_length``
:Aliases: al, aln_len
:Standalone executables: pk_alignment_length, pk_al, pk_aln_len
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/alignment_length.inc

Guidance, interpretation, and examples
--------------------------------------

.. image:: /_static/docs_img/aln_len.png 
   :align: center
   :width: 75%


Length of an input alignment is calculated using this function.

Longer alignments are associated with a strong phylogenetic signal.
   
The association between alignment length and phylogenetic signal
was determined by Shen et al., Genome Biology and Evolution (2016),
doi: 10.1093/gbe/evw179.

.. code-block:: shell

		phykit aln_len <alignment> [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--json*: optional argument to print results as JSON
