.. _cmd-alignment_length_no_gaps:
.. _command-alignment_length_no_gaps:

Alignment length no gaps
========================

Alignment length excluding gapped sites

Command identity
----------------

:Canonical command: ``alignment_length_no_gaps``
:Handler: ``alignment_length_no_gaps``
:Aliases: aln_len_no_gaps, alng
:Standalone executables: pk_alignment_length_no_gaps, pk_aln_len_no_gaps, pk_alng
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/alignment_length_no_gaps.inc

Guidance, interpretation, and examples
--------------------------------------

.. image:: /_static/docs_img/aln_len_no_gaps.png 
   :align: center
   :width: 75%


Calculate alignment length excluding sites with gaps.

Longer alignments when excluding sites with gaps are
associated with a strong phylogenetic signal.

PhyKIT reports three tab delimited values:
col1: number of sites without gaps
col2: total number of sites
col3: percentage of sites without gaps

The association between alignment length when excluding sites
with gaps and phylogenetic signal was determined by Shen 
et al., Genome Biology and Evolution (2016), 
doi: 10.1093/gbe/evw179.

.. code-block:: shell

		phykit aln_len_no_gaps <alignment> [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--json*: optional argument to print results as JSON
