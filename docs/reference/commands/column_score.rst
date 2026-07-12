.. _cmd-column_score:
.. _command-column_score:

Column score
============

Column score for alignment quality

Command identity
----------------

:Canonical command: ``column_score``
:Handler: ``column_score``
:Aliases: cs
:Standalone executables: pk_column_score, pk_cs
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/column_score.inc

Guidance, interpretation, and examples
--------------------------------------

.. image:: /_static/docs_img/column_score.png 
   :align: center
   :width: 75%


Calculates column score.

Column score is an accuracy metric for a multiple alignment relative
to a reference alignment. It is calculated by summing the correctly
aligned columns over all columns in an alignment. Thus, values range
from 0 to 1 and higher values indicate more accurate alignments.

Column score is calculated following Thompson et al., Nucleic
Acids Research (1999), doi: 10.1093/nar/27.13.2682.

.. code-block:: shell

	phykit column_score <alignment> --reference <reference_alignment> [--json]

Options: |br|
*<alignment>*: first argument after function name should be a query
fasta alignment file to be scored for accuracy |br|
*-r/--reference*: reference alignment to compare the query alignment
to |br|
*--json*: optional argument to print results as JSON
