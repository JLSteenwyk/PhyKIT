.. _cmd-sum_of_pairs_score:
.. _command-sum_of_pairs_score:

Sum-of-pairs score
==================

Sum-of-pairs alignment quality score

Command identity
----------------

:Canonical command: ``sum_of_pairs_score``
:Handler: ``sum_of_pairs_score``
:Aliases: sop, sops
:Standalone executables: pk_sum_of_pairs_score, pk_sop, pk_sops
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/sum_of_pairs_score.inc

Guidance, interpretation, and examples
--------------------------------------

Calculates sum-of-pairs score.

Sum-of-pairs is an accuracy metric for a multiple alignment relative
to a reference alignment. It is calculated by summing the correctly
aligned residue pairs over all pairs of sequences. Thus, values range
from 0 to 1 and higher values indicate more accurate alignments.

Sum-of-pairs score is calculated following Thompson et al., Nucleic
Acids Research (1999), doi: 10.1093/nar/27.13.2682.

.. code-block:: shell

	phykit sum_of_pairs_score <alignment> --reference <reference_alignment> [--json]

Options: |br|
*<alignment>*: first argument after function name should be a query
fasta alignment file to be scored for accuracy |br|
*-r/--reference*: reference alignment to compare the query alignment
to |br|
*--json*: optional argument to print results as JSON
