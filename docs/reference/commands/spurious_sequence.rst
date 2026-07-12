.. _cmd-spurious_sequence:
.. _command-spurious_sequence:

Spurious homolog identification
===============================

Identify spurious sequences in alignments

Command identity
----------------

:Canonical command: ``spurious_sequence``
:Handler: ``spurious_sequence``
:Aliases: spurious_seq, ss
:Standalone executables: pk_spurious_sequence, pk_spurious_seq, pk_ss
:Categories: Homology assessment

Runtime interface
-----------------

.. include:: /_generated/commands/spurious_sequence.inc

Guidance, interpretation, and examples
--------------------------------------

Determines potentially spurious homologs using branch lengths.

Identifies potentially spurious sequences and reports
tips in the phylogeny that could possibly be removed
from the associated multiple sequence alignment. PhyKIT
does so by identifying and reporting long terminal branches
defined as branches that are equal to or greater than 20 times the median
length of all branches.

PhyKIT reports the following information
col1: name of tip that is a putatively spurious sequence
col2: length of branch leading to putatively spurious sequence
col3: threshold used to identify putatively spurious sequences
col4: median branch length in the phylogeny

If there are no putatively spurious sequences, "None" is reported.

Using this method to identify potentially spurious sequences
was, to my knowledge, first introduced by Shen et al., (2018)
Cell doi: 10.1016/j.cell.2018.10.023. 

.. code-block:: shell

   phykit spurious_seq <file> -f/--factor [--json]

Options: |br|
*<file>*: first argument after function name should be a tree file |br|
*-f/--factor*: factor to multiply median branch length by to calculate
the threshold of long branches. (Default: 20) |br|
*--json*: optional argument to print results as JSON
