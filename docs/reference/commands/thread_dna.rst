.. _cmd-thread_dna:
.. _command-thread_dna:

Protein-to-nucleotide alignment
===============================

Thread nucleotide onto protein alignment

Command identity
----------------

:Canonical command: ``thread_dna``
:Handler: ``thread_dna``
:Aliases: p2n, pal2nal
:Standalone executables: pk_thread_dna, pk_p2n, pk_pal2nal
:Categories: Alignment & dataset utilities

Runtime interface
-----------------

.. include:: /_generated/commands/thread_dna.inc

Guidance, interpretation, and examples
--------------------------------------

Thread DNA sequence onto a protein alignment to create a
codon-based alignment. 

This function requires that input alignments be in fasta format.
Codon alignments are then printed to stdout. Note, paired
sequences are assumed to have the same name between the 
protein and nucleotide file. The order does not matter.

To thread nucleotide sequences over a trimmed amino acid
alignment, provide PhyKIT with a log file specifying which
sites have been trimmed and which have been kept. The log
file must be formatted the same as the log files outputted
by the alignment trimming toolkit ClipKIT (see -l in ClipKIT
documentation.) Details about ClipKIT can be seen here:
https://github.com/JLSteenwyk/ClipKIT.

If using a ClipKIT log file, the untrimmed protein alignment
should be provided in the -p/--protein argument.

.. code-block:: shell

   phykit thread_dna -p <file> -n <file> [-c/--clipkit-log-file <file>] [-s] [--json]

Options: |br|
*-p/--protein*: protein alignment file |br|
*-n/--nucleotide*: nucleotide sequence file |br|
*-c/--clipkit-log-file*: ClipKIT output log file. The legacy
``--clipkit_log_file`` spelling is also accepted. |br|
*-s/--stop*: if used, stop codons will be removed from the output |br|
*--json*: optional argument to print results as JSON
