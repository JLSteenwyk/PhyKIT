.. _cmd-faidx:
.. _command-faidx:

Faidx
=====

Extract entries from FASTA files

Command identity
----------------

:Canonical command: ``faidx``
:Handler: ``faidx``
:Aliases: ge, get_entry
:Standalone executables: pk_faidx, pk_ge, pk_get_entry
:Categories: Alignment & dataset utilities

Runtime interface
-----------------

.. include:: /_generated/commands/faidx.inc

Guidance, interpretation, and examples
--------------------------------------

.. image:: /_static/docs_img/faidx.png 
   :align: center
   :width: 75%


Extracts a sequence entry from a fasta file.

This function works similarly to the faidx function 
in samtools, but does not require an indexing step.

To obtain multiple entries, input multiple entries separated
by a comma (,). For example, if you want entries 
named "seq_0" and "seq_1", the string "seq_0,seq_1"
should be associated with the -e argument.

.. code-block:: shell

	phykit faidx <fasta> -e/--entry <fasta entry> [--json]

Options: |br|
*<fasta>*: first argument after function name should be a fasta file |br|
*-e/--entry*: entry name to be extracted from the inputted fasta file |br|
*--json*: optional argument to print results as JSON
