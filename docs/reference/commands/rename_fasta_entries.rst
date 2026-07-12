.. _cmd-rename_fasta_entries:
.. _command-rename_fasta_entries:

Rename FASTA entries
====================

Rename entries in a FASTA file

Command identity
----------------

:Canonical command: ``rename_fasta_entries``
:Handler: ``rename_fasta_entries``
:Aliases: rename_fasta
:Standalone executables: pk_rename_fasta_entries, pk_rename_fasta
:Categories: Alignment & dataset utilities

Runtime interface
-----------------

.. include:: /_generated/commands/rename_fasta_entries.inc

Guidance, interpretation, and examples
--------------------------------------

Renames fasta entries.

Renaming fasta entries will follow the scheme of a tab-delimited
file wherein the first column is the current fasta entry name and
the second column is the new fasta entry name in the resulting 
output alignment. Note, the input fasta file does not need to be
an alignment file.

.. code-block:: shell

	phykit rename_fasta_entries <fasta> -i/--idmap <idmap> [-o/--output <output_file>] [--json]

Options: |br|
*<fasta>*: first argument after function name should be a FASTA file |br|
*-i/--idmap*: identifier map of current FASTA names (col1) and desired FASTA names (col2) |br|
*--json*: optional argument to print results as JSON
