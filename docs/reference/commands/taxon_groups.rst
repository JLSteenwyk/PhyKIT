.. _cmd-taxon_groups:
.. _command-taxon_groups:

Taxon groups
============

Group alignment or tree files by shared taxon sets (works with both FASTA and Newick)

Command identity
----------------

:Canonical command: ``taxon_groups``
:Handler: ``taxon_groups``
:Aliases: shared_taxa, tgroups
:Standalone executables: pk_taxon_groups, pk_shared_taxa, pk_tgroups
:Categories: Alignment & dataset utilities

Runtime interface
-----------------

.. include:: /_generated/commands/taxon_groups.inc

Guidance, interpretation, and examples
--------------------------------------

Determine which tree or FASTA files share the same set of taxa.
Reads a file listing paths to gene trees or alignments and groups
them by their taxon set (exact match). Reports groups sorted by
size (largest first), with the taxa present in each group.

Useful for identifying subsets of genes with identical taxon
sampling for concatenation or comparative analysis.

.. code-block:: shell

   phykit taxon_groups -l <file> [-f trees|fasta] [--json]

Options: |br|
*-l/--list*: file listing paths to gene trees or FASTA files (one per line).
Blank lines and lines starting with # are skipped. Relative paths are resolved
relative to the list file's directory. |br|
*-f/--format*: input format: ``trees`` (Newick) or ``fasta`` (FASTA alignment).
Default: ``trees``. |br|
*--json*: optional argument to print results as JSON
