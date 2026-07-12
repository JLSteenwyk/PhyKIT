.. _cmd-occupancy_filter:
.. _command-occupancy_filter:

Occupancy filter
================

Filter alignments or trees by cross-file taxon occupancy (works with both FASTA and Newick)

Command identity
----------------

:Canonical command: ``occupancy_filter``
:Handler: ``occupancy_filter``
:Aliases: filter_occupancy, occ_filter
:Standalone executables: pk_occupancy_filter, pk_filter_occupancy, pk_occ_filter
:Categories: Alignment & dataset utilities

Runtime interface
-----------------

.. include:: /_generated/commands/occupancy_filter.inc

Guidance, interpretation, and examples
--------------------------------------

Filter alignments and/or trees by cross-file taxon occupancy. Given a
list of alignment or tree files, counts how many files each taxon appears
in and retains only taxa meeting a minimum threshold. Outputs filtered
copies of each input file.

This is useful for phylogenomics workflows where you want to ensure all
taxa in your dataset are present in at least N genes before concatenation
or downstream analysis. For FASTA files, sequences of removed taxa are
dropped. For tree files, tips of removed taxa are pruned.

**Example:** Given 10 alignment files and ``-t 0.5`` (the default), only
taxa present in at least 5 of the 10 alignments will be retained. New
filtered alignment files are written with the removed taxa excluded.

**Threshold interpretation:**

- Values **between 0 and 1** (inclusive) are treated as a **fraction** of the
  total number of files. For example, ``-t 0.5`` means 50% of files;
  ``-t 1.0`` means 100% (taxon must be in every file).
- Values **greater than 1** are treated as an **absolute count**. For example,
  ``-t 5`` means the taxon must appear in at least 5 files.
- The default is ``0.5`` (50% occupancy).

.. code-block:: shell

   # Keep taxa in at least 50% of files (default)
   phykit occupancy_filter -l alignment_list.txt

   # Keep taxa in all files (100% occupancy)
   phykit occupancy_filter -l alignment_list.txt -t 1.0

   # Keep taxa in at least 20 files
   phykit occupancy_filter -l alignment_list.txt -t 20

   # Filter trees instead of alignments
   phykit occ_filter -l tree_list.txt -f trees -t 0.5 -o filtered_trees/

.. code-block:: shell

   phykit occupancy_filter -l <file_list> [-f/--format fasta|trees]
       [-t/--threshold <float>] [-o/--output-dir <dir>] [--suffix <str>] [--json]

Options: |br|
*-l/--list*: file listing paths to alignment or tree files, one per line (required) |br|
*-f/--format*: input file format — ``fasta`` (default) or ``trees`` |br|
*-t/--threshold*: minimum occupancy to retain a taxon. Values between 0 and 1 (inclusive) are treated as a fraction (e.g., ``0.5`` = 50%, ``1.0`` = 100%); values > 1 are treated as an absolute count (default: ``0.5``) |br|
*-o/--output-dir*: directory for filtered output files (default: same directory as input) |br|
*--suffix*: suffix added to output filenames before the extension (default: ``.filtered``) |br|
*--json*: output results as JSON
