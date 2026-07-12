.. _cmd-alignment_subsample:
.. _command-alignment_subsample:

Alignment subsampling
=====================

Randomly subsample genes, partitions, or sites

Command identity
----------------

:Canonical command: ``alignment_subsample``
:Handler: ``alignment_subsample``
:Aliases: aln_subsample, subsample
:Standalone executables: pk_alignment_subsample, pk_aln_subsample, pk_subsample
:Categories: Alignment & dataset utilities

Runtime interface
-----------------

.. include:: /_generated/commands/alignment_subsample.inc

Guidance, interpretation, and examples
--------------------------------------

Randomly subsample genes, partitions, or sites from phylogenomic datasets.
Supports three modes:

- **genes**: Given a file listing alignment paths, randomly select N of them.
- **partitions**: Given a supermatrix + RAxML-style partition file, randomly
  select N partitions and extract their columns into a new alignment.
- **sites**: Given a single alignment, randomly select N columns.

Sampling can be without replacement (default, for jackknife-style tests) or
with replacement using ``--bootstrap`` (for gene/site bootstrapping).
Use ``--seed`` for reproducibility.

.. code-block:: shell

   # Subsample 50 genes from a list of 200
   phykit alignment_subsample --mode genes -l gene_list.txt --number 50 -o sub50

   # Subsample 50% of partitions from a supermatrix
   phykit alignment_subsample --mode partitions -a concat.fa -p concat.partition \
       --fraction 0.5 --seed 42 -o sub50pct

   # Bootstrap-resample sites from an alignment
   phykit alignment_subsample --mode sites -a alignment.fa \
       --number 500 --bootstrap --seed 1 -o boot1

Options: |br|
*--mode*: subsampling mode: genes, partitions, or sites (required) |br|
*-a/--alignment*: input alignment file in FASTA format (required for partitions and sites modes) |br|
*-l/--list*: file listing alignment paths, one per line (required for genes mode) |br|
*-p/--partition*: RAxML-style partition file (required for partitions mode) |br|
*--number*: exact number of items to select (mutually exclusive with --fraction) |br|
*--fraction*: fraction of items to select, 0.0 to 1.0 (mutually exclusive with --number) |br|
*--bootstrap*: sample with replacement (default: without replacement) |br|
*--seed*: random seed for reproducibility |br|
*-o/--output*: output file prefix (default: subsampled) |br|
*--json*: optional argument to print results as JSON
