.. _cmd-dstatistic:
.. _command-dstatistic:

D-statistic (ABBA-BABA test)
============================

Patterson's D-statistic for detecting introgression

Command identity
----------------

:Canonical command: ``dstatistic``
:Handler: ``dstatistic``
:Aliases: abba_baba, dstat
:Standalone executables: pk_dstatistic, pk_abba_baba, pk_dstat
:Categories: Introgression & gene flow

Runtime interface
-----------------

.. include:: /_generated/commands/dstatistic.inc

Guidance, interpretation, and examples
--------------------------------------

Compute Patterson's D-statistic (ABBA-BABA test) for detecting
introgression or gene flow between taxa. Given a four-taxon alignment
with topology ``(((P1, P2), P3), Outgroup)``, the test counts ABBA
and BABA site patterns:

- **ABBA**: P1 has the ancestral allele, P2 and P3 share the derived
  allele — suggests introgression between P2 and P3
- **BABA**: P2 has the ancestral allele, P1 and P3 share the derived
  allele — suggests introgression between P1 and P3
- **D = (ABBA - BABA) / (ABBA + BABA)**: D = 0 under ILS alone;
  D significantly different from 0 indicates introgression

Note: the D-statistic identifies which pair of lineages exchanged
genetic material but cannot determine the direction of gene flow
within that pair.

Significance is assessed via block jackknife (Green et al. 2010),
producing a Z-score and p-value.

Two input modes are supported:

1. **Site patterns** (``-a``): counts ABBA/BABA from a FASTA alignment.
   Only biallelic sites with no gaps or ambiguous characters are considered.
   Significance via block jackknife.

2. **Gene trees** (``-g``): counts discordant quartet topologies from gene
   trees. Gene trees can have any number of taxa; only the quartet induced
   by the four specified taxa is evaluated. Significance via chi-squared test.

.. code-block:: shell

   # Site-pattern mode
   phykit dstatistic -a <alignment> --p1 <taxon> --p2 <taxon> --p3 <taxon> \
       --outgroup <taxon> [--block-size 100] [--json]

   # Gene-tree mode
   phykit dstatistic -g <gene_trees> --p1 <taxon> --p2 <taxon> --p3 <taxon> \
       --outgroup <taxon> [--json]

Options: |br|
*-a/--alignment*: FASTA alignment file (site-pattern mode) |br|
*-g/--gene-trees*: gene trees file, one Newick per line (gene-tree mode; trees can have any number of taxa) |br|
*--p1*: taxon name for P1, sister to P2 (required) |br|
*--p2*: taxon name for P2, sister to P1, potential recipient of gene flow from P3 (required) |br|
*--p3*: taxon name for P3, donor lineage (required) |br|
*--outgroup*: outgroup taxon name, determines ancestral allele (required) |br|
*--block-size*: block size in sites for jackknife standard error estimation (default: 100; alignment mode only) |br|
*--support*: minimum branch support threshold; gene tree branches below this value are collapsed and treated as unresolved (gene-tree mode only) |br|
*--json*: optional argument to print results as JSON
