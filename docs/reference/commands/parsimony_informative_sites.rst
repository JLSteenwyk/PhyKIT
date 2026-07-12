.. _cmd-parsimony_informative_sites:
.. _command-parsimony_informative_sites:

Parsimony informative sites
===========================

Count parsimony informative sites

Command identity
----------------

:Canonical command: ``parsimony_informative_sites``
:Handler: ``parsimony_informative_sites``
:Aliases: pis
:Standalone executables: pk_parsimony_informative_sites, pk_pis
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/parsimony_informative_sites.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate the number and percentage of parsimony
informative sites in an alignment.

The number of parsimony informative sites in an alignment
is associated with a strong phylogenetic signal.

PhyKIT reports three tab delimited values:
col1: number of parsimony informative sites
col2: total number of sites
col3: percentage of parsimony informative sites

The association between the number of parsimony informative
sites and phylogenetic signal was determined by Shen 
et al., Genome Biology and Evolution (2016), 
doi: 10.1093/gbe/evw179 and Steenwyk et al., PLOS Biology
(2020), doi: 10.1371/journal.pbio.3001007.

.. code-block:: shell

		phykit parsimony_informative_sites <alignment> [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--json*: optional argument to print results as JSON
