.. _cmd-variable_sites:
.. _command-variable_sites:

Variable sites
==============

Count variable sites in an alignment

Command identity
----------------

:Canonical command: ``variable_sites``
:Handler: ``variable_sites``
:Aliases: vs
:Standalone executables: pk_variable_sites, pk_vs
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/variable_sites.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate the number of variable sites in an alignment.

The number of variable sites in an alignment is
associated with a strong phylogenetic signal.

PhyKIT reports three tab delimited values:
col1: number of variable sites
col2: total number of sites
col3: percentage of variable sites

The association between the number of variable sites and
phylogenetic signal was determined by Shen et al.,
Genome Biology and Evolution (2016), 
doi: 10.1093/gbe/evw179.

.. code-block:: shell

   phykit variable_sites <alignment> [--json]

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--json*: optional argument to print results as JSON
