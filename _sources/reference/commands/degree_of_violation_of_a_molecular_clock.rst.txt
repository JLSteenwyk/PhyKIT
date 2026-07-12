.. _cmd-degree_of_violation_of_a_molecular_clock:
.. _command-degree_of_violation_of_a_molecular_clock:

Degree of violation of the molecular clock
==========================================

Measure molecular clock violation

Command identity
----------------

:Canonical command: ``degree_of_violation_of_a_molecular_clock``
:Handler: ``dvmc``
:Aliases: dvmc
:Standalone executables: pk_degree_of_violation_of_a_molecular_clock, pk_dvmc
:Categories: Tree summary statistics

Runtime interface
-----------------

.. include:: /_generated/commands/degree_of_violation_of_a_molecular_clock.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate degree of violation of a molecular clock (or DVMC) in a phylogeny.

Lower DVMC values are thought to be desirable because they are indicative
of a lower degree of violation in the molecular clock assumption.

Typically, outgroup taxa are not included in molecular clock analysis. Thus,
prior to calculating DVMC from a single gene tree, users may want to prune
outgroup taxa from the phylogeny. To prune tips from a phylogeny, see the 
prune_tree function. 

Calculate DVMC in a tree following Liu et al., PNAS (2017), doi: 10.1073/pnas.1616744114.

.. code-block:: shell

   phykit degree_of_violation_of_a_molecular_clock <tree> [--json]

Options: |br|
*<tree>*: input file tree name |br|
*--json*: optional argument to print results as JSON
