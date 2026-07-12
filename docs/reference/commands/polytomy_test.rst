.. _cmd-polytomy_test:
.. _command-polytomy_test:

Polytomy testing
================

Test for polytomies in a tree

Command identity
----------------

:Canonical command: ``polytomy_test``
:Handler: ``polytomy_test``
:Aliases: polyt, polyt_test, ptt
:Standalone executables: pk_polytomy_test, pk_polyt, pk_polyt_test, pk_ptt
:Categories: Tree comparison & consensus

Runtime interface
-----------------

.. include:: /_generated/commands/polytomy_test.inc

Guidance, interpretation, and examples
--------------------------------------

Conduct a polytomy test for three clades in a phylogeny.

Polytomy tests can be used to identify putative radiations
as well as identify well supported alternative topologies.

The polytomy testing function takes as input a file with
the three groups of taxa to test the relationships for and
a single column file with the names of the desired tree files
to use for polytomy testing. Next, the script examines
support for the grouping of the three taxa using triplets
and gene support frequencies. 

This function can account for uncertainty in gene trees - 
that is, the input phylogenies can have collapsed bipartitions.

Thereafter, a chi-squared test is conducted to determine if there
is evidence to reject the null hypothesis wherein the null 
hypothesis is that the three possible topologies among the three
groups are equally supported. This test is done using gene support
frequencies.

.. code-block:: shell

   phykit polytomy_test -t/--trees <trees> -g/--groups <groups> [--json]

Options: |br|
*-t/--trees <trees>*: single column file with the names of 
phylogenies to use for polytomy testing |br|
*-g/--groups*: a tab-delimited file with the grouping designations
to test. Lines starting with comments are not considered. Names of
individual taxa should be separated by a semi-colon ';' |br|
*--json*: optional argument to print results as JSON

For example, the groups file could look like the following:

.. code-block:: shell

   #label group0  group1  group2
   name_of_test    tip_name_A;tip_name_B   tip_name_C  tip_name_D;tip_name_E
