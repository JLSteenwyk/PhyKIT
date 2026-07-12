.. _cmd-long_branch_score:
.. _command-long_branch_score:

Long branch score
=================

Identify long branches in a tree

Command identity
----------------

:Canonical command: ``long_branch_score``
:Handler: ``lb_score``
:Aliases: lb_score, lbs
:Standalone executables: pk_long_branch_score, pk_lb_score, pk_lbs
:Categories: Tree summary statistics

Runtime interface
-----------------

.. include:: /_generated/commands/long_branch_score.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate long branch (LB) scores in a phylogeny.

Lower LB scores are thought to be desirable because
they are indicative of taxa or trees that likely do
not have issues with long branch attraction.

LB score is the mean pairwise patristic distance of
taxon i compared to all other taxa divided by the average
pairwise patristic distance. 

PhyKIT reports summary statistics. To obtain LB scores
for each taxon, use the -v/--verbose option. 

LB scores are calculated following Struck, Evolutionary 
Bioinformatics (2014), doi: 10.4137/EBO.S14239.  

.. code-block:: shell

   phykit long_branch_score <tree> [-v/--verbose] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-v/--verbose*: optional argument to print all LB score values |br|
*--json*: optional argument to print results as JSON
