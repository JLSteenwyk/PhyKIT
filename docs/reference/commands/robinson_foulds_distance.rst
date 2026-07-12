.. _cmd-robinson_foulds_distance:
.. _command-robinson_foulds_distance:

Robinson-Foulds distance
========================

Topological distance between trees

Command identity
----------------

:Canonical command: ``robinson_foulds_distance``
:Handler: ``rf_distance``
:Aliases: rf, rf_dist, rf_distance
:Standalone executables: pk_robinson_foulds_distance, pk_rf, pk_rf_dist, pk_rf_distance
:Categories: Tree comparison & consensus

Runtime interface
-----------------

.. include:: /_generated/commands/robinson_foulds_distance.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate Robinson-Foulds (RF) distance between two trees.

Low RF distances reflect greater similarity between two phylogenies. 
This function prints out two values, the plain RF value and the
normalized RF value, which are separated by a tab. Normalized RF values
are calculated by taking the plain RF value and dividing it by 2(n-3)
where n is the number of tips in the phylogeny. Prior to calculating
an RF value, PhyKIT will first determine the number of shared tips
between the two input phylogenies and prune them to a common set of
tips. Thus, users can input trees with different topologies and 
infer an RF value among subtrees with shared tips.

PhyKIT will print out 
col 1: the plain RF distance and 
col 2: the normalized RF distance.

RF distances are calculated following Robinson & Foulds, Mathematical 
Biosciences (1981), doi: 10.1016/0025-5564(81)90043-2.

.. code-block:: shell

   phykit robinson_foulds_distance <tree_file_zero> <tree_file_one> [--json]

Options: |br|
*<tree_file_zero>*: first argument after function name should be a tree file |br|
*<tree_file_one>*: second argument after function name should be a tree file |br|
*--json*: optional argument to print results as JSON
