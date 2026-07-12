.. _cmd-hidden_paralogy_check:
.. _command-hidden_paralogy_check:

Hidden paralogy check
=====================

Check for hidden paralogy in gene trees

Command identity
----------------

:Canonical command: ``hidden_paralogy_check``
:Handler: ``hidden_paralogy_check``
:Aliases: clan_check
:Standalone executables: pk_hidden_paralogy_check, pk_clan_check
:Categories: Homology assessment

Runtime interface
-----------------

.. include:: /_generated/commands/hidden_paralogy_check.inc

Guidance, interpretation, and examples
--------------------------------------

Scan tree for evidence of hidden paralogy.

This analysis can be used to identify hidden paralogy. 
Specifically, this method will examine if a set of
well known monophyletic taxa are, in fact, monophyletic.
If they are not, the evolutionary history of the gene may
be subject to hidden paralogy. This analysis is typically
done with single-copy orthologous genes.

Requires a clade file, which specifies which monophyletic
lineages to check for. Multiple monophyletic
lineages can be specified. Each lineage should
be specified on a single line and each tip name 
(or taxon name) should be separated by a space.
For example, if it is anticipated that tips
"A", "B", and "C" are monophyletic and "D",
"E", and "F" are expected to be monophyletic, the
clade file should be formatted as follows: |br|
" |br|
A B C |br|
D E F |br|
"

The output will report if the specified taxa were monophyletic
or not. The number of rows will reflect how many groups of taxa
were checked for monophyly. For example,
if there were three rows of clades in the -c file, there will be
three rows in the output
where the first row in the output corresponds to the 
results of the first row in the clade file. |br|

The concept behind this analysis follows
Siu-Ting et al., Molecular Biology and Evolution (2019),
doi: 10.1093/molbev/msz067.

.. code-block:: shell

   phykit hidden_paralogy_check <tree> -c/--clade <clade_file> [--json]

Options: |br|
*-t/--tree*: input file tree name |br|
*-c/--clade*: clade file detailing which monophyletic lineages should
be scanned for |br|
*--json*: optional argument to print results as JSON
