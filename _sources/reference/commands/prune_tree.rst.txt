.. _cmd-prune_tree:
.. _command-prune_tree:

Prune tree
==========

Prune taxa from a tree

Command identity
----------------

:Canonical command: ``prune_tree``
:Handler: ``prune_tree``
:Aliases: prune
:Standalone executables: pk_prune_tree, pk_prune
:Categories: Tree manipulation & utilities

Runtime interface
-----------------

.. include:: /_generated/commands/prune_tree.inc

Guidance, interpretation, and examples
--------------------------------------

Prune tips from a phylogeny.

Provide a single column file with the names of the tips
in the input phylogeny you would like to prune from the
tree.

.. code-block:: shell

   phykit prune_tree <tree> <list_of_taxa> [-o/--output <output_file>] [-k/--keep]
       [--ignore-branch-labels] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*<list_of_taxa>*: single column file with the names of the tips to remove
from the phylogeny |br|
*-o/--output*: name of output file for the pruned phylogeny.
Default output will have the same name as the input file but with the suffix
".pruned" |br|
*-k/--keep*: optional argument. If used instead of pruning taxa in <list_of_taxa>,
keep them |br|
*--ignore-branch-labels*: optional argument. Strip HyPhy/aBSREL-style ``{...}``
branch labels (e.g., ``Hydlep{FG}``) from tip names when matching against
``<list_of_taxa>``. The labels are preserved on the surviving tips in the
output tree, so a list of bare names produced by tools like ``seqkit seq -n``
will keep (or remove) the foreground-tagged tips intact. |br|
*--json*: optional argument to print results as JSON

Tip: when preparing trees for branch-site selection tests such as aBSREL or
RELAX, use ``--ignore-branch-labels`` together with ``-k`` to subset a labeled
species tree to the taxa present in a given orthogroup alignment without
losing the ``{FG}`` (foreground) annotations on the labeled branches.
