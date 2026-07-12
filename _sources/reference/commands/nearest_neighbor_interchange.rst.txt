.. _cmd-nearest_neighbor_interchange:
.. _command-nearest_neighbor_interchange:

Nearest neighbor interchange
============================

Generate NNI tree rearrangements

Command identity
----------------

:Canonical command: ``nearest_neighbor_interchange``
:Handler: ``nearest_neighbor_interchange``
:Aliases: nni
:Standalone executables: pk_nearest_neighbor_interchange, pk_nni
:Categories: Tree manipulation & utilities

Runtime interface
-----------------

.. include:: /_generated/commands/nearest_neighbor_interchange.inc

Guidance, interpretation, and examples
--------------------------------------

Generate nearest neighbor interchange (NNI) moves for a binary
rooted tree.

By default, all NNI rearrangements are emitted (one tree per neighbor),
and the output file will have the same name as the input file but with
the suffix ``.nnis``. The input phylogeny is written at the top of the
output unless ``--no-input-tree`` is supplied.

Pass ``--branch`` (a comma-separated pair of taxa) or ``--branches``
(a file with one pair per line) to restrict output to the two NNI
rearrangements around a specific internal branch. The branch is
identified as the edge leading to the most recent common ancestor
(MRCA) of the supplied taxa. This is useful for branch-by-branch
likelihood comparison in IQ-TREE / RAxML: re-score the bundle of
output trees against an alignment and compare site-likelihoods at the
labelled split.

.. code-block:: shell

   phykit nearest_neighbor_interchange <tree> [-o/--output <output_file>]
       [--branch A,B] [--branches <file>] [--no-input-tree] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-o/--output*: optional argument to specify output file name |br|
*--branch*: optional argument. Two taxa separated by a comma
(e.g., ``--branch A,B``) whose MRCA defines the internal branch to
rearrange. Emits the two alternative NNI topologies for that branch. |br|
*--branches*: optional argument. Path to a file with one taxon pair per
line (comma- or tab-separated). Lines may optionally start with a label
followed by the two taxa (``label<TAB>A<TAB>B``), in which case the
label is echoed in the ``--json`` report. Lines that are blank or
start with ``#`` are ignored. Each pair yields two NNI trees, emitted
in the order the pairs appear. |br|
*--no-input-tree*: optional argument. Omit the input phylogeny from the
top of the output file. |br|
*--json*: optional argument to print summary metadata as JSON

Example: explore likelihoods around a single internal split

.. code-block:: shell

   # Emit the input tree + the 2 NNIs around the MRCA of taxa A, B
   phykit nni species.tre --branch A,B -o candidates.nwk

   # Then evaluate likelihoods in IQ-TREE
   iqtree -s alignment.fa -te candidates.nwk -pre nni_compare
