.. _cmd-neighbor_net:
.. _command-neighbor_net:

NeighborNet
===========

NeighborNet phylogenetic network from distance matrix (Bryant & Moulton 2004)

Command identity
----------------

:Canonical command: ``neighbor_net``
:Handler: ``neighbor_net``
:Aliases: nnet
:Standalone executables: pk_neighbor_net, pk_nnet
:Categories: Tree comparison & consensus

Runtime interface
-----------------

.. include:: /_generated/commands/neighbor_net.inc

Guidance, interpretation, and examples
--------------------------------------

Construct a NeighborNet phylogenetic network from pairwise distances
(Bryant & Moulton 2004). Unlike a consensus network (which summarizes
conflict across gene trees), NeighborNet infers a splits graph directly
from a distance matrix — analogous to Neighbor-Joining but producing a
network instead of a tree.

Input is either a FASTA alignment (distances computed internally) or a
pre-computed distance matrix CSV. The algorithm:

1. Computes pairwise distances (p-distance, identity, or Jukes-Cantor)
2. Builds a NJ tree to determine the circular taxon ordering
3. Enumerates all circular splits compatible with that ordering
4. Estimates split weights via non-negative least squares (NNLS)
5. Visualizes as a planar splits graph (Buneman graph)

.. code-block:: shell

   # From an alignment
   phykit neighbor_net -a alignment.fa -o network.pdf

   # From a pre-computed distance matrix
   phykit neighbor_net --distance-matrix distances.csv -o network.pdf

   # With Jukes-Cantor correction
   phykit neighbor_net -a alignment.fa -o network.pdf --metric jc

Options: |br|
*-a/--alignment*: FASTA alignment file (computes distances internally) |br|
*--distance-matrix*: pre-computed distance matrix CSV with taxon headers |br|
*-o/--output*: output figure path (required) |br|
*--metric*: distance metric when using -a: identity, p-distance (default), or jc (Jukes-Cantor) |br|
*--max-splits*: maximum splits for visualization (default: 30) |br|
*--json*: optional argument to print results as JSON
