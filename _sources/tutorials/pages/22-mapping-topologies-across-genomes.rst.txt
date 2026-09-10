.. _tutorial-22:

Tutorial 22: Mapping topologies across genomes
==================================================

Objectives
----------

Map six synthetic gene trees to two chromosomes and distinguish ctenophore-sister,
sponge-sister, and a ctenophore-plus-sponge clade. These examples demonstrate
classification; they are not empirical evidence about animal relationships.

Prerequisites and working directory
-----------------------------------

Install PhyKIT with ``topology_landscape`` support. Download
:download:`the complete example </data/topology_landscape_tutorial.tar.gz>`
into a new working directory, then unpack it:

.. code-block:: shell

   tar -xzf topology_landscape_tutorial.tar.gz
   cd topology_landscape

The manifest ``genes.tsv`` links gene IDs to relative Newick paths.
``reference.bed`` supplies zero-based half-open locations; ``groups.json``
defines Outgroup, Ctenophores, Sponges, and OtherAnimals. No sequence alignments
or externally inferred roots are required.

Related command references
--------------------------

- :doc:`Topology landscape </reference/commands/topology_landscape>`
- :doc:`Quartet pie </reference/commands/quartet_pie>`

Workflow
--------

Generate physical neighborhoods and chromosome plots:

.. code-block:: shell

   phykit topology_landscape --manifest genes.tsv --groups groups.json \
       --coordinates example bed reference.bed --outgroup Outgroup \
       --labels Ctenophore-sister Sponge-sister "Ctenophore + sponge" \
       --window-bp 500 --output-prefix animal_root --plot --json > animal_root.json

There are six genes overall: one of each of the three resolved topologies,
one unresolved, one missing a required group, and one with incompatible group
relationships. Thus the resolved denominator is three, not six.
The first chr1 window has four genes: three resolved and one unresolved.
Each topology has proportion 1/3 among those three resolved genes.
The second chr1 window is empty. The final chr1 window contains the
insufficiently sampled gene, and chr2 contains the incompatible gene.

Use fixed numbers of genes and PDF figures:

.. code-block:: shell

   phykit topomap --manifest genes.tsv --groups groups.json \
       --coordinates example bed reference.bed --window-genes 2 \
       --output-prefix gene_windows --plot --plot-output gene_windows.pdf

Select the reference edge separating outgroup/ctenophore from sponge/other:

.. code-block:: shell

   phykit topomap --manifest genes.tsv --reference-tree ctenophore.tre \
       --branch-taxa outgroup ctenophore \
       --coordinates example bed reference.bed --output-prefix reference_edge

Reference-edge mode derives its own four groups. Here ctenophore2 is outside
the reference tree, so this example does not have the same group membership
or incompatible-gene count as the explicit-groups example. Always inspect the
reported groups before comparing analyses.

Expected artifacts
------------------

- ``animal_root.genes.tsv``, ``animal_root.neighborhoods.tsv``, and
  ``animal_root.overall.tsv`` contain gene and neighborhood counts.
- ``animal_root.diagnostics.json`` reports matching and deduplication diagnostics.
- ``animal_root.json`` records groups, hypothesis labels, classifications, and paths.
- ``animal_root.1.png`` and ``animal_root.2.png`` are the two chromosome figures.

Troubleshooting
---------------

- Gene IDs must match exactly; transcript IDs are not automatically gene IDs.
- With ``--min-support``, unlabelled edges collapse unless another missing-support
  policy is explicitly selected. Use ``--support-scale 1`` for fractional labels.
- ``--interval`` selects gene midpoints, not all genes overlapping a region.
- Repeated coordinate rows merge to one span per gene; cross-chromosome mappings
  within a reference are rejected. Different references require different names.
- Concordance fractions are not statistical confidence, and local conflict does
  not by itself establish introgression or recombination breakpoints.

Related topology mapping includes Twisst, `Martin and Van Belleghem (2017)
<https://doi.org/10.1534/genetics.116.194720>`_. PhyKIT's strict whole-gene
classification is not Twisst's fractional subtree weighting.
