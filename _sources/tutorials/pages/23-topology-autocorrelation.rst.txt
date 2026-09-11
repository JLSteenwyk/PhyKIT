.. _tutorial-23:

23. Topology autocorrelation
============================

Objectives
----------

Compare the spatial organization of identical topology totals, interpret
excess agreement, and examine pointwise uncertainty across block sizes.

Prerequisites and working directory
-----------------------------------

Use a PhyKIT installation containing ``topology_autocorrelation``. This command
is under development on ``main`` and is not part of the 2.7.0 PyPI release.
Download :download:`the synthetic example archive </data/topology_autocorrelation_tutorial.tar.gz>`
into an empty working directory, then extract it:

.. code-block:: shell

   tar -xzf topology_autocorrelation_tutorial.tar.gz
   cd topology_autocorrelation

``clustered.tsv`` and ``interleaved.tsv`` contain the same 12 gene positions
and four genes per topology. In the first, topology labels occur in three
contiguous runs; in the second, they alternate. ``simulated.tsv`` contains
10,000 synthetic genes on two chromosomes, a spatial reset process with
250-bp persistence scale, and 20% random unresolved labels (seed 712).
These are illustrative labels, not empirical evidence about animal roots.

Related command references
--------------------------

- :doc:`/reference/commands/topology_autocorrelation`
- :doc:`/reference/commands/topology_landscape`

Workflow
--------

1. Analyze the clustered and interleaved examples with identical bins:

.. code-block:: shell

   phykit topology_autocorrelation --classified-genes clustered.tsv --distance-edges 0 101 201 401 801 --output-prefix clustered --plot --json
   phykit topology_autocorrelation --classified-genes interleaved.tsv --distance-edges 0 101 201 401 801 --output-prefix interleaved --plot --json

For the first bin, both have 11 pairs and a baseline of ``3/11``. The
clustered example has observed agreement ``9/11`` and excess ``6/11``.
The interleaved example has observed agreement zero and excess ``-3/11``.
Thus equal gene totals do not imply equal spatial organization.
No confidence intervals are calculated without block sizes.

2. Request uncertainty on the larger simulated example:

.. code-block:: shell

   phykit topo_ac --classified-genes simulated.tsv --distance-edges 0 201 501 --block-sizes 5000 10000 --replicates 999 --seed 7 --output-prefix simulated --plot --json

Inspect ``simulated.uncertainty.tsv`` for separate size-specific intervals,
status, withholding reasons, and sensitivity flags. The figure overlays
available pointwise intervals; it is not a simultaneous confidence band.
Inspect ``simulated.blocks.tsv`` for uneven sampling and unresolved counts.
A known simulation persistence parameter is not a clustering-range estimate
from this command.

3. For raw trees and coordinates, use the inputs prepared in
:doc:`22-mapping-topologies-across-genomes` and substitute this command for
``topology_landscape``. Both share the same classification rules. Alternatively,
pass that tutorial's ``animal_root.genes.tsv`` to ``--classified-genes``.

Expected artifacts
------------------

Each prefix produces ``.autocorrelation.tsv``, ``.chromosomes.tsv``,
``.genes.tsv``, ``.uncertainty.tsv``, ``.blocks.tsv``, ``.json``, and a PNG plot.
The first two runs produce header-only uncertainty/block TSVs; the simulated
run includes resampling results. JSON contains parameters and diagnostic counts.

Troubleshooting
---------------

- Unknown command: install the development version containing this function.
- Withheld intervals: inspect reasons; a short chromosome or sparse bin cannot
  supply enough spatial blocks. Reducing block size can violate the distance
  safeguard, and increasing it reduces replication. Do not force an interval.
- Null sensitivity: fewer than two sizes yielded admissible intervals.
- Empty bins: their estimates are undefined, not zero excess.
- Support options with classified input: these cannot re-filter a TSV; rerun
  using the original trees and the desired support policy.
- Apparent clustering: spatially changing label frequencies or biased missingness
  can also produce excess agreement. Clustering does not establish introgression,
  recombination breakpoints, or a number of independent loci.
