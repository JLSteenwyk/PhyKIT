.. _cmd-projected_covarying_rates:
.. _command-projected_covarying_rates:

Projected covarying rates
=========================

Topology-robust evolutionary-rate covariation by reference-edge projection

Command identity
----------------

:Canonical command: ``projected_covarying_rates``
:Handler: ``projected_covarying_rates``
:Aliases: pcover
:Standalone executables: pk_projected_covarying_rates, pk_pcover
:Categories: Evolutionary rate analysis
:Status: Experimental

Runtime interface
-----------------

.. include:: /_generated/commands/projected_covarying_rates.inc

Guidance, interpretation, and examples
--------------------------------------

Use ``projected_covarying_rates`` to estimate evolutionary-rate covariation
between two genes when their inferred tree topologies differ. Unlike
:doc:`covarying_evolutionary_rates`, this command does not require either gene
tree to match the reference topology after taxon pruning.

This is an experimental PhyKIT estimator. It is not the published CovER
branch-ratio statistic, and coefficients from the two methods are not
interchangeable.

Inputs and preprocessing
^^^^^^^^^^^^^^^^^^^^^^^^

Provide two gene trees and one reference tree in Newick format. Every tree
must have named tips and finite, nonnegative branch lengths. The reference is
typically a species tree, but it may be another biologically justified common
topology.

PhyKIT takes the intersection of the taxa in all three trees and prunes each
tree to that shared set. At least three shared taxa are required. Taxa present
in only one or two inputs do not contribute to the result.

Method
^^^^^^

For the shared taxa, PhyKIT:

1. Calculates all selected pairwise patristic distances in each gene tree.
2. Represents each identifiable unrooted reference-tree edge as a column in a
   path design matrix. A row contains one for every reference edge separating
   that taxon pair and zero otherwise.
3. Uses nonnegative least squares (NNLS) to estimate reference-edge lengths
   whose induced pairwise distances best approximate each gene's observed
   distances.
4. Divides each projected gene-edge length by the corresponding reference-edge
   length, removes undefined zero-reference rates and rates above
   ``--max-rate``, and Z-transforms each retained gene-rate vector.
5. Reports Pearson's correlation between the two standardized projected-rate
   vectors.

The default ``uniform`` weighting gives every taxon pair equal influence.
``--weighting inverse_reference`` weights a pair inversely by its reference
distance, increasing the relative influence of shorter paths.

Pairwise distances cannot separately identify the two complementary branches
below a degree-two root. PhyKIT therefore combines those branch lengths into
one unrooted reference edge. The root stem is excluded. Edge labels contain
the lexicographically canonical, smaller side of each split, with taxon names
separated by semicolons.

Projection diagnostics
^^^^^^^^^^^^^^^^^^^^^^

Each gene receives a normalized root mean squared projection error (NRMSE):

.. math::

   \mathrm{NRMSE} = \sqrt{\frac{\sum_i w_i(\hat{d}_i-d_i)^2}
   {\sum_i w_i d_i^2}}

An NRMSE of zero means that the gene distances are exactly additive on the
reference topology. Larger values mean that more of the gene tree's distance
signal cannot be represented by that topology. There is no universal NRMSE
cutoff; evaluate its empirical distribution for the trees and taxa in the
study.

Interpret a correlation together with both NRMSE values and the number of
retained edges. A high correlation accompanied by poor projection fits is not
strong evidence of shared branchwise rate history. Projection preserves a
common comparison basis, but it does not identify one-to-one homologous
branches between discordant topologies.

The reported p-value is the ordinary two-sided Pearson p-value with retained
reference edges treated as observations. Phylogenetic edges and projected
estimates are not statistically independent, so this p-value is descriptive.
For gene-network analyses, use a study-specific empirical null or threshold.

Scaling and reproducibility
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The full analysis uses all shared-taxon pairs until it reaches ``--max-pairs``
or the internal design-memory limit. If subsampling is needed, pairs are drawn
without replacement using ``--seed`` and then sorted into a deterministic
order. PhyKIT rejects a subsample that cannot uniquely identify every
reference edge; increase ``--max-pairs`` or change ``--seed`` in that case.

Record ``distance_pair_count``, ``distance_pairs_used``, ``--max-pairs``, and
``--seed`` when reporting a subsampled analysis. Check that conclusions are
stable across several seeds when pair subsampling is substantial.

Output
^^^^^^

Default text output reports the correlation, descriptive p-value, shared-taxon
count, used and available distance-pair counts, reference and retained edge
counts, and both projection NRMSE values. ``--verbose`` adds one tab-delimited
row per reference edge with projected lengths, relative rates, Z-scores, and a
filter status.

``--json`` returns those summary values plus:

- ``branches``: all reference-edge records, including excluded edges;
- ``tree_zero_projection`` and ``tree_one_projection``: NRMSE, weighted sum of
  squared errors, and maximum absolute residual;
- ``shared_taxa``: the exact sorted taxon intersection;
- ``weighting`` and ``max_rate``: settings needed to interpret the result.

Examples
^^^^^^^^

Run the uniform-weighted analysis and retain complete diagnostics:

.. code-block:: shell

   phykit projected_covarying_rates gene_a.tre gene_b.tre \
       --reference species.tre --json > projected_rates.json

Emphasize shorter reference paths and save a scatter plot:

.. code-block:: shell

   phykit pcover gene_a.tre gene_b.tre -r species.tre \
       --weighting inverse_reference --plot \
       --plot-output gene_a_gene_b_projected_rates.png

Use the published CovER calculation instead when constrained gene trees have
the same rooted topology as the reference:

.. code-block:: shell

   phykit cover constrained_gene_a.tre constrained_gene_b.tre \
       -r species.tre

Citations
^^^^^^^^^

The evolutionary-rate covariation context follows `Steenwyk et al. (2022)
<https://doi.org/10.1126/sciadv.abn0105>`_. The numerical projection uses the
nonnegative least-squares problem described by `Lawson and Hanson (1995)
<https://doi.org/10.1137/1.9781611971217>`_. The reference-edge projection and
its NRMSE diagnostic are experimental PhyKIT methodology and should be cited
as such rather than attributed to those publications.
