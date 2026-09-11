.. _cmd-topology_autocorrelation:
.. _command-topology_autocorrelation:

Topology autocorrelation
========================

Measure how excess agreement among focal gene-tree topologies varies with
physical genomic distance. Unlike :doc:`topology_landscape`, this command
calculates distance-dependent statistics rather than genomic window counts.

Command identity
----------------

:Canonical command: ``topology_autocorrelation``
:Handler: ``topology_autocorrelation``
:Aliases: topo_ac
:Standalone executables: pk_topology_autocorrelation, pk_topo_ac
:Categories: Tree comparison & consensus

Runtime interface
-----------------

.. include:: /_generated/commands/topology_autocorrelation.inc

Guidance, interpretation, and examples
--------------------------------------

Inputs and selection
^^^^^^^^^^^^^^^^^^^^

Supply either ``--classified-genes landscape.genes.tsv`` or ``--manifest``
with raw trees, ``--coordinates``, and ``--groups`` or ``--reference-tree``.
Raw mode uses the same strict classification, support filtering, reference-edge
selection, mapping diagnostics, and midpoint selection as :doc:`topology_landscape`.
The command does not first write an intermediate landscape or compute windows.

Classified TSV requires ``reference``, ``chromosome``, ``gene_id``, ``start``,
``end``, ``anchor``, and ``classification``. Coordinates are zero-based half-open;
the anchor must equal ``(start + end - 1) // 2``. Accepted classes are
``topology_1``, ``topology_2``, ``topology_3``, ``unresolved``,
``insufficient_sampling``, and ``incompatible_groups``. Other landscape columns
are retained as metadata. Duplicate reference/gene pairs, invalid anchors,
and conflicting classifications for the same gene across references are errors.
In raw mode, coordinate record merging follows the landscape rules instead.

Classification labels and support policies cannot be recovered from a standalone
TSV. Existing classifications are retained without re-filtering. Supply
``--labels`` to name the three alternatives; raw-tree and support options are
rejected in classified mode rather than silently ignored.

``--chromosome`` is repeatable. ``--interval START END`` retains anchors in
that half-open interval on each selected chromosome and reference. Selection
happens before estimating chromosome frequencies. Output ordering is lexical
by reference/chromosome, then increasing distance and metric order.

Distances and denominators
^^^^^^^^^^^^^^^^^^^^^^^^^^

Distance is the absolute difference between two gene anchors. Each unordered
pair of resolved genes contributes once; distinct genes at the same anchor
have distance zero. No pairs cross chromosomes, reference genomes, or selection
boundaries. Unresolved/unclassifiable genes do not contribute to pair
denominators, but remain in per-gene and diagnostic output.

Set ``--distance-edges 0 1000 5000 10000`` for custom bins, or combine
``--bin-width`` (default 10,000 bp) with ``--max-distance`` (default 1,000,000 bp).
Every bin is ``[lower, upper)``, including the last: a pair at exactly the
maximum distance is excluded. A final equal-width bin may be shorter.
At most 1,000 bins are supported. Empty bins remain visible with zero pairs
and null estimates, not apparent evidence for zero excess.

Let chromosome c have n resolved genes, including n_k in topology k.
Its finite-population baseline joint support is
``q_ck = n_k (n_k - 1) / (n (n - 1))``. This is the expected probability that
both members of a pair support k under random reassignment of observed labels
without replacement. Aggregate baseline agreement is the sum over k.

In each distance bin, ``observed`` is the fraction of resolved pairs supporting
the same topology (metric ``agreement``), or both supporting the specified
topology (metrics ``topology_1`` through ``topology_3``).
``excess = observed - baseline``. Topology-specific excesses sum to aggregate
excess. ``enrichment_ratio = observed / baseline`` is secondary and null when
the baseline is zero. These quantities are not Pearson correlations or
conditional probabilities given one endpoint's topology.

Reference-level estimates weight chromosome estimates and baselines by the
number of pairs in that distance bin. They never estimate a baseline from
pooled chromosome frequencies. The baseline can therefore change across bins
as contributing chromosomes change. ``contributing_genes`` counts distinct
pair endpoints, while ``contributing_chromosomes`` counts chromosomes with pairs.
Singleton chromosomes have undefined baselines and contribute no pairs.

Block resampling
^^^^^^^^^^^^^^^^

Without ``--block-sizes``, output is descriptive and interval/block TSV files
contain headers only. To request uncertainty, supply at least two distinct
positive sizes, for example ``--block-sizes 50000 100000 --replicates 999 --seed 7``.
At least 499 replicates are required; the default is 999. Seeds reproduce runs
with the same implementation, inputs, and numerical-library versions.

Each chromosome's observed span from its first to last selected anchor is
divided into ``max(1, floor(span / requested_size))`` near-equal physical blocks.
Widths differ by at most one base; there is no short terminal block. Empty
blocks remain eligible. This is an observed-gene span, not a known chromosome
length. Block diagnostics include actual bounds and all six class counts.

The bootstrap draws the original number of blocks with replacement independently
within each chromosome, using the same draws across all bins and metrics.
It resamples precomputed contributions, not genes joined into a synthetic
genome. Original pairs spanning a block boundary are retained; no new pairs,
wrapping, or adjacency are introduced.

Each pair denominator is split equally between its endpoints. For topology
indicator I and observed chromosome frequency p, its endpoint-i contribution is
``0.5 (I_i-p)(I_j-p) + p I_i - 0.5 p^2``. The two endpoint contributions sum
exactly to ``I_i I_j``. This algebraic centering assigns the linear frequency
term to its owning gene and reduces block-boundary artifacts. These signed
resampling contributions are separate from the actual reported pair counts.
Each replicate recomputes chromosome frequencies, finite baselines, and
pair-weighted excesses. Percentiles yield 95% pointwise intervals for excess.

Intervals are withheld for a bin if any contributing chromosome has fewer than
20 pair-occupied blocks or its actual block widths are less than five times
the bin's upper distance. They are also withheld if fewer than 95% of
replicates are valid or their distribution is degenerate. Reasons are explicit.
If admissible interval widths differ by more than a factor of two across
sizes, ``block_size_sensitive`` is true. With fewer than two admissible sizes,
sensitivity is null with a diagnostic, not a declaration of stability.

These safeguards do not verify stationarity or make neighboring blocks
independent. Intervals assume approximately stationary, short-range-dependent
marked gene processes within chromosomes and sufficient spatial replication.
They describe repeated spatial sampling, not a random-label null distribution.
They are not simultaneous confidence bands: scanning all bins for a pointwise
interval that excludes zero is not a multiple-testing-corrected analysis.
No p-values, significance calls, or clustering-range estimates are provided.

Outputs and limitations
^^^^^^^^^^^^^^^^^^^^^^^

Required ``--output-prefix PREFIX`` writes:

- ``.autocorrelation.tsv``: per-reference, per-bin, per-metric estimates and counts.
- ``.chromosomes.tsv``: corresponding chromosome-level descriptive estimates.
- ``.genes.tsv``: selected mapped genes, including unresolved/unclassifiable genes.
- ``.uncertainty.tsv``: intervals by reference, bin, metric, and requested block size.
- ``.blocks.tsv``: physical block bounds and class-count diagnostics for every size.
- ``.json``: complete results, parameters, input metadata, warnings, and diagnostics.

``--json`` also prints the complete payload. JSON null and empty TSV cells
mean undefined. TSV list/dictionary cells use JSON encoding. Existing parent
directories are required; outputs cannot overwrite input files or each other.

``--plot`` draws aggregate excess, topology-specific excess, and resolved-pair
counts. Dashed zero lines are the baseline for excess. Available 95% pointwise
intervals are shaded for requested sizes; exact size-specific values are in TSV.
Each reference gets a separate figure. PNG, PDF, SVG, JPEG, and TIFF are supported.
Multiple references add numbered filename suffixes. Configure labels, colors,
dimensions, fonts, and legend placement using standard plotting options.

Spatially changing topology frequencies, gene density, or nonrandom unresolved
genes can cause apparent clustering. Inspect class-count block diagnostics and
per-gene maps. Agreement concerns estimated gene trees, not error-free genealogies.
This analysis does not establish introgression, identify recombination breakpoints,
or estimate the number of independent loci. ``clustering_range`` is explicitly
``not_estimable``; no arbitrary first-zero-crossing rule is used.

Examples and related methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See :doc:`/tutorials/pages/23-topology-autocorrelation` for two datasets with
identical topology totals but different spatial structure, and a simulated
example for block-size sensitivity. The method specification and reproducible
calibration scripts are in ``docs/plans`` and ``benchmarks`` in the repository.

Spatial clustering of gene-tree support was studied by
`Pollard et al. (2006), PLoS Genetics, doi:10.1371/journal.pgen.0020173
<https://doi.org/10.1371/journal.pgen.0020173>`_. Resampling precomputed local
contributions is motivated by
`Loh and Stein (2004), Statistica Sinica 14:69-101
<https://web.njit.edu/~loh/Papers/Sinica.LohStein.2004.pdf>`_. Our categorical
estimator and centering are adaptations requiring separate calibration,
not their K-function estimator. This is not a claim of a novel spatial
statistical concept. Raw focal classification uses :doc:`topology_landscape`,
which explains its distinction from Twisst's fractional topology weighting.
