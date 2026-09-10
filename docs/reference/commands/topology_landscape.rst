.. _cmd-topology_landscape:
.. _command-topology_landscape:

Topology landscape
==================

Map support for a focal quartet relationship across reference genomes.

Command identity
----------------

:Canonical command: ``topology_landscape``
:Handler: ``topology_landscape``
:Aliases: topomap
:Standalone executables: pk_topology_landscape, pk_topomap
:Categories: Tree comparison & consensus

Runtime interface
-----------------

.. include:: /_generated/commands/topology_landscape.inc

Guidance, interpretation, and examples
--------------------------------------

Use this command to locate individual genes supporting competing relationships,
such as ctenophore-sister versus sponge-sister. It counts strict gene-tree
concordance, not confidence in a species tree. No likelihood tests, p-values,
introgression calls, or recombination breakpoints are calculated. Neighboring
genes are not necessarily independent evolutionary observations.

Inputs
^^^^^^

``--manifest`` is a tab-separated table with ``gene_id`` and ``tree`` headers.
Each identifier and tree path must occur once; each file contains exactly one
Newick tree with unique named tips. Relative paths resolve beside the manifest,
not relative to the working directory. Branch lengths are not used.

Repeat ``--coordinates REFERENCE FORMAT PATH`` for coordinate files. FORMAT is
``bed`` (at least four tab-separated columns: chromosome, start, end, gene ID)
or ``tsv`` (headers: ``gene_id``, ``chromosome``, ``start``, ``end``).
Both formats use zero-based, half-open intervals: ``0 <= start < end``.
BED track/browser lines, comments, and blank lines are ignored. Additional BED
columns or TSV columns do not change classification or weighting.

The reference identifier is mandatory. Files sharing an identifier are combined;
different identifiers remain separate even if chromosome names match.
Repeated gene records on the same chromosome merge to the span from their
minimum start to maximum end, with diagnostics. This supports repeated intervals
but does not infer transcript-to-gene identifiers. Map transcripts to gene IDs
before use. A gene on multiple chromosomes of one reference is an error.
GFF3 is not directly supported in this release; export gene-level BED or TSV
using an annotation-aware parser and the intended gene-ID attribute.

Manifest genes absent from all coordinate files cause an error by default.
``--unmapped skip`` omits them from maps and lists them in diagnostics.
Coordinate-only identifiers are reported and excluded. Diagnostics also list
manifest genes missing from each reference, and merged-record counts.

Selecting a relationship
^^^^^^^^^^^^^^^^^^^^^^^^

``--groups`` reads a JSON object containing exactly four named, nonempty,
disjoint taxon lists. Object order defines groups A, B, C, D and the three
resolutions: ``topology_1 = AB|CD``, ``topology_2 = AC|BD``,
``topology_3 = AD|BC``. JSON key order in result output is not significant;
use the reported labels to interpret classification IDs.

``--outgroup NAME`` moves that group first without reordering the other groups.
Its partner in each split is the inferred sister group among the three ingroups.
``--labels LABEL1 LABEL2 LABEL3`` customizes the displayed hypothesis names.
An outgroup grouping is explicit; arbitrary Newick roots are never trusted.

Alternatively, supply ``--reference-tree TREE --branch-taxa TAXON ...``. The
taxa must exactly identify either complete side of an existing internal edge.
After suppressing artificial degree-two vertices, both endpoints must have
degree three. The two components on the selected side become A/B, and the
opposite components become C/D, sorted lexicographically within each side.
The first resolution therefore matches the selected reference edge.
Reference mode reports the resulting groups so the choice can be audited.

Strict classification
^^^^^^^^^^^^^^^^^^^^^

All sampled representatives are retained; taxa outside the groups are ignored.
Each group needs at least one sampled member, not necessarily all listed members.
``sampled_counts`` reports actual representation; ``sampling_status=complete``
means all four groups are represented.

Each multi-tip group must form an unrooted clan: an edge separates its sampled
members from all other sampled groups. If a retained edge contradicts this
separation, classify the gene as ``incompatible_groups``. If no contradictory
edge exists but a clan-defining edge is missing, classify it as ``unresolved``.
Otherwise, classify by the group-pair split, or ``unresolved`` if absent.
Missing any group gives ``insufficient_sampling`` before topology evaluation.
Singleton groups automatically satisfy the clan criterion. This strict rule
does not select representatives, count quartet votes, or assign fractional weights.

Support filtering
^^^^^^^^^^^^^^^^^

Without ``--min-support``, every existing split is retained regardless of its
support. With a threshold, splits below it are collapsed. ``--support-scale``
must explicitly be 100 (default) or 1; numeric labels outside that scale are
rejected rather than silently converted. Nonnumeric internal names are not
support values. No support is inferred from branch length.

If multiple branches represent the same induced split after root suppression
or excluding extra taxa, use the minimum available numeric support. This is a
conservative path rule, not a probability for the induced split. When any label
on that path is absent, ``--missing-support collapse`` removes it by default;
``keep`` retains it subject to known numeric labels, and ``error`` aborts.
The policy applies only when a threshold is supplied. Output records both the
minimum numeric resolution support and whether any contributing label was absent.

Coordinates, windows, and output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A merged gene span is assigned to integer position ``(start + end - 1) // 2``.
This midpoint, not interval overlap, determines chromosome-interval selection
and window membership. Every gene contributes once per reference to an overall
summary and to exactly one nonoverlapping neighborhood.

``--window-bp`` defaults to 1,000,000. Physical windows begin at zero and stop
at the last retained gene end, not an inferred chromosome length. With
``--interval START END``, they begin at START and stop at END. Empty windows
between retained genes are included, with zero counts and undefined proportions.
An entirely empty selection is an error; chromosomes without mapped genes are
not synthesized. The interval applies to each selected chromosome/reference.

``--window-genes`` instead partitions sorted genes into equal-count bins,
including unclassifiable genes. Ties in midpoint are broken by gene ID. The
last bin may be smaller. Its coordinate extent covers the first through last
anchor; different rank bins may have overlapping coordinate extents when tied.
Gene-count summary plots therefore use window index, while the gene track
retains reference coordinates. Chromosomes and references sort lexicographically.

The required output prefix produces:

- ``.genes.tsv``: coordinates, anchor, classification, sampling counts, ignored
  taxa, and resolution support. Dictionary/list cells contain JSON.
- ``.neighborhoods.tsv``: bounds, index, mode, counts, and proportions per window.
- ``.overall.tsv``: per-reference counts and proportions after selection.
- ``.diagnostics.json``: mapping omissions and merged records.

``total_genes`` includes all six classifications. ``informative_genes`` includes
only the three resolved topologies and is the denominator of topology proportions.
When it is zero, proportions are JSON null or empty TSV cells, never zero support.
``--json`` prints the full result including groups, labels, every manifest gene's
classification, mapped genes, summaries, diagnostics, and output paths.

``--plot`` writes three-panel figures: categorical gene marks, resolved topology
fractions, and all-versus-resolved gene counts. Each reference/chromosome gets a
separate figure; multiple figures add numbered suffixes to ``--plot-output``.
PNG, PDF, SVG, JPEG, and TIFF are supported. Parent output directories must exist.

Example and related methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^

See :doc:`/tutorials/pages/22-mapping-topologies-across-genomes` for a runnable
synthetic animal-root example with known counts. For concordance on the species
tree rather than genomic coordinates, use :doc:`quartet_pie`.

Topology mapping is established methodology, not a new statistical test.
Related topology weighting is described by `Martin and Van Belleghem (2017),
Genetics, doi:10.1534/genetics.116.194720
<https://doi.org/10.1534/genetics.116.194720>`_. This command differs from Twisst:
it assigns strictly classifiable genes whole votes rather than weighting sampled
subtree resolutions. Cite PhyKIT for this implementation and describe the strict
classification and coordinate-window rules used in the analysis.
