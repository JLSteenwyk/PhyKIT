.. _cmd-alignment_outlier_taxa:
.. _command-alignment_outlier_taxa:

Alignment outlier taxa
======================

Identify outlier taxa in alignments

Command identity
----------------

:Canonical command: ``alignment_outlier_taxa``
:Handler: ``alignment_outlier_taxa``
:Aliases: aot, outlier_taxa
:Standalone executables: pk_alignment_outlier_taxa, pk_aot, pk_outlier_taxa
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/alignment_outlier_taxa.inc

Guidance, interpretation, and examples
--------------------------------------

Identify potential outlier taxa in an alignment and explicitly report why each taxon was flagged.

The following features are evaluated per taxon:
1) ``gap_rate``: fraction of gap/ambiguous symbols |br|
2) ``occupancy``: fraction of valid symbols |br|
3) ``composition_distance``: Euclidean distance from the median composition profile |br|
4) ``long_branch_proxy``: mean pairwise sequence distance to other taxa |br|
5) ``rcvt``: relative composition variability per taxon |br|
6) ``entropy_burden``: average site entropy over this taxon's valid positions

If a taxon exceeds one or more feature-specific thresholds, PhyKIT reports the exact
feature(s), observed value(s), threshold(s), and explanation(s) for the flag.

.. code-block:: shell

	phykit alignment_outlier_taxa <alignment> [--gap-z <float>] [--composition-z <float>] [--distance-z <float>] [--rcvt-z <float>] [--occupancy-z <float>] [--entropy-z <float>] [--json]

Example output:

.. code-block:: shell

	features_evaluated	gap_rate,occupancy,composition_distance,long_branch_proxy,rcvt,entropy_burden
	thresholds	gap_rate>0.0;occupancy<0.4;composition_distance>0.1;long_branch_proxy>0.4181;rcvt>0.0095;entropy_burden>0.35
	taxon_d	composition_distance=1.3454>0.1;long_branch_proxy=1.0>0.4181	Unusual sequence composition profile relative to other taxa. | High mean pairwise sequence distance to other taxa.
	taxon_e	gap_rate=0.6>0.0;occupancy=0.4<0.4	High fraction of gap/ambiguous symbols compared to other taxa. | Low fraction of valid symbols compared to other taxa.

Options: |br|
*<alignment>*: first argument after function name should be an alignment file |br|
*--gap-z*: z-threshold used for high-gap outlier detection (default: 3.0) |br|
*--composition-z*: z-threshold used for composition-distance outlier detection (default: 3.0) |br|
*--distance-z*: z-threshold used for long-branch-proxy outlier detection (default: 3.0) |br|
*--rcvt-z*: z-threshold used for RCVT outlier detection (default: 3.0) |br|
*--occupancy-z*: z-threshold used for low-occupancy outlier detection (default: 3.0) |br|
*--entropy-z*: z-threshold used for entropy-burden outlier detection (default: 3.0) |br|
*--json*: optional argument to print results as JSON
