.. _cmd-identity_matrix:
.. _command-identity_matrix:

Identity matrix
===============

Pairwise sequence identity heatmap

Command identity
----------------

:Canonical command: ``identity_matrix``
:Handler: ``identity_matrix``
:Aliases: id_matrix, seqid
:Standalone executables: pk_identity_matrix, pk_id_matrix, pk_seqid
:Categories: Alignment quality & statistics

Runtime interface
-----------------

.. include:: /_generated/commands/identity_matrix.inc

Guidance, interpretation, and examples
--------------------------------------

Compute a pairwise sequence identity matrix from an alignment and plot it as a
clustered heatmap.

For each pair of taxa, identity is defined as the fraction of non-gap,
non-ambiguous columns that are identical. Gaps, '?', 'N', 'X', and '*' in
either sequence cause a column to be skipped.

The matrix can be displayed as identity (default) or p-distance (1 - identity)
using --metric. Ordering can be by hierarchical clustering (default), tree tip
order (--sort tree --tree <file>), or alphabetical (--sort alpha).

.. code-block:: shell

	phykit identity_matrix -a <alignment> -o <output>
		[--metric identity|p-distance] [--tree <file>]
		[--sort alpha|cluster|tree] [--partition <file>] [--json]
		[--fig-width <float>] [--fig-height <float>] [--dpi <int>]
		[--no-title] [--title <str>]
		[--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
		[--title-fontsize <float>] [--axis-fontsize <float>]

Options: |br|
*-a/--alignment*: alignment file (FASTA or other supported format) |br|
*-o/--output*: output figure path (.png, .pdf, .svg) |br|
*--metric*: 'identity' (fraction matching) or 'p-distance' (1 - identity); default: identity |br|
*--tree*: tree file for tree-guided ordering (Newick format) |br|
*--sort*: ordering method: 'cluster' (hierarchical), 'tree' (requires --tree), or 'alpha' (alphabetical); default: cluster |br|
*--partition*: RAxML-style partition file (reserved for future use) |br|
*--json*: output structured JSON instead of plain text |br|
*--fig-width*: figure width in inches (auto-scaled if omitted) |br|
*--fig-height*: figure height in inches (auto-scaled if omitted) |br|
*--dpi*: resolution in DPI (default: 300) |br|
*--no-title*: hide the plot title |br|
*--title*: custom title text |br|
*--ylabel-fontsize*: font size for y-axis labels; 0 to hide |br|
*--xlabel-fontsize*: font size for x-axis labels; 0 to hide |br|
*--title-fontsize*: font size for the title |br|
*--axis-fontsize*: font size for axis labels
