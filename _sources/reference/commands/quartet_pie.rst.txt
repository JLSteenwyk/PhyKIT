.. _cmd-quartet_pie:
.. _command-quartet_pie:

Quartet pie chart
=================

Phylogram with quartet concordance pie charts at internal nodes

Command identity
----------------

:Canonical command: ``quartet_pie``
:Handler: ``quartet_pie``
:Aliases: qpie, quartet_pie_chart
:Standalone executables: pk_quartet_pie, pk_qpie, pk_quartet_pie_chart
:Categories: Tree comparison & consensus

Runtime interface
-----------------

.. include:: /_generated/commands/quartet_pie.inc

Guidance, interpretation, and examples
--------------------------------------

Draw a phylogram with pie charts at internal nodes showing quartet
concordance proportions. In native mode (-g provided), computes gene
concordance factors (gCF, gDF1, gDF2) from a species tree and gene trees
via bipartition matching. In ASTRAL mode (no -g), parses q1/q2/q3
annotations from ASTRAL ``-t 2`` or wASTRAL ``--support 3`` output.

Pie slices show concordant (blue), discordant alt 1 (red), and
discordant alt 2 (gray) proportions. Use ``--annotate`` to add numeric
values near each pie.

.. image:: /_static/quartet_pie_example.png
   :alt: PhyKIT quartet pie example figure
   :align: center
   :width: 90%

.. code-block:: shell

	phykit quartet_pie -t <species_tree> [-g <gene_trees>] -o <output>
		[--annotate] [--branch-labels] [--csv <file>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-t/--tree*: species tree file (required) |br|
*-g/--gene-trees*: gene trees file, one Newick tree per line (optional; if omitted, ASTRAL ``-t 2`` or wASTRAL ``--support 3`` annotations are parsed) |br|
*-o/--output*: output figure path (required; supports .png, .pdf, .svg) |br|
*--annotate*: show gCF/gDF values as text near each pie chart |br|
*--branch-labels*: show the number of concordant gene trees (blue, above branch) and LPP support (red, below branch) on each internal branch. In ASTRAL mode, values come from ``f1`` and ``pp1`` annotations; in native mode, the concordant count is computed from gene trees |br|
*--csv*: output per-branch concordance values (gCF, gDF1, gDF2, counts) as a CSV file |br|
*--pie-size*: scale factor for pie chart size (default: 1.0; use 2.0 for double, 0.5 for half) |br|
*--colors*: comma-separated colors for concordant, disc1, disc2 (default: blue, red, gray) |br|
*--ladderize*: ladderize (sort) the tree before plotting |br|
*--cladogram*: draw cladogram (equal branch lengths, tips aligned) instead of phylogram |br|
*--circular*: draw circular (radial/fan) phylogram instead of rectangular |br|
*--color-file*: color annotation file for tip labels, clade ranges, and branch colors (iTOL-inspired TSV format) |br|
*--json*: optional argument to output per-node concordance values as JSON
