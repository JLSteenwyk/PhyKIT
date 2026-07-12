.. _cmd-trait_rate_map:
.. _command-trait_rate_map:

Trait rate map
==============

Per-branch evolutionary rate map for a continuous trait

Command identity
----------------

:Canonical command: ``trait_rate_map``
:Handler: ``trait_rate_map``
:Aliases: branch_rates, rate_map
:Standalone executables: pk_trait_rate_map, pk_branch_rates, pk_rate_map
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/trait_rate_map.inc

Guidance, interpretation, and examples
--------------------------------------

Estimate per-branch evolutionary rates for a continuous trait and
display them as a branch-colored phylogram. Ancestral states are
reconstructed via Felsenstein's weighted-average method. Per-branch
rate is the squared standardized contrast (proportional to local
Brownian motion rate): ``rate = (child_val - parent_val)^2 / branch_length``.

.. code-block:: shell

   phykit trait_rate_map -t <tree> -d <trait_data> -o <output>
       [--trait <column>] [--json]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]

Options: |br|
*-t/--tree*: a tree file in Newick format |br|
*-d/--trait_data*: tab-delimited trait file (two-column: taxon<tab>value with no header; or multi-column with header when --trait is used) |br|
*-o/--output*: output plot file path (required) |br|
*--trait*: column name to use from a multi-column trait file (if omitted, two-column format is expected) |br|
*--fig-width*: figure width in inches (auto-scaled if omitted) |br|
*--fig-height*: figure height in inches (auto-scaled if omitted) |br|
*--dpi*: resolution in DPI (default: 300) |br|
*--no-title*: hide the plot title |br|
*--title*: custom title text |br|
*--legend-position*: legend location (e.g., "upper right", "none" to hide) |br|
*--ylabel-fontsize*: font size for y-axis labels; 0 to hide |br|
*--xlabel-fontsize*: font size for x-axis labels; 0 to hide |br|
*--title-fontsize*: font size for the title |br|
*--axis-fontsize*: font size for axis labels |br|
*--colors*: comma-separated colors (hex or named) |br|
*--ladderize*: ladderize (sort) the tree before plotting |br|
*--cladogram*: draw cladogram (equal branch lengths, tips aligned) instead of phylogram |br|
*--circular*: draw circular (radial/fan) phylogram instead of rectangular |br|
*--color-file*: color annotation file for tip labels, clade ranges, and branch colors (iTOL-inspired TSV format) |br|
*--json*: optional argument to print results as JSON
