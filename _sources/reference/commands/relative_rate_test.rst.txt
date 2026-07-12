.. _cmd-relative_rate_test:
.. _command-relative_rate_test:

Relative rate test
==================

Relative rate test between lineages

Command identity
----------------

:Canonical command: ``relative_rate_test``
:Handler: ``relative_rate_test``
:Aliases: rrt, tajima_rrt
:Standalone executables: pk_relative_rate_test, pk_rrt, pk_tajima_rrt
:Categories: Evolutionary rate analysis

Runtime interface
-----------------

.. include:: /_generated/commands/relative_rate_test.inc

Guidance, interpretation, and examples
--------------------------------------

Tajima's relative rate test (Tajima, *Genetics*, 1993).

Tests whether two ingroup lineages have evolved at equal rates
since diverging from their common ancestor. The outgroup is
automatically inferred from the rooted tree as the earliest-diverging
taxon (the single taxon on the smaller side of the root split).
All pairwise ingroup combinations are tested with Bonferroni
and Benjamini-Hochberg FDR multiple testing correction.

At each alignment column, the test classifies informative sites:

- **m1**: the first ingroup taxon differs from the outgroup, but the second matches
- **m2**: the second ingroup taxon differs from the outgroup, but the first matches
- Sites where both differ or both match the outgroup are uninformative and skipped
- Sites with gaps or ambiguous characters are skipped

Test statistic: ``chi2 = (m1 - m2)^2 / (m1 + m2)``, with 1 degree of freedom.

**Single alignment mode:**

.. code-block:: shell

   phykit relative_rate_test -a <alignment> -t <rooted_tree> [-v/--verbose]
       [--plot-output <path>] [--fig-width <float>] [--fig-height <float>] [--dpi <int>]
       [--no-title] [--title <str>] [--legend-position <str>] [--ylabel-fontsize <float>]
       [--xlabel-fontsize <float>] [--title-fontsize <float>] [--axis-fontsize <float>]
       [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

**Batch mode (multiple alignments, one shared tree):**

.. code-block:: shell

   phykit relative_rate_test -l <alignment_list> -t <rooted_tree> [-v/--verbose]
       [--plot-output <path>] [--fig-width <float>] [--fig-height <float>] [--dpi <int>]
       [--no-title] [--title <str>] [--legend-position <str>] [--ylabel-fontsize <float>]
       [--xlabel-fontsize <float>] [--title-fontsize <float>] [--axis-fontsize <float>]
       [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-a/--alignment*: a single alignment file |br|
*-l/--alignment-list*: a file with one alignment path per line (batch mode) |br|
*-t/--tree*: a rooted tree file |br|
*-v/--verbose*: print additional detail |br|
*--plot-output*: save a pairwise p-value heatmap to this path |br|
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

Single mode outputs one row per ingroup pair with m1, m2, chi-squared,
raw p-value, Bonferroni-corrected p-value, FDR-corrected p-value, and a
significance indicator. Batch mode aggregates across genes, reporting the
number and percentage of genes rejecting equal rates for each pair.

When ``--plot-output`` is provided, a symmetric heatmap of
-log10(FDR-corrected p-values) is saved. Cells with significant
pairs (FDR < 0.05) are marked with an asterisk. Darker colors
indicate stronger evidence for unequal evolutionary rates.
The dashed line on the colorbar marks the significance threshold
(-log10(0.05) ≈ 1.3); cells above this line are significant.

**Equal rates** — all taxa evolve at similar rates, so no pairwise
comparison reaches significance. The heatmap is uniformly pale with no
asterisks:

.. image:: /_static/img/rrt_equal_rates.png
   :align: center
   :width: 80%

|

**Unequal rates** — taxa B and C have accumulated many more substitutions
than A and D. Significant pairs (dark red, marked with ``*``) indicate
lineages evolving at detectably different rates:

.. image:: /_static/img/rrt_unequal_rates.png
   :align: center
   :width: 80%

|

Validated against R's ``pegas::rr.test()`` — chi-squared and p-values
match to machine precision.
