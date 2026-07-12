.. _cmd-discordance_asymmetry:
.. _command-discordance_asymmetry:

Discordance asymmetry
=====================

Test for asymmetric discordance (gene flow detection)

Command identity
----------------

:Canonical command: ``discordance_asymmetry``
:Handler: ``discordance_asymmetry``
:Aliases: da, disc_asym
:Standalone executables: pk_discordance_asymmetry, pk_da, pk_disc_asym
:Categories: Introgression & gene flow

Runtime interface
-----------------

.. include:: /_generated/commands/discordance_asymmetry.inc

Guidance, interpretation, and examples
--------------------------------------

Test whether the two discordant NNI alternative topologies at each species tree
branch are equally frequent. Under incomplete lineage sorting (ILS) alone, the
two minor NNI alternatives (gDF1 and gDF2) should appear at equal frequency.
When they are significantly asymmetric, it suggests introgression or gene flow
between specific lineages.

For each internal branch of the species tree, a two-sided binomial test (H0:
P(alt1) = 0.5) is applied, and p-values are corrected for multiple testing
using Benjamini-Hochberg FDR.

**Interpreting the plot:** When ``--plot`` is used, branches are colored by
the **asymmetry ratio** = max(gDF1, gDF2) / (gDF1 + gDF2):

- **Blue/cool colors (ratio ~ 0.5):** The two discordant topologies are
  roughly equal, consistent with ILS alone — no evidence of gene flow.
- **Red/warm colors (ratio → 1.0):** One discordant topology is much more
  frequent than the other, suggesting introgression or gene flow.
- **Red stars:** Branches where the asymmetry is statistically significant
  (FDR < 0.05) — these are introgression candidates.

Use ``--annotate`` to display gCF values on each branch. Use
``--ylabel-fontsize 0`` to hide tip labels for large trees, and
``--legend-position none`` to hide both the legend and colorbar.

.. code-block:: shell

   phykit discordance_asymmetry -t <species_tree> -g <gene_trees> [--plot <output>] [-v]
       [--annotate]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-t/--tree*: a species tree file |br|
*-g/--gene-trees*: multi-Newick file of gene trees (branch lengths not required) |br|
*--plot*: optional output path for asymmetry phylogram (PNG) |br|
*-v/--verbose*: print per-branch details |br|
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

Example output:

.. code-block:: text

   branch                          n_conc  n_alt1  n_alt2  asym_ratio     binom_p       fdr_p   gene_flow
   ------------------------------------------------------------------------------------------------------
   bear,dog,raccoon                     6       0       1       1.000      1.0000      1.0000           -
   bear,raccoon                         7       1       2       0.667      1.0000      1.0000           -
   cat,monkey                          10       0       0          NA          NA          NA           -
   cat,monkey,weasel                    9       1       0       1.000      1.0000      1.0000           -
   sea_lion,seal                        9       1       0       1.000      1.0000      1.0000           -
   ---
   Summary: 4 branches tested, 0 significant (FDR<0.05)

Each row corresponds to an internal branch of the species tree identified by the
smaller partition of taxa. The ``n_conc`` column shows concordant gene trees, while
``n_alt1`` and ``n_alt2`` show the counts for the two NNI alternative topologies.
The ``asym_ratio`` is max(n_alt1, n_alt2) / (n_alt1 + n_alt2), ranging from 0.5
(perfectly symmetric) to 1.0 (maximally asymmetric). ``NA`` indicates no discordant
gene trees were observed. The ``binom_p`` is the two-sided binomial test p-value,
``fdr_p`` is the Benjamini-Hochberg corrected p-value, and ``gene_flow`` shows
which NNI alternative is favored when the result is significant (FDR < 0.05).

**Interpretation:**

- **Symmetric discordance** (asym_ratio near 0.5, not significant): Consistent with
  ILS alone — both NNI alternatives arise with equal frequency from random coalescent
  sorting.
- **Asymmetric discordance** (asym_ratio near 1.0, significant): Suggests gene flow
  or introgression. The favored NNI alternative indicates which lineages are exchanging
  genetic material. For alt1, gene flow is between the C1 lineage and the S (sibling)
  lineage; for alt2, gene flow is between C2 and S.

To generate a visualization:

.. code-block:: shell

   phykit discordance_asymmetry -t <species_tree> -g <gene_trees> --plot asymmetry.png

.. image:: /_static/img/discordance_asymmetry_plot.png
   :alt: PhyKIT discordance asymmetry plot figure
   :align: center
   :width: 80%

The plot shows a phylogram with internal branches colored by asymmetry ratio using
a diverging colormap (blue = symmetric/0.5, red = highly asymmetric/1.0). Significant
branches (FDR < 0.05) are marked with a red star. Internal nodes are annotated with
gCF values.
