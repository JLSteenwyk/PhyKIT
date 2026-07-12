.. _cmd-evo_tempo_map:
.. _command-evo_tempo_map:

Evolutionary tempo mapping
==========================

Detect rate-topology associations in gene trees

Command identity
----------------

:Canonical command: ``evo_tempo_map``
:Handler: ``evo_tempo_map``
:Aliases: etm
:Standalone executables: pk_evo_tempo_map, pk_etm
:Categories: Tree comparison & consensus

Runtime interface
-----------------

.. include:: /_generated/commands/evo_tempo_map.inc

Guidance, interpretation, and examples
--------------------------------------

Detect rate-topology associations by comparing branch length distributions
between concordant and discordant gene trees at each species tree branch.

Under the multispecies coalescent, discordant gene trees should have shorter
internal branches near the discordant node (because the coalescence happened
deeper, in the ancestral population). Deviations from this expectation suggest
substitution rate heterogeneity correlated with topology, which could indicate
adaptive evolution, different selective pressures in hybridizing lineages, or
systematic error from model misspecification.

For each internal branch of the species tree, gene trees are classified as
concordant or discordant via bipartition matching (same as gCF). The homologous
branch length is extracted from each gene tree and the two groups are compared
using a Mann-Whitney U test and a permutation test (1000 permutations). P-values
are corrected for multiple testing using Benjamini-Hochberg FDR.

A global treeness (internal/total branch length ratio) comparison between
concordant and discordant gene trees is also reported.

.. code-block:: shell

   phykit evo_tempo_map -t <species_tree> -g <gene_trees> [--plot <output>] [-v]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-t/--tree*: a species tree file |br|
*-g/--gene-trees*: multi-Newick file of gene trees with branch lengths |br|
*--plot*: optional output path for box/strip plot (PNG) |br|
*-v/--verbose*: print per-gene-tree classification details |br|
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

   branch                          n_conc  n_disc    med_conc    med_disc      U_pval   perm_pval       fdr_p
   ----------------------------------------------------------------------------------------------------------
   bear,dog,raccoon                     6       1    3.875000    3.600000          NA          NA          NA
   bear,raccoon                         7       3    0.880000    0.700000    0.516667    0.077000    0.516667
   cat,monkey                          10       0   20.450000          NA          NA          NA          NA
   cat,monkey,weasel                    9       1    2.120000    2.800000          NA          NA          NA
   sea_lion,seal                        9       1    7.500000    7.200000          NA          NA          NA
   ---
   Global treeness: concordant=0.126489 (n=6), discordant=0.119014 (n=4)
   Branches tested: 1, significant (FDR<0.05): 0

Each row corresponds to an internal branch of the species tree identified by the
smaller partition of taxa. The ``n_conc`` and ``n_disc`` columns show how many gene
trees are concordant or discordant at that branch. The ``med_conc`` and ``med_disc``
columns show the median branch length (in substitutions/site) for each group. The
``U_pval`` is the two-sided Mann-Whitney U test p-value, ``perm_pval`` is the
permutation test p-value (1000 permutations), and ``fdr_p`` is the
Benjamini-Hochberg corrected p-value. Branches with fewer than 2 gene trees
in either group show ``NA`` for p-values.

The global treeness comparison tests whether concordant gene trees have
systematically different ratios of internal to total branch lengths.

To generate a visualization:

.. code-block:: shell

   phykit evo_tempo_map -t <species_tree> -g <gene_trees> --plot tempo_map.png

.. image:: /_static/img/tutorial_etm_plot.png
   :align: center
   :width: 80%

The plot shows grouped box plots with jittered data points for each species tree
branch, comparing branch lengths between concordant (blue) and discordant (orange)
gene trees. Branches where the FDR-corrected p-value is below 0.05 are marked
with an asterisk.
