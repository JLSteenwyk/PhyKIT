.. _cmd-rate_heterogeneity:
.. _command-rate_heterogeneity:

Rate heterogeneity test (multi-rate Brownian motion)
====================================================

Test for rate heterogeneity in trait evolution

Command identity
----------------

:Canonical command: ``rate_heterogeneity``
:Handler: ``rate_heterogeneity``
:Aliases: brownie, rh
:Standalone executables: pk_rate_heterogeneity, pk_brownie, pk_rh
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/rate_heterogeneity.inc

Guidance, interpretation, and examples
--------------------------------------

Test for rate heterogeneity across phylogenetic regimes using multi-rate
Brownian motion (O'Meara et al. 2006), analogous to R's
``phytools::brownie.lite()``.

Fits single-rate vs. multi-rate BM models and performs a likelihood ratio
test. Users specify a tree, continuous trait data, and a regime file
mapping tips to regimes (tab-delimited, ``taxon<tab>regime_label``).

Regime assignments to internal branches are inferred via Fitch parsimony.
Per-regime VCV matrices are decomposed and per-regime sigma-squared values
are estimated via maximum likelihood.

Optionally, a parametric bootstrap can be run to compute a simulated
p-value (``-n/--nsim``). A horizontal phylogram with branches colored by
regime can be generated using ``--plot``.

An effect size metric R²_regime is also reported, measuring the variance
reduction from regime-specific rates vs. a single rate, weighted by the
number of tips per regime.

.. code-block:: shell

   phykit rate_heterogeneity -t <tree> -d <trait_data> -r <regime_data> [-n <nsim>] [--seed <seed>] [--plot <output.png>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-t/--tree*: a tree file in Newick format |br|
*-d/--trait_data*: tab-delimited trait file (taxon<tab>value) |br|
*-r/--regime_data*: tab-delimited regime file (taxon<tab>regime_label) |br|
*-n/--nsim*: number of parametric bootstrap simulations (default: 0) |br|
*--seed*: random seed for reproducibility |br|
*--plot*: output plot file path for phylogram with colored branches |br|
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

**R validation:** Validated against ``phytools`` in R
(see ``tests/r_validation/validate_brownie_r2.R``).
