.. _cmd-threshold_model:
.. _command-threshold_model:

Threshold model
===============

Felsenstein threshold model for trait correlation

Command identity
----------------

:Canonical command: ``threshold_model``
:Handler: ``threshold_model``
:Aliases: thresh, thresh_bayes, threshbayes, threshold
:Standalone executables: pk_threshold_model, pk_thresh, pk_thresh_bayes, pk_threshbayes, pk_threshold
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/threshold_model.inc

Guidance, interpretation, and examples
--------------------------------------

Estimate the evolutionary correlation between two traits using the
Felsenstein (2012) threshold model via MCMC. Binary discrete characters
are modelled as arising from continuous latent "liabilities" that evolve
under Brownian motion and cross a threshold at 0. This lets you estimate
correlations between binary traits (or between a binary and a continuous
trait) using BM rather than Mk transition rates.

This is the Python equivalent of ``phytools::threshBayes`` in R.

The sampler uses a Gibbs / Metropolis-Hastings hybrid:

- **Gibbs step**: sample each discrete tip's liability from a truncated
  normal conditioned on all other values
- **Metropolis step**: update sigma2_1, sigma2_2 (log-normal proposal),
  ancestral values a1, a2 (normal proposal), and the correlation r
  (normal proposal with reflection on [-1, 1])
- **Adaptive tuning**: during burn-in, proposal variances are adjusted
  to target ~23% acceptance

Trait types:

- ``discrete``: binary (0/1). Liabilities < 0 map to state 0,
  liabilities > 0 map to state 1.
- ``continuous``: observed values used directly (no liability needed).

Any combination of two traits is supported: discrete+continuous,
discrete+discrete, or continuous+continuous.

.. code-block:: shell

   phykit threshold_model -t <tree> -d <trait_data> --traits <t1,t2> --types <type1,type2> [--ngen 100000] [--sample 100] [--burnin 0.2] [--seed <int>] [--plot <file>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-t/--tree*: a rooted phylogeny file with branch lengths (required) |br|
*-d/--trait-data*: tab-delimited trait file with header row (required) |br|
*--traits*: comma-separated pair of trait column names (required) |br|
*--types*: comma-separated pair of trait types, each ``discrete`` or ``continuous`` (required) |br|
*--ngen*: number of MCMC generations (default: 100000) |br|
*--sample*: sample frequency (default: 100) |br|
*--burnin*: burn-in fraction (default: 0.2) |br|
*--seed*: random seed for reproducibility |br|
*--plot*: output filename for trace and posterior density plot (3 rows x 2 columns: |br|
left = MCMC trace, right = posterior histogram with 95% HPD shading) |br|
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

Output (text mode): |br|
``Trait 1: habitat (discrete, 2 states: 0, 1)`` |br|
``Trait 2: body_mass (continuous)`` |br|
``MCMC: 100000 generations, sampled every 100, burn-in 20%`` |br|
``---`` |br|
``Posterior correlation (r): 0.6234 (95% HPD: 0.312, 0.891)`` |br|
``Posterior sigma2_1: 1.234 (95% HPD: 0.456, 2.345)`` |br|
``Posterior sigma2_2: 0.567 (95% HPD: 0.234, 1.012)`` |br|
``Acceptance rates: r=0.234, sigma2_1=0.312, sigma2_2=0.287, a1=0.241, a2=0.228``

|

**Tutorial: habitat type and body mass in carnivores**

This example uses the classic 8-taxon carnivore tree to test whether
habitat type (a binary trait: 0 = non-arboreal, 1 = arboreal) is
correlated with body mass on the latent liability scale.

*Step 1: Prepare the trait file.*

Create a tab-delimited file with a header row. The first column is
the taxon name, followed by columns for each trait:

.. code-block:: none

   taxon	habitat	body_mass
   raccoon	0	1.04
   bear	0	2.39
   sea_lion	0	2.30
   seal	0	1.88
   monkey	1	0.60
   cat	1	0.56
   weasel	1	-0.30
   dog	0	1.18

*Step 2: Run the threshold model.*

.. code-block:: shell

   phykit threshold_model \
     -t carnivore.tre \
     -d traits.tsv \
     --traits habitat,body_mass \
     --types discrete,continuous \
     --ngen 100000 \
     --seed 42

*Step 3: Examine the trace plots for convergence.*

.. code-block:: shell

   phykit threshold_model \
     -t carnivore.tre \
     -d traits.tsv \
     --traits habitat,body_mass \
     --types discrete,continuous \
     --ngen 100000 \
     --seed 42 \
     --plot trace.png

This generates a 3-row x 2-column figure. The left column shows
MCMC trace plots for r, sigma2_1, and sigma2_2 (check for
stationarity and good mixing). The right column shows posterior
density histograms with red shading for the 95% HPD interval
and a dashed line at the posterior mean.

*Step 4: Get full posterior samples as JSON for custom analysis.*

.. code-block:: shell

   phykit threshold_model \
     -t carnivore.tre \
     -d traits.tsv \
     --traits habitat,body_mass \
     --types discrete,continuous \
     --ngen 100000 \
     --seed 42 \
     --json > posterior.json

The JSON output includes full posterior sample arrays, summary
statistics (mean, median, 95% HPD), and MCMC metadata.
