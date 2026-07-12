.. _cmd-dtt:
.. _command-dtt:

Disparity through time (DTT)
============================

Disparity-through-time analysis with MDI statistic (Harmon et al. 2003)

Command identity
----------------

:Canonical command: ``dtt``
:Handler: ``dtt``
:Aliases: disparity_through_time
:Standalone executables: pk_dtt, pk_disparity_through_time
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/dtt.inc

Guidance, interpretation, and examples
--------------------------------------

Compute disparity-through-time (DTT) curves following Harmon et al. (2003).
At each branching time the mean relative subclade disparity is calculated.
Under Brownian motion the DTT curve declines linearly from 1 to 0; deviations
from this expectation reveal tempo of morphological diversification.

The Morphological Disparity Index (MDI) measures the area between the
observed DTT and the Brownian motion null median. Positive MDI indicates
late disparity accumulation; negative MDI indicates early disparity
accumulation (consistent with adaptive radiation).

.. code-block:: shell

   phykit dtt -t <tree> --traits <traits_file>
       [--trait <column>] [--index avg_sq|avg_manhattan]
       [--nsim <int>] [--seed <int>]
       [--plot-output <file>] [--json]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>]

Options: |br|
*-t/--tree*: ultrametric tree file (required) |br|
*--traits*: TSV file with header row; first column ``taxon``, remaining columns are trait values (required) |br|
*--trait*: specific trait column name; if omitted all trait columns are used (multivariate DTT) |br|
*--index*: disparity index — ``avg_sq`` (average squared Euclidean distance, default) or ``avg_manhattan`` |br|
*--nsim*: number of Brownian motion simulations for null DTT envelope and MDI p-value (default: 0, no simulations) |br|
*--seed*: random seed for reproducibility |br|
*--plot-output*: output filename for the DTT plot (PNG, PDF, SVG) |br|
*--json*: output results as JSON

**Default output**: text summary with DTT curve data points.

.. code-block:: shell

   # Single trait DTT
   phykit dtt -t species.tre --traits traits.tsv --trait body_mass

   # Multivariate DTT (all traits)
   phykit dtt -t species.tre --traits traits.tsv

   # With BM null envelope and MDI p-value
   phykit dtt -t species.tre --traits traits.tsv --trait body_mass --nsim 1000 --seed 42

   # With plot output
   phykit dtt -t species.tre --traits traits.tsv --trait body_mass --nsim 1000 --plot-output dtt.png

   # JSON output
   phykit dtt -t species.tre --traits traits.tsv --trait body_mass --json

**R validation:** Validated against ``geiger::dtt()`` in R
(see ``tests/r_validation/validate_dtt.R``).
