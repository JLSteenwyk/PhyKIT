.. _cmd-ancestral_state_reconstruction:
.. _command-ancestral_state_reconstruction:

Ancestral state reconstruction
==============================

Reconstruct ancestral character states

Command identity
----------------

:Canonical command: ``ancestral_state_reconstruction``
:Handler: ``ancestral_state_reconstruction``
:Aliases: anc_recon, asr
:Standalone executables: pk_ancestral_state_reconstruction, pk_anc_recon, pk_asr
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/ancestral_state_reconstruction.inc

Guidance, interpretation, and examples
--------------------------------------

Estimate ancestral states using maximum likelihood. Supports both
continuous and discrete traits.

**Continuous traits** (``--type continuous``, default): Brownian Motion
model, analogous to R's ``phytools::fastAnc()`` and
``ape::ace(type="ML")``. Optionally produce a contMap plot.

Two methods are available for continuous traits:

- **fast** (default): Felsenstein's pruning/contrasts shortcut, O(n) time
- **ml**: full VCV-based ML with exact conditional CIs, O(n^3)

Both methods produce identical point estimates; ``ml`` gives exact
conditional confidence intervals.

**Discrete traits** (``--type discrete``): Mk model with marginal
posterior probabilities at each internal node, analogous to
``ape::ace(type="discrete")``. Optionally produce a pie-chart phylogeny
plot.

Three models are available for discrete traits:

- **ER** (default): equal rates
- **SYM**: symmetric rates
- **ARD**: all rates different

Input trait data can be either a two-column file (``taxon<tab>value``)
when ``-c`` is omitted, or a multi-trait file with header row when ``-c``
specifies which column to use.

.. code-block:: shell

   phykit ancestral_state_reconstruction -t <tree> -d <trait_data> [-c <trait>] [--type <type>] [-m <method>] [--model <model>] [--ci] [--plot <output>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*-t/--tree*: a phylogenetic tree file |br|
*-d/--trait_data*: trait data file (two-column or multi-trait with header) |br|
*-c/--trait*: trait column name (required for multi-trait files) |br|
*--type*: trait type: ``continuous`` or ``discrete`` (default: ``continuous``) |br|
*-m/--method*: method to use: ``fast`` or ``ml`` (continuous only; default: ``fast``) |br|
*--model*: Mk model: ``ER``, ``SYM``, or ``ARD`` (discrete only; default: ``ER``) |br|
*--ci*: include 95% confidence intervals (continuous only) |br|
*--plot*: output path for plot (requires matplotlib) |br|
*--plot-ci*: draw confidence interval bars at internal nodes on the contMap plot (requires --ci and --plot) |br|
*--ci-size*: scale factor for CI bar size (default: 1.0; use 2.0 for larger, 0.5 for smaller) |br|
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
*--json*: output results as JSON

Example contMap plot generated with the ``--plot`` option. Branches are colored
by interpolated ancestral trait values:

.. image:: /_static/img/asr_example.png
   :alt: PhyKIT asr example figure
   :align: center
   :width: 80%
