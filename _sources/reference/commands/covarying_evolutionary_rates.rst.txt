.. _cmd-covarying_evolutionary_rates:
.. _command-covarying_evolutionary_rates:

Covarying evolutionary rates
============================

Detect covariation in evolutionary rates

Command identity
----------------

:Canonical command: ``covarying_evolutionary_rates``
:Handler: ``covarying_evolutionary_rates``
:Aliases: cover
:Standalone executables: pk_covarying_evolutionary_rates, pk_cover
:Categories: Evolutionary rate analysis

Runtime interface
-----------------

.. include:: /_generated/commands/covarying_evolutionary_rates.inc

Guidance, interpretation, and examples
--------------------------------------

Determine if two genes have a signature of covariation with one another.
Genes that have covarying evolutionary histories tend to have 
similar functions and expression levels.

Input two gene phylogenies and calculate the Pearson correlation among
their relative evolutionary rates. The two gene trees and the reference
tree may initially contain different taxa. PhyKIT prunes all three trees
to their shared taxa, after which they must have the same rooted topology.

For each corresponding terminal and internal branch, PhyKIT calculates
``gene branch length / reference branch length``. A root stem is not a
phylogenetic branch and is excluded. Branches with a zero or missing
reference length are also excluded because their ratio is undefined;
zero-length gene branches are retained. A branch pair is removed if either
relative rate is greater than five. The remaining values are Z-transformed,
and their Pearson correlation coefficient is reported.

PhyKIT reports two tab delimited values:
col1: correlation coefficient
col2: p-value

This implementation is the CovER procedure described by `Steenwyk et al.
(2022) <https://doi.org/10.1126/sciadv.abn0105>`_. It is based on the
mirror-tree principle developed by `Sato et al. (2005)
<https://doi.org/10.1093/bioinformatics/bti564>`_, and the biological basis
of evolutionary rate covariation was evaluated by `Clark et al. (2012)
<https://doi.org/10.1101/gr.132647.111>`_.

Comparison with other ERC software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An independent implementation of this branch-ratio procedure using R
``ape`` and ``stats::cor.test`` reproduces PhyKIT's coefficient and p-value.
The validation is maintained in ``tests/validation/validate_covarying_rates.R``.

Do not expect the same numerical coefficient from `ERC2.0
<https://github.com/nclark-lab/erc>`_ or `RERconverge
<https://github.com/nclark-lab/RERconverge>`_. Those programs estimate
relative rates as regression residuals and include additional normalization
such as heteroskedasticity correction and winsorization. These are related
ERC methods, not independent implementations of the PhyKIT/CovER statistic.

The reported p-value is the standard two-sided Pearson correlation p-value,
using the retained branches as observations. Because branches in a phylogeny
are not strictly independent, interpret this parametric p-value cautiously;
for gene-network analyses, establish an empirical coefficient threshold or
null distribution appropriate for the data set.

.. code-block:: shell

   phykit covarying_evolutionary_rates <tree_file_zero> <tree_file_one> -r/--reference <reference_tree_file> [-v/--verbose] [--plot] [--plot-output <path>]
       [--fig-width <float>] [--fig-height <float>] [--dpi <int>] [--no-title] [--title <str>]
       [--legend-position <str>] [--ylabel-fontsize <float>] [--xlabel-fontsize <float>]
       [--title-fontsize <float>] [--axis-fontsize <float>] [--colors <str>] [--ladderize] [--cladogram] [--circular] [--color-file <file>] [--json]

Options: |br|
*<tree_file_zero>*: first argument after function name should be a tree file |br|
*<tree_file_one>*: second argument after function name should be a tree file |br| 
*-r/--reference*: a tree whose corresponding branch lengths normalize the two
gene trees. Typically, this is a putative species tree. |br|
*-v/--verbose*: print out corrected branch lengths shared between tree 0 and tree 1 |br|
*--plot*: save a covarying-rates scatter plot |br|
*--plot-output*: output path for plot (default: ``covarying_rates_plot.png``) |br|
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
