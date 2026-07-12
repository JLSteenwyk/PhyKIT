.. _cmd-chronogram:
.. _command-chronogram:

Chronogram
==========

Time-calibrated tree with geological timescale (rectangular or circular)

Command identity
----------------

:Canonical command: ``chronogram``
:Handler: ``chronogram``
:Aliases: chrono, time_tree
:Standalone executables: pk_chronogram, pk_chrono, pk_time_tree
:Categories: Tree manipulation & utilities

Runtime interface
-----------------

.. include:: /_generated/commands/chronogram.inc

Guidance, interpretation, and examples
--------------------------------------

Plot a chronogram (time-calibrated phylogeny) with geological timescale
bands. Requires an ultrametric (or approximately ultrametric) tree and
the root age in millions of years (Ma).

Geological epoch, period, or era bands are drawn behind the tree as
colored stripes based on the ICS 2024 International Chronostratigraphic
Chart. A labeled timescale bar is displayed below the tree. The time
axis runs from past (left) to present (right).

**95% HPD confidence intervals** are automatically detected and drawn
when the input tree contains BEAST (``height_95%_HPD``) or MCMCTree
(``95%HPD``) node annotations. The intervals appear as translucent
blue bars at each internal node — no extra flags needed. Trees without
annotations are plotted without bars.

The ``--timescale`` option controls the level of detail:

- **auto** (default): selects epochs for trees <= 66 Ma, periods for
  <= 252 Ma, eras for deeper timescales
- **epoch**: Cenozoic and Mesozoic epochs (Holocene through Early Triassic)
- **period**: geological periods (Quaternary through Cambrian)
- **era**: Cenozoic, Mesozoic, Paleozoic

.. image:: /_static/chronogram_ultrametric.png
   :align: center
   :width: 90%

*Rectangular chronogram with epoch bands, node age labels, and a
geological timescale bar.*

.. image:: /_static/chronogram_ultrametric_circular.png
   :align: center
   :width: 60%

*Circular chronogram with concentric geological epoch rings and
radial time tick marks.*

.. code-block:: shell

   phykit chronogram -t <tree> --root-age <float> --plot-output <file>
       [--timescale auto|epoch|period|era] [--node-ages]
       [--circular] [--ladderize] [--color-file <file>] [--json]

Options: |br|
*-t/--tree*: ultrametric tree file (required) |br|
*--root-age*: age of the root in millions of years (Ma; required) |br|
*--plot-output*: output figure path (.png, .pdf, .svg; required) |br|
*--timescale*: timescale level — ``auto`` (default), ``epoch``, ``period``, or ``era`` |br|
*--node-ages*: label internal nodes with divergence times (Ma) |br|
*--circular*: draw circular chronogram with concentric geological rings |br|
*--ladderize*: ladderize (sort) the tree before plotting |br|
*--color-file*: color annotation file for tip labels, clade ranges, and branch colors (iTOL-inspired TSV format) |br|
*--json*: output node ages as JSON
