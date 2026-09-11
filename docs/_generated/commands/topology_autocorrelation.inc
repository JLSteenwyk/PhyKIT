Synopsis
^^^^^^^^

.. code-block:: text

   phykit topology_autocorrelation [--classified-genes <classified_genes>] [--manifest <manifest>] [--coordinates <coordinates>] [--groups <groups>] [--reference-tree <reference_tree>] [--branch-taxa <branch_taxa>] [--outgroup <outgroup>] [--min-support <min_support>] [--support-scale <support_scale>] [--missing-support <missing_support>] [--unmapped <unmapped>] [--labels <labels>] [--chromosome <chromosome>] [--interval <interval>] [--distance-edges <distance_edges>] [--bin-width <bin_width>] [--max-distance <max_distance>] [--block-sizes <block_sizes>] [--replicates <replicates>] [--seed <seed>] --output-prefix <output_prefix> [--plot] [--plot-output <plot_output>] [--fig-width <fig_width>] [--fig-height <fig_height>] [--dpi <dpi>] [--title <title>] [--no-title] [--legend-position <legend_position>] [--colors <colors>] [--xlabel-fontsize <xlabel_fontsize>] [--ylabel-fontsize <ylabel_fontsize>] [--axis-fontsize <axis_fontsize>] [--title-fontsize <title_fontsize>] [--json]

Arguments
^^^^^^^^^

This table is generated from the live command parser. It is the authoritative
source for accepted spellings, required arguments, types, defaults, and choices.

.. list-table::
   :header-rows: 1
   :widths: 28 10 12 18 20

   * - Argument
     - Required
     - Type
     - Default
     - Choices
   * - ``--classified-genes``
     - false
     - string
     - none
     - any
   * - ``--manifest``
     - false
     - string
     - none
     - any
   * - ``--coordinates``
     - false
     - string
     - none
     - any
   * - ``--groups``
     - false
     - string
     - none
     - any
   * - ``--reference-tree``
     - false
     - string
     - none
     - any
   * - ``--branch-taxa``
     - false
     - string
     - none
     - any
   * - ``--outgroup``
     - false
     - string
     - none
     - any
   * - ``--min-support``
     - false
     - float
     - none
     - any
   * - ``--support-scale``
     - false
     - int
     - none
     - 1, 100
   * - ``--missing-support``
     - false
     - string
     - none
     - collapse, keep, error
   * - ``--unmapped``
     - false
     - string
     - none
     - error, skip
   * - ``--labels``
     - false
     - string
     - none
     - any
   * - ``--chromosome``
     - false
     - string
     - none
     - any
   * - ``--interval``
     - false
     - int
     - none
     - any
   * - ``--distance-edges``
     - false
     - int
     - none
     - any
   * - ``--bin-width``
     - false
     - int
     - none
     - any
   * - ``--max-distance``
     - false
     - int
     - none
     - any
   * - ``--block-sizes``
     - false
     - int
     - none
     - any
   * - ``--replicates``
     - false
     - int
     - 999
     - any
   * - ``--seed``
     - false
     - int
     - 0
     - any
   * - ``--output-prefix, -o``
     - true
     - string
     - required
     - any
   * - ``--plot``
     - false
     - boolean
     - false
     - any
   * - ``--plot-output``
     - false
     - string
     - none
     - any
   * - ``--fig-width``
     - false
     - float
     - none
     - any
   * - ``--fig-height``
     - false
     - float
     - none
     - any
   * - ``--dpi``
     - false
     - int
     - 300
     - any
   * - ``--title``
     - false
     - string
     - none
     - any
   * - ``--no-title``
     - false
     - boolean
     - false
     - any
   * - ``--legend-position``
     - false
     - string
     - none
     - any
   * - ``--colors``
     - false
     - string
     - none
     - any
   * - ``--xlabel-fontsize``
     - false
     - float
     - none
     - any
   * - ``--ylabel-fontsize``
     - false
     - float
     - none
     - any
   * - ``--axis-fontsize``
     - false
     - float
     - none
     - any
   * - ``--title-fontsize``
     - false
     - float
     - none
     - any
   * - ``--json``
     - false
     - boolean
     - false
     - any

Output and errors
^^^^^^^^^^^^^^^^^

``--json`` provides the command's structured JSON representation. Output-file options: ``--output-prefix``, ``--plot-output``. Invalid command syntax exits with status 2. Input
validation and scientific limitations are described in the guidance below.
