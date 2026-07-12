Synopsis
^^^^^^^^

.. code-block:: text

   phykit phylo_anova --tree <tree> --traits <traits> [--group-column <group_column>] [--method <method>] [--permutations <permutations>] [--pairwise] [--plot-output <plot_output>] [--plot-type <plot_type>] [--seed <seed>] [--fig-width <fig_width>] [--fig-height <fig_height>] [--dpi <dpi>] [--no-title] [--title <title>] [--legend-position <legend_position>] [--ylabel-fontsize <ylabel_fontsize>] [--xlabel-fontsize <xlabel_fontsize>] [--title-fontsize <title_fontsize>] [--axis-fontsize <axis_fontsize>] [--colors <colors>] [--ladderize] [--cladogram] [--circular] [--color-file <color_file>] [--json]

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
   * - ``-t, --tree``
     - true
     - str
     - required
     - any
   * - ``--traits``
     - true
     - str
     - required
     - any
   * - ``--group-column``
     - false
     - str
     - none
     - any
   * - ``--method``
     - false
     - str
     - auto
     - auto, anova, manova
   * - ``--permutations``
     - false
     - int
     - 1000
     - any
   * - ``--pairwise``
     - false
     - boolean
     - false
     - any
   * - ``--plot-output``
     - false
     - str
     - none
     - any
   * - ``--plot-type``
     - false
     - str
     - auto
     - auto, boxplot, phylomorphospace
   * - ``--seed``
     - false
     - int
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
   * - ``--no-title``
     - false
     - boolean
     - false
     - any
   * - ``--title``
     - false
     - str
     - none
     - any
   * - ``--legend-position``
     - false
     - str
     - none
     - any
   * - ``--ylabel-fontsize``
     - false
     - float
     - none
     - any
   * - ``--xlabel-fontsize``
     - false
     - float
     - none
     - any
   * - ``--title-fontsize``
     - false
     - float
     - none
     - any
   * - ``--axis-fontsize``
     - false
     - float
     - none
     - any
   * - ``--colors``
     - false
     - str
     - none
     - any
   * - ``--ladderize``
     - false
     - boolean
     - false
     - any
   * - ``--cladogram``
     - false
     - boolean
     - false
     - any
   * - ``--circular``
     - false
     - boolean
     - false
     - any
   * - ``--color-file``
     - false
     - str
     - none
     - any
   * - ``--json``
     - false
     - boolean
     - false
     - any

Output and errors
^^^^^^^^^^^^^^^^^

``--json`` provides the command's structured JSON representation. Output-file options: ``--plot-output``. Invalid command syntax exits with status 2. Input
validation and scientific limitations are described in the guidance below.
