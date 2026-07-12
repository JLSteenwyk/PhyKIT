Synopsis
^^^^^^^^

.. code-block:: text

   phykit evolutionary_rate_per_site <alignment> [--plot] [--plot-output <plot_output>] [--fig-width <fig_width>] [--fig-height <fig_height>] [--dpi <dpi>] [--no-title] [--title <title>] [--legend-position <legend_position>] [--ylabel-fontsize <ylabel_fontsize>] [--xlabel-fontsize <xlabel_fontsize>] [--title-fontsize <title_fontsize>] [--axis-fontsize <axis_fontsize>] [--colors <colors>] [--ladderize] [--cladogram] [--circular] [--color-file <color_file>] [--json]

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
   * - ``alignment``
     - true
     - str
     - required
     - any
   * - ``--plot``
     - false
     - boolean
     - false
     - any
   * - ``--plot-output``
     - false
     - str
     - evolutionary_rate_per_site_plot.png
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
