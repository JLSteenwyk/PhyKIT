Synopsis
^^^^^^^^

.. code-block:: text

   phykit alignment_outlier_taxa <alignment> [--gap-z <gap_z>] [--composition-z <composition_z>] [--distance-z <distance_z>] [--rcvt-z <rcvt_z>] [--occupancy-z <occupancy_z>] [--entropy-z <entropy_z>] [--json]

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
   * - ``--gap-z``
     - false
     - float
     - 3.0
     - any
   * - ``--composition-z``
     - false
     - float
     - 3.0
     - any
   * - ``--distance-z``
     - false
     - float
     - 3.0
     - any
   * - ``--rcvt-z``
     - false
     - float
     - 3.0
     - any
   * - ``--occupancy-z``
     - false
     - float
     - 3.0
     - any
   * - ``--entropy-z``
     - false
     - float
     - 3.0
     - any
   * - ``--json``
     - false
     - boolean
     - false
     - any

Output and errors
^^^^^^^^^^^^^^^^^

``--json`` provides the command's structured JSON representation. Unless the guidance below states otherwise, results are emitted as command output. Invalid command syntax exits with status 2. Input
validation and scientific limitations are described in the guidance below.
