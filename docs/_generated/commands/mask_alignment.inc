Synopsis
^^^^^^^^

.. code-block:: text

   phykit mask_alignment <alignment> [--max_gap <max_gap>] [--min_occupancy <min_occupancy>] [--max_entropy <max_entropy>] [--json]

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
   * - ``-g, --max_gap``
     - false
     - float
     - 1.0
     - any
   * - ``-o, --min_occupancy``
     - false
     - float
     - 0.0
     - any
   * - ``-e, --max_entropy``
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

``--json`` provides the command's structured JSON representation. Unless the guidance below states otherwise, results are emitted as command output. Invalid command syntax exits with status 2. Input
validation and scientific limitations are described in the guidance below.
