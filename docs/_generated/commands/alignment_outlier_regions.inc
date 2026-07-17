Synopsis
^^^^^^^^

.. code-block:: text

   phykit alignment_outlier_regions <alignment> [--cutoff <cutoff>] [--report <report>] [--mask-output <mask_output>] [--mask-character <mask_character>] [--json]

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
   * - ``--cutoff``
     - false
     - float
     - 3.0
     - any
   * - ``--report``
     - false
     - str
     - none
     - any
   * - ``--mask-output``
     - false
     - str
     - none
     - any
   * - ``--mask-character``
     - false
     - str
     - ?
     - any
   * - ``--json``
     - false
     - boolean
     - false
     - any

Output and errors
^^^^^^^^^^^^^^^^^

``--json`` provides the command's structured JSON representation. Output-file options: ``--mask-output``. Invalid command syntax exits with status 2. Input
validation and scientific limitations are described in the guidance below.
