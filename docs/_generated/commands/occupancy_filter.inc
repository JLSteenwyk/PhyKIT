Synopsis
^^^^^^^^

.. code-block:: text

   phykit occupancy_filter --list <list> [--format <format>] [--threshold <threshold>] [--output-dir <output_dir>] [--suffix <suffix>] [--json]

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
   * - ``-l, --list``
     - true
     - str
     - required
     - any
   * - ``-f, --format``
     - false
     - str
     - fasta
     - fasta, trees
   * - ``-t, --threshold``
     - false
     - float
     - 0.5
     - any
   * - ``-o, --output-dir``
     - false
     - str
     - none
     - any
   * - ``--suffix``
     - false
     - str
     - .filtered
     - any
   * - ``--json``
     - false
     - boolean
     - false
     - any

Output and errors
^^^^^^^^^^^^^^^^^

``--json`` provides the command's structured JSON representation. Output-file options: ``--output-dir``. Invalid command syntax exits with status 2. Input
validation and scientific limitations are described in the guidance below.
