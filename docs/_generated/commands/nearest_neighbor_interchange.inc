Synopsis
^^^^^^^^

.. code-block:: text

   phykit nearest_neighbor_interchange <tree> [--output <output>] [--branch <branch>] [--branches <branches>] [--no-input-tree] [--json]

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
   * - ``tree``
     - true
     - str
     - required
     - any
   * - ``-o, --output``
     - false
     - str
     - none
     - any
   * - ``--branch``
     - false
     - str
     - none
     - any
   * - ``--branches``
     - false
     - str
     - none
     - any
   * - ``--no-input-tree``
     - false
     - boolean
     - false
     - any
   * - ``--json``
     - false
     - boolean
     - false
     - any

Output and errors
^^^^^^^^^^^^^^^^^

``--json`` provides the command's structured JSON representation. Output-file options: ``--output``. Invalid command syntax exits with status 2. Input
validation and scientific limitations are described in the guidance below.
