Synopsis
^^^^^^^^

.. code-block:: text

   phykit prune_tree <tree> <list_of_taxa> [--output <output>] [--keep <keep>] [--ignore-branch-labels] [--json]

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
   * - ``list_of_taxa``
     - true
     - str
     - required
     - any
   * - ``-o, --output``
     - false
     - str
     - none
     - any
   * - ``-k, --keep``
     - false
     - str2bool
     - false
     - any
   * - ``--ignore-branch-labels``
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
