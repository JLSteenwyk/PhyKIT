Synopsis
^^^^^^^^

.. code-block:: text

   phykit phylo_impute --tree <tree> --trait-data <trait_data> --output <output> [--gene-trees <gene_trees>] [--json]

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
   * - ``-d, --trait-data``
     - true
     - str
     - required
     - any
   * - ``-o, --output``
     - true
     - str
     - required
     - any
   * - ``-g, --gene-trees``
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

``--json`` provides the command's structured JSON representation. Output-file options: ``--output``. Invalid command syntax exits with status 2. Input
validation and scientific limitations are described in the guidance below.
