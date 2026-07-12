Synopsis
^^^^^^^^

.. code-block:: text

   phykit dstatistic [--alignment <alignment>] [--gene-trees <gene_trees>] --p1 <p1> --p2 <p2> --p3 <p3> --outgroup <outgroup> [--block-size <block_size>] [--support <support>] [--json]

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
   * - ``-a, --alignment``
     - false
     - str
     - none
     - any
   * - ``-g, --gene-trees``
     - false
     - str
     - none
     - any
   * - ``--p1``
     - true
     - str
     - required
     - any
   * - ``--p2``
     - true
     - str
     - required
     - any
   * - ``--p3``
     - true
     - str
     - required
     - any
   * - ``--outgroup``
     - true
     - str
     - required
     - any
   * - ``--block-size``
     - false
     - int
     - 100
     - any
   * - ``--support``
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
