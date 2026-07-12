Synopsis
^^^^^^^^

.. code-block:: text

   phykit phylogenetic_signal --tree <tree> --trait_data <trait_data> [--method <method>] [--permutations <permutations>] [--gene-trees <gene_trees>] [--multivariate] [--json]

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
   * - ``-d, --trait_data``
     - true
     - str
     - required
     - any
   * - ``-m, --method``
     - false
     - str
     - blombergs_k
     - blombergs_k, lambda
   * - ``-p, --permutations``
     - false
     - int
     - 1000
     - any
   * - ``-g, --gene-trees``
     - false
     - str
     - none
     - any
   * - ``--multivariate``
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

``--json`` provides the command's structured JSON representation. Unless the guidance below states otherwise, results are emitted as command output. Invalid command syntax exits with status 2. Input
validation and scientific limitations are described in the guidance below.
