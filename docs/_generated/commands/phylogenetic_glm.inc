Synopsis
^^^^^^^^

.. code-block:: text

   phykit phylogenetic_glm --tree <tree> --trait_data <trait_data> --response <response> --predictors <predictors> --family <family> [--method <method>] [--btol <btol>] [--log-alpha-bound <log_alpha_bound>] [--gene-trees <gene_trees>] [--json]

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
   * - ``-y, --response``
     - true
     - str
     - required
     - any
   * - ``-x, --predictors``
     - true
     - str
     - required
     - any
   * - ``--family``
     - true
     - str
     - required
     - binomial, poisson
   * - ``--method``
     - false
     - str
     - none
     - logistic_MPLE, poisson_GEE
   * - ``--btol``
     - false
     - float
     - 10
     - any
   * - ``--log-alpha-bound``
     - false
     - float
     - 4
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

``--json`` provides the command's structured JSON representation. Unless the guidance below states otherwise, results are emitted as command output. Invalid command syntax exits with status 2. Input
validation and scientific limitations are described in the guidance below.
