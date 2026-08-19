Synopsis
^^^^^^^^

.. code-block:: text

   phykit episodic_rate_covariation <tree_zero> <tree_one> --reference <reference> [--verbose] [--weighting <weighting>] [--max-rate <max_rate>] [--max-pairs <max_pairs>] [--seed <seed>] [--permutations <permutations>] [--depth-bins <depth_bins>] [--min-edges <min_edges>] [--alpha <alpha>] [--output <output>] [--annotated-tree <annotated_tree>] [--json]

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
   * - ``tree_zero``
     - true
     - str
     - required
     - any
   * - ``tree_one``
     - true
     - str
     - required
     - any
   * - ``-r, --reference``
     - true
     - str
     - required
     - any
   * - ``-v, --verbose``
     - false
     - boolean
     - false
     - any
   * - ``--weighting``
     - false
     - string
     - uniform
     - uniform, inverse_reference
   * - ``--max-rate``
     - false
     - float
     - 5.0
     - any
   * - ``--max-pairs``
     - false
     - int
     - 50000
     - any
   * - ``--seed``
     - false
     - int
     - 0
     - any
   * - ``--permutations``
     - false
     - int
     - 999
     - any
   * - ``--depth-bins``
     - false
     - int
     - 4
     - any
   * - ``--min-edges``
     - false
     - int
     - 3
     - any
   * - ``--alpha``
     - false
     - float
     - 0.05
     - any
   * - ``-o, --output``
     - false
     - str
     - none
     - any
   * - ``--annotated-tree``
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

``--json`` provides the command's structured JSON representation. Output-file options: ``--output``, ``--annotated-tree``. Invalid command syntax exits with status 2. Input
validation and scientific limitations are described in the guidance below.
