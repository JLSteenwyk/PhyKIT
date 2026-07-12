Synopsis
^^^^^^^^

.. code-block:: text

   phykit consensus_tree --trees <trees> [--method <method>] [--missing-taxa <missing_taxa>] [--json]

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
   * - ``-t, --trees``
     - true
     - str
     - required
     - any
   * - ``-m, --method``
     - false
     - str
     - majority
     - strict, majority
   * - ``--missing-taxa``
     - false
     - str
     - error
     - error, shared
   * - ``--json``
     - false
     - boolean
     - false
     - any

Output and errors
^^^^^^^^^^^^^^^^^

``--json`` provides the command's structured JSON representation. Unless the guidance below states otherwise, results are emitted as command output. Invalid command syntax exits with status 2. Input
validation and scientific limitations are described in the guidance below.
