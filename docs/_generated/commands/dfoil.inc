Synopsis
^^^^^^^^

.. code-block:: text

   phykit dfoil --alignment <alignment> --p1 <p1> --p2 <p2> --p3 <p3> --p4 <p4> --outgroup <outgroup> [--json]

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
     - true
     - str
     - required
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
   * - ``--p4``
     - true
     - str
     - required
     - any
   * - ``--outgroup``
     - true
     - str
     - required
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
