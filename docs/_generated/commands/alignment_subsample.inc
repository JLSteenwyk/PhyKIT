Synopsis
^^^^^^^^

.. code-block:: text

   phykit alignment_subsample --mode <mode> [--alignment <alignment>] [--list <list>] [--partition <partition>] [--number <number>] [--fraction <fraction>] [--seed <seed>] [--bootstrap] [--output <output>] [--json]

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
   * - ``--mode``
     - true
     - str
     - required
     - genes, partitions, sites
   * - ``-a, --alignment``
     - false
     - str
     - none
     - any
   * - ``-l, --list``
     - false
     - str
     - none
     - any
   * - ``-p, --partition``
     - false
     - str
     - none
     - any
   * - ``--number``
     - false
     - int
     - none
     - any
   * - ``--fraction``
     - false
     - float
     - none
     - any
   * - ``--seed``
     - false
     - int
     - none
     - any
   * - ``--bootstrap``
     - false
     - boolean
     - false
     - any
   * - ``-o, --output``
     - false
     - str
     - subsampled
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
