Synopsis
^^^^^^^^

.. code-block:: text

   phykit codon_dnds <alignment> [--method <method>] [--genetic-code <genetic_code>] [--kappa <kappa>] [--codon-frequency <codon_frequency>] [--stop-policy <stop_policy>] [--reference <reference>] [--verbose] [--json]

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
   * - ``alignment``
     - true
     - str
     - required
     - any
   * - ``--method``
     - false
     - upper
     - NG86
     - NG86, LWL85, YN00, ML
   * - ``--genetic-code``
     - false
     - int
     - 1
     - any
   * - ``--kappa``
     - false
     - float
     - 1.0
     - any
   * - ``--codon-frequency``
     - false
     - upper
     - F3X4
     - F1X4, F3X4, F61
   * - ``--stop-policy``
     - false
     - string
     - error
     - error, skip
   * - ``--reference``
     - false
     - str
     - none
     - any
   * - ``-v, --verbose``
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
