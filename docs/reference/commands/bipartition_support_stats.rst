.. _cmd-bipartition_support_stats:
.. _command-bipartition_support_stats:

Bipartition support statistics
==============================

Summary statistics of bipartition support values

Command identity
----------------

:Canonical command: ``bipartition_support_stats``
:Handler: ``bipartition_support_stats``
:Aliases: bss
:Standalone executables: pk_bipartition_support_stats, pk_bss
:Categories: Tree summary statistics

Runtime interface
-----------------

.. include:: /_generated/commands/bipartition_support_stats.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate summary statistics for bipartition support.

High bipartition support values are thought to be desirable because
they are indicative of greater certainty in tree topology.

To obtain all bipartition support values, use the -v/--verbose option.
In addition to support values for each node, the names of all terminal
branch tips are also included. Each terminal branch name is separated
with a semi-colon (;).

.. code-block:: shell

   phykit bipartition_support_stats <tree> [-v/--verbose]
       [--thresholds <comma-separated-floats>] [--json]

Options: |br|
*<tree>*: first argument after function name should be a tree file |br|
*-v/--verbose*: optional argument to print all bipartition support values |br|
*--thresholds*: optional comma-separated support cutoffs; prints count and
fraction of bipartitions below each cutoff |br|
*--json*: optional argument to print results as JSON

Example JSON output (summary mode):

.. code-block:: shell

   phykit bipartition_support_stats test.tre --thresholds 70,90 --json
   {"summary": {"maximum": 100, "mean": 95.71428571428571, "median": 100, "minimum": 85, "seventy_fifth": 100.0, "standard_deviation": 7.319250547113999, "twenty_fifth": 92.5, "variance": 53.57142857142857}, "thresholds": [{"count_below": 0, "fraction_below": 0.0, "threshold": 70.0}, {"count_below": 2, "fraction_below": 0.2857142857142857, "threshold": 90.0}], "verbose": false}

Example JSON output (verbose mode):

.. code-block:: shell

   phykit bipartition_support_stats test.tre -v --json
   {"bipartitions": [{"support": 85, "terminals": ["taxon_a", "taxon_b"]}, {"support": 100, "terminals": ["taxon_c", "taxon_d"]}], "thresholds": [], "verbose": true}
