.. _cmd-evolutionary_rate:
.. _command-evolutionary_rate:

Evolutionary rate
=================

Calculate tree-based evolutionary rate

Command identity
----------------

:Canonical command: ``evolutionary_rate``
:Handler: ``evolutionary_rate``
:Aliases: evo_rate
:Standalone executables: pk_evolutionary_rate, pk_evo_rate
:Categories: Tree summary statistics

Runtime interface
-----------------

.. include:: /_generated/commands/evolutionary_rate.inc

Guidance, interpretation, and examples
--------------------------------------

Calculate a tree-based estimation of the evolutionary rate of a gene.

Evolutionary rate is the total tree length divided by the number
of terminals.

Calculate evolutionary rate following Telford et al., Proceedings
of the Royal Society B (2014). 

.. code-block:: shell

   phykit evolutionary_rate <tree> [--json]

Options: |br|
*<tree>*: input file tree name |br|
*--json*: optional argument to print results as JSON
