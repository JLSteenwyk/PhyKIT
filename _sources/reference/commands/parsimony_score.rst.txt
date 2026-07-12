.. _cmd-parsimony_score:
.. _command-parsimony_score:

Parsimony score
===============

Fitch parsimony score of a tree given an alignment

Command identity
----------------

:Canonical command: ``parsimony_score``
:Handler: ``parsimony_score``
:Aliases: pars, parsimony
:Standalone executables: pk_parsimony_score, pk_pars, pk_parsimony
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/parsimony_score.inc

Guidance, interpretation, and examples
--------------------------------------

Compute the Fitch (1971) maximum parsimony score of a tree given an
alignment. The parsimony score is the minimum number of character state
changes required to explain the alignment on the given tree topology.
Gap characters (-, N, X, ?) are treated as wildcards.

Cross-validated against R's phangorn::parsimony(method="fitch").

.. code-block:: shell

	phykit parsimony_score -t <tree> -a <alignment> [-v/--verbose] [--json]

Options: |br|
*-t/--tree*: tree file (required) |br|
*-a/--alignment*: alignment file in FASTA format (required) |br|
*-v/--verbose*: print per-site parsimony scores |br|
*--json*: optional argument to print results as JSON

**R validation:** Validated against ``phangorn`` in R
(see ``tests/r_validation/validate_parsimony.R``).
