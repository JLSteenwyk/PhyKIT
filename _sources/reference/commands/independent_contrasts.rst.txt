.. _cmd-independent_contrasts:
.. _command-independent_contrasts:

Independent contrasts (PIC)
===========================

Felsenstein's phylogenetically independent contrasts

Command identity
----------------

:Canonical command: ``independent_contrasts``
:Handler: ``independent_contrasts``
:Aliases: phylo_contrasts, pic
:Standalone executables: pk_independent_contrasts, pk_phylo_contrasts, pk_pic
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/independent_contrasts.inc

Guidance, interpretation, and examples
--------------------------------------

Compute Felsenstein's (1985) phylogenetically independent contrasts for a
continuous trait on a phylogeny. Each internal node yields one standardized
contrast, producing n-1 contrasts for n tips. Multifurcations are
automatically resolved.

Cross-validated against R's ape::pic(). The sum of squared contrasts
matches R exactly.

.. code-block:: shell

	phykit independent_contrasts -t <tree> -d <trait_data> [--json]

Options: |br|
*-t/--tree*: tree file (required) |br|
*-d/--trait_data*: trait data file, two columns: taxon<tab>value (required) |br|
*--json*: optional argument to print results as JSON

**R validation:** Validated against ``ape`` in R
(see ``tests/r_validation/validate_pic.R``).
