.. _cmd-fit_discrete:
.. _command-fit_discrete:

Discrete trait evolution model comparison (fitDiscrete)
=======================================================

Compare ER, SYM, ARD Mk models

Command identity
----------------

:Canonical command: ``fit_discrete``
:Handler: ``fit_discrete``
:Aliases: fd, fitdiscrete
:Standalone executables: pk_fit_discrete, pk_fd, pk_fitdiscrete
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/fit_discrete.inc

Guidance, interpretation, and examples
--------------------------------------

Compare models of discrete trait evolution on a phylogeny. Fits ER
(Equal Rates), SYM (Symmetric), and ARD (All Rates Different) Mk
models of discrete character evolution via maximum likelihood.
Compares models using AIC and BIC.

Analogous to R's geiger::fitDiscrete(). Cross-validated against R's
geiger package.

.. code-block:: shell

	phykit fit_discrete -t <tree> -d <trait_data> -c <trait>
		[--models ER,SYM,ARD] [--json]

Options: |br|
*-t/--tree*: tree file (required) |br|
*-d/--trait_data*: trait data file in TSV format (required) |br|
*-c/--trait*: column name for the discrete trait in the data file (required) |br|
*--models*: comma-separated list of models to fit (default: ``ER,SYM,ARD``) |br|
*--json*: optional argument to print results as JSON

**R validation:** Validated against ``geiger`` in R
(see ``tests/r_validation/validate_fit_discrete.R``).
