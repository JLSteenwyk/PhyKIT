.. _cmd-phylo_impute:
.. _command-phylo_impute:

Phylogenetic imputation
=======================

Impute missing trait values using phylogenetic relationships

Command identity
----------------

:Canonical command: ``phylo_impute``
:Handler: ``phylo_impute``
:Aliases: impute, phylo_imp
:Standalone executables: pk_phylo_impute, pk_impute, pk_phylo_imp
:Categories: Trait evolution

Runtime interface
-----------------

.. include:: /_generated/commands/phylo_impute.inc

Guidance, interpretation, and examples
--------------------------------------

Impute missing continuous trait values using phylogenetic relationships
and between-trait correlations under Brownian motion. Missing values
(``NA``, ``?``, or empty) in a multi-trait TSV are predicted from:

1. **Phylogenetic neighbors**: closely related species with observed
   values contribute more to the imputation via the phylogenetic VCV
2. **Trait correlations**: if a taxon has observed values for other
   traits, between-trait covariance improves the prediction

Reports imputed values with standard errors and 95% confidence
intervals. The output TSV is a drop-in replacement for the input with
all missing values filled in.

.. code-block:: shell

   phykit phylo_impute -t <tree> -d <trait_data> -o <output> [--json]
   phykit phylo_impute -t <tree> -d <trait_data> -o <output> -g <gene_trees>

Options: |br|
*-t/--tree*: tree file in Newick format (required) |br|
*-d/--trait-data*: multi-trait TSV with header; missing values as NA, ?, or empty (required) |br|
*-o/--output*: output TSV file with imputed values (required) |br|
*-g/--gene-trees*: gene trees for discordance-aware VCV |br|
*--json*: optional argument to print results as JSON

**R validation:** Validated against ``Rphylopars`` in R
(see ``tests/r_validation/validate_phylo_impute.R``).
