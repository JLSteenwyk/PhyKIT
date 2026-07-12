.. _cmd-phylo_gwas:
.. _command-phylo_gwas:

Phylo GWAS
==========

Phylogenetic genome-wide association study

Command identity
----------------

:Canonical command: ``phylo_gwas``
:Handler: ``phylo_gwas``
:Aliases: pgwas
:Standalone executables: pk_phylo_gwas, pk_pgwas
:Categories: Phylogenetic comparative methods

Runtime interface
-----------------

.. include:: /_generated/commands/phylo_gwas.inc

Guidance, interpretation, and examples
--------------------------------------

Phylogenetic genome-wide association study following the Pease et al.
(2016) approach. Performs per-site association tests between alignment
columns and a phenotype, applies Benjamini-Hochberg FDR correction,
optionally classifies significant associations as monophyletic or
polyphyletic using a phylogenetic tree, and produces a Manhattan plot.

Categorical phenotypes use Fisher's exact test (2 groups) or chi-squared
test (>2 groups). Continuous phenotypes use point-biserial correlation.
Only biallelic sites are tested; invariant and multiallelic sites are
skipped. Sites with gaps or ambiguous characters are also skipped.

.. code-block:: shell

   phykit phylo_gwas -a <alignment> -d <phenotype> -o <output>
     [-t <tree>] [-p <partition>] [--alpha 0.05]
     [--exclude-monophyletic] [--csv <file>] [--json]

Options: |br|
*-a/--alignment*: FASTA alignment file |br|
*-d/--phenotype*: two-column TSV file (taxon<tab>phenotype) |br|
*-o/--output*: output Manhattan plot path |br|
*-t/--tree*: optional Newick tree for monophyletic/polyphyletic classification |br|
*-p/--partition*: optional RAxML-style partition file for gene annotations |br|
*--alpha*: FDR significance threshold (default: 0.05) |br|
*--exclude-monophyletic*: exclude monophyletic associations from results |br|
*--dot-size*: scale factor for dot size in the Manhattan plot (default: 1.0; use 2.0 for double, 0.5 for half) |br|
*--csv*: output per-site results as CSV to the specified file |br|
*--fig-width*: figure width in inches (auto-scaled if omitted) |br|
*--fig-height*: figure height in inches (auto-scaled if omitted) |br|
*--dpi*: resolution in DPI (default: 300) |br|
*--no-title*: hide the plot title |br|
*--title*: custom title text |br|
*--legend-position*: legend location (e.g., "upper right", "none") |br|
*--colors*: comma-separated colors (hex or named) |br|
*--json*: optional argument to print results as JSON
