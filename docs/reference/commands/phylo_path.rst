.. _cmd-phylo_path:
.. _command-phylo_path:

Phylogenetic path analysis
==========================

Compare causal DAGs via d-separation + PGLS (von Hardenberg & Gonzalez-Voyer 2013)

Command identity
----------------

:Canonical command: ``phylo_path``
:Handler: ``phylo_path``
:Aliases: phylopath, ppath
:Standalone executables: pk_phylo_path, pk_phylopath, pk_ppath
:Categories: Phylogenetic comparative methods

Runtime interface
-----------------

.. include:: /_generated/commands/phylo_path.inc

Guidance, interpretation, and examples
--------------------------------------

Phylogenetic path analysis (von Hardenberg & Gonzalez-Voyer 2013).
Compare competing causal hypotheses (DAGs) about trait evolution while
accounting for phylogenetic non-independence, using Shipley's
d-separation test combined with PGLS.

**How it works:**

1. Define a set of candidate causal models as directed acyclic graphs (DAGs).
   Each DAG specifies which traits causally influence which other traits.
2. For each DAG, derive the **basis set** — the minimal set of conditional
   independence claims implied by the graph (Shipley 2000).
3. Test each independence claim using PGLS with Pagel's lambda. If the
   predictor is non-significant after conditioning, the independence holds.
4. Combine p-values into **Fisher's C** statistic. A model with a high p-value
   (C-test) is consistent with the data; a low p-value means the model is rejected.
5. Rank models by **CICc** (C-statistic Information Criterion, corrected for
   sample size). Lower CICc = better fit.
6. Estimate path coefficients (standardized PGLS regression slopes) and
   **conditionally average** them across well-supported models weighted by CICc.

**Interpreting the output:**

- **CICc and weights**: models with delta CICc < 2 have substantial support.
  Weights sum to 1 and represent relative evidence.
- **Model p-value**: p > 0.05 means the data are consistent with the DAG's
  independence claims. A rejected model (p < 0.05) has causal structure
  contradicted by the data.
- **Path coefficients**: standardized (z-scored) regression slopes. Positive
  = positive causal effect; magnitude indicates strength. CIs not
  overlapping zero indicate significance.

.. image:: /_static/phylo_path_dag.png
   :align: center
   :width: 60%

*DAG plot showing model-averaged path coefficients. Blue = positive
effect, line thickness proportional to coefficient magnitude.*

**Model file format:**

.. code-block:: text

   # model_name: edge1, edge2, ...
   direct: body_mass->brain_size, body_mass->longevity
   mediated: body_mass->brain_size, brain_size->longevity
   full: body_mass->brain_size, body_mass->longevity, brain_size->longevity

.. code-block:: shell

   phykit phylo_path -t <tree> --traits <traits_file> --models <models_file>
       [--best-only] [--plot-output <file>] [--csv <file>] [--json]

Options: |br|
*-t/--tree*: species tree file (required) |br|
*--traits*: TSV file with taxon and continuous trait columns (required) |br|
*--models*: model definition file with candidate DAGs (required) |br|
*--best-only*: report only best model coefficients instead of model averaging |br|
*--plot-output*: output DAG plot file (.png, .pdf, .svg) |br|
*--csv*: output CSV with model comparison table and path coefficients |br|
*--json*: output results as JSON

**R validation:** Fisher's C, CICc, weights, and path coefficients validated
against ``phylopath::phylo_path()`` in R
(see ``tests/r_validation/validate_phylo_path.R``).
