.. _tutorial-10:
.. _tutorial-phylogenetic-glm-for-binary-and-count-data:

Tutorial 10: Phylogenetic GLM for binary and count data
=======================================================

Objectives
----------

- Complete the phylogenetic glm for binary and count data workflow.
- Interpret the reported values and generated artifacts in their scientific context.
- Identify the canonical command references for each analysis step.

Prerequisites and working directory
-----------------------------------

Install the current PhyKIT release and create a dedicated working directory.
Download the data linked in this tutorial into that directory before running
the commands. All paths below are relative to this directory.

.. code-block:: shell

   mkdir phykit-tutorial-10
   cd phykit-tutorial-10

Related command references
--------------------------

- :doc:`Phylogenetic Glm </reference/commands/phylogenetic_glm>`

Workflow
--------

Standard PGLS handles continuous response variables. When the response is binary
(e.g., presence/absence) or count data (e.g., number of offspring), a Generalized
Linear Model is needed. Phylogenetic GLM extends GLM to account for phylogenetic
non-independence among species.

**Hypothetical study question.** Is body mass a significant predictor of a binary
trait (e.g., dietary specialization) or a count trait (e.g., litter size) after
accounting for phylogenetic relationships?

PhyKIT's ``phylogenetic_glm`` command (aliases: ``phylo_glm``, ``pglm``) supports
two families: binomial (logistic MPLE) and Poisson (GEE). |br|

.. centered::
   Download test data:
   :download:`Mammal phylogeny </data/tree_simple.tre>`;
   :download:`Binary and count trait data </data/tree_simple_glm_traits.tsv>`

|

Step 1: Fit a binomial (logistic) model
****************************************

Test whether body mass predicts a binary trait:

.. code-block:: shell

   phykit phylogenetic_glm \
       -t tree_simple.tre \
       -d tree_simple_glm_traits.tsv \
       -y binary_trait \
       -x body_mass \
       --family binomial

.. code-block:: text

   Phylogenetic GLM (Logistic MPLE)

   Formula: binary_trait ~ body_mass
   Family: binomial, Method: logistic_MPLE

   Estimated alpha: 0.0183

   Coefficients:
                           Estimate   Std.Error     z-value     p-value
   (Intercept)             -10.0000     23.6347     -0.4231    0.672218
   body_mass                10.0000     22.0222      0.4541    0.649766
   ---
   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1

   Log-likelihood: -5.3694    AIC: 16.7387
   Number of observations: 8

The logistic MPLE method jointly estimates regression coefficients and the
phylogenetic signal parameter alpha. The coefficients hit the btol boundary
in this example due to quasi-complete separation in the small dataset.

Step 2: Fit a Poisson model for count data
*******************************************

Test whether body mass predicts count data:

.. code-block:: shell

   phykit phylogenetic_glm \
       -t tree_simple.tre \
       -d tree_simple_glm_traits.tsv \
       -y count_trait \
       -x body_mass \
       --family poisson

.. code-block:: text

   Phylogenetic GLM (Poisson GEE)

   Formula: count_trait ~ body_mass
   Family: poisson, Method: poisson_GEE

   Overdispersion (phi): 0.1730

   Coefficients:
                           Estimate   Std.Error     z-value     p-value
   (Intercept)               0.6741      0.1678      4.0176    0.000059    ***
   body_mass                 0.5968      0.0877      6.8082    0.000000    ***
   ---
   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1

   Log-likelihood: -13.4024    AIC: 30.8048
   Number of observations: 8

The Poisson GEE uses phylogenetic correlations derived from the tree and
reports the overdispersion parameter phi. Both predictors are highly significant.

Step 3: Export results as JSON
******************************

.. code-block:: shell

   phykit phylogenetic_glm \
       -t tree_simple.tre \
       -d tree_simple_glm_traits.tsv \
       -y count_trait \
       -x body_mass \
       --family poisson \
       --json

Summary
*******

In this tutorial, we used phylogenetic GLMs to model binary and count response
variables while accounting for phylogenetic non-independence. The key steps were:
(1) fitting a logistic model for binary data with the binomial family,
(2) fitting a Poisson model for count data, and (3) exporting results as JSON.
Phylogenetic GLM complements PGLS by handling non-continuous response variables.

For methodological details, see
`Ives and Garland (2010) <https://doi.org/10.1086/656006>`_ for logistic MPLE
and `Paradis and Claude (2002) <https://doi.org/10.1046/j.1420-9101.2002.00405.x>`_
for Poisson GEE. The R equivalent is ``phylolm::phyloglm()``.

|

Expected artifacts
------------------

Each step identifies its expected terminal output or generated files. Confirm
that those artifacts exist before continuing to the next step; filenames are
relative to the tutorial working directory unless an absolute path is shown.

Troubleshooting
---------------

- Run ``phykit <command> --help`` to compare an invocation with the live interface.
- Confirm that downloaded files are in the current working directory and retain
  the filenames shown in the tutorial.
- For parsing errors, compare taxon names exactly across alignments, trees, and
  trait tables, including capitalization and underscores.
- See :doc:`Troubleshooting </troubleshooting/index>` for installation, format,
  and error-reporting guidance.
