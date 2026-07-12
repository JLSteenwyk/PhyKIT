.. _tutorial-09:
.. _tutorial-phylogenetic-regression-pgls:

Tutorial 9: Phylogenetic regression (PGLS)
==========================================

Objectives
----------

- Complete the phylogenetic regression (pgls) workflow.
- Interpret the reported values and generated artifacts in their scientific context.
- Identify the canonical command references for each analysis step.

Prerequisites and working directory
-----------------------------------

Install the current PhyKIT release and create a dedicated working directory.
Download the data linked in this tutorial into that directory before running
the commands. All paths below are relative to this directory.

.. code-block:: shell

   mkdir phykit-tutorial-09
   cd phykit-tutorial-09

Related command references
--------------------------

- :doc:`Phylogenetic Regression </reference/commands/phylogenetic_regression>`

Workflow
--------

Standard linear regression assumes that data points are independent and identically
distributed. In comparative biology, species are not independent because they share
evolutionary history. Phylogenetic Generalized Least Squares (PGLS) addresses this by
incorporating the phylogenetic variance-covariance matrix into the regression model,
properly accounting for the expected covariance among species due to shared ancestry
(`Grafen 1989 <https://doi.org/10.1098/rstb.1989.0106>`_;
`Martins and Hansen 1997 <https://doi.org/10.1086/286013>`_;
`Freckleton et al. 2002 <https://doi.org/10.1086/343873>`_).

**Hypothetical study question.** A classic question in comparative biology is the
relationship between brain size and body mass across species. Specifically: does brain
size predict body mass after accounting for the phylogenetic non-independence among
species? Is the relationship significant, and how much variance does it explain?

PhyKIT's ``phylogenetic_regression`` command (aliases: ``phylo_regression``, ``pgls``)
fits PGLS regressions under a Brownian motion model. |br|

.. centered::
   Download test data:
   :download:`Mammal phylogeny </data/tree_simple.tre>`;
   :download:`Multi-trait data </data/tree_simple_multi_traits.tsv>`

|

Step 0: Prepare data
********************

Three input files/specifications are needed: a phylogenetic tree, a tab-delimited
multi-trait file with a header row, and the specification of response (``-y``) and
predictor (``-x``) variables. The trait file is the same format used for phylogenetic
PCA and phylomorphospace.

|

Step 1: Run a simple PGLS regression
*************************************

Test whether brain size predicts body mass:

.. code-block:: shell

   phykit phylogenetic_regression \
       -t tree_simple.tre \
       -d tree_simple_multi_traits.tsv \
       -y body_mass \
       -x brain_size

.. code-block:: text

   Phylogenetic Generalized Least Squares (PGLS)

   Formula: body_mass ~ brain_size

   Coefficients:
                           Estimate   Std.Error     t-value     p-value
   (Intercept)              -1.3350      0.2000     -6.6733    0.000548    ***
   brain_size                1.3778      0.0877     15.7140    0.000004    ***
   ---
   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1

   Residual standard error: 0.0349 on 6 degrees of freedom
   Multiple R-squared: 0.9763    Adjusted R-squared: 0.9723
   F-statistic: 246.93 on 1 and 6 DF    p-value: 0.000004
   Log-likelihood: 3.3957    AIC: -0.7915

**Interpretation.** After accounting for phylogenetic relatedness, brain size is a highly
significant predictor of body mass (t = 15.71, p < 0.0001). The coefficient estimate of
1.38 indicates that for each unit increase in log brain size, log body mass increases by
approximately 1.38 units. The model explains 97.6% of the variance (R-squared = 0.9763),
with an F-statistic of 246.93.

The negative intercept (-1.34) means that at the reference level of brain size (log brain
size = 0), the predicted log body mass is negative. The significance codes (``***``)
indicate p < 0.001.

|

Step 2: Run a multiple regression
***********************************

Test whether adding longevity as a second predictor improves the model:

.. code-block:: shell

   phykit phylogenetic_regression \
       -t tree_simple.tre \
       -d tree_simple_multi_traits.tsv \
       -y body_mass \
       -x brain_size longevity

.. code-block:: text

   Phylogenetic Generalized Least Squares (PGLS)

   Formula: body_mass ~ brain_size + longevity

   Coefficients:
                           Estimate   Std.Error     t-value     p-value
   (Intercept)              -1.2194      0.3985     -3.0601    0.028099    *
   brain_size                1.4466      0.2205      6.5608    0.001233    **
   longevity                -0.0879      0.2545     -0.3455    0.743755
   ---
   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1

   Residual standard error: 0.0377 on 5 degrees of freedom
   Multiple R-squared: 0.9768    Adjusted R-squared: 0.9676
   F-statistic: 105.40 on 2 and 5 DF    p-value: 0.000082
   Log-likelihood: 3.4901    AIC: 1.0197

**Interpretation.** Longevity is not a significant predictor of body mass (t = -0.35,
p = 0.74) after accounting for brain size and phylogenetic relatedness. The adjusted
R-squared decreases slightly (0.9676 vs 0.9723), and the AIC increases (1.02 vs -0.79),
confirming that adding longevity does not improve the model. Brain size remains the
dominant predictor.

|

Step 3: Export results as JSON
******************************

For downstream scripting, results can be exported as JSON:

.. code-block:: shell

   phykit phylogenetic_regression \
       -t tree_simple.tre \
       -d tree_simple_multi_traits.tsv \
       -y body_mass \
       -x brain_size \
       --json

The JSON output includes all regression statistics: coefficients, standard errors, t-values,
p-values, R-squared, adjusted R-squared, F-statistic, log-likelihood, AIC, fitted values,
and residuals for each taxon.

|

Summary
*******

In this tutorial, we used PGLS to test the relationship between body mass and brain size
while accounting for phylogenetic non-independence. The key steps were: (1) fitting a simple
regression with one predictor, (2) extending to multiple predictors to evaluate model
improvement, and (3) exporting results as JSON. PGLS is essential whenever comparing traits
across species because ignoring phylogenetic structure inflates degrees of freedom and can
produce spurious correlations.

For methodological details, see
`Freckleton et al. (2002) <https://doi.org/10.1086/343873>`_ and
`Symonds and Blomberg (2014) <https://doi.org/10.1007/978-3-662-43550-2_5>`_.
The R equivalent is ``caper::pgls()`` or ``nlme::gls()`` with ``ape::corBrownian()``.

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
