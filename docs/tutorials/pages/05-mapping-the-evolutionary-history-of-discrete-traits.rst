.. _tutorial-05:
.. _tutorial-mapping-the-evolutionary-history-of-discrete-traits:

Tutorial 5: Mapping the evolutionary history of discrete traits
===============================================================

Objectives
----------

- Complete the mapping the evolutionary history of discrete traits workflow.
- Interpret the reported values and generated artifacts in their scientific context.
- Identify the canonical command references for each analysis step.

Prerequisites and working directory
-----------------------------------

Install the current PhyKIT release and create a dedicated working directory.
Download the data linked in this tutorial into that directory before running
the commands. All paths below are relative to this directory.

.. code-block:: shell

   mkdir phykit-tutorial-05
   cd phykit-tutorial-05

Related command references
--------------------------

- :doc:`Stochastic Character Map </reference/commands/stochastic_character_map>`

Workflow
--------

A common question in comparative biology is: how did a discrete trait, such as diet, habitat,
or reproductive strategy, evolve across a phylogeny? Simply labeling tips on a tree does not
tell us *when* or *how often* transitions between states occurred. Stochastic character mapping
(`Huelsenbeck et al. 2003 <https://doi.org/10.1080/10635150390192780>`_;
`Bollback 2006 <https://doi.org/10.1186/1471-2148-6-88>`_)
addresses this by simulating plausible evolutionary histories of a discrete trait along each
branch, conditioned on the observed tip states and a fitted substitution model. This approach
is widely used in macroevolutionary studies to quantify the tempo and mode of trait evolution
(`O'Meara 2012 <https://doi.org/10.1146/annurev-ecolsys-110411-160331>`_).

**Hypothetical study question.** Suppose we are studying a group of eight mammal species and
want to understand how diet (carnivore, herbivore, or omnivore) has evolved across the
phylogeny. Specifically, we want to know: (1) Is the transition rate between all dietary states
equal, or are some transitions more frequent than others? (2) How much evolutionary time has
been spent in each dietary state? (3) How many transitions between states occurred on average?

PhyKIT's ``stochastic_character_map`` command (alias: ``simmap``) lets us answer these
questions directly from the command line.

In this tutorial, we will use test data included with PhyKIT: an eight-taxon mammal phylogeny
and a tab-delimited file assigning each species to a dietary category. |br|

.. centered::
   Download test data:
   :download:`Mammal phylogeny </data/tree_simple.tre>`;
   :download:`Discrete trait data </data/tree_simple_discrete_traits.tsv>`

|

Step 0: Prepare data
********************

Two input files are needed: a phylogenetic tree in Newick format and a tab-delimited trait
file with a header row. The trait file looks like this:

.. code-block:: text

   taxon	diet
   raccoon	carnivore
   bear	carnivore
   sea_lion	carnivore
   seal	herbivore
   monkey	herbivore
   cat	omnivore
   weasel	omnivore
   dog	omnivore

The first column is the taxon name (must match the tree tip labels) and subsequent columns
contain discrete trait values. You specify which column to use with the ``-c/--trait`` flag.
Lines starting with ``#`` are treated as comments and blank lines are ignored.

Download both files above into the tutorial working directory.

|

Step 1: Fit a substitution model and run stochastic mapping
***********************************************************

First, run stochastic character mapping with the default equal-rates (ER) model and 100
simulations:

.. code-block:: shell

   phykit simmap \
       -t tree_simple.tre \
       -d tree_simple_discrete_traits.tsv \
       -c diet \
       -n 100 \
       --seed 42

This produces the following output:

.. code-block:: text

   Stochastic Character Mapping (SIMMAP)

   Model: ER
   Number of simulations: 100

   Fitted Q matrix:
                    carnivore   herbivore    omnivore
   carnivore          -0.1138      0.0569      0.0569
   herbivore           0.0569     -0.1138      0.0569
   omnivore            0.0569      0.0569     -0.1138

   Log-likelihood: -8.7874

   Mean dwelling times:
     carnivore         99.08 (35.7%)
     herbivore         89.18 (32.2%)
     omnivore          89.02 (32.1%)

   Mean transitions:
     carnivore -> herbivore:  5.51
     carnivore -> omnivore:  5.69
     herbivore -> carnivore:  5.48
     herbivore -> omnivore:  5.37
     omnivore -> carnivore:  4.94
     omnivore -> herbivore:  4.89
     Total:                      31.88

**Interpreting the output.** The fitted Q matrix shows the instantaneous rates of transition
between states. Under the ER model, all off-diagonal rates are equal (0.0569 per unit branch
length). The log-likelihood of -8.79 is the maximized log-likelihood of the data given the
model.

The mean dwelling times tell us how much total evolutionary time (summed across all branches)
was spent in each state, averaged over all 100 simulated histories. Here, the three dietary
states have roughly equal dwelling times, reflecting the ER model's symmetry and the
distribution of states across the tree.

The mean transition counts show how many times each type of state change occurred on average.
These are averaged across 100 stochastic maps.

|

Step 2: Compare substitution models
************************************

Is the equal-rates model adequate, or do different transitions have different rates? We can
compare the ER model with the all-rates-different (ARD) model:

.. code-block:: shell

   phykit simmap \
       -t tree_simple.tre \
       -d tree_simple_discrete_traits.tsv \
       -c diet \
       -m ARD \
       -n 100 \
       --seed 42

You can also try the symmetric (SYM) model, which assumes forward and reverse rates between
any pair of states are equal but allows different pairs to differ:

.. code-block:: shell

   phykit simmap \
       -t tree_simple.tre \
       -d tree_simple_discrete_traits.tsv \
       -c diet \
       -m SYM \
       -n 100 \
       --seed 42

Compare the log-likelihoods across models to assess fit. Because ARD has more parameters than
SYM, which has more than ER, a likelihood ratio test or AIC comparison can be used to
determine whether the additional parameters are justified.

|

Step 3: Generate a stochastic character map plot
************************************************

To visualize one of the simulated character histories on the phylogeny, use the ``--plot``
flag:

.. code-block:: shell

   phykit simmap \
       -t tree_simple.tre \
       -d tree_simple_discrete_traits.tsv \
       -c diet \
       -n 100 \
       --seed 42 \
       --plot simmap_diet.png

This generates a horizontal phylogram with branches colored by the mapped character state.
Each branch segment is colored according to the state occupied during that interval,
reflecting one of the simulated character histories. A legend maps colors to states.

|

Step 4: Export results as JSON
******************************

For downstream analysis or scripting, results can be exported as JSON:

.. code-block:: shell

   phykit simmap \
       -t tree_simple.tre \
       -d tree_simple_discrete_traits.tsv \
       -c diet \
       -n 100 \
       --seed 42 \
       --json

The JSON output includes the fitted Q matrix, log-likelihood, state list, mean dwelling times
and proportions, and mean transition counts. This can be parsed with standard JSON tools
(e.g., ``jq``, Python's ``json`` module) for further analysis.

|

Step 5: Ensure reproducibility
******************************

The ``--seed`` flag ensures that the stochastic simulations are reproducible. Running the same
command with the same seed will produce identical results. This is important for
reproducibility in publications. Omitting ``--seed`` produces different results each time due
to the stochastic nature of the simulations.

|

Summary
*******

In this tutorial, we used stochastic character mapping to reconstruct the evolutionary history
of diet across a mammal phylogeny. The key steps were: (1) fitting a continuous-time Markov
chain rate matrix to estimate transition rates, (2) comparing equal-rates and
all-rates-different models to assess whether transition rates vary, (3) simulating character
histories to estimate dwelling times and transition counts, and (4) plotting a stochastic
character map for visualization.

For methodological details, see
`Huelsenbeck et al. (2003) <https://doi.org/10.1080/10635150390192780>`_ and
`Bollback (2006) <https://doi.org/10.1186/1471-2148-6-88>`_.
The R equivalent is ``phytools::make.simmap()``
(`Revell 2012 <https://doi.org/10.1111/j.2041-210X.2011.00169.x>`_).

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
