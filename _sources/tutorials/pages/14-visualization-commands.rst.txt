.. _tutorial-14:
.. _tutorial-visualization-commands:

Tutorial 14: Visualization commands
===================================

Objectives
----------

- Complete the visualization commands workflow.
- Interpret the reported values and generated artifacts in their scientific context.
- Identify the canonical command references for each analysis step.

Prerequisites and working directory
-----------------------------------

Install the current PhyKIT release and create a dedicated working directory.
Download the data linked in this tutorial into that directory before running
the commands. All paths below are relative to this directory.

.. code-block:: shell

   mkdir phykit-tutorial-14
   cd phykit-tutorial-14

Related command references
--------------------------

- :doc:`Cont Map </reference/commands/cont_map>`
- :doc:`Density Map </reference/commands/density_map>`
- :doc:`Phenogram </reference/commands/phenogram>`
- :doc:`Cophylo </reference/commands/cophylo>`

Workflow
--------

PhyKIT provides several phylogenetic visualization commands analogous to
R's ``phytools`` plotting functions. All produce publication-quality plots
saved at 300 DPI.

.. centered::
   Download test data:
   :download:`Mammal phylogeny </data/tree_simple.tre>`;
   :download:`Continuous trait data </data/tree_simple_traits.tsv>`;
   :download:`Discrete trait data </data/tree_simple_discrete_traits.tsv>`

**contMap** — continuous trait mapping onto a phylogram:

.. code-block:: shell

   phykit cont_map -t tree_simple.tre -d tree_simple_traits.tsv -o contmap.png

**densityMap** — posterior state probabilities from stochastic character mapping:

.. code-block:: shell

   phykit density_map -t tree_simple.tre -d tree_simple_discrete_traits.tsv \
       -c diet -o densitymap.png -n 100 --seed 42

**phenogram (traitgram)** — trait values vs. distance from root:

.. code-block:: shell

   phykit phenogram -t tree_simple.tre -d tree_simple_traits.tsv -o phenogram.png

**cophylo (tanglegram)** — two trees facing each other with connecting lines:

.. code-block:: shell

   phykit cophylo -t tree_simple.tre -t2 tree_simple.tre -o tanglegram.png

All visualization commands support ``--json`` output. The ``density_map`` command
also supports ``--seed`` for reproducibility and ``-n`` to set the number of
stochastic mapping simulations. The ``cophylo`` command supports ``-m/--mapping``
for a tab-delimited file matching taxa between trees with different tip names.

For methodological details, see
`Revell (2012) <https://doi.org/10.1111/j.2041-210X.2011.00169.x>`_.

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
