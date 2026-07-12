.. _getting-started:

Installation and first analysis
===============================

Requirements
------------

PhyKIT supports Python 3.10, 3.11, 3.12, and 3.13. Use a virtual environment
to keep its dependencies separate from other projects.

Install from PyPI
-----------------

.. code-block:: shell

   python -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install phykit
   phykit version

On Windows PowerShell, activate the environment with
``.venv\Scripts\Activate.ps1``.

Install with conda
------------------

.. code-block:: shell

   conda create -n phykit -c conda-forge -c bioconda phykit
   conda activate phykit
   phykit version

Run a first analysis
--------------------

Create a working directory and download a versioned example alignment.

.. code-block:: shell

   mkdir phykit-first-analysis
   cd phykit-first-analysis
   curl -LO https://jlsteenwyk.github.io/PhyKIT/data/Steenwyk_etal_mBio_2019_EOG091N44MS.aln.fa
   phykit alignment_length Steenwyk_etal_mBio_2019_EOG091N44MS.aln.fa

Expected output:

.. code-block:: text

   624

PhyKIT read an aligned FASTA file and reported its 624 aligned sites. Request
structured output with ``--json``:

.. code-block:: shell

   phykit alignment_length Steenwyk_etal_mBio_2019_EOG091N44MS.aln.fa --json

Find the next command
---------------------

- Alignment quality and composition: :doc:`Alignment quality and statistics </reference/categories/alignment-quality-statistics>`
- Tree summaries: :doc:`Tree summary statistics </reference/categories/tree-summary-statistics>`
- Tree comparison: :doc:`Tree comparison and consensus </reference/categories/tree-comparison-consensus>`
- Trait evolution: :doc:`Trait evolution </reference/categories/trait-evolution>`
- Comparative methods: :doc:`Phylogenetic comparative methods </reference/categories/phylogenetic-comparative-methods>`
- Introgression and gene flow: :doc:`Introgression and gene flow </reference/categories/introgression-gene-flow>`

For any command, display the interface installed in the active environment::

   phykit <command> --help
