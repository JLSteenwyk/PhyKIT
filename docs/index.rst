.. image:: _static/img/logo.png
   :alt: PhyKIT logo
   :width: 55%
   :align: center
   :target: https://jlsteenwyk.com/PhyKIT

PhyKIT
======

PhyKIT is a command-line toolkit for processing and analyzing multiple
sequence alignments, phylogenetic trees, and comparative trait data.

Start here
----------

- :doc:`Installation and first analysis <getting_started/index>`
- :doc:`Choose a command by analytical task <reference/index>`
- :doc:`Follow a complete workflow <tutorials/index>`
- :doc:`Prepare input files <formats/index>`
- :doc:`Troubleshoot an error <troubleshooting/index>`

First analysis
--------------

Install PhyKIT with Python 3.10 or newer, download the example alignment,
and calculate its length::

   python -m pip install phykit
   curl -LO https://jlsteenwyk.github.io/PhyKIT/data/Steenwyk_etal_mBio_2019_EOG091N44MS.aln.fa
   phykit alignment_length Steenwyk_etal_mBio_2019_EOG091N44MS.aln.fa

Expected output::

   624

The value is the number of aligned sites. Continue with
:doc:`Tutorial 1 <tutorials/pages/01-summarizing-information-content>` to
measure additional information-content and tree-quality statistics.

Citation
--------

If you use PhyKIT, cite *PhyKIT: a broadly applicable UNIX shell toolkit for
processing and analyzing phylogenomic data*, Bioinformatics,
`doi:10.1093/bioinformatics/btab096 <https://doi.org/10.1093/bioinformatics/btab096>`_.

.. toctree::
   :maxdepth: 4

   getting_started/index
   reference/index
   tutorials/index
   formats/index
   glossary/index
   troubleshooting/index
   frequently_asked_questions/index
   usage/index
   about/index
   change_log/index
   other_software/index
