.. image:: _static/img/logo.png
   :width: 55%
   :align: center
   :target: https://jlsteenwyk.com/PhyKIT

^^^^^


PhyKIT, a toolkit for the UNIX shell environment with 30 functions that process multiple
sequence alignments and phylogenies for broad applications

If you found clipkit useful, please cite *MANUSCRIPT TITLE*. bioRxiv. doi: |doiLink|_.

.. _doiLink: BIORXIV LINK
.. |doiLink| replace:: TEXT

Quick Start
-----------
**1) Installation**

To install, use the following commands:

.. code-block:: shell

	pip install phykit

|

To install from source, use the following commands:

.. code-block:: shell

	git clone https://github.com/JLSteenwyk/PhyKIT.git
	cd PhyKIT/
	make install

If you run into permission errors when executing *make install*, create a 
virtual environemnt for your installation:

.. code-block:: shell

	git clone https://github.com/JLSteenwyk/PhyKIT.git
	cd PhyKIT/
	python -m venv .venv
	source .venv/bin/activate
	make install

Note, the virtual environment must be activated to use phykit.

|

**2) Usage**

Get the help message from PhyKIT:

.. code-block:: shell

	phykit -h

|


^^^^

.. toctree::
	:maxdepth: 4

	about/index
	analytic_functions/index
	helper_functions/index
	frequently_asked_questions/index

^^^^

