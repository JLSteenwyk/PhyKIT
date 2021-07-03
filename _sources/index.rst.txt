.. image:: _static/img/logo.png
   :width: 55%
   :align: center
   :target: https://jlsteenwyk.com/PhyKIT

^^^^^


PhyKIT, a toolkit for the UNIX shell environment with numerous functions that process multiple
sequence alignments and phylogenies for broad applications

If you found PhyKIT useful, please cite *PhyKIT: a broadly applicable UNIX shell toolkit for processing
and analyzing phylogenomic data*. Bioinformatics. doi: |doiLink|_.

.. _doiLink: https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btab096/6131675?redirectedFrom=fulltext
.. |doiLink| replace:: 10.1093/bioinformatics/btab096

Quick Start
-----------
**1) Installation**

To install using *pip*, we strongly recommend building a virtual environment to avoid 
software dependency issues. To do so, execute the following commands:

.. code-block:: shell

	# create virtual environment
	python -m venv .venv
	# activate virtual environment
	source .venv/bin/activate
	# install phykit
	pip install phykit

**Note, the virtual environment must be activated to use phykit.**

After using PhyKIT, you may wish to deactivate your virtual environment and can do so using the following command:

.. code-block:: shell

	# deactivate virtual environment
	deactivate

|

Similarly, to install from source, we strongly recommend using a virtual environment. To do so, use the
following commands:

.. code-block:: shell

	# download
	git clone https://github.com/JLSteenwyk/PhyKIT.git
	cd PhyKIT/
	# create virtual environment
	python -m venv .venv
	# activate virtual environment
	source .venv/bin/activate
	# install
	make install

To deactivate your virtual environment, use the following command:

.. code-block:: shell

	# deactivate virtual environment
	deactivate

**Note, the virtual environment must be activated to use phykit.**

|

To install via anaconda, execute the following command:

.. code-block:: shell

	conda install -c jlsteenwyk phykit

Visit here for more information:
https://anaconda.org/JLSteenwyk/phykit

**2) Usage**

Get the help message from PhyKIT:

.. code-block:: shell

	phykit -h

|

^^^^


.. toctree::
	:maxdepth: 4

	about/index
	usage/index
	tutorials/index
	change_log/index
	other_software/index
	frequently_asked_questions/index

^^^^
