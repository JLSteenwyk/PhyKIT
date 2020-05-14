run.treeness:
	python3 -m phykit-runner treeness ./tests/sample_files/Yeasts_2832_eMRC_reference_renamed.tree

run.internode_labeler:
	python3 -m phykit-runner internode_labeler ./tests/sample_files/tree_simple.tre

install:
	# install so clipkit command is available in terminal
	python setup.py install

develop:
	# https://setuptools.readthedocs.io/en/latest/setuptools.html#development-mode
	python setup.py develop

test: test.unit test.integration

test.unit:
	python -m pytest -m "not integration"

test.integration:
	rm -rf output/
	mkdir output/
	python -m pytest --basetemp=output -m "integration"
