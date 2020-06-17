## Run alignment-based functions
run.alignment_length:
	python3 -m phykit-runner alignment_length ./tests/sample_files/simple.fa 

run.alignment_length_no_gaps:
	python3 -m phykit-runner alignment_length_no_gaps ./tests/sample_files/simple.fa

run.gc_content:
	python3 phykit-runner.py gc_content ./tests/sample_files/simple.fa

run.rename_fasta_entries:
	python3 phykit-runner.py rename_fasta_entries ./tests/sample_files/simple.fa -i tests/sample_files/simple_fasta_idmap.txt 

run.parsimony_informative_sites:
	python3 -m phykit-runner parsimony_informative_sites ./tests/sample_files/simple.fa

run.pairwise_identity:
	python3 phykit-runner.py pairwise_identity ./tests/sample_files/simple.fa

run.rcv:
	python3 -m phykit-runner rcv ./tests/sample_files/simple.fa

run.variable_sites:
	python3 -m phykit-runner variable_sites ./tests/sample_files/simple.fa


## Run tree-based functions
run.bipartition_support_stats:
	python3 -m phykit-runner bipartition_support_stats ./tests/sample_files/small_Aspergillus_tree.tre 

run.branch_length_multiplier:
	python3 -m phykit-runner.py branch_length_multiplier ./tests/sample_files/tree_simple.tre -f 2

run.covarying_evolutionary_rates:
	python3 phykit-runner.py covarying_evolutionary_rates ./tests/sample_files/tree_simple.tre ./tests/sample_files/tree_simple_1.tre -r ./tests/sample_files/tree_simple_2.tre  

run.dvmc:
	python3 -m phykit-runner dvmc -t ./tests/sample_files/tree_simple.tre -r ./tests/sample_files/tree_simple.outgroup.txt

run.internal_branch_stats:
	python3 -m phykit-runner internal_branch_stats ./tests/sample_files/tree_simple.tre

run.internode_labeler:
	python3 -m phykit-runner internode_labeler ./tests/sample_files/tree_simple.tre

run.lb_score:
	python3 -m phykit-runner lb_score ./tests/sample_files/tree_simple.tre

run.patristic_distances:
	python3 -m phykit-runner patristic_distances ./tests/sample_files/tree_simple.tre

run.print_tree:
	python3 -m phykit-runner print_tree ./tests/sample_files/tree_simple.tre

run.rename_tree_tips:
	python3 phykit-runner.py rename_tree_tips ./tests/sample_files/tree_simple.tre -i ./tests/sample_files/tree_simple_idmap.txt

run.rf_distance:
	python3 -m phykit-runner rf_distance ./tests/sample_files/tree_simple.tre ./tests/sample_files/tree_simple_other_topology.tre

run.saturation:
	python3 phykit-runner.py saturation -t ./tests/sample_files/12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit.treefile -a ./tests/sample_files/12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit

run.spurious_sequence:
	python3 -m phykit-runner spurious_sequence ./tests/sample_files/tree_simple.tre  

run.tip_labels:
	python3 -m phykit-runner tip_labels ./tests/sample_files/tree_simple.tre

run.total_tree_length:
	python3 -m phykit-runner total_tree_length ./tests/sample_files/tree_simple.tre

run.treeness:
	python3 -m phykit-runner treeness ./tests/sample_files/Yeasts_2832_eMRC_reference_renamed.tree

run.treeness_over_rcv:
	python3 -m phykit-runner treeness_over_rcv -t ./tests/sample_files/tree_simple.tre -a ./tests/sample_files/simple.fa 


## Helper commands
run.create_concatenation_matrix:
	python3 phykit-runner.py create_concatenation_matrix -a tests/sample_files/alignment_list_for_create_concat_matrix.txt -p test

run.thread_dna:
	python3 -m phykit-runner thread_dna -p ./tests/sample_files/EOG091N44MS.fa.mafft -n ./tests/sample_files/EOG091N44MS.fa


## Install, develop, and testing make commands
install:
	# install so phykit command is available in terminal
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
