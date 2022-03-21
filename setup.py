from os import path
from setuptools import setup, find_packages

from phykit.version import __version__

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md")) as f:
    long_description = f.read()

CLASSIFIERS = [
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Topic :: Scientific/Engineering',
]

REQUIRES = [
    "biopython>=1.79",
    "numpy>=1.18.2",
    "scipy>=1.4.1",
    #"ete3>=3.1.2",
    "cython"
]

setup(
    name="phykit",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Jacob L. Steenwyk",
    author_email="jlsteenwyk@gmail.com",
    url="https://github.com/jlsteenwyk/phykit",
    packages=find_packages(),
    classifiers=CLASSIFIERS,
    entry_points={
        "console_scripts": [
            "phykit = phykit.phykit:main",
            "pk_alignment_length = phykit.phykit:alignment_length", # Alignment-based functions
            "pk_aln_len = phykit.phykit:alignment_length",
            "pk_al = phykit.phykit:alignment_length",
            "pk_alignment_length_no_gaps = phykit.phykit:alignment_length_no_gaps",
            "pk_aln_len_no_gaps = phykit.phykit:alignment_length_no_gaps",
            "pk_alng = phykit.phykit:alignment_length_no_gaps",
            "pk_column_score = phykit.phykit:column_score",
            "pk_cs = phykit.phykit:column_score",
            "pk_faidx = phykit.phykit:faidx",
            "pk_get_entry = phykit.phykit:faidx",
            "pk_ge = phykit.phykit:faidx",
            "pk_gc_content = phykit.phykit:gc_content",
            "pk_gc = phykit.phykit:gc_content",
            "pk_pairwise_identity = phykit.phykit:pairwise_identity",
            "pk_pairwise_id = phykit.phykit:pairwise_identity",
            "pk_pi = phykit.phykit:pairwise_identity",
            "pk_parsimony_informative_sites = phykit.phykit:parsimony_informative_sites",
            "pk_pis = phykit.phykit:parsimony_informative_sites",
            "pk_relative_composition_variability = phykit.phykit:rcv",
            "pk_rel_comp_var = phykit.phykit:rcv",
            "pk_rcv = phykit.phykit:rcv",
            "pk_rename_fasta_entries = phykit.phykit:rename_fasta_entries",
            "pk_rename_fasta = phykit.phykit:rename_fasta_entries",
            "pk_sum_of_pairs_score = phykit.phykit:sum_of_pairs_score",
            "pk_sops = phykit.phykit:sum_of_pairs_score",
            "pk_sop = phykit.phykit:sum_of_pairs_score",
            "pk_variable_sites = phykit.phykit:variable_sites",
            "pk_vs = phykit.phykit:variable_sites",
            "pk_bipartition_support_stats = phykit.phykit:bipartition_support_stats", # Tree-based functions
            "pk_bss = phykit.phykit:bipartition_support_stats",
            "pk_branch_length_multiplier = phykit.phykit:branch_length_multiplier",
            "pk_blm = phykit.phykit:branch_length_multiplier",
            "pk_collapse_branches = phykit.phykit:collapse_branches",
            "pk_collapse = phykit.phykit:collapse_branches",
            "pk_cb = phykit.phykit:collapse_branches",
            "pk_covarying_evolutionary_rates = phykit.phykit:covarying_evolutionary_rates",
            "pk_cover = phykit.phykit:covarying_evolutionary_rates",
            "pk_dvmc = phykit.phykit:dvmc",
            "pk_degree_of_violation_of_a_molecular_clock = phykit.phykit:dvmc",
            "pk_evo_rate = phykit.phykit:evolutionary_rate",
            "pk_evolutionary_rate = phykit.phykit:evolutionary_rate",
            "pk_hidden_paralogy_check = phykit.phykit:hidden_paralogy_check",
            "pk_clan_check = phykit.phykit:hidden_paralogy_check",
            "pk_internal_branch_stats = phykit.phykit:internal_branch_stats",
            "pk_ibs = phykit.phykit:internal_branch_stats",
            "pk_internode_labeler = phykit.phykit:internode_labeler",
            "pk_il = phykit.phykit:internode_labeler",
            "pk_last_common_ancestor_subtree = phykit.phykit:last_common_ancestor_subtree",
            "pk_lca_subtree = phykit.phykit:last_common_ancestor_subtree",
            "pk_lb_score = phykit.phykit:lb_score",
            "pk_long_branch_score = phykit.phykit:lb_score",
            "pk_lbs = phykit.phykit:lb_score",
            "pk_monophyly_check = phykit.phykit:monophyly_check",
            "pk_is_monophyletic = phykit.phykit:monophyly_check",
            "pk_nearest_neighbor_interchange = phykit.phykit:nearest_neighbor_interchange",
            "pk_nni = phykit.phykit:nearest_neighbor_interchange",
            "pk_patristic_distances = phykit.phykit:patristic_distances",
            "pk_pd = phykit.phykit:patristic_distances",
            "pk_polytomy_test = phykit.phykit:polytomy_test",
            "pk_polyt_test = phykit.phykit:polytomy_test",
            "pk_ptt = phykit.phykit:polytomy_test",
            "pk_polyt = phykit.phykit:polytomy_test",
            "pk_print_tree = phykit.phykit:print_tree",
            "pk_print = phykit.phykit:print_tree",
            "pk_pt = phykit.phykit:print_tree",
            "pk_prune_tree = phykit.phykit:prune_tree",
            "pk_prune = phykit.phykit:prune_tree",
            "pk_rename_tree_tips = phykit.phykit:rename_tree_tips",
            "pk_rename_tree = phykit.phykit:rename_tree_tips",
            "pk_rename_tips = phykit.phykit:rename_tree_tips",
            "pk_rf_distance = phykit.phykit:rf_distance",
            "pk_robinson_foulds_distance = phykit.phykit:rf_distance",
            "pk_rf_dist = phykit.phykit:rf_distance",
            "pk_rf = phykit.phykit:rf_distance",
            "pk_root_tree = phykit.phykit:root_tree",
            "pk_root = phykit.phykit:root_tree",
            "pk_rt = phykit.phykit:root_tree",
            "pk_spurious_sequence = phykit.phykit:spurious_sequence",
            "pk_spurious_seq = phykit.phykit:spurious_sequence",
            "pk_ss = phykit.phykit:spurious_sequence",
            "pk_terminal_branch_stats = phykit.phykit:terminal_branch_stats",
            "pk_tbs = phykit.phykit:terminal_branch_stats",
            "pk_tip_labels = phykit.phykit:tip_labels",
            "pk_labels = phykit.phykit:tip_labels",
            "pk_tree_labels = phykit.phykit:tip_labels",
            "pk_tl = phykit.phykit:tip_labels",
            "pk_tip_to_tip_distance = phykit.phykit:tip_to_tip_distance",
            "pk_t2t_dist = phykit.phykit:tip_to_tip_distance",
            "pk_t2t = phykit.phykit:tip_to_tip_distance",
            "pk_tip_to_tip_node_distance = phykit.phykit:tip_to_tip_node_distance",
            "pk_t2t_node_dist = phykit.phykit:tip_to_tip_node_distance",
            "pk_t2t_nd = phykit.phykit:tip_to_tip_node_distance",
            "pk_total_tree_length = phykit.phykit:total_tree_length",
            "pk_tree_len = phykit.phykit:total_tree_length",
            "pk_treeness = phykit.phykit:treeness",
            "pk_tness = phykit.phykit:treeness",
            "pk_saturation = phykit.phykit:saturation", # Alignment- and tree-based functions
            "pk_sat = phykit.phykit:saturation",
            "pk_treeness_over_rcv = phykit.phykit:treeness_over_rcv",
            "pk_toverr = phykit.phykit:treeness_over_rcv",
            "pk_tor = phykit.phykit:treeness_over_rcv",
            "pk_create_concatenation_matrix = phykit.phykit:create_concatenation_matrix", # Helper functions
            "pk_create_concat = phykit.phykit:create_concatenation_matrix",
            "pk_cc = phykit.phykit:create_concatenation_matrix",
            "pk_thread_dna = phykit.phykit:thread_dna",
            "pk_pal2nal = phykit.phykit:thread_dna",
            "pk_p2n = phykit.phykit:thread_dna",
        ]
    },
    version=__version__,
    include_package_data=True,
    install_requires=REQUIRES,
)

## push new version to pypi
# rm -rf dist
# python3 setup.py sdist bdist_wheel --universal
# twine upload dist/* -r pypi
# then push to anaconda
# 