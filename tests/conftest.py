import pytest
from argparse import Namespace
from pathlib import Path

from Bio import AlignIO
from Bio import Phylo


here = Path(__file__)


def pytest_configure(config):
    config.addinivalue_line("markers", "integration: mark as integration test")

# alignment fixtures
@pytest.fixture
def alignment_simple(mocker):
    return AlignIO.read(open(f"{here.parent}/sample_files/simple.fa"), "fasta")

@pytest.fixture
def alignment_complex(mocker):
    return AlignIO.read(open(f"{here.parent}/sample_files/12_YPR191W_Anc_7.548_codon_aln.fasta.clipkit"), "fasta")

@pytest.fixture
def alignments(mocker):
    alignment_list = f"{here.parent}/sample_files/alignment_list_for_create_concat_matrix.txt"
    alignments = [line.rstrip('\n') for line in open(alignment_list)]
    return alignments


# tree fixtures
@pytest.fixture
def tree_zero_branch_length(mocker):
    return Phylo.read(
        f"{here.parent}/sample_files/tree_zero_branch_length.tre", "newick",
    )


@pytest.fixture
def tree_simple(mocker):
    return Phylo.read(f"{here.parent}/sample_files/tree_simple.tre", "newick",)


@pytest.fixture
def tree_simple_other(mocker):
    return Phylo.read(
        f"{here.parent}/sample_files/tree_simple_other_topology.tre", "newick",
    )


@pytest.fixture
def tree_simple_outgroup(mocker):
    return [
        line.rstrip("\n")
        for line in open(f"{here.parent}/sample_files/tree_simple.outgroup.txt")
    ]


@pytest.fixture
def small_aspergillus_tree(mocker):
    return Phylo.read(
        f"{here.parent}/sample_files/small_Aspergillus_tree.tre", "newick",
    )

