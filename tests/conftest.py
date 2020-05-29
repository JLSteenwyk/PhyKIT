import pytest
from argparse import Namespace
from pathlib import Path

from Bio import AlignIO
from Bio import Phylo


here = Path(__file__)


def pytest_configure(config):
    config.addinivalue_line("markers", "integration: mark as integration test")


@pytest.fixture
def args():
    kwargs = dict(
        alignment="/some/path/to/file.fa", 
        tree="/some/path/to/file.tre", 
        root="/home/path/to/file.txt"
        )
    return Namespace(**kwargs)

# tree fixtures
@pytest.fixture
def tree_zero_branch_length(mocker):
    return Phylo.read(
        f"{here.parent}/sample_files/tree_zero_branch_length.tre",
        "newick",
    )

@pytest.fixture
def tree_simple(mocker):
    return Phylo.read(
        f"{here.parent}/sample_files/tree_simple.tre", "newick",
    )

@pytest.fixture
def tree_simple_outgroup(mocker):
    return [
        line.rstrip('\n') for line in open(f"{here.parent}/sample_files/tree_simple.outgroup.txt")
    ]

@pytest.fixture
def small_aspergillus_tree(mocker):
    return Phylo.read(
        f"{here.parent}/sample_files/small_Aspergillus_tree.tre", "newick",
    )

# alignment fixtures
@pytest.fixture
def alignment_simple(mocker):
    return AlignIO.read(open(
        f"{here.parent}/sample_files/simple.fa"), "fasta"
    )