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
    kwargs = dict(tree="/some/path/to/file.tre", alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)

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
def alignment_simple(mocker):
    return AlignIO.read(open(
        f"{here.parent}/sample_files/simple.fa"), "fasta"
    )