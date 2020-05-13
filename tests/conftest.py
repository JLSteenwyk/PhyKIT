import pytest
from pathlib import Path

from Bio import Phylo

here = Path(__file__)


def pytest_configure(config):
    config.addinivalue_line("markers", "integration: mark as integration test")
