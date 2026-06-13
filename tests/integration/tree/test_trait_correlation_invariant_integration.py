import os
import sys
from pathlib import Path

import pytest
from mock import patch

from phykit.phykit import Phykit

here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
MULTI_TRAITS_INVARIANT = str(SAMPLE_FILES / "tree_simple_multi_traits_invariant.tsv")


@pytest.mark.integration
class TestTraitCorrelationInvariant:
    """A multi-trait file containing an invariant (zero-variance) column must
    not crash clustering; the invariant column is dropped with a warning."""

    def test_invariant_column_dropped(self, tmp_path, capsys):
        out = str(tmp_path / "corr.png")
        testargs = [
            "phykit",
            "trait_correlation",
            "-t", TREE_SIMPLE,
            "-d", MULTI_TRAITS_INVARIANT,
            "-o", out,
            "--cluster",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        # plot is produced despite the invariant column
        assert os.path.exists(out)
        captured = capsys.readouterr()
        # invariant column reported on stderr and excluded from the matrix
        assert "constant_trait" in captured.err
        assert "dropping 1 invariant" in captured.err
        assert "constant_trait" not in captured.out
