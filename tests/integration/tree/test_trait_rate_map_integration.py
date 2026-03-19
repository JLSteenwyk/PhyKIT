import json
import os
import sys
import tempfile
from pathlib import Path

import pytest
from mock import patch

from phykit.phykit import Phykit


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


@pytest.mark.integration
class TestTraitRateMapIntegration:
    def test_basic_invocation(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            testargs = [
                "phykit",
                "trait_rate_map",
                "-t", TREE_SIMPLE,
                "-d", TRAITS_FILE,
                "-o", tmppath,
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()

            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_alias_rate_map(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            testargs = [
                "phykit",
                "rate_map",
                "-t", TREE_SIMPLE,
                "-d", TRAITS_FILE,
                "-o", tmppath,
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()

            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_alias_branch_rates(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            testargs = [
                "phykit",
                "branch_rates",
                "-t", TREE_SIMPLE,
                "-d", TRAITS_FILE,
                "-o", tmppath,
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()

            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_circular(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            testargs = [
                "phykit",
                "trait_rate_map",
                "-t", TREE_SIMPLE,
                "-d", TRAITS_FILE,
                "-o", tmppath,
                "--circular",
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()

            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_cladogram(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            testargs = [
                "phykit",
                "trait_rate_map",
                "-t", TREE_SIMPLE,
                "-d", TRAITS_FILE,
                "-o", tmppath,
                "--cladogram",
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()

            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            testargs = [
                "phykit",
                "trait_rate_map",
                "-t", TREE_SIMPLE,
                "-d", TRAITS_FILE,
                "-o", tmppath,
                "--json",
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()

            # The last print call should be JSON
            payload = json.loads(mocked_print.call_args.args[0])
            assert payload["n_taxa"] == 8
            assert "mean_rate" in payload
            assert "min_rate" in payload
            assert "max_rate" in payload
            assert "branches" in payload
            assert len(payload["branches"]) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)
