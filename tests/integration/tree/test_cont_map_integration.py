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


@pytest.mark.integration
class TestContMapIntegration:
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
                "cont_map",
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

    def test_alias_contmap(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            testargs = [
                "phykit",
                "contmap",
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

    def test_alias_cmap(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            testargs = [
                "phykit",
                "cmap",
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
                "cont_map",
                "-t", TREE_SIMPLE,
                "-d", TRAITS_FILE,
                "-o", tmppath,
                "--json",
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()

            # The last print call should be JSON
            payload = json.loads(mocked_print.call_args.args[0])
            assert payload["n_tips"] == 8
            assert "sigma2" in payload
            assert "tip_values" in payload
            assert "ancestral_estimates" in payload
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)
