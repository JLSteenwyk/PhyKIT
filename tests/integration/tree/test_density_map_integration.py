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
DISCRETE_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_discrete_traits.tsv")


@pytest.mark.integration
class TestDensityMapIntegration:
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
                "density_map",
                "-t", TREE_SIMPLE,
                "-d", DISCRETE_TRAITS_FILE,
                "-c", "diet",
                "-n", "10",
                "--seed", "42",
                "-o", tmppath,
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()

            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_alias_densitymap(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            testargs = [
                "phykit",
                "densitymap",
                "-t", TREE_SIMPLE,
                "-d", DISCRETE_TRAITS_FILE,
                "-c", "diet",
                "-n", "10",
                "--seed", "42",
                "-o", tmppath,
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()

            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_alias_dmap(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            testargs = [
                "phykit",
                "dmap",
                "-t", TREE_SIMPLE,
                "-d", DISCRETE_TRAITS_FILE,
                "-c", "diet",
                "-n", "10",
                "--seed", "42",
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
                "density_map",
                "-t", TREE_SIMPLE,
                "-d", DISCRETE_TRAITS_FILE,
                "-c", "diet",
                "-n", "10",
                "--seed", "42",
                "-o", tmppath,
                "--json",
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()

            payload = json.loads(mocked_print.call_args.args[0])
            assert payload["n_tips"] == 8
            assert "states" in payload
            assert payload["n_sim"] == 10
            assert "plot_output" in payload
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)
