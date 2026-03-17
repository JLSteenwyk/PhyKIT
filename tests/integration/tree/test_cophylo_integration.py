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
TREE1 = str(SAMPLE_FILES / "tree_simple.tre")
TREE2 = str(SAMPLE_FILES / "tree_simple_other_topology.tre")


@pytest.mark.integration
class TestCophyloIntegration:
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
                "cophylo",
                "-t", TREE1,
                "-t2", TREE2,
                "-o", tmppath,
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()

            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_alias_tanglegram(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            testargs = [
                "phykit",
                "tanglegram",
                "-t", TREE1,
                "-t2", TREE2,
                "-o", tmppath,
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()

            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_alias_tangle(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            testargs = [
                "phykit",
                "tangle",
                "-t", TREE1,
                "-t2", TREE2,
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
                "cophylo",
                "-t", TREE1,
                "-t2", TREE2,
                "-o", tmppath,
                "--json",
            ]
            with patch.object(sys, "argv", testargs):
                Phykit()

            payload = json.loads(mocked_print.call_args.args[0])
            assert payload["n_tips_tree1"] == 8
            assert payload["n_tips_tree2"] == 8
            assert payload["n_matched"] == 8
            assert "matched_taxa" in payload
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_cophylo_circular(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            testargs = [
                "phykit",
                "cophylo",
                "-t", TREE1,
                "-t2", TREE2,
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
