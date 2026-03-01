import json
import os
import tempfile

import pytest
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.cont_map import ContMap


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            output="out.png",
        )
        svc = ContMap.__new__(ContMap)
        parsed = svc.process_args(args)
        assert parsed["tree_file_path"] == "t.tre"
        assert parsed["trait_data_path"] == "d.tsv"
        assert parsed["output_path"] == "out.png"
        assert parsed["json_output"] is False

    def test_overrides(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            output="out.png",
            json=True,
        )
        svc = ContMap.__new__(ContMap)
        parsed = svc.process_args(args)
        assert parsed["output_path"] == "out.png"
        assert parsed["json_output"] is True


class TestRun:
    def test_plot_created(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = Namespace(
                tree=TREE_SIMPLE,
                trait_data=TRAITS_FILE,
                output=tmppath,
                json=False,
            )
            svc = ContMap(args)
            svc.run()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_json_output(self, capsys):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = Namespace(
                tree=TREE_SIMPLE,
                trait_data=TRAITS_FILE,
                output=tmppath,
                json=True,
            )
            svc = ContMap(args)
            svc.run()
            captured = capsys.readouterr()
            # The last line of output should be valid JSON
            lines = captured.out.strip().split("\n")
            payload = json.loads(lines[-1])
            assert "n_tips" in payload
            assert payload["n_tips"] == 8
            assert "sigma2" in payload
            assert payload["sigma2"] > 0
            assert "tip_values" in payload
            assert "ancestral_estimates" in payload
            assert "plot_output" in payload
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)
