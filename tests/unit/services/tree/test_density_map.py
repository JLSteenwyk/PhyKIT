import json
import os
import tempfile

import pytest
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.density_map import DensityMap


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
DISCRETE_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_discrete_traits.tsv")


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            trait="diet",
            output="out.png",
        )
        svc = DensityMap.__new__(DensityMap)
        parsed = svc.process_args(args)
        assert parsed["tree_file_path"] == "t.tre"
        assert parsed["trait_data_path"] == "d.tsv"
        assert parsed["trait_column"] == "diet"
        assert parsed["n_sim"] == 100
        assert parsed["seed"] is None
        assert parsed["output_path"] == "out.png"
        assert parsed["json_output"] is False

    def test_overrides(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            trait="diet",
            nsim=200,
            seed=99,
            output="out.png",
            json=True,
        )
        svc = DensityMap.__new__(DensityMap)
        parsed = svc.process_args(args)
        assert parsed["n_sim"] == 200
        assert parsed["seed"] == 99
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
                trait_data=DISCRETE_TRAITS_FILE,
                trait="diet",
                nsim=10,
                seed=42,
                output=tmppath,
                json=False,
            )
            svc = DensityMap(args)
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
                trait_data=DISCRETE_TRAITS_FILE,
                trait="diet",
                nsim=10,
                seed=42,
                output=tmppath,
                json=True,
            )
            svc = DensityMap(args)
            svc.run()
            captured = capsys.readouterr()
            lines = captured.out.strip().split("\n")
            payload = json.loads(lines[-1])
            assert "n_tips" in payload
            assert payload["n_tips"] == 8
            assert "states" in payload
            assert payload["n_sim"] == 10
            assert "plot_output" in payload
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)
