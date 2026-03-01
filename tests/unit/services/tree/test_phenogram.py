import json
import os
import tempfile
import pytest
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from phykit.services.tree.phenogram import Phenogram

here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(tree="t.tre", trait_data="d.tsv", output="out.png", json=False)
        svc = Phenogram.__new__(Phenogram)
        parsed = svc.process_args(args)
        assert parsed["tree_file_path"] == "t.tre"
        assert parsed["output_path"] == "out.png"
        assert parsed["json_output"] is False

    def test_json_override(self):
        args = Namespace(tree="t.tre", trait_data="d.tsv", output="out.png", json=True)
        svc = Phenogram.__new__(Phenogram)
        parsed = svc.process_args(args)
        assert parsed["json_output"] is True


class TestRun:
    def test_plot_created(self):
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name
        try:
            args = Namespace(tree=TREE_SIMPLE, trait_data=TRAITS_FILE, output=tmppath, json=False)
            svc = Phenogram(args)
            svc.run()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name
        try:
            args = Namespace(tree=TREE_SIMPLE, trait_data=TRAITS_FILE, output=tmppath, json=True)
            svc = Phenogram(args)
            svc.run()
            payload = json.loads(mocked_print.call_args.args[0])
            assert "n_tips" in payload
            assert payload["n_tips"] == 8
            assert "sigma2" in payload
            assert "tip_values" in payload
            assert "ancestral_estimates" in payload
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)
