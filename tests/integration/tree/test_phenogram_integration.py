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
class TestPhenogramIntegration:
    def test_basic_invocation(self):
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name
        try:
            testargs = ["phykit", "phenogram", "-t", TREE_SIMPLE, "-d", TRAITS_FILE, "-o", tmppath]
            with patch.object(sys, "argv", testargs):
                Phykit()
            assert os.path.exists(tmppath)
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_alias_traitgram(self):
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name
        try:
            testargs = ["phykit", "traitgram", "-t", TREE_SIMPLE, "-d", TRAITS_FILE, "-o", tmppath]
            with patch.object(sys, "argv", testargs):
                Phykit()
            assert os.path.exists(tmppath)
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_alias_tg(self):
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name
        try:
            testargs = ["phykit", "tg", "-t", TREE_SIMPLE, "-d", TRAITS_FILE, "-o", tmppath]
            with patch.object(sys, "argv", testargs):
                Phykit()
            assert os.path.exists(tmppath)
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name
        try:
            testargs = ["phykit", "phenogram", "-t", TREE_SIMPLE, "-d", TRAITS_FILE, "-o", tmppath, "--json"]
            with patch.object(sys, "argv", testargs):
                Phykit()
            payload = json.loads(mocked_print.call_args.args[0])
            assert payload["n_tips"] == 8
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)
