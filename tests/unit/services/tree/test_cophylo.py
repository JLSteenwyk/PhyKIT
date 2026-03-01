import json
import os
import tempfile

import pytest
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.cophylo import Cophylo


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE1 = str(SAMPLE_FILES / "tree_simple.tre")
TREE2 = str(SAMPLE_FILES / "tree_simple_other_topology.tre")


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree1="t1.tre",
            tree2="t2.tre",
            output="out.png",
            mapping=None,
        )
        svc = Cophylo.__new__(Cophylo)
        parsed = svc.process_args(args)
        assert parsed["tree1_file_path"] == "t1.tre"
        assert parsed["tree2_file_path"] == "t2.tre"
        assert parsed["output_path"] == "out.png"
        assert parsed["mapping_file_path"] is None
        assert parsed["json_output"] is False

    def test_overrides(self):
        args = Namespace(
            tree1="t1.tre",
            tree2="t2.tre",
            output="out.png",
            mapping="map.tsv",
            json=True,
        )
        svc = Cophylo.__new__(Cophylo)
        parsed = svc.process_args(args)
        assert parsed["mapping_file_path"] == "map.tsv"
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
                tree1=TREE1,
                tree2=TREE2,
                output=tmppath,
                mapping=None,
                json=False,
            )
            svc = Cophylo(args)
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
                tree1=TREE1,
                tree2=TREE2,
                output=tmppath,
                mapping=None,
                json=True,
            )
            svc = Cophylo(args)
            svc.run()
            captured = capsys.readouterr()
            lines = captured.out.strip().split("\n")
            payload = json.loads(lines[-1])
            assert "n_tips_tree1" in payload
            assert "n_tips_tree2" in payload
            assert "n_matched" in payload
            assert payload["n_tips_tree1"] == 8
            assert payload["n_tips_tree2"] == 8
            assert payload["n_matched"] == 8
            assert "matched_taxa" in payload
            assert "plot_output" in payload
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_with_mapping_file(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        # Create a mapping file that maps a subset of taxa
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as mf:
            mf.write("raccoon\traccoon\n")
            mf.write("bear\tbear\n")
            mf.write("dog\tdog\n")
            mf.write("cat\tcat\n")
            mappath = mf.name

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = Namespace(
                tree1=TREE1,
                tree2=TREE2,
                output=tmppath,
                mapping=mappath,
                json=True,
            )
            svc = Cophylo(args)
            svc.run()
            captured_lines = []  # json printed via print_json
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)
            if os.path.exists(mappath):
                os.unlink(mappath)


class TestValidation:
    def test_too_few_shared_taxa(self):
        """Two trees with < 2 shared taxa should raise an error."""
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        # Create two tiny trees with no overlap
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tre", delete=False
        ) as f1:
            f1.write("(A:1,B:1);")
            tree1_path = f1.name

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tre", delete=False
        ) as f2:
            f2.write("(C:1,D:1);")
            tree2_path = f2.name

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            from phykit.errors import PhykitUserError

            args = Namespace(
                tree1=tree1_path,
                tree2=tree2_path,
                output=tmppath,
                mapping=None,
                json=False,
            )
            with pytest.raises((PhykitUserError, SystemExit)):
                svc = Cophylo(args)
                svc.run()
        finally:
            for p in [tree1_path, tree2_path, tmppath]:
                if os.path.exists(p):
                    os.unlink(p)
