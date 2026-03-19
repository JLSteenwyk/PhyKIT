import json
import os
import tempfile

import pytest
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.trait_rate_map import TraitRateMap


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


def _make_args(**kwargs):
    defaults = dict(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        output="/tmp/test_rate_map.png",
        trait=None,
        json=False,
    )
    defaults.update(kwargs)
    return Namespace(**defaults)


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            output="out.png",
        )
        svc = TraitRateMap.__new__(TraitRateMap)
        parsed = svc.process_args(args)
        assert parsed["tree_file_path"] == "t.tre"
        assert parsed["trait_data_path"] == "d.tsv"
        assert parsed["output_path"] == "out.png"
        assert parsed["trait_column"] is None
        assert parsed["json_output"] is False

    def test_with_trait_column(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            output="out.png",
            trait="body_mass",
            json=True,
        )
        svc = TraitRateMap.__new__(TraitRateMap)
        parsed = svc.process_args(args)
        assert parsed["trait_column"] == "body_mass"
        assert parsed["json_output"] is True


class TestAncestralReconstruction:
    def test_node_values_computed(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(output=tmppath)
            svc = TraitRateMap(args)
            tree = svc.read_tree_file()
            tree_tips = svc.get_tip_names_from_tree(tree)
            trait_values = svc._parse_single_trait_data(TRAITS_FILE, tree_tips)
            node_values = svc._ancestral_reconstruction(tree, trait_values)

            # Should have values for all nodes (tips + internal)
            all_nodes = list(tree.find_clades())
            for node in all_nodes:
                assert id(node) in node_values, f"Missing value for node {node.name}"
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)


class TestBranchRates:
    def test_rates_nonnegative(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(output=tmppath)
            svc = TraitRateMap(args)
            tree = svc.read_tree_file()
            tree_tips = svc.get_tip_names_from_tree(tree)
            trait_values = svc._parse_single_trait_data(TRAITS_FILE, tree_tips)
            node_values = svc._ancestral_reconstruction(tree, trait_values)
            parent_map = svc._build_parent_map(tree)
            node_labels = svc._label_internal_nodes(tree)
            branch_rates = svc._compute_branch_rates(tree, node_values, parent_map, node_labels)

            for entry in branch_rates:
                assert entry["rate"] >= 0, f"Negative rate for branch {entry['parent']} -> {entry['child']}"
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_rate_zero_for_no_change(self):
        """If parent and child have the same value, rate should be 0."""
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        # Create a trait file where all taxa have the same value
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            equal_trait_path = f.name
            f.write("raccoon\t1.0\n")
            f.write("bear\t1.0\n")
            f.write("sea_lion\t1.0\n")
            f.write("seal\t1.0\n")
            f.write("monkey\t1.0\n")
            f.write("cat\t1.0\n")
            f.write("weasel\t1.0\n")
            f.write("dog\t1.0\n")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(trait_data=equal_trait_path, output=tmppath)
            svc = TraitRateMap(args)
            tree = svc.read_tree_file()
            tree_tips = svc.get_tip_names_from_tree(tree)
            trait_values = svc._parse_single_trait_data(equal_trait_path, tree_tips)
            node_values = svc._ancestral_reconstruction(tree, trait_values)
            parent_map = svc._build_parent_map(tree)
            node_labels = svc._label_internal_nodes(tree)
            branch_rates = svc._compute_branch_rates(tree, node_values, parent_map, node_labels)

            for entry in branch_rates:
                assert entry["rate"] == pytest.approx(0.0, abs=1e-10), \
                    f"Rate should be 0 for branch {entry['parent']} -> {entry['child']}, got {entry['rate']}"
        finally:
            if os.path.exists(equal_trait_path):
                os.unlink(equal_trait_path)
            if os.path.exists(tmppath):
                os.unlink(tmppath)


class TestRun:
    def test_creates_png(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(output=tmppath)
            svc = TraitRateMap(args)
            svc.run()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_cladogram_mode(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(output=tmppath, cladogram=True)
            svc = TraitRateMap(args)
            svc.run()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_circular_mode(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(output=tmppath, circular=True)
            svc = TraitRateMap(args)
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
            args = _make_args(output=tmppath, json=True)
            svc = TraitRateMap(args)
            svc.run()
            captured = capsys.readouterr()
            lines = captured.out.strip().split("\n")
            payload = json.loads(lines[-1])
            assert "n_taxa" in payload
            assert payload["n_taxa"] == 8
            assert "mean_rate" in payload
            assert "min_rate" in payload
            assert "max_rate" in payload
            assert "branches" in payload
            assert len(payload["branches"]) > 0
            assert "output_file" in payload
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_multi_trait_with_column(self):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            args = _make_args(
                trait_data=MULTI_TRAITS_FILE,
                output=tmppath,
                trait="body_mass",
                json=True,
            )
            svc = TraitRateMap(args)
            svc.run()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)
