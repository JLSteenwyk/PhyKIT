import os
import subprocess
import sys

import pytest
import json
import numpy as np
from argparse import Namespace
from io import StringIO
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

import phykit.services.tree.phylomorphospace as phylomorphospace_module
from phykit.services.tree.phylomorphospace import (
    Phylomorphospace,
    _root_distance_max,
)
from phykit.helpers.trait_parsing import parse_multi_trait_file
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


def test_module_import_does_not_import_numpy_or_matplotlib():
    code = """
import sys
import phykit.services.tree.phylomorphospace as module

assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert callable(module.parse_multi_trait_file)
assert callable(module.trait_matrix_from_rows)
assert "pickle" not in sys.modules
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.trait_parsing" not in sys.modules
assert "numpy" not in sys.modules
assert "matplotlib" not in sys.modules
assert "matplotlib.pyplot" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=MULTI_TRAITS_FILE,
        trait_x="body_mass",
        trait_y="brain_size",
        color_by=None,
        plot_output="phylomorphospace_plot.png",
        json=False,
    )


@pytest.fixture
def auto_select_two_trait_file(tmp_path):
    """Create a trait file with exactly 2 traits for auto-selection."""
    trait_file = tmp_path / "two_traits.tsv"
    trait_file.write_text(
        "taxon\ttrait_a\ttrait_b\n"
        "raccoon\t1.0\t2.0\n"
        "bear\t2.0\t3.0\n"
        "sea_lion\t3.0\t4.0\n"
        "seal\t4.0\t5.0\n"
        "monkey\t5.0\t6.0\n"
        "cat\t6.0\t7.0\n"
        "weasel\t7.0\t8.0\n"
        "dog\t8.0\t9.0\n"
    )
    return str(trait_file)


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(tree="t.tre", trait_data="d.tsv")
        svc = Phylomorphospace(args)
        assert svc.tree_file_path == "t.tre"
        assert svc.trait_data_path == "d.tsv"
        assert svc.trait_x is None
        assert svc.trait_y is None
        assert svc.color_by is None
        assert svc.plot_output == "phylomorphospace_plot.png"
        assert svc.json_output is False

    def test_overrides(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            trait_x="body_mass",
            trait_y="brain_size",
            color_by="longevity",
            plot_output="my_plot.png",
            json=True,
        )
        svc = Phylomorphospace(args)
        assert svc.trait_x == "body_mass"
        assert svc.trait_y == "brain_size"
        assert svc.color_by == "longevity"
        assert svc.plot_output == "my_plot.png"
        assert svc.json_output is True

    def test_auto_select_two_trait_file(self, auto_select_two_trait_file, tmp_path):
        plot_path = str(tmp_path / "auto.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=auto_select_two_trait_file,
            trait_x=None,
            trait_y=None,
            color_by=None,
            plot_output=plot_path,
            json=False,
        )
        svc = Phylomorphospace(args)
        svc.run()
        assert Path(plot_path).exists()


class TestRun:
    def test_run_uses_unmodified_tree_read(self, mocker):
        args = Namespace(
            tree="/some/path/to/file.tre",
            trait_data="/some/path/to/traits.tsv",
            trait_x="x",
            trait_y="y",
            color_by=None,
            plot_output="plot.png",
            json=False,
        )
        svc = Phylomorphospace(args)
        tree = object()
        mocked_read = mocker.patch.object(
            Phylomorphospace,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(
            Phylomorphospace,
            "read_tree_file",
            side_effect=AssertionError("run should not copy cached trees"),
        )
        mocked_validate = mocker.patch.object(
            Phylomorphospace, "validate_tree"
        )
        mocker.patch.object(
            Phylomorphospace,
            "get_tip_names_from_tree",
            return_value=["a", "b", "c"],
        )
        mocker.patch(
            "phykit.services.tree.phylomorphospace.parse_multi_trait_file",
            return_value=(
                ["x", "y"],
                {
                    "a": [1.0, 2.0],
                    "b": [2.0, 3.0],
                    "c": [3.0, 4.0],
                },
            ),
        )
        tree_pruned = object()
        mocked_reconstruct = mocker.patch.object(
            Phylomorphospace,
            "_reconstruct_ancestral_scores",
            return_value=({}, {}, tree_pruned),
        )
        mocked_plot = mocker.patch.object(
            Phylomorphospace, "_plot_phylomorphospace"
        )
        mocker.patch("builtins.print")

        svc.run()

        mocked_read.assert_called_once_with()
        mocked_validate.assert_called_once_with(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="phylomorphospace",
        )
        assert mocked_reconstruct.call_args.args[0] is tree
        assert mocked_reconstruct.call_args.args[1].tolist() == [
            [1.0, 2.0],
            [2.0, 3.0],
            [3.0, 4.0],
        ]
        assert mocked_plot.call_args.args[1].tolist() == [
            [1.0, 2.0],
            [2.0, 3.0],
            [3.0, 4.0],
        ]
        assert mocked_plot.call_args.args[8].tolist() == [
            [1.0, 2.0],
            [2.0, 3.0],
            [3.0, 4.0],
        ]

    def test_basic_text_output(self, default_args, tmp_path, capsys):
        default_args.plot_output = str(tmp_path / "test.png")
        svc = Phylomorphospace(default_args)
        svc.run()
        captured = capsys.readouterr()
        assert "Saved phylomorphospace plot" in captured.out
        assert Path(default_args.plot_output).exists()

    def test_json_output(self, default_args, tmp_path, capsys):
        default_args.json = True
        default_args.plot_output = str(tmp_path / "test.png")
        svc = Phylomorphospace(default_args)
        svc.run()
        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        assert payload["trait_x"] == "body_mass"
        assert payload["trait_y"] == "brain_size"
        assert "tip_data" in payload
        assert "raccoon" in payload["tip_data"]
        assert "plot_output" in payload

    def test_trait_not_found_error(self, tmp_path):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            trait_x="nonexistent",
            trait_y="brain_size",
            color_by=None,
            plot_output=str(tmp_path / "test.png"),
            json=False,
        )
        svc = Phylomorphospace(args)
        with pytest.raises(PhykitUserError):
            svc.run()

    def test_only_one_trait_flag_error(self, tmp_path):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            trait_x="body_mass",
            trait_y=None,
            color_by=None,
            plot_output=str(tmp_path / "test.png"),
            json=False,
        )
        svc = Phylomorphospace(args)
        with pytest.raises(PhykitUserError):
            svc.run()

    def test_plot_creates_file(self, default_args, tmp_path):
        default_args.plot_output = str(tmp_path / "test.png")
        svc = Phylomorphospace(default_args)
        svc.run()
        assert Path(default_args.plot_output).exists()
        assert Path(default_args.plot_output).stat().st_size > 0


class TestPlot:
    def test_parse_numeric_color_file_uses_fromiter(self, tmp_path, monkeypatch):
        color_file = tmp_path / "colors.tsv"
        ordered_names = ["a", "b", "c"]
        color_file.write_text(
            "   # ignored\n"
            "a\t1.5\textra\n"
            "b\t2.25\textra\n"
            "c\t3.75\textra\n"
        )
        svc = Phylomorphospace(
            Namespace(
                tree="unused",
                trait_data="unused",
                trait_x="x",
                trait_y="y",
                color_by=str(color_file),
                plot_output=str(tmp_path / "plot.png"),
                json=False,
            )
        )

        def fail_array(*_args, **_kwargs):
            raise AssertionError("numeric color files should use np.fromiter")

        monkeypatch.setattr(phylomorphospace_module.np, "array", fail_array)

        values, categories, kind = svc._parse_color_by(
            str(color_file),
            [],
            np.empty((3, 0)),
            ordered_names,
        )

        assert kind == "continuous"
        assert categories == []
        np.testing.assert_allclose(values, [1.5, 2.25, 3.75])

    def test_with_color_by_column(self, tmp_path):
        plot_path = str(tmp_path / "color_col.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            trait_x="body_mass",
            trait_y="brain_size",
            color_by="longevity",
            plot_output=plot_path,
            json=False,
        )
        svc = Phylomorphospace(args)
        svc.run()
        assert Path(plot_path).exists()
        assert Path(plot_path).stat().st_size > 0

    def test_with_color_by_discrete_file(self, tmp_path):
        color_file = tmp_path / "groups.tsv"
        color_file.write_text(
            "raccoon\tgroup_A\n"
            "bear\tgroup_A\n"
            "sea_lion\tgroup_B\n"
            "seal\tgroup_B\n"
            "monkey\tgroup_A\n"
            "cat\tgroup_B\n"
            "weasel\tgroup_A\n"
            "dog\tgroup_B\n"
        )
        plot_path = str(tmp_path / "color_disc.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            trait_x="body_mass",
            trait_y="brain_size",
            color_by=str(color_file),
            plot_output=plot_path,
            json=False,
        )
        svc = Phylomorphospace(args)
        svc.run()
        assert Path(plot_path).exists()
        assert Path(plot_path).stat().st_size > 0


class TestReconstructAncestralScores:
    @pytest.fixture(autouse=True)
    def _make_trait_data(self, tmp_path):
        """Create local trait data so tests don't depend on sample_files."""
        self._trait_file = str(tmp_path / "traits.tsv")
        Path(self._trait_file).write_text(
            "taxon\tbody_mass\tbrain_size\n"
            "raccoon\t1.04\t1.60\n"
            "bear\t2.39\t2.66\n"
            "sea_lion\t2.30\t2.74\n"
            "seal\t1.88\t2.45\n"
            "monkey\t0.60\t1.85\n"
            "cat\t0.56\t1.30\n"
            "weasel\t-0.30\t0.85\n"
            "dog\t1.18\t1.87\n"
        )

    def test_prune_setup_uses_fast_tip_names_and_shared_prune_helper(self, mocker):
        class ContainsFailList(list):
            def __contains__(self, item):
                raise AssertionError("ordered_names should be converted to a set")

        svc = self._make_svc()
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree_copy = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ContainsFailList(["A", "C", "D"])
        data = np.array([[0.0, 0.0], [2.0, 2.0], [4.0, 4.0]])
        spy = mocker.spy(svc, "get_tip_names_from_tree")
        prune_helper = mocker.spy(svc, "_tips_to_prune_for_ordered_names")
        fast_copy = mocker.patch.object(svc, "_fast_copy", return_value=tree_copy)
        prune = mocker.spy(svc, "prune_tree_using_taxa_list")

        _, _, tree_pruned = svc._reconstruct_ancestral_scores(
            tree, data, ordered_names
        )

        assert spy.call_count == 1
        prune_helper.assert_called_once_with(["A", "B", "C", "D"], ordered_names)
        fast_copy.assert_called_once_with(tree)
        prune.assert_called_once_with(tree_copy, ["B"])
        assert {tip.name for tip in tree_pruned.get_terminals()} == {"A", "C", "D"}

    def test_reconstruct_skips_copy_when_all_tips_are_shared(self, mocker):
        svc = self._make_svc()
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        data = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
        fast_copy = mocker.patch.object(
            svc,
            "_fast_copy",
            side_effect=AssertionError("all-shared reconstruction should not copy tree"),
        )

        _, _, tree_pruned = svc._reconstruct_ancestral_scores(
            tree, data, ordered_names
        )

        fast_copy.assert_not_called()
        assert tree_pruned is tree

    def _make_svc(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=self._trait_file,
            trait_x="body_mass",
            trait_y="brain_size",
            color_by=None,
            plot_output="unused.png",
            json=False,
        )
        return Phylomorphospace(args)

    def _get_data(self, svc):
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits = parse_multi_trait_file(
            svc.trait_data_path, tree_tips
        )
        ordered_names = sorted(traits.keys())
        x_idx = trait_names.index("body_mass")
        y_idx = trait_names.index("brain_size")
        data = np.array(
            [[traits[name][x_idx], traits[name][y_idx]] for name in ordered_names]
        )
        return tree, ordered_names, data

    def test_all_nodes_scored(self):
        svc = self._make_svc()
        tree, ordered_names, data = self._get_data(svc)
        node_estimates, node_distances, tree_pruned = (
            svc._reconstruct_ancestral_scores(tree, data, ordered_names)
        )

        # All clades should have estimates
        all_clades = list(tree_pruned.find_clades())
        for clade in all_clades:
            assert id(clade) in node_estimates

    def test_tips_match_input(self):
        svc = self._make_svc()
        tree, ordered_names, data = self._get_data(svc)
        node_estimates, _, tree_pruned = svc._reconstruct_ancestral_scores(
            tree, data, ordered_names
        )

        # Check that tip estimates match input data
        for tip in tree_pruned.get_terminals():
            if tip.name in ordered_names:
                idx = ordered_names.index(tip.name)
                est = node_estimates[id(tip)]
                np.testing.assert_array_almost_equal(est, data[idx])

    def test_root_in_range(self):
        svc = self._make_svc()
        tree, ordered_names, data = self._get_data(svc)
        node_estimates, _, tree_pruned = svc._reconstruct_ancestral_scores(
            tree, data, ordered_names
        )

        root_est = node_estimates[id(tree_pruned.root)]
        # Root estimate should be within the range of tip values
        for dim in range(2):
            assert root_est[dim] >= data[:, dim].min() - 1.0
            assert root_est[dim] <= data[:, dim].max() + 1.0

    def test_node_distances_fast_path_does_not_call_distance(self, monkeypatch):
        svc = self._make_svc()
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:3):0.5;"), "newick")
        ordered_names = ["A", "B", "C"]
        data = np.array([[0.0, 0.0], [2.0, 2.0], [4.0, 4.0]])

        def fail_distance(self, *args, **kwargs):
            raise AssertionError("distance fallback should not be called")

        monkeypatch.setattr(type(tree), "distance", fail_distance)
        _, node_distances, tree_pruned = svc._reconstruct_ancestral_scores(
            tree, data, ordered_names
        )

        assert node_distances[id(tree_pruned.root)] == pytest.approx(0.0)
        tip_distances = {
            tip.name: node_distances[id(tip)]
            for tip in tree_pruned.get_terminals()
        }
        assert tip_distances == pytest.approx({"A": 2.0, "B": 2.0, "C": 3.0})

    def test_root_distance_max_scans_values_once(self):
        class SinglePassValues:
            def __init__(self, values):
                self.values = values
                self.iterations = 0

            def __iter__(self):
                self.iterations += 1
                if self.iterations > 1:
                    raise AssertionError("distances should be scanned once")
                return iter(self.values)

        values = SinglePassValues([0.0, 2.0, 1.5, 3.0])

        assert _root_distance_max(values) == 3.0
        assert values.iterations == 1
        assert _root_distance_max([]) == 1.0

    def test_reconstruct_ancestral_scores_uses_direct_traversal(self, monkeypatch):
        svc = self._make_svc()
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:3):0.5;"), "newick")
        ordered_names = ["A", "B", "C"]
        data = np.array([[0.0, 0.0], [2.0, 0.0], [4.0, 0.0]])

        def fail_generic_traversal(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_generic_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_generic_traversal)
        monkeypatch.setattr(TreeMixin, "depths", fail_generic_traversal)
        monkeypatch.setattr(type(tree), "distance", fail_generic_traversal)

        node_estimates, node_distances, tree_pruned = (
            svc._reconstruct_ancestral_scores(tree, data, ordered_names)
        )

        root = tree_pruned.root
        left_internal = root.clades[0]
        tips = {
            root.clades[0].clades[0].name: root.clades[0].clades[0],
            root.clades[0].clades[1].name: root.clades[0].clades[1],
            root.clades[1].name: root.clades[1],
        }

        np.testing.assert_allclose(node_estimates[id(left_internal)], [1.0, 0.0])
        np.testing.assert_allclose(node_estimates[id(root)], [2.0, 0.0])
        assert node_distances[id(root)] == pytest.approx(0.0)
        assert node_distances[id(tips["A"])] == pytest.approx(2.0)
        assert node_distances[id(tips["C"])] == pytest.approx(3.0)

    def test_preorder_clades_direct_preserves_order_with_mixed_child_counts(self):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("direct preorder should push children explicitly")

        class Clade:
            def __init__(self, name, clades=None):
                self.name = name
                self.clades = NoReversedList(clades or [])

        tree = type(
            "Tree",
            (),
            {
                "root": Clade(
                    "root",
                    [
                        Clade("A"),
                        Clade("inner", [Clade("B")]),
                        Clade("poly", [Clade("C"), Clade("D"), Clade("E")]),
                    ],
                )
            },
        )()

        clades = Phylomorphospace._preorder_clades_direct(tree)

        assert [clade.name for clade in clades] == [
            "root",
            "A",
            "inner",
            "B",
            "poly",
            "C",
            "D",
            "E",
        ]

    def test_plot_phylomorphospace_uses_direct_preorder(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        svc = self._make_svc()
        svc.plot_output = str(tmp_path / "phylomorphospace.png")
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:3):0.5;"), "newick")
        ordered_names = ["A", "B", "C"]
        data = np.array([[0.0, 0.0], [2.0, 0.0], [4.0, 0.0]])
        node_estimates = {}
        node_distances = {}
        stack = [(tree.root, 0.0)]
        while stack:
            clade, distance = stack.pop()
            node_estimates[id(clade)] = np.array([distance, distance * 0.5])
            node_distances[id(clade)] = distance
            for child in clade.clades:
                stack.append((child, distance + (child.branch_length or 0.0)))

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("plot branch setup should use direct preorder")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        svc._plot_phylomorphospace(
            tree,
            data,
            ordered_names,
            node_estimates,
            node_distances,
            "body_mass",
            "brain_size",
            ["body_mass", "brain_size"],
            data,
        )

        assert Path(svc.plot_output).exists()
