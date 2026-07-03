import json
import os
import subprocess
import sys
import tempfile
from io import StringIO

import numpy as np
import pytest
from argparse import Namespace
from pathlib import Path
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

import phykit.services.tree.density_map as module
from phykit.services.tree.density_map import DensityMap
from phykit.services.tree.stochastic_character_map import StochasticCharacterMap


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
DISCRETE_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_discrete_traits.tsv")


def test_module_import_does_not_import_numpy_or_matplotlib():
    code = """
import sys
import phykit.services.tree.density_map as module

assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.circular_layout" not in sys.modules
assert "phykit.helpers.color_annotations" not in sys.modules
assert "phykit.services.tree.stochastic_character_map" not in sys.modules
assert "numpy" not in sys.modules
assert "matplotlib" not in sys.modules
assert "matplotlib.pyplot" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = module._LazyNumpy()

    first_array = lazy_np.array
    second_array = lazy_np.array

    assert first_array is second_array
    assert lazy_np.__dict__["array"] is first_array
    assert lazy_np._module is not None


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


class TestTextOutput:
    def test_print_text_output_batches_summary(self, mocker):
        svc = DensityMap.__new__(DensityMap)
        svc.n_sim = 250
        svc.output_path = "density.png"
        printed = mocker.patch("builtins.print")

        svc._print_text_output(["A", "B", "C"], 12)

        printed.assert_called_once_with(
            "Density Map\n"
            "\nStates: 3 (A, B, C)\n"
            "Number of tips: 12\n"
            "Number of simulations: 250\n"
            "Saved density map plot: density.png"
        )


class TestBranchPosteriors:
    def test_compute_branch_posteriors_scans_histories_once(self, monkeypatch):
        svc = DensityMap.__new__(DensityMap)
        clade_id = 101
        states = ["A", "B", "C"]
        mappings = [
            {"branch_histories": {clade_id: [(0.0, 0), (0.35, 1), (0.75, 2)]}},
            {"branch_histories": {clade_id: [(0.0, 1), (0.55, 2)]}},
            {"branch_histories": {}},
        ]

        def fail_state_lookup(*_args, **_kwargs):
            raise AssertionError("branch posterior path should scan histories once")

        monkeypatch.setattr(
            DensityMap,
            "_get_state_at_position",
            staticmethod(fail_state_lookup),
        )

        posteriors = svc._compute_branch_posteriors(
            mappings,
            clade_id,
            branch_length=1.0,
            states=states,
            n_segments=4,
        )

        expected = np.array(
            [
                [1 / 3, 1 / 3, 0.0],
                [0.0, 2 / 3, 0.0],
                [0.0, 1 / 3, 1 / 3],
                [0.0, 0.0, 2 / 3],
            ]
        )
        np.testing.assert_allclose(posteriors, expected)

    def test_zero_length_branch_posteriors_preserve_missing_history_denominator(self):
        svc = DensityMap.__new__(DensityMap)
        clade_id = 202
        mappings = [
            {"branch_histories": {clade_id: [(0.0, 0)]}},
            {"branch_histories": {clade_id: [(0.0, 1)]}},
            {"branch_histories": {}},
        ]

        posteriors = svc._compute_branch_posteriors(
            mappings,
            clade_id,
            branch_length=0.0,
            states=["A", "B"],
            n_segments=3,
        )

        expected = np.array(
            [
                [1 / 3, 1 / 3],
                [1 / 3, 1 / 3],
                [1 / 3, 1 / 3],
            ]
        )
        np.testing.assert_allclose(posteriors, expected)


class TestRun:
    def test_terminal_clades_fast_path_does_not_call_get_terminals(self, monkeypatch):
        tree = Phylo.read(StringIO("((a:1,b:2):3,c:4);"), "newick")

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("standard tree should use direct traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        tips = DensityMap._terminal_clades(tree)

        assert [tip.name for tip in tips] == ["a", "b", "c"]

    def test_direct_traversals_preserve_mixed_child_order(self, monkeypatch):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1,G:1);"),
            "newick",
        )
        expected_preorder = list(tree.find_clades(order="preorder"))
        expected_terminals = tree.get_terminals()

        def fail_traversal(*args, **kwargs):
            raise AssertionError("standard tree should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        assert list(DensityMap._iter_preorder(tree.root)) == expected_preorder
        assert DensityMap._terminal_clades(tree) == expected_terminals

    def _stub_expensive_run_tail(self, svc, monkeypatch, captured_tree):
        def fit_q_matrix(_self, tree, *_args, **_kwargs):
            captured_tree["tree"] = tree
            import numpy as np

            return np.array([[-1.0, 1.0], [1.0, -1.0]]), -1.25

        monkeypatch.setattr(
            StochasticCharacterMap,
            "_fit_q_matrix",
            fit_q_matrix,
        )
        monkeypatch.setattr(
            StochasticCharacterMap,
            "_felsenstein_pruning",
            lambda _self, tree, *_args, **_kwargs: ({id(tree.root): [1.0, 1.0]}, None),
        )
        monkeypatch.setattr(
            StochasticCharacterMap,
            "_build_parent_map",
            lambda *_args, **_kwargs: {},
        )
        monkeypatch.setattr(
            StochasticCharacterMap,
            "_build_simulation_metadata",
            lambda *_args, **_kwargs: ([], []),
        )
        monkeypatch.setattr(svc, "_plot_density_map", lambda *_args, **_kwargs: None)
        monkeypatch.setattr(svc, "_terminal_clades", lambda tree: tree.get_terminals())

    def test_all_tips_present_uses_read_only_tree_without_copy_or_prune(
        self, monkeypatch, capsys
    ):
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        tip_states = {"A": "x", "B": "x", "C": "y", "D": "y"}
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=DISCRETE_TRAITS_FILE,
            trait="diet",
            nsim=0,
            seed=42,
            output="density.png",
            json=False,
        )
        svc = DensityMap(args)
        captured_tree = {}

        monkeypatch.setattr(svc, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(
            svc,
            "read_tree_file",
            lambda: (_ for _ in ()).throw(
                AssertionError("mutable tree read should not be used")
            ),
        )
        monkeypatch.setattr(
            StochasticCharacterMap,
            "_parse_discrete_trait_file",
            lambda *_args, **_kwargs: tip_states,
        )
        monkeypatch.setattr(
            svc,
            "_fast_copy",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("tree should not be copied")
            ),
        )
        monkeypatch.setattr(
            StochasticCharacterMap,
            "_prune_tree_to_tip_states",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("no-prune path should skip pruning")
            ),
        )
        self._stub_expensive_run_tail(svc, monkeypatch, captured_tree)

        svc.run()

        capsys.readouterr()
        assert captured_tree["tree"] is tree
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}

    def test_missing_tip_states_copy_before_pruning(self, monkeypatch, capsys):
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        tip_states = {"A": "x", "B": "x", "C": "y"}
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=DISCRETE_TRAITS_FILE,
            trait="diet",
            nsim=0,
            seed=42,
            output="density.png",
            json=False,
        )
        svc = DensityMap(args)
        original_fast_copy = svc._fast_copy
        copied_trees = []
        captured_tree = {}

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        monkeypatch.setattr(svc, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(
            StochasticCharacterMap,
            "_parse_discrete_trait_file",
            lambda *_args, **_kwargs: tip_states,
        )
        monkeypatch.setattr(svc, "_fast_copy", copy_spy)
        self._stub_expensive_run_tail(svc, monkeypatch, captured_tree)

        svc.run()

        capsys.readouterr()
        assert len(copied_trees) == 1
        assert captured_tree["tree"] is copied_trees[0]
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}
        assert {tip.name for tip in copied_trees[0].get_terminals()} == {
            "A",
            "B",
            "C",
        }

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

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_density_map_uses_direct_tree_traversal(self, monkeypatch, circular):
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
                circular=circular,
            )
            svc = DensityMap(args)
            tree = Phylo.read(StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")
            scm = StochasticCharacterMap.__new__(StochasticCharacterMap)
            parent_map = scm._build_parent_map(tree)

            def fail_traversal(*args, **kwargs):
                raise AssertionError("plot setup should use direct traversal")

            monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
            monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

            svc._plot_density_map(
                tree, [], ["0", "1"], tmppath,
                parent_map=parent_map, scm=scm,
            )

            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_plot_density_map_reuses_preorder_for_node_positions(
        self, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import phykit.helpers.plot_config as plot_config

        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=DISCRETE_TRAITS_FILE,
            trait="diet",
            nsim=10,
            seed=42,
            output=str(tmp_path / "density_preorder_positions.png"),
            json=False,
            circular=False,
        )
        svc = DensityMap(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")
        scm = StochasticCharacterMap.__new__(StochasticCharacterMap)
        parent_map = scm._build_parent_map(tree)
        expected_ids = [id(clade) for clade in svc._iter_preorder(tree.root)]
        original_compute_node_positions = plot_config.compute_node_positions
        calls = []

        def assert_preorder_reused(
            tree_arg, parent_map_arg, cladogram=False, preorder_clades=None
        ):
            assert preorder_clades is not None
            calls.append([id(clade) for clade in preorder_clades])
            return original_compute_node_positions(
                tree_arg,
                parent_map_arg,
                cladogram=cladogram,
                preorder_clades=preorder_clades,
            )

        monkeypatch.setattr(
            plot_config, "compute_node_positions", assert_preorder_reused
        )

        svc._plot_density_map(
            tree, [], ["0", "1"], args.output,
            parent_map=parent_map, scm=scm,
        )

        assert calls == [expected_ids]
        assert os.path.exists(args.output)

    def test_plot_density_map_reuses_clade_lists_for_circular_coords(
        self, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import phykit.helpers.circular_layout as circular_layout

        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=DISCRETE_TRAITS_FILE,
            trait="diet",
            nsim=10,
            seed=42,
            output=str(tmp_path / "density_circular_precomputed_coords.png"),
            json=False,
            circular=True,
            no_title=True,
            ylabel_fontsize=0,
        )
        svc = DensityMap(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")
        scm = StochasticCharacterMap.__new__(StochasticCharacterMap)
        parent_map = scm._build_parent_map(tree)
        preorder_clades = list(svc._iter_preorder(tree.root))
        tips = [clade for clade in preorder_clades if not clade.clades]
        expected_preorder_ids = [id(clade) for clade in preorder_clades]
        expected_terminal_ids = [id(tip) for tip in tips]
        original_compute_circular_coords = circular_layout.compute_circular_coords
        calls = []

        def assert_clade_lists_reused(
            tree_arg,
            node_x_arg,
            parent_map_arg,
            preorder_clades=None,
            terminal_clades=None,
        ):
            assert preorder_clades is not None
            assert terminal_clades is not None
            calls.append(
                (
                    [id(clade) for clade in preorder_clades],
                    [id(clade) for clade in terminal_clades],
                )
            )
            return original_compute_circular_coords(
                tree_arg,
                node_x_arg,
                parent_map_arg,
                preorder_clades=preorder_clades,
                terminal_clades=terminal_clades,
            )

        monkeypatch.setattr(
            circular_layout, "compute_circular_coords", assert_clade_lists_reused
        )

        svc._plot_density_map(
            tree, [], ["0", "1"], args.output,
            parent_map=parent_map, scm=scm,
        )

        assert calls == [(expected_preorder_ids, expected_terminal_ids)]
        assert os.path.exists(args.output)

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_density_map_batches_branch_segments(
        self, monkeypatch, tmp_path, circular
    ):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=DISCRETE_TRAITS_FILE,
            trait="diet",
            nsim=10,
            seed=42,
            output=str(tmp_path / f"density_batched_{circular}.png"),
            json=False,
            circular=circular,
        )
        svc = DensityMap(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")
        scm = StochasticCharacterMap.__new__(StochasticCharacterMap)
        parent_map = scm._build_parent_map(tree)
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("density-map branches should be batched")

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        svc._plot_density_map(
            tree, [], ["0", "1", "2"], args.output,
            parent_map=parent_map, scm=scm,
        )

        assert len(line_collections) >= 2
        assert os.path.exists(args.output)
        assert os.path.getsize(args.output) > 0
