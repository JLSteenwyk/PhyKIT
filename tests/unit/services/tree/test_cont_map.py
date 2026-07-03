import json
import os
import subprocess
import sys
import tempfile

import pytest
from argparse import Namespace
from io import StringIO
from pathlib import Path

import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, TreeMixin

from phykit.errors import PhykitUserError
import phykit.services.tree.cont_map as module
from phykit.services.tree.cont_map import ContMap, _value_range


def test_module_import_does_not_import_numpy():
    code = """
import sys
import phykit.services.tree.cont_map as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "pickle" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.circular_layout" not in sys.modules
assert "phykit.helpers.color_annotations" not in sys.modules
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = module._LazyNumpy()

    first_array = lazy_np.array
    second_array = lazy_np.array

    assert first_array is second_array
    assert lazy_np.__dict__["array"] is first_array
    assert lazy_np._module is not None


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


class TestTraitParsing:
    def test_comments_and_blanks(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "# comment\n\nraccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\n"
            "seal\t4.0\nmonkey\t5.0\ncat\t6.0\nweasel\t7.0\ndog\t8.0\n"
        )
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=str(trait_file), output="plot.png", json=False
        )
        svc = ContMap(args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        traits = svc._parse_single_trait_data(str(trait_file), tree_tips)

        assert len(traits) == 8
        assert traits["raccoon"] == pytest.approx(1.0)

    def test_all_shared_trait_file_emits_no_warnings(self, tmp_path, capsys):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "raccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\nseal\t4.0\n"
            "monkey\t5.0\ncat\t6.0\nweasel\t7.0\ndog\t8.0\n"
        )
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=str(trait_file), output="plot.png", json=False
        )
        svc = ContMap(args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        traits = svc._parse_single_trait_data(str(trait_file), tree_tips)

        stderr = capsys.readouterr().err
        assert set(traits) == set(tree_tips)
        assert stderr == ""

    def test_extra_columns_error(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "raccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\ndog\t4.0\nextra\t1.0\t2.0\n"
        )
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=str(trait_file), output="plot.png", json=False
        )
        svc = ContMap(args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_single_trait_data(str(trait_file), tree_tips)

        assert (
            "Line 5 in trait file has 3 columns; expected 2."
            in exc_info.value.messages
        )

    def test_non_numeric_trait_error(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("raccoon\t1.0\nbear\tbad\nsea_lion\t3.0\ndog\t4.0\n")
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=str(trait_file), output="plot.png", json=False
        )
        svc = ContMap(args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_single_trait_data(str(trait_file), tree_tips)

        assert (
            "Non-numeric trait value 'bad' for taxon 'bear' on line 2."
            in exc_info.value.messages
        )


class TestRun:
    def test_run_all_tips_present_uses_read_only_tree_without_copy_or_prune(
        self, mocker
    ):
        args = Namespace(
            tree="/some/path/to/file.tre",
            trait_data="/some/path/to/traits.tsv",
            output="plot.png",
            json=False,
        )
        svc = ContMap(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,c:2);"), "newick")
        mocked_read = mocker.patch.object(
            ContMap,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(
            ContMap,
            "read_tree_file",
            side_effect=AssertionError("run should not copy cached trees"),
        )
        mocked_validate = mocker.patch.object(ContMap, "validate_tree")
        mocker.patch.object(
            ContMap,
            "get_tip_names_from_tree",
            return_value=["a", "b", "c"],
        )
        mocker.patch.object(
            ContMap,
            "_parse_single_trait_data",
            return_value={"a": 1.0, "b": 2.0, "c": 3.0},
        )
        mocker.patch.object(
            ContMap,
            "_fast_copy",
            side_effect=AssertionError("all-shared tree should not be copied"),
        )
        mocker.patch.object(
            ContMap,
            "prune_tree_using_taxa_list",
            side_effect=AssertionError("all-shared tree should not be pruned"),
        )
        mocker.patch.object(ContMap, "_label_internal_nodes", return_value={})
        mocked_anc = mocker.patch.object(
            ContMap,
            "_fast_anc",
            return_value=({}, 1.0),
        )
        mocker.patch.object(ContMap, "_plot_contmap")
        mocker.patch("builtins.print")

        svc.run()

        mocked_read.assert_called_once_with()
        mocked_validate.assert_called_once_with(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="contMap visualization",
        )
        assert mocked_anc.call_args.args[0] is tree

    def test_run_missing_trait_taxa_copies_before_pruning(self, monkeypatch):
        args = Namespace(
            tree="/some/path/to/file.tre",
            trait_data="/some/path/to/traits.tsv",
            output="plot.png",
            json=False,
        )
        svc = ContMap(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")
        original_fast_copy = svc._fast_copy
        copied_trees = []
        captured = {}

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        def fast_anc(tree_arg, *_args):
            captured["tree"] = tree_arg
            return {}, 1.0

        monkeypatch.setattr(svc, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(svc, "_fast_copy", copy_spy)
        monkeypatch.setattr(
            svc,
            "_parse_single_trait_data",
            lambda *_args: {"a": 1.0, "b": 2.0, "c": 3.0},
        )
        monkeypatch.setattr(svc, "_label_internal_nodes", lambda *_args: {})
        monkeypatch.setattr(svc, "_fast_anc", fast_anc)
        monkeypatch.setattr(svc, "_plot_contmap", lambda *_args: None)
        monkeypatch.setattr("builtins.print", lambda *_args, **_kwargs: None)

        svc.run()

        assert len(copied_trees) == 1
        assert captured["tree"] is copied_trees[0]
        assert {tip.name for tip in tree.get_terminals()} == {"a", "b", "c", "d"}
        assert {tip.name for tip in copied_trees[0].get_terminals()} == {
            "a",
            "b",
            "c",
        }

    def test_run_ladderize_copies_before_mutating_tree(self, monkeypatch):
        args = Namespace(
            tree="/some/path/to/file.tre",
            trait_data="/some/path/to/traits.tsv",
            output="plot.png",
            json=False,
            ladderize=True,
        )
        svc = ContMap(args)
        tree = Phylo.read(StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")
        original_fast_copy = svc._fast_copy
        copied_trees = []
        captured = {}

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        def fast_anc(tree_arg, *_args):
            captured["tree"] = tree_arg
            return {}, 1.0

        monkeypatch.setattr(svc, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(svc, "_fast_copy", copy_spy)
        monkeypatch.setattr(
            svc,
            "_parse_single_trait_data",
            lambda *_args: {"a": 1.0, "b": 2.0, "c": 3.0, "d": 4.0},
        )
        monkeypatch.setattr(svc, "_label_internal_nodes", lambda *_args: {})
        monkeypatch.setattr(svc, "_fast_anc", fast_anc)
        monkeypatch.setattr(svc, "_plot_contmap", lambda *_args: None)
        monkeypatch.setattr("builtins.print", lambda *_args, **_kwargs: None)

        svc.run()

        assert len(copied_trees) == 1
        assert captured["tree"] is copied_trees[0]
        assert {tip.name for tip in tree.get_terminals()} == {"a", "b", "c", "d"}

    def test_print_text_output_batches_summary(self, mocker):
        svc = ContMap.__new__(ContMap)
        svc.output_path = "contmap.png"
        printed = mocker.patch("builtins.print")

        svc._print_text_output(12, 1.234567)

        printed.assert_called_once_with(
            "Continuous Trait Map (contMap)\n"
            "\nNumber of tips: 12\n"
            "Sigma-squared (BM rate): 1.2346\n"
            "Saved contMap plot: contmap.png"
        )

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

    def test_run_reuses_initial_tip_names_for_prune_setup(self, mocker):
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
            spy = mocker.spy(svc, "get_tip_names_from_tree")
            mocker.patch.object(ContMap, "_plot_contmap")
            mocker.patch("phykit.services.tree.cont_map.print_json")

            svc.run()

            assert spy.call_count == 1
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)


class TestFastAnc:
    def test_iter_postorder_matches_biopython_order(self):
        svc = ContMap.__new__(ContMap)
        tree = Phylo.read(
            StringIO("((A:1,B:2,C:3):4,(D:5,E:6):7,F:8):0;"),
            "newick",
        )

        direct = list(svc._iter_postorder(tree.root))
        reference = list(tree.find_clades(order="postorder"))

        assert [id(clade) for clade in direct] == [
            id(clade) for clade in reference
        ]

    def test_fast_anc_uses_direct_tree_traversal(self, monkeypatch):
        svc = ContMap.__new__(ContMap)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        x = np.array([1.0, 2.0, 3.0, 4.0])

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("optimized path should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        labels = svc._label_internal_nodes(tree)
        estimates, sigma2 = svc._fast_anc(tree, x, ordered_names, labels)

        assert set(estimates) == {"N1", "N2", "N3"}
        assert estimates["N1"] == pytest.approx(2.5)
        assert sigma2 > 0

    def test_build_parent_map_handles_mixed_child_counts(self, monkeypatch):
        svc = ContMap.__new__(ContMap)
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_iter_preorder(*args, **kwargs):
            raise AssertionError("parent map should build directly")

        monkeypatch.setattr(svc, "_iter_preorder", fail_iter_preorder)

        parent_map = svc._build_parent_map(tree)

        terminal, binary, trifurcating = tree.root.clades
        assert id(tree.root) not in parent_map
        assert parent_map[id(terminal)] is tree.root
        assert parent_map[id(binary)] is tree.root
        assert parent_map[id(trifurcating)] is tree.root
        assert all(parent_map[id(child)] is binary for child in binary.clades)
        assert all(
            parent_map[id(child)] is trifurcating for child in trifurcating.clades
        )

    def test_fast_anc_builds_labels_from_direct_child_lists(self, monkeypatch):
        svc = ContMap.__new__(ContMap)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        x = np.array([1.0, 2.0, 3.0, 4.0])
        labels = svc._label_internal_nodes(tree)

        def fail_is_terminal(self):
            raise AssertionError("optimized path should use clade child lists")

        monkeypatch.setattr(Clade, "is_terminal", fail_is_terminal)

        estimates, sigma2 = svc._fast_anc(tree, x, ordered_names, labels)

        assert set(estimates) == {"N1", "N2", "N3"}
        assert sigma2 > 0

    def test_sigma2_from_contrasts_handles_polytomy(self):
        svc = ContMap.__new__(ContMap)
        tree = Phylo.read(StringIO("(A:1,B:1,C:1);"), "newick")
        values = {
            tip.name: float(index + 1)
            for index, tip in enumerate(tree.root.clades)
        }
        node_estimates = {id(tip): values[tip.name] for tip in tree.root.clades}
        node_variances = {id(tip): 0.0 for tip in tree.root.clades}

        sigma2 = svc._compute_sigma2_from_contrasts(
            tree,
            np.array([1.0, 2.0, 3.0]),
            node_estimates,
            node_variances,
            ["A", "B", "C"],
        )

        expected_first = ((1.0 - 2.0) ** 2) / 2.0
        combined_est = (1.0 + 2.0) / 2.0
        combined_var = 0.5
        expected_second = ((3.0 - combined_est) ** 2) / (1.0 + combined_var)
        assert sigma2 == pytest.approx((expected_first + expected_second) / 2.0)


class TestPlotContMap:
    def test_value_range_scans_values_once(self):
        class SinglePassValues:
            def __init__(self, values):
                self.values = values
                self.iterations = 0

            def __iter__(self):
                self.iterations += 1
                if self.iterations > 1:
                    raise AssertionError("values should be scanned once")
                return iter(self.values)

        values = SinglePassValues([3.0, 1.0, 4.0, 2.0])

        assert _value_range(values) == (1.0, 4.0)
        assert values.iterations == 1
        assert _value_range([]) is None
        assert _value_range([2.0, 2.0]) == (2.0, 2.0)

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_contmap_uses_direct_tree_traversal(self, monkeypatch, circular):
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
                circular=circular,
            )
            svc = ContMap(args)
            tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
            ordered_names = ["A", "B", "C", "D"]
            trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
            x = np.array([trait_values[name] for name in ordered_names])
            node_labels = svc._label_internal_nodes(tree)
            node_estimates, _sigma2 = svc._fast_anc(
                tree, x, ordered_names, node_labels
            )

            def fail_traversal(*args, **kwargs):
                raise AssertionError("plot setup should use direct traversal")

            monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
            monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

            svc._plot_contmap(
                tree, node_estimates, node_labels,
                trait_values, "trait", tmppath,
            )

            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_prepare_contmap_plot_data_single_pass(self, monkeypatch):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            output="unused.png",
            json=False,
        )
        svc = ContMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        node_labels = svc._label_internal_nodes(tree)
        node_estimates = {
            label: float(index)
            for index, label in enumerate(node_labels.values(), start=1)
        }
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}

        def fail_parent_map(*_args, **_kwargs):
            raise AssertionError("plot setup should build parent map once")

        monkeypatch.setattr(svc, "_build_parent_map", fail_parent_map)

        parent_map, preorder, all_estimates, tips = (
            svc._prepare_contmap_plot_data(
                tree, node_estimates, node_labels, trait_values
            )
        )

        assert len(parent_map) == len(preorder) - 1
        assert [tip.name for tip in tips] == ["A", "B", "C", "D"]
        assert all_estimates[id(tips[0])] == pytest.approx(1.0)
        assert all_estimates[id(tree.root)] == pytest.approx(
            node_estimates[node_labels[id(tree.root)]]
        )

    def test_contmap_scalar_arcs_binary_children_use_indexed_path(self):
        class IndexedOnlyList(list):
            def __iter__(self):
                raise AssertionError("binary arc setup should index children")

        root = Clade(name="root")
        left = Clade(name="left")
        right = Clade(name="right")
        root.clades = IndexedOnlyList([left, right])
        left.clades = []
        right.clades = []

        coords = {
            id(root): {"radius": 5.0, "angle": 0.0},
            id(left): {"radius": 6.0, "angle": 2.0},
            id(right): {"radius": 6.0, "angle": 1.0},
        }
        all_estimates = {id(root): 3.0}

        arcs = ContMap._contmap_scalar_arcs(
            [root, left, right],
            coords,
            all_estimates,
        )

        assert arcs == [(0, 0, 5.0, 1.0, 2.0, 3.0)]

    def test_iter_preorder_preserves_order_without_reversed(self):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("preorder traversal should push children directly")

        svc = ContMap.__new__(ContMap)
        root = Clade(name="root")
        left = Clade(name="left")
        middle = Clade(name="middle")
        right = Clade(name="right")
        left_a = Clade(name="left_a")
        left_b = Clade(name="left_b")
        middle_a = Clade(name="middle_a")
        left.clades = NoReversedList([left_a, left_b])
        middle.clades = NoReversedList([middle_a])
        right.clades = NoReversedList()
        root.clades = NoReversedList([left, middle, right])

        preorder_names = [clade.name for clade in svc._iter_preorder(root)]

        assert preorder_names == [
            "root",
            "left",
            "left_a",
            "left_b",
            "middle",
            "middle_a",
            "right",
        ]

    def test_plot_contmap_reuses_preorder_for_node_positions(
        self, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import phykit.helpers.plot_config as plot_config

        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            output=str(tmp_path / "contmap_preorder_positions.png"),
            json=False,
            circular=False,
        )
        svc = ContMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        x = np.array([trait_values[name] for name in ordered_names])
        node_labels = svc._label_internal_nodes(tree)
        node_estimates, _sigma2 = svc._fast_anc(
            tree, x, ordered_names, node_labels
        )
        expected_preorder_ids = [
            id(clade)
            for clade in svc._prepare_contmap_plot_data(
                tree, node_estimates, node_labels, trait_values
            )[1]
        ]
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

        svc._plot_contmap(
            tree, node_estimates, node_labels, trait_values, "trait", args.output
        )

        assert calls == [expected_preorder_ids]
        assert os.path.exists(args.output)

    def test_plot_contmap_reuses_clade_lists_for_circular_coords(
        self, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import phykit.helpers.circular_layout as circular_layout

        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            output=str(tmp_path / "contmap_circular_precomputed_coords.png"),
            json=False,
            circular=True,
            no_title=True,
            ylabel_fontsize=0,
        )
        svc = ContMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        x = np.array([trait_values[name] for name in ordered_names])
        node_labels = svc._label_internal_nodes(tree)
        node_estimates, _sigma2 = svc._fast_anc(
            tree, x, ordered_names, node_labels
        )
        _, preorder_clades, _, tips = svc._prepare_contmap_plot_data(
            tree, node_estimates, node_labels, trait_values
        )
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

        svc._plot_contmap(
            tree, node_estimates, node_labels, trait_values, "trait", args.output
        )

        assert calls == [(expected_preorder_ids, expected_terminal_ids)]
        assert os.path.exists(args.output)

    def test_rectangular_plot_batches_gradient_branches(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes

        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            output=str(tmp_path / "contmap_batched.png"),
            json=False,
            circular=False,
        )
        svc = ContMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        x = np.array([trait_values[name] for name in ordered_names])
        node_labels = svc._label_internal_nodes(tree)
        node_estimates, _sigma2 = svc._fast_anc(
            tree, x, ordered_names, node_labels
        )

        def fail_plot(*args, **kwargs):
            raise AssertionError("rectangular contMap branches should be batched")

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)

        svc._plot_contmap(
            tree, node_estimates, node_labels,
            trait_values, "trait", args.output,
        )

        assert os.path.exists(args.output)

    def test_rectangular_plot_uses_collection_scalar_colormap(
        self, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            output=str(tmp_path / "contmap_scalar_colormap.png"),
            json=False,
            circular=False,
        )
        svc = ContMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        x = np.array([trait_values[name] for name in ordered_names])
        node_labels = svc._label_internal_nodes(tree)
        node_estimates, _sigma2 = svc._fast_anc(
            tree, x, ordered_names, node_labels
        )
        collection_arrays = []
        original_add_collection = matplotlib.axes.Axes.add_collection

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                array = collection.get_array()
                if array is not None:
                    collection_arrays.append(array)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "add_collection", capture_collection)

        svc._plot_contmap(
            tree,
            node_estimates,
            node_labels,
            trait_values,
            "trait",
            args.output,
        )

        array_lengths = {len(array) for array in collection_arrays}
        assert 300 in array_lengths
        assert 6 in array_lengths
        assert os.path.exists(args.output)

    def test_circular_plot_batches_gradient_branches(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            output=str(tmp_path / "contmap_circular_batched.png"),
            json=False,
            circular=True,
            no_title=True,
            ylabel_fontsize=0,
        )
        svc = ContMap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        x = np.array([trait_values[name] for name in ordered_names])
        node_labels = svc._label_internal_nodes(tree)
        node_estimates, _sigma2 = svc._fast_anc(
            tree, x, ordered_names, node_labels
        )
        line_collections = []
        original_add_collection = matplotlib.axes.Axes.add_collection

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "add_collection", capture_collection)

        svc._plot_contmap(
            tree, node_estimates, node_labels,
            trait_values, "trait", args.output,
        )

        segment_counts = {len(collection.get_segments()) for collection in line_collections}
        assert 180 in segment_counts
        assert 3 in segment_counts
        assert os.path.exists(args.output)
