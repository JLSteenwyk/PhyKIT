import json
import os
import subprocess
import sys
import tempfile
from unittest.mock import call

import pytest
from argparse import Namespace
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
from io import StringIO

import phykit.services.tree.cophylo as module
from phykit.services.tree.cophylo import Cophylo
from phykit.errors import PhykitUserError


def test_module_import_does_not_import_numpy_or_biophylo():
    code = """
import sys
import phykit.services.tree.cophylo as module
assert callable(module.print_json)
assert callable(module.compute_circular_coords)
assert callable(module.draw_circular_branches)
assert callable(module.draw_circular_tip_labels)
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.circular_layout" not in sys.modules
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "phykit.helpers.color_annotations" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


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


class TestMappingParsing:
    def test_mapping_file_skips_comments_and_blanks(self, tmp_path):
        svc = Cophylo.__new__(Cophylo)
        mapping_file = tmp_path / "mapping.tsv"
        mapping_file.write_text(
            "# comment\n"
            "\n"
            "A\tA_prime\n"
            "B\tB_prime\n"
            "C\tC_prime\n"
        )

        mapping = svc._parse_mapping_file(str(mapping_file))

        assert mapping == {
            "A": "A_prime",
            "B": "B_prime",
            "C": "C_prime",
        }

    def test_mapping_file_skips_whitespace_prefixed_comments(self, tmp_path):
        svc = Cophylo.__new__(Cophylo)
        mapping_file = tmp_path / "mapping.tsv"
        mapping_file.write_text(
            "   # comment\n"
            "\n"
            "A\tA_prime\n"
            "B\tB_prime\n"
        )

        mapping = svc._parse_mapping_file(str(mapping_file))

        assert mapping == {
            "A": "A_prime",
            "B": "B_prime",
        }

    def test_mapping_file_rejects_wrong_column_count(self, tmp_path):
        svc = Cophylo.__new__(Cophylo)
        mapping_file = tmp_path / "mapping.tsv"
        mapping_file.write_text("A\tA_prime\nB\tB_prime\textra\n")

        with pytest.raises(PhykitUserError):
            svc._parse_mapping_file(str(mapping_file))

    def test_mapping_file_reports_missing_path(self, tmp_path):
        svc = Cophylo.__new__(Cophylo)
        missing = tmp_path / "missing.tsv"

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_mapping_file(str(missing))

        assert "no such file or directory" in exc_info.value.messages[0]


class TestRun:
    def test_run_reports_second_tree_read_failure(self, mocker):
        svc = Cophylo(
            Namespace(
                tree1=TREE1,
                tree2="missing.tre",
                output="out.png",
                mapping=None,
                json=False,
            )
        )
        mocker.patch.object(
            svc,
            "_read_prepared_tree",
            side_effect=[object(), ValueError("invalid tree")],
        )

        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()

        assert "missing.tre" in exc_info.value.messages[0]

    def test_run_filters_mapping_and_ladderizes_both_trees(self, mocker):
        svc = Cophylo(
            Namespace(
                tree1=TREE1,
                tree2=TREE2,
                output="out.png",
                mapping="mapping.tsv",
                json=False,
                ladderize=True,
            )
        )
        tree1 = mocker.Mock()
        tree2 = mocker.Mock()
        mocker.patch.object(
            svc,
            "_read_prepared_tree",
            side_effect=[tree1, tree2],
        )
        mocker.patch.object(
            svc,
            "get_tip_names_from_tree",
            side_effect=[["A", "B"], ["X", "Y"]],
        )
        mocker.patch.object(
            svc,
            "_parse_mapping_file",
            return_value={"A": "X", "B": "Y", "missing": "Y"},
        )
        optimize = mocker.patch.object(
            svc,
            "_optimize_tip_order",
            return_value=({"A": 0, "B": 1}, {"X": 0, "Y": 1}),
        )
        plot = mocker.patch.object(svc, "_plot_cophylo")
        mocker.patch.object(svc, "_print_text_output")

        svc.run()

        tree1.ladderize.assert_called_once_with()
        tree2.ladderize.assert_called_once_with()
        optimize.assert_called_once_with(
            tree1,
            tree2,
            {"A": "X", "B": "Y"},
            mutate_tree2=True,
        )
        plot.assert_called_once()

    def test_run_rejects_mapping_with_too_few_valid_taxa(self, mocker):
        svc = Cophylo(
            Namespace(
                tree1=TREE1,
                tree2=TREE2,
                output="out.png",
                mapping="mapping.tsv",
                json=False,
            )
        )
        mocker.patch.object(
            svc,
            "_read_prepared_tree",
            side_effect=[object(), object()],
        )
        mocker.patch.object(
            svc,
            "get_tip_names_from_tree",
            side_effect=[["A", "B"], ["X", "Y"]],
        )
        mocker.patch.object(
            svc,
            "_parse_mapping_file",
            return_value={"A": "X", "B": "missing"},
        )

        with pytest.raises(PhykitUserError) as exc_info:
            svc.run()

        assert exc_info.value.messages[0].startswith("Only 1 valid mapped taxa")

    def test_run_prepares_both_trees_without_forcing_rectangular_copies(self, mocker):
        args = Namespace(
            tree1=TREE1,
            tree2=TREE2,
            output="out.png",
            mapping=None,
            json=False,
        )
        svc = Cophylo(args)
        tree1 = object()
        tree2 = object()
        tips = ["A", "B"]

        prepare_tree = mocker.patch.object(
            svc,
            "_read_prepared_tree",
            side_effect=[tree1, tree2],
        )
        mocker.patch.object(
            svc, "get_tip_names_from_tree", side_effect=[tips, tips]
        )
        mocker.patch.object(
            svc,
            "_optimize_tip_order",
            return_value=({"A": 0, "B": 1}, {"A": 0, "B": 1}),
        )
        mocker.patch.object(svc, "_plot_cophylo")
        mocker.patch("builtins.print")

        svc.run()

        assert prepare_tree.call_args_list == [
            call(TREE1, "tree_file_path", "tree1", force_copy=False),
            call(TREE2, "tree2_file_path", "tree2", force_copy=False),
        ]

    def test_run_default_mapping_scans_tip_names_once(self, mocker):
        class CountingList(list):
            def __init__(self, values):
                super().__init__(values)
                self.iterations = 0

            def __iter__(self):
                self.iterations += 1
                return super().__iter__()

        args = Namespace(
            tree1=TREE1,
            tree2=TREE2,
            output="out.png",
            mapping=None,
            json=False,
        )
        svc = Cophylo(args)
        tree1 = object()
        tree2 = object()
        tips1 = CountingList(["A", "B", "C"])
        tips2 = CountingList(["A", "B", "C"])

        mocker.patch.object(
            svc,
            "_read_prepared_tree",
            side_effect=[tree1, tree2],
        )
        mocker.patch.object(
            svc,
            "get_tip_names_from_tree",
            side_effect=[tips1, tips2],
        )
        mocker.patch.object(
            svc,
            "_optimize_tip_order",
            return_value=({"A": 0, "B": 1, "C": 2}, {"A": 0, "B": 1, "C": 2}),
        )
        mocker.patch.object(svc, "_plot_cophylo")
        mocker.patch("builtins.print")

        svc.run()

        assert tips1.iterations == 1
        assert tips2.iterations == 1

    def test_run_rectangular_branch_length_trees_skip_fast_copy(
        self, mocker, tmp_path
    ):
        tree1_path = tmp_path / "tree1.tre"
        tree2_path = tmp_path / "tree2.tre"
        tree1_path.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
        tree2_path.write_text("((C:1,D:1):1,(A:1,B:1):1);\n")
        args = Namespace(
            tree1=str(tree1_path),
            tree2=str(tree2_path),
            output=str(tmp_path / "out.png"),
            mapping=None,
            json=False,
        )
        svc = Cophylo(args)
        mocker.patch.object(
            svc,
            "_fast_copy",
            side_effect=AssertionError(
                "rectangular cophylo should use read-only trees when possible"
            ),
        )
        mocker.patch.object(svc, "_plot_cophylo")
        mocker.patch("builtins.print")

        svc.run()

    def test_read_prepared_tree_fills_missing_lengths_on_copy(self, tmp_path):
        tree_path = tmp_path / "tree.tre"
        tree_path.write_text("((A,B),(C,D));\n")
        svc = Cophylo(
            Namespace(
                tree1=str(tree_path),
                tree2=str(tree_path),
                output=str(tmp_path / "out.png"),
                mapping=None,
                json=False,
            )
        )

        prepared = svc._read_prepared_tree(
            str(tree_path),
            "tree_file_path",
            "tree1",
        )
        cached = svc._read_tree_with_error(
            str(tree_path),
            "tree_file_path",
            copy_tree=False,
        )

        assert all(
            clade.branch_length == 1.0
            for clade in prepared.root.clades
        )
        assert all(
            clade.branch_length is None
            for clade in cached.root.clades
        )

    def test_nonmutating_rotated_tip_order_matches_mutating_rotation(self):
        svc = Cophylo.__new__(Cophylo)
        tree_mutating = Phylo.read(StringIO("((C,D),(A,B));"), "newick")
        tree_readonly = Phylo.read(StringIO("((C,D),(A,B));"), "newick")
        original_child_order = [
            child.clades[0].name for child in tree_readonly.root.clades
        ]
        target_order = {"A": 0, "B": 1, "C": 2, "D": 3}
        reverse_mapping = {name: name for name in target_order}

        svc._rotate_tree(tree_mutating, target_order, reverse_mapping)
        mutating_order = svc._get_tip_order(tree_mutating)
        readonly_order = svc._get_rotated_tip_order(
            tree_readonly,
            target_order,
            reverse_mapping,
        )

        assert readonly_order == mutating_order
        assert [
            child.clades[0].name for child in tree_readonly.root.clades
        ] == original_child_order

    def test_validate_tree_uses_direct_standard_tree_traversal(self, monkeypatch):
        svc = Cophylo.__new__(Cophylo)
        tree = Phylo.read(StringIO("((A,B),(C,D));"), "newick")

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("standard tree validation should not collect terminals")

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("standard tree validation should not use find_clades")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(tree, "find_clades", fail_find_clades)

        svc._validate_tree(tree, "tree")

        assert tree.root.branch_length is None
        assert all(
            clade.branch_length == 1.0
            for clade in tree.root.clades
        )

    def test_validate_standard_tree_handles_mixed_child_counts(self, monkeypatch):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("standard tree validation should not collect terminals")

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("standard tree validation should not use find_clades")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(tree, "find_clades", fail_find_clades)

        assert Cophylo._validate_standard_tree(tree) == 6
        assert tree.root.branch_length is None
        assert tree.root.clades[1].clades[1].branch_length == 1.0

    def test_validate_tree_supports_generic_tree_protocol(self):
        class GenericClade:
            def __init__(self, name, branch_length=None, children=()):
                self.name = name
                self.branch_length = branch_length
                self.clades = list(children)

        tips = [GenericClade("A"), GenericClade("B", branch_length=2.0)]

        class GenericTree:
            root = object()

            def get_terminals(self):
                return tips

            def find_clades(self, order=None):
                return iter(tips)

        tree = GenericTree()
        svc = Cophylo.__new__(Cophylo)

        assert svc._validate_tree(tree, "generic tree") is True
        assert tips[0].branch_length == 1.0
        assert tips[1].branch_length == 2.0
        assert Cophylo._validate_standard_tree(tree) is None

    def test_validate_tree_rejects_generic_single_tip_tree(self):
        tip = Namespace(name="A", branch_length=1.0, clades=[])

        class GenericTree:
            root = object()

            def get_terminals(self):
                return [tip]

            def find_clades(self, order=None):
                return iter([tip])

        with pytest.raises(PhykitUserError) as exc_info:
            Cophylo.__new__(Cophylo)._validate_tree(GenericTree(), "tree")

        assert "at least 2 tips" in exc_info.value.messages[0]

    def test_standard_tree_validation_handles_malformed_child(self):
        malformed = Namespace(root=Namespace(clades=[object()]))

        assert Cophylo._validate_standard_tree(malformed) is None

    def test_traversal_helpers_fall_back_for_generic_tree(self):
        tip_a = Namespace(name="A", clades=[])
        tip_b = Namespace(name="B", clades=[])
        internal = Namespace(name=None, clades=[tip_a, tip_b])

        class GenericTree:
            root = object()

            def find_clades(self, order=None):
                if order == "postorder":
                    return iter([tip_a, tip_b, internal])
                return iter([internal, tip_a, tip_b])

            def get_terminals(self):
                return [tip_a, tip_b]

        tree = GenericTree()
        svc = Cophylo.__new__(Cophylo)

        assert svc._build_parent_map(tree) == {
            id(tip_a): internal,
            id(tip_b): internal,
        }
        assert Cophylo._preorder_clades(tree) == [internal, tip_a, tip_b]
        assert Cophylo._postorder_clades(tree) == [tip_a, tip_b, internal]
        assert Cophylo._terminal_clades(tree) == [tip_a, tip_b]

    @pytest.mark.parametrize(
        ("helper", "expected"),
        [
            (Cophylo._preorder_clades, "preorder"),
            (Cophylo._postorder_clades, "postorder"),
            (Cophylo._terminal_clades, "terminals"),
        ],
    )
    def test_traversal_helpers_recover_from_malformed_child(self, helper, expected):
        root = Namespace(clades=[object()])

        class MalformedTree:
            def __init__(self):
                self.root = root

            def find_clades(self, order=None):
                return [expected]

            def get_terminals(self):
                return [expected]

        assert helper(MalformedTree()) == [expected]

    def test_get_tip_order_uses_fast_terminal_names_for_parsed_tree(
        self, monkeypatch
    ):
        svc = Cophylo.__new__(Cophylo)
        tree = Phylo.read(StringIO("((A,B),(C,D));"), "newick")

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("parsed trees should use fast terminal names")

        monkeypatch.setattr(tree, "get_terminals", fail_get_terminals)

        assert svc._get_tip_order(tree) == ["A", "B", "C", "D"]

    def test_print_text_output_batches_summary(self, mocker):
        svc = Cophylo.__new__(Cophylo)
        svc.output_path = "out.png"
        printed = mocker.patch("builtins.print")

        svc._print_text_output(12, 10, 8)

        printed.assert_called_once_with(
            "Cophylogenetic Plot (Tanglegram)\n"
            "\nTree 1 tips: 12\n"
            "Tree 2 tips: 10\n"
            "Matched taxa: 8\n"
            "Saved cophylo plot: out.png"
        )

    def test_draw_phylogram_uses_direct_standard_tree_traversal(self, monkeypatch):
        svc = Cophylo.__new__(Cophylo)
        svc.plot_config = Namespace(cladogram=False, color_file=None)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tip_order = {"A": 0, "B": 1, "C": 2, "D": 3}

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("standard phylogram drawing should not collect terminals")

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("standard phylogram drawing should not use find_clades")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        class Ax:
            spines = {
                key: type("Spine", (), {"set_visible": lambda self, value: None})()
                for key in ("top", "right", "left")
            }

            def plot(self, *args, **kwargs):
                return None

            def text(self, *args, **kwargs):
                return None

            def set_ylim(self, *args, **kwargs):
                return None

            def set_yticks(self, *args, **kwargs):
                return None

            def set_xlabel(self, *args, **kwargs):
                return None

        svc._draw_phylogram(Ax(), tree, tip_order, direction="right")

    def test_draw_phylogram_batches_real_matplotlib_branches(self, monkeypatch):
        pytest.importorskip("matplotlib")
        import matplotlib
        import matplotlib.axes

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        svc = Cophylo.__new__(Cophylo)
        svc.plot_config = Namespace(cladogram=False, color_file=None)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tip_order = {"A": 0, "B": 1, "C": 2, "D": 3}
        fig, ax = plt.subplots()

        def fail_plot(*args, **kwargs):
            raise AssertionError("real phylogram branches should be batched")

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)

        try:
            svc._draw_phylogram(ax, tree, tip_order, direction="right")
            assert len(ax.collections) == 2
        finally:
            plt.close(fig)

    def test_plot_cophylo_rect_batches_association_lines(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib
        import matplotlib.axes

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection

        args = Namespace(
            tree1=TREE1,
            tree2=TREE2,
            output=str(tmp_path / "cophylo.png"),
            mapping=None,
            json=False,
            circular=False,
            color_file=None,
            no_title=True,
            ylabel_fontsize=0,
        )
        svc = Cophylo(args)
        tree1 = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree2 = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        order = {"A": 0, "B": 1, "C": 2, "D": 3}
        mapping = {name: name for name in order}
        output_path = tmp_path / "cophylo_associations.png"

        original_plot = matplotlib.axes.Axes.plot
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def fail_association_plot(self, *args, **kwargs):
            if kwargs.get("color") == "gray" and kwargs.get("zorder") == 1:
                raise AssertionError("association lines should use LineCollection")
            return original_plot(self, *args, **kwargs)

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_association_plot)
        monkeypatch.setattr(matplotlib.axes.Axes, "add_collection", capture_collection)

        try:
            svc._plot_cophylo(tree1, tree2, mapping, order, order, str(output_path))
            assert any(
                len(collection.get_segments()) == len(mapping)
                for collection in line_collections
            )
            assert output_path.exists()
        finally:
            plt.close("all")

    def test_plot_cophylo_uses_tip_order_lengths_for_size(
        self, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        args = Namespace(
            tree1=TREE1,
            tree2=TREE2,
            output=str(tmp_path / "cophylo.png"),
            mapping=None,
            json=False,
            circular=False,
            color_file=None,
            no_title=True,
            ylabel_fontsize=0,
        )
        svc = Cophylo(args)
        tree1 = Phylo.read(StringIO("((A:1,B:1):1,C:1);"), "newick")
        tree2 = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree1_order = {"A": 0, "B": 1, "C": 2}
        tree2_order = {"A": 0, "B": 1, "C": 2, "D": 3}
        mapping = {"A": "A", "B": "B", "C": "C"}
        observed = {}

        def fail_terminal_scan(*_args, **_kwargs):
            raise AssertionError("plot dispatch should reuse tip-order lengths")

        def capture_rect(
            _tree1,
            _tree2,
            _mapping,
            _tree1_order,
            _tree2_order,
            _output_path,
            _config,
            n_max,
            _plt,
        ):
            observed["n_max"] = n_max

        monkeypatch.setattr(Cophylo, "_terminal_clades", fail_terminal_scan)
        monkeypatch.setattr(svc, "_plot_cophylo_rect", capture_rect)

        try:
            svc._plot_cophylo(
                tree1,
                tree2,
                mapping,
                tree1_order,
                tree2_order,
                str(tmp_path / "cophylo_dispatch.png"),
            )
        finally:
            plt.close("all")

        assert observed["n_max"] == 4

    def test_plot_cophylo_rect_skips_redundant_tight_layout(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib
        import matplotlib.figure

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        args = Namespace(
            tree1=TREE1,
            tree2=TREE2,
            output=str(tmp_path / "cophylo.png"),
            mapping=None,
            json=False,
            circular=False,
            color_file=None,
            no_title=True,
            ylabel_fontsize=0,
        )
        svc = Cophylo(args)
        tree1 = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree2 = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        order = {"A": 0, "B": 1, "C": 2, "D": 3}
        mapping = {name: name for name in order}
        output_path = tmp_path / "cophylo_no_tight_layout.png"

        def fail_tight_layout(self, *args, **kwargs):
            raise AssertionError("rectangular cophylo should rely on tight savefig")

        monkeypatch.setattr(
            matplotlib.figure.Figure,
            "tight_layout",
            fail_tight_layout,
        )

        try:
            svc._plot_cophylo(tree1, tree2, mapping, order, order, str(output_path))
        finally:
            plt.close("all")

        assert output_path.exists()

    def test_circular_plot_reuses_clade_lists_for_layout_helpers(
        self, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        args = Namespace(
            tree1=TREE1,
            tree2=TREE2,
            output=str(tmp_path / "cophylo.png"),
            mapping=None,
            json=False,
            circular=True,
            color_file=None,
            no_title=True,
            ylabel_fontsize=0,
        )
        svc = Cophylo(args)
        tree1 = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree2 = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        order = {"A": 0, "B": 1, "C": 2, "D": 3}
        mapping = {name: name for name in order}
        output_path = tmp_path / "cophylo_layout_lists.png"
        layout_calls = []
        original_circular_coords = module.compute_circular_coords

        def spy_circular_coords(*args, **kwargs):
            layout_calls.append(
                (
                    args[0],
                    kwargs.get("preorder_clades"),
                    kwargs.get("terminal_clades"),
                )
            )
            return original_circular_coords(*args, **kwargs)

        monkeypatch.setattr(module, "compute_circular_coords", spy_circular_coords)

        try:
            svc._plot_cophylo(tree1, tree2, mapping, order, order, str(output_path))
        finally:
            plt.close("all")

        assert output_path.exists()
        assert len(layout_calls) == 2
        for tree, preorder_clades, terminal_clades in layout_calls:
            expected_preorder = svc._preorder_clades(tree)
            expected_terminals = [
                clade for clade in expected_preorder if not clade.clades
            ]
            assert len(preorder_clades) == len(expected_preorder)
            assert len(terminal_clades) == len(expected_terminals)
            assert all(
                actual is expected
                for actual, expected in zip(preorder_clades, expected_preorder)
            )
            assert all(
                actual is expected
                for actual, expected in zip(terminal_clades, expected_terminals)
            )

    def test_circular_plot_batches_association_lines(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib
        import matplotlib.figure
        import matplotlib.patches
        from matplotlib.collections import LineCollection

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        args = Namespace(
            tree1=TREE1,
            tree2=TREE2,
            output=str(tmp_path / "cophylo.png"),
            mapping=None,
            json=False,
            circular=True,
            color_file=None,
            no_title=True,
            ylabel_fontsize=0,
        )
        svc = Cophylo(args)
        tree1 = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree2 = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        order = {"A": 0, "B": 1, "C": 2, "D": 3}
        mapping = {name: name for name in order}
        output_path = tmp_path / "cophylo_associations_circular.png"
        original_add_artist = matplotlib.figure.Figure.add_artist
        line_collections = []

        def fail_connection_patch(*_args, **_kwargs):
            raise AssertionError("circular association lines should be batched")

        def count_artist(self, artist, *args, **kwargs):
            if isinstance(artist, LineCollection):
                line_collections.append(artist)
            return original_add_artist(self, artist, *args, **kwargs)

        monkeypatch.setattr(
            matplotlib.patches,
            "ConnectionPatch",
            fail_connection_patch,
        )
        monkeypatch.setattr(matplotlib.figure.Figure, "add_artist", count_artist)

        try:
            svc._plot_cophylo(tree1, tree2, mapping, order, order, str(output_path))
        finally:
            plt.close("all")

        assert any(
            len(collection.get_segments()) == len(mapping)
            and collection.get_linewidths()[0] == pytest.approx(0.5)
            for collection in line_collections
        )
        assert output_path.exists()

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_batches_color_clade_overlay(self, monkeypatch, tmp_path, circular):
        pytest.importorskip("matplotlib")
        import matplotlib
        import matplotlib.axes
        import phykit.helpers.color_annotations as color_annotations

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection

        color_file = tmp_path / "colors.tsv"
        color_file.write_text("A,B\tclade\t#ff0000\tAB\n")
        args = Namespace(
            tree1=TREE1,
            tree2=TREE2,
            output=str(tmp_path / "cophylo.png"),
            mapping=None,
            json=False,
            circular=circular,
            color_file=str(color_file),
            no_title=True,
            ylabel_fontsize=0,
        )
        svc = Cophylo(args)
        tree1 = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree2 = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        order = {"A": 0, "B": 1, "C": 2, "D": 3}
        mapping = {name: name for name in order}

        original_plot = matplotlib.axes.Axes.plot
        original_add_collection = matplotlib.axes.Axes.add_collection
        original_parse_color_file = color_annotations.parse_color_file
        line_collections = []
        parse_calls = 0

        def fail_overlay_plot(self, *args, **kwargs):
            if kwargs.get("zorder") == 2:
                raise AssertionError("clade overlay branches should use LineCollection")
            return original_plot(self, *args, **kwargs)

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        def counting_parse_color_file(*args, **kwargs):
            nonlocal parse_calls
            parse_calls += 1
            return original_parse_color_file(*args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_overlay_plot)
        monkeypatch.setattr(matplotlib.axes.Axes, "add_collection", capture_collection)
        monkeypatch.setattr(
            color_annotations, "parse_color_file", counting_parse_color_file
        )

        output_path = tmp_path / f"cophylo_color_{circular}.png"
        try:
            svc._plot_cophylo(tree1, tree2, mapping, order, order, str(output_path))
            assert len(line_collections) >= 4
            assert output_path.exists()
            assert parse_calls == 1
        finally:
            plt.close("all")

    def test_rotate_tree_matches_terminal_scan_reference(self):
        svc = Cophylo.__new__(Cophylo)
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);"),
            "newick",
        )
        target_order = {"a": 4, "b": 3, "c": 2, "d": 1, "e": 0}
        reverse_mapping = {
            "A": "a",
            "B": "b",
            "C": "c",
            "D": "d",
            "E": "e",
        }

        def reference_rotate(clade):
            if clade.is_terminal():
                return
            for child in clade.clades:
                reference_rotate(child)
            if len(clade.clades) != 2:
                return

            means = []
            for child in clade.clades:
                positions = [
                    target_order[reverse_mapping[tip.name]]
                    for tip in child.get_terminals()
                    if tip.name in reverse_mapping
                ]
                means.append(sum(positions) / len(positions) if positions else 0.0)

            if means[0] > means[1]:
                clade.clades[0], clade.clades[1] = clade.clades[1], clade.clades[0]

        expected = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);"),
            "newick",
        )
        reference_rotate(expected.root)

        svc._rotate_tree(tree, target_order, reverse_mapping)

        assert [tip.name for tip in tree.get_terminals()] == [
            tip.name for tip in expected.get_terminals()
        ]

    def test_rotate_tree_handles_unmapped_descendants(self):
        svc = Cophylo.__new__(Cophylo)
        tree = Phylo.read(
            StringIO("((A:1,X:1):1,(B:1,(C:1,Y:1):1):1);"),
            "newick",
        )
        target_order = {"a": 2, "b": 1, "c": 0}
        reverse_mapping = {"A": "a", "B": "b", "C": "c"}

        svc._rotate_tree(tree, target_order, reverse_mapping)

        assert [tip.name for tip in tree.get_terminals()] == [
            "C", "Y", "B", "X", "A",
        ]

    def test_rotate_tree_uses_direct_postorder_traversal(self, monkeypatch):
        svc = Cophylo.__new__(Cophylo)
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        target_order = {"A": 3, "B": 2, "C": 1, "D": 0}
        reverse_mapping = {name: name for name in target_order}

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree rotation should not use find_clades")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        svc._rotate_tree(tree, target_order, reverse_mapping)

        assert svc._get_tip_order(tree) == ["D", "C", "B", "A"]

    def test_postorder_clades_preserves_biophylo_order(self, monkeypatch):
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);"),
            "newick",
        )
        expected_ids = [
            id(clade) for clade in tree.find_clades(order="postorder")
        ]

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree helper should not use find_clades")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        observed_ids = [id(clade) for clade in Cophylo._postorder_clades(tree)]

        assert observed_ids == expected_ids

    def test_postorder_clades_does_not_use_generic_traversal(self, monkeypatch):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree helper should not use find_clades")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        assert Cophylo._postorder_clades(tree)

    def test_preorder_clades_preserves_biophylo_order(self, monkeypatch):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )
        expected_ids = [
            id(clade) for clade in tree.find_clades(order="preorder")
        ]

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree helper should not use find_clades")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        observed_ids = [id(clade) for clade in Cophylo._preorder_clades(tree)]

        assert observed_ids == expected_ids

    def test_terminal_clades_preserves_biophylo_order(self, monkeypatch):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )
        expected_ids = [id(clade) for clade in tree.get_terminals()]

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("standard tree helper should not use get_terminals")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        observed_ids = [id(clade) for clade in Cophylo._terminal_clades(tree)]

        assert observed_ids == expected_ids

    def test_assign_internal_y_positions_uses_binary_fast_path(self):
        from Bio.Phylo.BaseTree import Clade

        class IndexedOnlyList(list):
            def __iter__(self):
                raise AssertionError("binary y-position helper should not iterate")

        left = Clade(name="A")
        right = Clade(name="B")
        binary = Clade(clades=IndexedOnlyList([left, right]))
        singleton = Clade(name="C")
        unary = Clade(clades=[singleton])
        extra = Clade(name="D")
        root = Clade(clades=[binary, unary, extra])
        preorder_clades = [root, binary, left, right, unary, singleton, extra]
        node_y = {
            id(left): 0.0,
            id(right): 4.0,
            id(singleton): 8.0,
            id(extra): 12.0,
        }

        Cophylo._assign_internal_y_positions(preorder_clades, node_y)

        assert node_y[id(binary)] == 2.0
        assert node_y[id(unary)] == 8.0
        assert node_y[id(root)] == pytest.approx((2.0 + 8.0 + 12.0) / 3.0)

    def test_build_parent_map_handles_mixed_child_counts(self, monkeypatch):
        svc = Cophylo.__new__(Cophylo)
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("parent map should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

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
