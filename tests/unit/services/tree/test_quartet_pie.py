"""
Unit tests for quartet_pie (quartet concordance pie chart visualization).

Tests initialization, process_args, plot creation, and JSON output.
Ground truth gCF/gDF values validated against manual bipartition
matching computation on sample data.
"""
import pytest
import builtins
import subprocess
import sys
from argparse import Namespace
from io import StringIO
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

import phykit.services.tree.quartet_pie as quartet_pie_module
from phykit.services.tree.quartet_pie import QuartetPie


def test_module_import_does_not_import_numpy_matplotlib_or_biophylo():
    code = """
import sys
import phykit.services.tree.quartet_pie as module

assert callable(module.print_json)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "numpy" not in sys.modules
assert "matplotlib" not in sys.modules
assert "matplotlib.pyplot" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.quartet_utils" not in sys.modules
assert "phykit.helpers.circular_layout" not in sys.modules
assert "phykit.helpers.color_annotations" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def args(tmp_path):
    return Namespace(
        tree="tests/sample_files/tree_simple.tre",
        gene_trees="tests/sample_files/gene_trees_simple.nwk",
        output=str(tmp_path / "qpie.png"),
        annotate=False,
    )


@pytest.fixture(autouse=True)
def clear_quartet_pie_plot_cache():
    previous_cache = quartet_pie_module._QUARTET_PIE_PLOT_CACHE.copy()
    quartet_pie_module._QUARTET_PIE_PLOT_CACHE.clear()
    try:
        yield
    finally:
        quartet_pie_module._QUARTET_PIE_PLOT_CACHE.clear()
        quartet_pie_module._QUARTET_PIE_PLOT_CACHE.update(previous_cache)


class TestQuartetPieInit:
    def test_init_sets_fields(self, args):
        qp = QuartetPie(args)
        assert qp.tree_file_path == args.tree
        assert qp.gene_trees_path == args.gene_trees
        assert qp.output_path == args.output
        assert qp.annotate is False
        assert qp.json_output is False

    def test_process_args_no_gene_trees(self, tmp_path):
        args = Namespace(
            tree="t.tre", output=str(tmp_path / "out.png"),
        )
        qp = QuartetPie(args)
        assert qp.gene_trees_path is None


class TestQuartetPiePlot:
    def test_creates_png(self, args):
        qp = QuartetPie(args)
        qp.run()
        assert Path(args.output).exists()
        assert Path(args.output).stat().st_size > 0

    def test_plot_without_color_file_skips_color_annotation_import(
        self, args, monkeypatch
    ):
        original_import = builtins.__import__

        def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "phykit.helpers.color_annotations":
                raise AssertionError(
                    "color annotation helpers should only load for --color-file"
                )
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", guarded_import)

        qp = QuartetPie(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        proportions = {
            id(clade): (0.7, 0.2, 0.1, 7, 2, 1)
            for clade in qp._iter_preorder(tree.root)
            if clade.clades and clade != tree.root
        }

        qp._plot_quartet_pie(tree, proportions, args.output, {})

        assert Path(args.output).exists()

    def test_repeated_plot_reuses_cached_rendered_bytes(
        self, args, monkeypatch, tmp_path
    ):
        previous_cache = quartet_pie_module._QUARTET_PIE_PLOT_CACHE.copy()
        quartet_pie_module._QUARTET_PIE_PLOT_CACHE.clear()
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        first = QuartetPie(args)
        proportions = {
            id(clade): (0.7, 0.2, 0.1, 7, 2, 1)
            for clade in first._iter_preorder(tree.root)
            if clade.clades and clade != tree.root
        }

        first._plot_quartet_pie(tree, proportions, args.output, {})
        first_bytes = Path(args.output).read_bytes()

        second_args = Namespace(
            tree=args.tree,
            gene_trees=args.gene_trees,
            output=str(tmp_path / "qpie_cached.png"),
            annotate=False,
        )
        second = QuartetPie(second_args)
        second_tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        second_proportions = {
            id(clade): (0.7, 0.2, 0.1, 7, 2, 1)
            for clade in second._iter_preorder(second_tree.root)
            if clade.clades and clade != second_tree.root
        }
        original_import = builtins.__import__

        def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "matplotlib" or name.startswith("matplotlib."):
                raise AssertionError("cached plot should not import matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", guarded_import)
        try:
            second._plot_quartet_pie(
                second_tree, second_proportions, second_args.output, {}
            )
            assert Path(second_args.output).read_bytes() == first_bytes
        finally:
            quartet_pie_module._QUARTET_PIE_PLOT_CACHE.clear()
            quartet_pie_module._QUARTET_PIE_PLOT_CACHE.update(previous_cache)

    def test_creates_pdf(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees="tests/sample_files/gene_trees_simple.nwk",
            output=str(tmp_path / "qpie.pdf"),
            annotate=False,
        )
        qp = QuartetPie(args)
        qp.run()
        assert Path(args.output).exists()

    def test_annotate_flag(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees="tests/sample_files/gene_trees_simple.nwk",
            output=str(tmp_path / "qpie_annot.png"),
            annotate=True,
        )
        qp = QuartetPie(args)
        qp.run()
        assert Path(args.output).exists()

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_uses_direct_tree_traversal(self, monkeypatch, tmp_path, circular):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees=None,
            output=str(tmp_path / "qpie_direct.png"),
            annotate=False,
            branch_labels=True,
            circular=circular,
        )
        qp = QuartetPie(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clades = list(qp._iter_preorder(tree.root))
        proportions = {
            id(clade): (0.7, 0.2, 0.1, 7, 2, 1)
            for clade in clades
            if clade.clades and clade != tree.root
        }
        branch_info = {
            id(clade): {"f1": 7, "pp1": 0.9}
            for clade in clades
            if clade.clades and clade != tree.root
        }

        def fail_traversal(*args, **kwargs):
            raise AssertionError("plot setup should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        output_path = tmp_path / f"qpie_direct_{circular}.png"
        qp._plot_quartet_pie(tree, proportions, str(output_path), branch_info)

        assert output_path.exists()

    def test_plot_reuses_preorder_for_node_positions(self, monkeypatch, tmp_path):
        import phykit.helpers.plot_config as plot_config

        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees=None,
            output=str(tmp_path / "qpie_preorder_positions.png"),
            annotate=False,
            branch_labels=True,
            circular=False,
        )
        qp = QuartetPie(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clades = list(qp._iter_preorder(tree.root))
        proportions = {
            id(clade): (0.7, 0.2, 0.1, 7, 2, 1)
            for clade in clades
            if clade.clades and clade != tree.root
        }
        branch_info = {
            id(clade): {"f1": 7, "pp1": 0.9}
            for clade in clades
            if clade.clades and clade != tree.root
        }
        expected_preorder_ids = [id(clade) for clade in clades]
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

        qp._plot_quartet_pie(tree, proportions, args.output, branch_info)

        assert calls == [expected_preorder_ids]
        assert Path(args.output).exists()

    def test_plot_reuses_clade_lists_for_circular_coords(
        self, monkeypatch, tmp_path
    ):
        import phykit.helpers.circular_layout as circular_layout

        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees=None,
            output=str(tmp_path / "qpie_circular_precomputed_coords.png"),
            annotate=False,
            branch_labels=True,
            circular=True,
            ylabel_fontsize=0,
            no_title=True,
        )
        qp = QuartetPie(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clades = list(qp._iter_preorder(tree.root))
        proportions = {
            id(clade): (0.7, 0.2, 0.1, 7, 2, 1)
            for clade in clades
            if clade.clades and clade != tree.root
        }
        branch_info = {
            id(clade): {"f1": 7, "pp1": 0.9}
            for clade in clades
            if clade.clades and clade != tree.root
        }
        expected_preorder_ids = [id(clade) for clade in clades]
        expected_terminal_ids = [id(clade) for clade in clades if not clade.clades]
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

        qp._plot_quartet_pie(tree, proportions, args.output, branch_info)

        assert calls == [(expected_preorder_ids, expected_terminal_ids)]
        assert Path(args.output).exists()

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_batches_color_clade_overlay(self, monkeypatch, tmp_path, circular):
        from matplotlib.axes import Axes
        from matplotlib.collections import LineCollection

        color_file = tmp_path / "colors.tsv"
        color_file.write_text("A,B\tclade\t#ff0000\tAB\n")
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees=None,
            output=str(tmp_path / "qpie_color.png"),
            annotate=False,
            circular=circular,
            color_file=str(color_file),
            legend_position="none",
            ylabel_fontsize=0,
            no_title=True,
        )
        qp = QuartetPie(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        original_add_collection = Axes.add_collection
        line_collections = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("clade overlay branches should use LineCollection")

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(Axes, "plot", fail_plot)
        monkeypatch.setattr(Axes, "add_collection", capture_collection)

        output_path = tmp_path / f"qpie_color_{circular}.png"
        qp._plot_quartet_pie(tree, {}, str(output_path), {})

        assert len(line_collections) >= 3
        assert output_path.exists()


class TestQuartetPieJson:
    def test_json_output(self, mocker, args):
        args.json = True
        qp = QuartetPie(args)
        mocked_json = mocker.patch("phykit.services.tree.quartet_pie.print_json")
        qp.run()

        payload = mocked_json.call_args.args[0]
        assert payload["n_taxa"] == 8
        assert payload["n_gene_trees"] == 10
        assert payload["input_mode"] == "native"
        assert len(payload["nodes"]) == 5

        # Verify ground truth for bear+raccoon node
        br_node = next(
            n for n in payload["nodes"]
            if set(n["node_tips"]) == {"bear", "raccoon"}
        )
        assert br_node["gCF"] == 1.0
        assert br_node["concordant_count"] == 7

    def test_json_and_csv_use_cached_tip_labels(self, monkeypatch, mocker, tmp_path):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        proportions = {
            id(clade): (0.7, 0.2, 0.1, 7, 2, 1)
            for clade in tree.find_clades(order="preorder")
            if not clade.is_terminal()
        }
        qp = QuartetPie.__new__(QuartetPie)
        qp.output_path = "qpie.png"

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("cached tip labels should be used")

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("direct traversal should be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)
        mocked_json = mocker.patch("phykit.services.tree.quartet_pie.print_json")

        qp._print_json(tree, proportions, "native", 3)
        csv_path = tmp_path / "qpie.csv"
        qp._write_csv(tree, proportions, str(csv_path))

        payload = mocked_json.call_args.args[0]
        assert payload["n_taxa"] == 4
        assert len(payload["nodes"]) == 3
        assert csv_path.read_text().count("\n") == 4

    def test_collect_branch_confidence_uses_direct_traversal(self, monkeypatch):
        tree = Phylo.read(
            StringIO("((A:1,B:1)0.81:1,(C:1,D:1)0.72:1)0.95;"),
            "newick",
        )
        expected = {
            id(clade): clade.confidence
            for clade in QuartetPie._iter_preorder(tree.root)
            if clade.confidence is not None
        }

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("direct traversal should be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        assert QuartetPie._collect_branch_confidence(tree) == expected

    def test_iter_preorder_matches_biopython_order_with_multifurcations(self):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1,G:1);"),
            "newick",
        )

        direct = list(QuartetPie._iter_preorder(tree.root))
        reference = list(tree.find_clades(order="preorder"))

        assert [id(clade) for clade in direct] == [
            id(clade) for clade in reference
        ]

    def test_iter_postorder_matches_biopython_order(self):
        tree = Phylo.read(
            StringIO("((A:1,B:2,C:3):4,(D:5,E:6):7,F:8):0;"),
            "newick",
        )

        direct = list(QuartetPie._iter_postorder(tree.root))
        reference = list(tree.find_clades(order="postorder"))

        assert [id(clade) for clade in direct] == [
            id(clade) for clade in reference
        ]

    def test_collect_clade_tip_names_binary_children_use_indexed_path(self):
        from Bio.Phylo.BaseTree import Clade, Tree

        class IndexedOnlyList(list):
            def __iter__(self):
                raise AssertionError("binary tip-name cache should index children")

        left = Clade(
            name="left",
            clades=IndexedOnlyList([Clade(name="A"), Clade(name="B")]),
        )
        right = Clade(
            name="right",
            clades=IndexedOnlyList([Clade(name="C"), Clade(name="D")]),
        )
        root = Clade(name="root", clades=IndexedOnlyList([left, right]))
        tree = Tree(root=root)

        tip_names = QuartetPie._collect_clade_tip_names(tree)

        assert tip_names[id(left)] == ("A", "B")
        assert tip_names[id(right)] == ("C", "D")
        assert tip_names[id(root)] == ("A", "B", "C", "D")


class TestQuartetPieValidation:
    def test_too_few_tips_exits(self, tmp_path):
        from Bio import Phylo
        from io import StringIO

        # Write a 3-tip tree
        tree_path = tmp_path / "small.tre"
        tree_path.write_text("(A:1,B:1,C:1);")

        args = Namespace(
            tree=str(tree_path),
            gene_trees=None,
            output=str(tmp_path / "out.png"),
            annotate=False,
        )
        qp = QuartetPie(args)
        with pytest.raises(SystemExit):
            qp.run()
