"""
Unit tests for phylo_heatmap (phylogenetic heatmap visualization).

Tests initialization, trait matrix parsing, plot creation, and
standardization. Analogous to R's phytools::phylo.heatmap().
"""
import builtins
import subprocess
import sys

import pytest
from argparse import Namespace
from io import StringIO
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, TreeMixin

import phykit.services.tree.phylo_heatmap as module
from phykit.services.tree.phylo_heatmap import PhyloHeatmap


def test_module_import_does_not_import_numpy():
    code = """
import sys
import phykit.services.tree.phylo_heatmap as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert callable(module.build_parent_map)
assert callable(module.compute_circular_coords)
assert callable(module.draw_circular_branches)
assert callable(module.draw_circular_tip_labels)
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "matplotlib" not in sys.modules
assert "matplotlib.pyplot" not in sys.modules
assert "scipy" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.circular_layout" not in sys.modules
assert "phykit.helpers.color_annotations" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = module._LazyNumpy()

    first_asarray = lazy_np.asarray
    second_asarray = lazy_np.asarray

    assert first_asarray is second_asarray
    assert lazy_np.__dict__["asarray"] is first_asarray
    assert lazy_np._module is not None


@pytest.fixture
def args(tmp_path):
    return Namespace(
        tree="tests/sample_files/tree_simple.tre",
        data="tests/sample_files/tree_simple_multi_traits.tsv",
        output=str(tmp_path / "heatmap.png"),
        split=0.3,
        standardize=False,
        cmap="viridis",
    )


class TestPhyloHeatmapInit:
    def test_init_sets_fields(self, args):
        ph = PhyloHeatmap(args)
        assert ph.tree_file_path == args.tree
        assert ph.data_path == args.data
        assert ph.split == 0.3
        assert ph.standardize is False
        assert ph.cmap_name == "viridis"
        assert ph.json_output is False

    def test_invalid_split_exits(self, tmp_path):
        args = Namespace(
            tree="t.tre", data="d.tsv",
            output=str(tmp_path / "out.png"),
            split=0.0,
        )
        with pytest.raises(SystemExit):
            PhyloHeatmap(args)

    def test_invalid_split_above_one_exits(self, tmp_path):
        args = Namespace(
            tree="t.tre", data="d.tsv",
            output=str(tmp_path / "out.png"),
            split=1.0,
        )
        with pytest.raises(SystemExit):
            PhyloHeatmap(args)


class TestTraitParsing:
    @pytest.mark.parametrize(
        "contents",
        [
            "# comment\n\n",
            "taxon\nA\nB\nC\n",
            "taxon\ttrait\n# no data\n",
        ],
    )
    def test_structurally_empty_matrices_are_rejected(
        self, args, tmp_path, contents
    ):
        ph = PhyloHeatmap(args)
        trait_file = tmp_path / "invalid.tsv"
        trait_file.write_text(contents)

        with pytest.raises(SystemExit):
            ph._parse_trait_matrix(str(trait_file), ["A", "B", "C"])

    def test_parse_multi_column_tsv(self, args):
        ph = PhyloHeatmap(args)
        tree = ph.read_tree_file()
        tree_tips = [t.name for t in tree.get_terminals()]
        trait_names, trait_data = ph._parse_trait_matrix(args.data, tree_tips)
        assert trait_names == ["body_mass", "brain_size", "longevity"]
        assert len(trait_data) == 8
        assert "raccoon" in trait_data
        assert len(trait_data["raccoon"]) == 3

    def test_missing_file_exits(self, args):
        ph = PhyloHeatmap(args)
        with pytest.raises(SystemExit):
            ph._parse_trait_matrix("/nonexistent.tsv", ["A", "B", "C"])

    def test_too_few_shared_taxa_exits(self, args, tmp_path):
        ph = PhyloHeatmap(args)
        f = tmp_path / "traits.tsv"
        f.write_text("taxon\tval\nX\t1.0\nY\t2.0\n")
        with pytest.raises(SystemExit):
            ph._parse_trait_matrix(str(f), ["A", "B", "C"])

    def test_non_numeric_exits(self, args, tmp_path):
        ph = PhyloHeatmap(args)
        f = tmp_path / "traits.tsv"
        f.write_text("taxon\tval\nA\tabc\nB\t2.0\nC\t3.0\n")
        with pytest.raises(SystemExit):
            ph._parse_trait_matrix(str(f), ["A", "B", "C"])

    def test_comments_and_blanks_around_header(self, args, tmp_path):
        ph = PhyloHeatmap(args)
        f = tmp_path / "traits.tsv"
        f.write_text(
            "# ignored before header\n"
            "\n"
            "taxon\ttrait1\ttrait2\n"
            "# ignored before data\n"
            "\n"
            "A\t1.0\t2.0\n"
            "B\t3.0\t4.0\n"
            "C\t5.0\t6.0\n"
        )

        trait_names, trait_data = ph._parse_trait_matrix(str(f), ["A", "B", "C"])

        assert trait_names == ["trait1", "trait2"]
        assert trait_data == {
            "A": [1.0, 2.0],
            "B": [3.0, 4.0],
            "C": [5.0, 6.0],
        }

    def test_parse_trait_matrix_skips_whitespace_prefixed_comments(self, args, tmp_path):
        ph = PhyloHeatmap(args)
        f = tmp_path / "traits.tsv"
        f.write_text(
            "   # ignored before header\n"
            "taxon\ttrait1\ttrait2\n"
            "A\t1.0\t2.0\n"
            "   # ignored between rows\n"
            "B\t3.0\t4.0\n"
            "C\t5.0\t6.0\n"
        )

        trait_names, trait_data = ph._parse_trait_matrix(str(f), ["A", "B", "C"])

        assert trait_names == ["trait1", "trait2"]
        assert trait_data == {
            "A": [1.0, 2.0],
            "B": [3.0, 4.0],
            "C": [5.0, 6.0],
        }

    def test_all_shared_trait_matrix_emits_no_warnings(self, args, tmp_path, capsys):
        ph = PhyloHeatmap(args)
        f = tmp_path / "traits.tsv"
        f.write_text(
            "taxon\ttrait1\ttrait2\n"
            "A\t1.0\t2.0\n"
            "B\t3.0\t4.0\n"
            "C\t5.0\t6.0\n"
        )

        trait_names, trait_data = ph._parse_trait_matrix(str(f), ["A", "B", "C"])

        assert trait_names == ["trait1", "trait2"]
        assert trait_data == {
            "A": [1.0, 2.0],
            "B": [3.0, 4.0],
            "C": [5.0, 6.0],
        }
        assert capsys.readouterr().err == ""

    def test_ordered_exact_trait_matrix_skips_set_validation(
        self, args, tmp_path, monkeypatch
    ):
        ph = PhyloHeatmap(args)
        f = tmp_path / "traits.tsv"
        f.write_text(
            "taxon\ttrait1\ttrait2\n"
            "A\t1.0\t2.0\n"
            "B\t3.0\t4.0\n"
            "C\t5.0\t6.0\n"
        )

        def fail_set(*_args, **_kwargs):
            raise AssertionError("exact ordered matrices should not build taxon sets")

        monkeypatch.setattr(builtins, "set", fail_set)

        trait_names, trait_data = ph._parse_trait_matrix(str(f), ["A", "B", "C"])

        assert trait_names == ["trait1", "trait2"]
        assert trait_data == {
            "A": [1.0, 2.0],
            "B": [3.0, 4.0],
            "C": [5.0, 6.0],
        }

    def test_wrong_column_count_exits(self, args, tmp_path):
        ph = PhyloHeatmap(args)
        f = tmp_path / "traits.tsv"
        f.write_text("taxon\ttrait1\ttrait2\nA\t1.0\t2.0\nB\t3.0\nC\t5.0\t6.0\n")

        with pytest.raises(SystemExit):
            ph._parse_trait_matrix(str(f), ["A", "B", "C"])

    def test_unordered_exact_taxa_emit_no_warnings(self, args, tmp_path, capsys):
        ph = PhyloHeatmap(args)
        trait_file = tmp_path / "unordered.tsv"
        trait_file.write_text(
            "taxon\ttrait\nC\t3.0\nA\t1.0\nB\t2.0\n"
        )

        names, data = ph._parse_trait_matrix(
            str(trait_file), ["A", "B", "C"]
        )

        assert names == ["trait"]
        assert data == {"C": [3.0], "A": [1.0], "B": [2.0]}
        assert capsys.readouterr().err == ""

    def test_taxon_mismatches_warn_and_return_shared_data(
        self, args, tmp_path, capsys
    ):
        ph = PhyloHeatmap(args)
        trait_file = tmp_path / "mismatch.tsv"
        trait_file.write_text(
            "taxon\ttrait\nA\t1.0\nB\t2.0\nC\t3.0\nextra\t4.0\n"
        )

        names, data = ph._parse_trait_matrix(
            str(trait_file), ["A", "B", "C", "missing"]
        )

        assert names == ["trait"]
        assert data == {"A": [1.0], "B": [2.0], "C": [3.0]}
        stderr = capsys.readouterr().err
        assert "missing" in stderr
        assert "extra" in stderr


class TestPhyloHeatmapPlot:
    def test_circular_constant_matrix_with_color_annotations(self, args, tmp_path):
        pytest.importorskip("matplotlib")
        tree = Phylo.read(args.tree, "newick")
        taxa = [tip.name for tip in tree.get_terminals()]
        data_file = tmp_path / "constant.tsv"
        data_file.write_text(
            "taxon\tconstant\n"
            + "".join(f"{taxon}\t1.0\n" for taxon in taxa)
        )
        color_file = tmp_path / "colors.tsv"
        color_file.write_text(
            "bear,raccoon\trange\t#ffcc00\tRange\n"
            "bear,raccoon\tclade\t#cc0033\tClade\n"
            "bear\tlabel\t#0033cc\n"
        )
        args.data = str(data_file)
        args.output = str(tmp_path / "circular-constant.png")
        args.circular = True
        args.standardize = True
        args.color_file = str(color_file)
        ph = PhyloHeatmap(args)

        ph.run()

        assert Path(args.output).exists()
        assert Path(args.output).stat().st_size > 0

    def test_clustered_columns_render_dendrogram(self, args):
        pytest.importorskip("matplotlib")
        pytest.importorskip("scipy")
        args.cluster_columns = True
        ph = PhyloHeatmap(args)

        ph.run()

        assert Path(args.output).exists()
        assert Path(args.output).stat().st_size > 0

    @pytest.mark.parametrize(
        "tree",
        [
            Namespace(root=object()),
            Namespace(root=Namespace(clades=[object()])),
        ],
    )
    def test_default_branch_length_scan_handles_generic_trees(self, tree):
        assert PhyloHeatmap._needs_default_branch_lengths(tree) is True

    def test_build_heatmap_matrix_uses_typed_asarray(self, monkeypatch):
        trait_data = {
            "B": [3, 4],
            "A": [1, 2],
            "C": [5, 6],
        }

        def fail_array(*_args, **_kwargs):
            raise AssertionError("heatmap matrix should use np.asarray")

        monkeypatch.setattr(module.np, "array", fail_array, raising=False)

        matrix = PhyloHeatmap._build_heatmap_matrix(["A", "B", "C"], trait_data)

        assert str(matrix.dtype) == "float64"
        assert matrix.tolist() == [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]

    def test_standardize_heatmap_matrix_finite_path_avoids_nan_reductions(
        self, monkeypatch
    ):
        import numpy as np

        matrix = np.array(
            [
                [1.0, 2.0],
                [2.0, 4.0],
                [3.0, 6.0],
            ],
        )

        def fail_nan_reduction(*_args, **_kwargs):
            raise AssertionError("finite heatmap matrix should use plain reductions")

        monkeypatch.setattr(module.np, "nanmean", fail_nan_reduction)
        monkeypatch.setattr(module.np, "nanstd", fail_nan_reduction)

        standardized = PhyloHeatmap._standardize_heatmap_matrix(matrix)

        expected = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0)
        np.testing.assert_allclose(standardized, expected)

    def test_standardize_heatmap_matrix_preserves_nan_fallback(self):
        import numpy as np

        matrix = np.array(
            [
                [1.0, 2.0],
                [np.nan, 4.0],
                [3.0, 6.0],
            ],
        )

        standardized = PhyloHeatmap._standardize_heatmap_matrix(matrix)
        col_means = np.nanmean(matrix, axis=0)
        col_stds = np.nanstd(matrix, axis=0)
        col_stds[col_stds == 0] = 1.0
        expected = (matrix - col_means) / col_stds

        np.testing.assert_allclose(standardized, expected, equal_nan=True)

    def test_run_uses_fast_tip_name_helper_for_setup(self, mocker, args):
        args.json = True
        ph = PhyloHeatmap(args)
        tip_name_spy = mocker.spy(ph, "get_tip_names_from_tree")
        mocker.patch.object(PhyloHeatmap, "_plot_phylo_heatmap")
        mocker.patch("phykit.services.tree.phylo_heatmap.print_json")

        ph.run()

        assert tip_name_spy.call_count == 1

    def test_run_all_tips_present_uses_read_only_tree_without_copy(
        self, monkeypatch, args, capsys
    ):
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        trait_data = {
            "A": [1.0],
            "B": [2.0],
            "C": [3.0],
            "D": [4.0],
        }
        ph = PhyloHeatmap(args)
        captured = {}

        monkeypatch.setattr(ph, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(
            ph,
            "read_tree_file",
            lambda: (_ for _ in ()).throw(
                AssertionError("mutable tree read should not be used")
            ),
        )
        monkeypatch.setattr(
            ph,
            "_fast_copy",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("tree should not be copied")
            ),
        )
        monkeypatch.setattr(
            ph,
            "_parse_trait_matrix",
            lambda *_args, **_kwargs: (["trait"], trait_data),
        )
        monkeypatch.setattr(
            ph,
            "_plot_phylo_heatmap",
            lambda tree_arg, tip_order, *_args, **_kwargs: captured.update(
                tree=tree_arg,
                tip_order=tip_order,
            ),
        )

        ph.run()

        capsys.readouterr()
        assert captured["tree"] is tree
        assert captured["tip_order"] == ["A", "B", "C", "D"]
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}

    def test_run_missing_traits_copies_before_pruning(
        self, monkeypatch, args, capsys
    ):
        class OrderedTraitData(dict):
            def __contains__(self, key):
                raise AssertionError("ordered prune path should not scan membership")

        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        trait_data = OrderedTraitData({
            "A": [1.0],
            "B": [2.0],
            "C": [3.0],
        })
        ph = PhyloHeatmap(args)
        original_fast_copy = ph._fast_copy
        copied_trees = []
        captured = {}

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        monkeypatch.setattr(ph, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(module.Tree, "_ORDERED_MAPPING_PRUNE_MIN_SIZE", 0)
        monkeypatch.setattr(ph, "_fast_copy", copy_spy)
        monkeypatch.setattr(
            ph,
            "_parse_trait_matrix",
            lambda *_args, **_kwargs: (["trait"], trait_data),
        )
        monkeypatch.setattr(
            ph,
            "_plot_phylo_heatmap",
            lambda tree_arg, tip_order, *_args, **_kwargs: captured.update(
                tree=tree_arg,
                tip_order=tip_order,
            ),
        )

        ph.run()

        capsys.readouterr()
        assert len(copied_trees) == 1
        assert captured["tree"] is copied_trees[0]
        assert captured["tip_order"] == ["A", "B", "C"]
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}
        assert {tip.name for tip in copied_trees[0].get_terminals()} == {
            "A",
            "B",
            "C",
        }

    def test_default_branch_length_scan_handles_mixed_child_counts(
        self, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree scan should not use find_clades")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        assert PhyloHeatmap._needs_default_branch_lengths(tree) is False
        tree.root.clades[2].clades[1].branch_length = None
        assert PhyloHeatmap._needs_default_branch_lengths(tree) is True

    def test_run_missing_branch_lengths_copies_before_validation(
        self, monkeypatch, args, capsys
    ):
        tree = Phylo.read(
            StringIO("((A,B):1,(C:1,D:1):1);"),
            "newick",
        )
        trait_data = {
            "A": [1.0],
            "B": [2.0],
            "C": [3.0],
            "D": [4.0],
        }
        ph = PhyloHeatmap(args)
        original_fast_copy = ph._fast_copy
        copied_trees = []
        captured = {}

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        monkeypatch.setattr(ph, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(ph, "_fast_copy", copy_spy)
        monkeypatch.setattr(
            ph,
            "_parse_trait_matrix",
            lambda *_args, **_kwargs: (["trait"], trait_data),
        )
        monkeypatch.setattr(
            ph,
            "_plot_phylo_heatmap",
            lambda tree_arg, tip_order, *_args, **_kwargs: captured.update(
                tree=tree_arg,
                tip_order=tip_order,
            ),
        )

        ph.run()

        capsys.readouterr()
        assert len(copied_trees) == 1
        assert captured["tree"] is copied_trees[0]
        assert tree.root.clades[0].clades[0].branch_length is None
        assert copied_trees[0].root.clades[0].clades[0].branch_length == 1e-8

    def test_creates_png(self, args):
        ph = PhyloHeatmap(args)
        ph.run()
        assert Path(args.output).exists()
        assert Path(args.output).stat().st_size > 0

    def test_nonclustered_rectangular_plot_does_not_import_scipy_cluster(
        self, tmp_path
    ):
        code = f"""
import sys
from argparse import Namespace
from phykit.services.tree.phylo_heatmap import PhyloHeatmap

args = Namespace(
    tree="tests/sample_files/tree_simple.tre",
    data="tests/sample_files/tree_simple_multi_traits.tsv",
    output={str(tmp_path / "heatmap_no_cluster.png")!r},
    split=0.3,
    standardize=False,
    cmap="viridis",
    cluster_columns=False,
    circular=False,
    json=False,
    no_title=True,
    ylabel_fontsize=0,
)
PhyloHeatmap(args).run()
assert "scipy.cluster.hierarchy" not in sys.modules
assert "scipy.spatial.distance" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)
        assert (tmp_path / "heatmap_no_cluster.png").exists()

    def test_creates_pdf(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            data="tests/sample_files/tree_simple_multi_traits.tsv",
            output=str(tmp_path / "heatmap.pdf"),
            split=0.3, standardize=False, cmap="viridis",
        )
        ph = PhyloHeatmap(args)
        ph.run()
        assert Path(args.output).exists()

    def test_standardize_flag(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            data="tests/sample_files/tree_simple_multi_traits.tsv",
            output=str(tmp_path / "heatmap_std.png"),
            split=0.3, standardize=True, cmap="viridis",
        )
        ph = PhyloHeatmap(args)
        ph.run()
        assert Path(args.output).exists()

    def test_custom_cmap(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            data="tests/sample_files/tree_simple_multi_traits.tsv",
            output=str(tmp_path / "heatmap_hot.png"),
            split=0.3, standardize=False, cmap="hot",
        )
        ph = PhyloHeatmap(args)
        ph.run()
        assert Path(args.output).exists()

    def test_custom_split(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            data="tests/sample_files/tree_simple_multi_traits.tsv",
            output=str(tmp_path / "heatmap_split.png"),
            split=0.5, standardize=False, cmap="viridis",
        )
        ph = PhyloHeatmap(args)
        ph.run()
        assert Path(args.output).exists()

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_uses_direct_tree_traversal(self, monkeypatch, tmp_path, circular):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            data="tests/sample_files/tree_simple_multi_traits.tsv",
            output=str(tmp_path / "heatmap_direct.png"),
            split=0.3,
            standardize=False,
            cmap="viridis",
            circular=circular,
        )
        ph = PhyloHeatmap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tip_order = ["A", "B", "C", "D"]
        trait_names = ["x", "y"]
        trait_data = {
            "A": [1.0, 2.0],
            "B": [2.0, 3.0],
            "C": [3.0, 4.0],
            "D": [4.0, 5.0],
        }

        def fail_traversal(*args, **kwargs):
            raise AssertionError("plot setup should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        ph._plot_phylo_heatmap(
            tree,
            tip_order,
            trait_names,
            trait_data,
            str(tmp_path / f"heatmap_direct_{circular}.png"),
        )

        assert Path(tmp_path / f"heatmap_direct_{circular}.png").exists()

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_reuses_clade_lists_for_layout_helpers(
        self, monkeypatch, tmp_path, circular
    ):
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            data="tests/sample_files/tree_simple_multi_traits.tsv",
            output=str(tmp_path / "heatmap_layout_lists.png"),
            split=0.3,
            standardize=False,
            cmap="viridis",
            circular=circular,
            ylabel_fontsize=0,
            no_title=True,
        )
        ph = PhyloHeatmap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tip_order = ["A", "B", "C", "D"]
        trait_names = ["x", "y"]
        trait_data = {
            "A": [1.0, 2.0],
            "B": [2.0, 3.0],
            "C": [3.0, 4.0],
            "D": [4.0, 5.0],
        }
        layout_calls = {}
        original_node_positions = module.compute_node_positions
        original_circular_coords = module.compute_circular_coords

        def spy_node_positions(*args, **kwargs):
            layout_calls["node_preorder"] = kwargs.get("preorder_clades")
            return original_node_positions(*args, **kwargs)

        def spy_circular_coords(*args, **kwargs):
            layout_calls["circular_preorder"] = kwargs.get("preorder_clades")
            layout_calls["circular_tips"] = kwargs.get("terminal_clades")
            return original_circular_coords(*args, **kwargs)

        monkeypatch.setattr(module, "compute_node_positions", spy_node_positions)
        monkeypatch.setattr(module, "compute_circular_coords", spy_circular_coords)

        ph._plot_phylo_heatmap(
            tree,
            tip_order,
            trait_names,
            trait_data,
            str(tmp_path / f"heatmap_layout_lists_{circular}.png"),
        )

        expected_preorder = list(ph._iter_preorder(tree.root))
        assert len(layout_calls["node_preorder"]) == len(expected_preorder)
        assert all(
            actual is expected
            for actual, expected in zip(
                layout_calls["node_preorder"], expected_preorder
            )
        )
        if circular:
            expected_tips = [
                clade for clade in expected_preorder if not clade.clades
            ]
            assert layout_calls["circular_preorder"] is layout_calls["node_preorder"]
            assert len(layout_calls["circular_tips"]) == len(expected_tips)
            assert all(
                actual is expected
                for actual, expected in zip(
                    layout_calls["circular_tips"], expected_tips
                )
            )
        else:
            assert "circular_preorder" not in layout_calls
        assert Path(tmp_path / f"heatmap_layout_lists_{circular}.png").exists()

    def test_iter_preorder_preserves_order_without_reversed(self):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("_iter_preorder should push children directly")

        ph = PhyloHeatmap.__new__(PhyloHeatmap)
        root = Clade(name="root")
        left = Clade(name="left")
        middle = Clade(name="middle")
        right = Clade(name="right")
        left.clades = NoReversedList([Clade(name="left_a"), Clade(name="left_b")])
        middle.clades = NoReversedList([Clade(name="middle_a")])
        right.clades = NoReversedList()
        root.clades = NoReversedList([left, middle, right])

        order = [clade.name for clade in ph._iter_preorder(root)]

        assert order == [
            "root",
            "left",
            "left_a",
            "left_b",
            "middle",
            "middle_a",
            "right",
        ]

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_batches_color_clade_overlay(self, monkeypatch, tmp_path, circular):
        from matplotlib.axes import Axes
        from matplotlib.collections import LineCollection

        color_file = tmp_path / "colors.tsv"
        color_file.write_text("A,B\tclade\t#ff0000\tAB\n")
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            data="tests/sample_files/tree_simple_multi_traits.tsv",
            output=str(tmp_path / f"heatmap_color_{circular}.png"),
            split=0.3,
            standardize=False,
            cmap="viridis",
            circular=circular,
            color_file=str(color_file),
            legend_position="none",
            ylabel_fontsize=0,
            no_title=True,
        )
        ph = PhyloHeatmap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tip_order = ["A", "B", "C", "D"]
        trait_names = ["x", "y"]
        trait_data = {
            "A": [1.0, 2.0],
            "B": [2.0, 3.0],
            "C": [3.0, 4.0],
            "D": [4.0, 5.0],
        }

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

        ph._plot_phylo_heatmap(
            tree,
            tip_order,
            trait_names,
            trait_data,
            str(tmp_path / f"heatmap_color_{circular}.png"),
        )

        assert len(line_collections) >= 3
        assert Path(tmp_path / f"heatmap_color_{circular}.png").exists()

    def test_circular_plot_batches_heatmap_wedges(self, monkeypatch, tmp_path):
        from matplotlib.axes import Axes
        from matplotlib.collections import PatchCollection

        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            data="tests/sample_files/tree_simple_multi_traits.tsv",
            output=str(tmp_path / "heatmap_wedges.png"),
            split=0.3,
            standardize=False,
            cmap="viridis",
            circular=True,
            color_file=None,
            legend_position="none",
            ylabel_fontsize=0,
            no_title=True,
        )
        ph = PhyloHeatmap(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tip_order = ["A", "B", "C", "D"]
        trait_names = ["x", "y", "z"]
        trait_data = {
            "A": [1.0, 2.0, 3.0],
            "B": [2.0, 3.0, 4.0],
            "C": [3.0, 4.0, 5.0],
            "D": [4.0, 5.0, 6.0],
        }

        original_add_collection = Axes.add_collection
        patch_collections = []

        def fail_add_patch(*args, **kwargs):
            raise AssertionError("circular heatmap wedges should use PatchCollection")

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, PatchCollection):
                patch_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(Axes, "add_patch", fail_add_patch)
        monkeypatch.setattr(Axes, "add_collection", capture_collection)

        ph._plot_phylo_heatmap(
            tree,
            tip_order,
            trait_names,
            trait_data,
            str(tmp_path / "heatmap_wedges.png"),
        )

        assert len(patch_collections) == 1
        expected_cells = len(tip_order) * len(trait_names)
        assert len(patch_collections[0].get_paths()) == expected_cells
        assert len(patch_collections[0].get_array()) == expected_cells
        assert Path(tmp_path / "heatmap_wedges.png").exists()


class TestPhyloHeatmapJson:
    def test_json_output(self, mocker, args):
        args.json = True
        ph = PhyloHeatmap(args)
        mocked_json = mocker.patch("phykit.services.tree.phylo_heatmap.print_json")
        ph.run()
        payload = mocked_json.call_args.args[0]
        assert payload["n_taxa"] == 8
        assert payload["n_traits"] == 3
        assert payload["trait_names"] == ["body_mass", "brain_size", "longevity"]
        assert len(payload["tip_order"]) == 8
