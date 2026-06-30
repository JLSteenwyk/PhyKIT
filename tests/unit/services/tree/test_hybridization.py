import json
import os
import subprocess
import sys
import tempfile

import pytest
from argparse import Namespace
from unittest.mock import patch

from phykit.errors import PhykitUserError
import phykit.services.tree.hybridization as hybridization_module

TREE_SIMPLE = "tests/sample_files/tree_simple.tre"
GENE_TREES = "tests/sample_files/gene_trees_simple.nwk"


def test_module_import_does_not_import_numpy_biophylo_or_matplotlib():
    code = """
import sys
import phykit.services.tree.hybridization as module

assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.circular_layout" not in sys.modules
assert "phykit.helpers.color_annotations" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "matplotlib" not in sys.modules
assert "matplotlib.pyplot" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def _make_args(**kwargs):
    defaults = dict(
        tree=TREE_SIMPLE,
        gene_trees=GENE_TREES,
        support=None,
        alpha=0.05,
        json=False,
        plot_output=None,
    )
    defaults.update(kwargs)
    return Namespace(**defaults)


def _make_svc(**kwargs):
    from phykit.services.tree.hybridization import Hybridization
    return Hybridization(_make_args(**kwargs))


def test_binomial_two_sided_p_value_matches_expected_values():
    assert hybridization_module._binomial_two_sided_p_value(5, 10) == pytest.approx(1.0)
    assert hybridization_module._binomial_two_sided_p_value(9, 10) == pytest.approx(0.021484375)
    assert hybridization_module._binomial_two_sided_p_value(1, 1) == pytest.approx(1.0)


def test_asymmetry_test_does_not_import_scipy_stats(monkeypatch):
    original_import = __import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.stats" or name.startswith("scipy.stats."):
            raise AssertionError("hybridization binomial p-values should not import scipy.stats")
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr("builtins.__import__", fake_import)

    result = _make_svc()._test_asymmetry(9, 1)
    assert result["p_value"] == pytest.approx(0.021484375)
    assert result["favored_alt"] == "alt1"


def test_count_split_matches_scans_gene_trees_once():
    from phykit.services.tree.hybridization import Hybridization

    concordant = frozenset({"a", "b"})
    alt1 = frozenset({"a", "c"})
    alt2 = frozenset({"a", "d"})

    class CountingSplits:
        def __init__(self, rows):
            self.rows = rows
            self.iterations = 0

        def __iter__(self):
            self.iterations += 1
            return iter(self.rows)

    gene_tree_splits = CountingSplits(
        [
            {concordant, alt1},
            {alt1},
            {alt2},
            {concordant, alt2},
            set(),
        ]
    )

    assert Hybridization._count_split_matches(
        gene_tree_splits, concordant, alt1, alt2
    ) == (2, 2, 2)
    assert gene_tree_splits.iterations == 1


class TestCanonicalSplit:
    def test_equal_size_returns_lexicographically_smaller_side(self):
        from phykit.services.tree.hybridization import Hybridization

        all_taxa = frozenset({"A", "B", "C", "D"})

        assert (
            Hybridization._canonical_split(frozenset({"C", "D"}), all_taxa)
            == frozenset({"A", "B"})
        )
        assert (
            Hybridization._canonical_split(frozenset({"A", "B"}), all_taxa)
            == frozenset({"A", "B"})
        )

    def test_empty_split_still_canonicalizes_to_empty_set(self):
        from phykit.services.tree.hybridization import Hybridization

        assert Hybridization._canonical_split(frozenset(), frozenset()) == frozenset()


class TestPerBranchCounts:
    def test_parse_gene_trees_skips_comments_blanks_and_whitespace(self, tmp_path):
        gene_trees = tmp_path / "gene_trees.nwk"
        gene_trees.write_text(
            "# ignored\n\n  (A:1.0,B:2.0,C:3.0);  \n# also ignored\n"
        )
        svc = _make_svc(gene_trees=str(gene_trees))

        trees = svc._parse_gene_trees(str(gene_trees))

        assert len(trees) == 1

    def test_parse_gene_tree_path_list_avoids_per_row_path_objects(
        self, tmp_path, monkeypatch
    ):
        from pathlib import Path
        import phykit.services.tree.hybridization as module

        (tmp_path / "one.nwk").write_text("(A:1,B:1);\n")
        (tmp_path / "two.nwk").write_text("(A:1,C:1);\n")
        gene_trees = tmp_path / "gene_tree_paths.txt"
        gene_trees.write_text("one.nwk\ntwo.nwk\n")
        parent_joins = 0

        class CountingParent:
            def __init__(self, path):
                self._path = path

            def __str__(self):
                return str(self._path)

            def __truediv__(self, other):
                nonlocal parent_joins
                parent_joins += 1
                return self._path / other

        class CountingPath:
            def __init__(self, path):
                self._path = Path(path)

            @property
            def parent(self):
                return CountingParent(self._path.parent)

            def open(self, *args, **kwargs):
                return self._path.open(*args, **kwargs)

        monkeypatch.setattr(module, "Path", CountingPath)
        svc = _make_svc(gene_trees=str(gene_trees))

        trees = svc._parse_gene_trees(str(gene_trees))

        assert len(trees) == 2
        assert parent_joins == 0

    def test_count_topologies_uses_fast_tip_name_helper_for_species_taxa(self, mocker):
        svc = _make_svc()
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        spy = mocker.spy(svc, "get_tip_names_from_tree")

        result, _ = svc._count_topologies(species_tree, gene_trees)

        assert result
        assert spy.call_count == 1

    def test_per_branch_counts(self):
        svc = _make_svc()
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result, _ = svc._count_topologies(species_tree, gene_trees)
        assert isinstance(result, dict)
        assert len(result) > 0
        for branch_key, data in result.items():
            assert "split" in data
            assert "n_concordant" in data
            assert "n_alt1" in data
            assert "n_alt2" in data
            total = data["n_concordant"] + data["n_alt1"] + data["n_alt2"]
            assert total <= 10

    def test_cached_clade_taxa_paths_do_not_call_get_terminals(self, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        svc = _make_svc()
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        all_taxa = frozenset({"A", "B", "C", "D"})
        parent_map = svc._build_parent_map(tree)
        clade_taxa = svc._collect_clade_taxa(tree)
        node = tree.root.clades[0]

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("cached clade taxa should be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        groups = svc._get_four_groups(tree, node, parent_map, all_taxa, clade_taxa)
        assert groups == (
            frozenset({"A"}),
            frozenset({"B"}),
            frozenset({"C"}),
            frozenset({"D"}),
        )

    def test_build_parent_map_handles_mixed_child_counts(self, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        svc = _make_svc()
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_find_clades(*args, **kwargs):
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

    def test_collect_clade_taxa_handles_mixed_child_counts(self, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        svc = _make_svc()
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("standard clade taxa helper should build directly")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        clade_taxa = svc._collect_clade_taxa(tree)
        terminal, binary, trifurcating = tree.root.clades

        assert clade_taxa[id(terminal)] == frozenset({"A"})
        assert clade_taxa[id(binary)] == frozenset({"B", "C"})
        assert clade_taxa[id(trifurcating)] == frozenset({"D", "E", "F"})
        assert clade_taxa[id(tree.root)] == frozenset({"A", "B", "C", "D", "E", "F"})

    def test_get_terminal_clades_preserves_order_with_mixed_child_counts(
        self, monkeypatch
    ):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        svc = _make_svc()
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("standard terminal helper should build directly")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        tips = svc._get_terminal_clades(tree)

        assert [tip.name for tip in tips] == ["A", "B", "C", "D", "E", "F"]


class TestAsymmetryRatioRange:
    def test_asymmetry_ratio_range(self):
        svc = _make_svc()
        # Test with various inputs
        for n1, n2 in [(5, 5), (9, 1), (1, 9), (3, 7), (0, 10), (10, 0)]:
            result = svc._test_asymmetry(n1, n2)
            if result["asymmetry_ratio"] is not None:
                assert 0.5 <= result["asymmetry_ratio"] <= 1.0

    def test_asymmetry_ratio_none_for_zero(self):
        svc = _make_svc()
        result = svc._test_asymmetry(0, 0)
        assert result["asymmetry_ratio"] is None


class TestReticulationCountNonneg:
    def test_run_reuses_unmodified_species_tree(self):
        svc = _make_svc()
        species_tree = object()
        gene_trees = [object()] * 5
        topology_counts = {
            "A,B": {
                "split": ["A", "B"],
                "n_concordant": 3,
                "n_alt1": 1,
                "n_alt2": 2,
            }
        }

        with patch.object(
            svc,
            "read_tree_file",
            side_effect=AssertionError("run should not copy the cached species tree"),
        ), patch.object(
            svc, "read_tree_file_unmodified", return_value=species_tree
        ) as read_unmodified, patch.object(
            svc, "_parse_gene_trees", return_value=gene_trees
        ), patch.object(
            svc, "_count_topologies",
            return_value=(topology_counts, frozenset({"A", "B", "C", "D"})),
        ) as count_topologies, patch.object(
            svc, "_output_text"
        ):
            svc.run()

        read_unmodified.assert_called_once_with()
        count_topologies.assert_called_once_with(species_tree, gene_trees, None)

    def test_ladderized_plot_uses_copied_species_tree(self):
        args = _make_args(plot_output="hybrid.png")
        args.ladderize = True
        from phykit.services.tree.hybridization import Hybridization

        svc = Hybridization(args)
        analysis_tree = object()
        plot_tree = type(
            "PlotTree",
            (),
            {"ladderize": lambda self: setattr(self, "ladderized", True)},
        )()
        gene_trees = [object()] * 5
        topology_counts = {
            "A,B": {
                "split": ["A", "B"],
                "n_concordant": 3,
                "n_alt1": 1,
                "n_alt2": 2,
            }
        }

        with patch.object(
            svc, "read_tree_file_unmodified", return_value=analysis_tree
        ), patch.object(
            svc, "read_tree_file", return_value=plot_tree
        ) as read_copy, patch.object(
            svc, "_parse_gene_trees", return_value=gene_trees
        ), patch.object(
            svc, "_count_topologies",
            return_value=(topology_counts, frozenset({"A", "B", "C", "D"})),
        ), patch.object(
            svc, "_output_text"
        ), patch.object(
            svc, "_plot"
        ) as plot:
            svc.run()

        read_copy.assert_called_once_with()
        assert getattr(plot_tree, "ladderized", False) is True
        plot.assert_called_once()
        assert plot.call_args.args[0] is plot_tree

    def test_reticulation_count_nonneg(self, capsys):
        svc = _make_svc()
        svc.run()
        captured = capsys.readouterr()
        # Parse the reticulation count from text output
        for line in captured.out.splitlines():
            if "Estimated reticulation events:" in line:
                count = int(line.split(":")[-1].strip())
                assert count >= 0
                break


class TestSupportThreshold:
    def test_support_threshold(self):
        svc_no_thresh = _make_svc()
        svc_with_thresh = _make_svc(support=99999.0)

        species_tree = svc_no_thresh.read_tree_file()
        gene_trees = svc_no_thresh._parse_gene_trees(GENE_TREES)

        result_no, _ = svc_no_thresh._count_topologies(species_tree, gene_trees)
        result_with, _ = svc_with_thresh._count_topologies(species_tree, gene_trees,
                                                           support_threshold=99999.0)

        # With a very high threshold, all gene tree bipartitions should be
        # collapsed, so concordant counts should differ or be zero
        # (gene trees in sample don't have support values, so confidence is None
        #  and the threshold check is skipped; counts should be equal)
        # This verifies the threshold logic doesn't crash.
        assert isinstance(result_with, dict)

    def test_support_threshold_with_supported_trees(self, tmp_path):
        """Test with gene trees that have support values."""
        from Bio import Phylo
        from io import StringIO

        sp_file = tmp_path / "species.tre"
        sp_file.write_text("((a,b),(c,d));")

        # Gene trees with support values
        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = [
            "((a,b)90,(c,d)90);",
            "((a,b)50,(c,d)50);",
            "((a,c)90,(b,d)90);",
        ]
        gt_file.write_text("\n".join(gt_lines))

        svc_low = _make_svc(tree=str(sp_file), gene_trees=str(gt_file), support=40.0)
        svc_high = _make_svc(tree=str(sp_file), gene_trees=str(gt_file), support=80.0)

        species_tree = svc_low.read_tree_file()
        gt_low = svc_low._parse_gene_trees(str(gt_file))
        gt_high = svc_high._parse_gene_trees(str(gt_file))

        result_low, _ = svc_low._count_topologies(species_tree, gt_low, support_threshold=40.0)
        result_high, _ = svc_high._count_topologies(species_tree, gt_high, support_threshold=80.0)

        # Both should produce results without error
        assert isinstance(result_low, dict)
        assert isinstance(result_high, dict)


class TestTextOutput:
    def test_text_output(self, capsys):
        svc = _make_svc()
        svc.run()
        captured = capsys.readouterr()
        assert "Hybridization Analysis" in captured.out
        assert "======================" in captured.out
        assert "Species tree branches:" in captured.out
        assert "Gene trees:" in captured.out
        assert "Estimated reticulation events:" in captured.out
        assert "Non-significant branches:" in captured.out

    def test_output_text_batches_report(self, monkeypatch):
        from phykit.services.tree.hybridization import Hybridization

        svc = Hybridization.__new__(Hybridization)
        branch_results = [
            {
                "taxa": ("a", "b"),
                "n_concordant": 20,
                "n_alt1": 15,
                "n_alt2": 5,
                "gcf": 0.5,
                "asymmetry_ratio": 0.75,
                "fdr_p": 0.012345,
                "favored_alt": "alt1",
                "significant": True,
            },
            {
                "taxa": ("c", "d"),
                "n_concordant": 0,
                "n_alt1": 0,
                "n_alt2": 0,
                "gcf": 1.0,
                "asymmetry_ratio": None,
                "fdr_p": None,
                "favored_alt": None,
                "significant": True,
            },
            {
                "taxa": ("e", "f"),
                "n_concordant": 3,
                "n_alt1": 1,
                "n_alt2": 1,
                "gcf": 0.6,
                "asymmetry_ratio": 0.5,
                "fdr_p": 0.5,
                "favored_alt": "-",
                "significant": False,
            },
        ]
        summary = {
            "n_branches": 3,
            "n_gene_trees": 10,
            "support_threshold": 80,
            "n_reticulations": 2,
            "alpha": 0.05,
        }
        printed = []

        def fake_print(*args, **kwargs):
            printed.append((args, kwargs))

        monkeypatch.setattr("builtins.print", fake_print)

        svc._output_text(branch_results, summary)

        header = (
            f"  {'Branch (taxa)':<35}"
            f"{'gCF':>8}"
            f"{'gDF1':>8}"
            f"{'gDF2':>8}"
            f"{'Ratio':>8}"
            f"{'FDR_p':>10}"
            f"{'Favored':>12}"
        )
        expected = "\n".join([
            "Hybridization Analysis",
            "======================",
            "Species tree branches: 3",
            "Gene trees: 10",
            "Support threshold: 80",
            "",
            "Estimated reticulation events: 2",
            "",
            "Significant branches (FDR < 0.05):",
            header,
            (
                f"  {'a,b':<35}"
                f"{0.5:>8.2f}"
                f"{0.375:>8.2f}"
                f"{0.125:>8.2f}"
                f"{'0.75':>8}"
                f"{'0.0123':>10}"
                f"{'alt1':>12}"
            ),
            (
                f"  {'c,d':<35}"
                f"{1.0:>8.2f}"
                f"{0.0:>8.2f}"
                f"{0.0:>8.2f}"
                f"{'NA':>8}"
                f"{'NA':>10}"
                f"{'-':>12}"
            ),
            "",
            "Non-significant branches: 1",
        ])
        assert printed == [((expected,), {})]


class TestJsonOutput:
    def test_json_output(self, capsys):
        svc = _make_svc(json=True)
        svc.run()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert "n_branches" in data
        assert "n_gene_trees" in data
        assert "alpha" in data
        assert "n_reticulations" in data
        assert "branches" in data
        assert isinstance(data["branches"], list)
        assert len(data["branches"]) > 0

        # Check branch entry has expected keys
        branch = data["branches"][0]
        for key in ["taxa", "n_concordant", "n_alt1", "n_alt2",
                     "gcf", "asymmetry_ratio", "p_value", "fdr_p",
                     "significant", "hybrid_score", "favored_alt"]:
            assert key in branch, f"Missing key: {key}"

    def test_json_n_gene_trees(self, capsys):
        svc = _make_svc(json=True)
        svc.run()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert data["n_gene_trees"] == 10


class TestPlotCreated:
    def test_plot_created(self, tmp_path):
        output = str(tmp_path / "test_hybrid.png")
        svc = _make_svc(plot_output=output)
        svc.run()
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    def test_plot_uses_cached_clade_taxa(self, tmp_path, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        output = str(tmp_path / "test_hybrid_cached.png")
        svc = _make_svc(plot_output=output)
        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        branch_results = [
            dict(
                taxa=["A", "B"],
                n_concordant=3,
                n_alt1=1,
                n_alt2=0,
                gcf=0.75,
                asymmetry_ratio=1.0,
                p_value=1.0,
                fdr_p=1.0,
                significant=False,
                hybrid_score=0.0,
                favored_alt=None,
            )
        ]

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("plot setup should use cached clade taxa")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        svc._plot(
            species_tree,
            branch_results,
            output,
            shared_taxa=frozenset({"A", "B", "C", "D"}),
        )

        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    def test_plot_reuses_direct_traversal_lists(self, tmp_path, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        output = str(tmp_path / "test_hybrid_direct_traversal.png")
        args = _make_args(plot_output=output)
        args.legend_position = "none"
        args.ylabel_fontsize = 0
        args.no_title = True
        from phykit.services.tree.hybridization import Hybridization
        svc = Hybridization(args)
        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        branch_results = [
            dict(
                taxa=["A", "B"],
                n_concordant=3,
                n_alt1=1,
                n_alt2=0,
                gcf=0.75,
                asymmetry_ratio=1.0,
                p_value=1.0,
                fdr_p=1.0,
                significant=False,
                hybrid_score=0.75,
                favored_alt=None,
            )
        ]

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("standard plotting should reuse direct traversals")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        svc._plot(
            species_tree,
            branch_results,
            output,
            shared_taxa=frozenset({"A", "B", "C", "D"}),
        )

        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    def test_rectangular_plot_batches_hybridization_branches(self, tmp_path, monkeypatch):
        from Bio import Phylo
        from io import StringIO
        from matplotlib.axes import Axes
        from matplotlib.collections import LineCollection

        output = str(tmp_path / "test_hybrid_batched_rect.png")
        args = _make_args(plot_output=output)
        args.legend_position = "none"
        args.ylabel_fontsize = 0
        args.no_title = True
        from phykit.services.tree.hybridization import Hybridization
        svc = Hybridization(args)
        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        branch_results = [
            dict(
                taxa=["A", "B"],
                n_concordant=3,
                n_alt1=1,
                n_alt2=0,
                gcf=0.75,
                asymmetry_ratio=1.0,
                p_value=1.0,
                fdr_p=1.0,
                significant=False,
                hybrid_score=0.75,
                favored_alt=None,
            ),
            dict(
                taxa=["C", "D"],
                n_concordant=2,
                n_alt1=2,
                n_alt2=0,
                gcf=0.5,
                asymmetry_ratio=1.0,
                p_value=1.0,
                fdr_p=1.0,
                significant=False,
                hybrid_score=0.0,
                favored_alt=None,
            ),
        ]

        original_add_collection = Axes.add_collection
        line_collections = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("branch rendering should use LineCollection")

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(Axes, "plot", fail_plot)
        monkeypatch.setattr(Axes, "add_collection", capture_collection)

        svc._plot(
            species_tree,
            branch_results,
            output,
            shared_taxa=frozenset({"A", "B", "C", "D"}),
        )

        assert len(line_collections) >= 2
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    def test_rectangular_plot_reuses_repeated_hybrid_score_colors(self):
        from Bio import Phylo
        from io import StringIO

        class Spine:
            def set_visible(self, value):
                self.visible = value

        class AxisLabel:
            def set_fontsize(self, value):
                self.fontsize = value

        class Axis:
            def __init__(self):
                self.label = AxisLabel()

        class FakeAxes:
            def __init__(self):
                self.collections = []
                self.spines = {
                    "top": Spine(),
                    "right": Spine(),
                    "left": Spine(),
                }
                self.xaxis = Axis()

            def add_collection(self, collection):
                self.collections.append(collection)

            def autoscale_view(self):
                pass

            def scatter(self, *args, **kwargs):
                pass

            def text(self, *args, **kwargs):
                pass

            def set_xlabel(self, *args, **kwargs):
                pass

            def set_yticks(self, *args, **kwargs):
                pass

            def set_title(self, *args, **kwargs):
                pass

        class CountingCmap:
            def __init__(self):
                self.calls = 0

            def __call__(self, value):
                self.calls += 1
                return "red"

        args = _make_args(plot_output="hybrid.png")
        from phykit.services.tree.hybridization import Hybridization
        svc = Hybridization(args)
        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        parent_map = svc._build_parent_map(species_tree)
        preorder_clades = svc._preorder_clades(species_tree)
        node_x = {id(species_tree.root): 0.0}
        node_y = {}
        tips = []
        for idx, clade in enumerate(preorder_clades):
            if clade.clades:
                continue
            tips.append(clade)
            node_y[id(clade)] = float(idx)
        for clade in preorder_clades:
            if clade is not species_tree.root:
                parent = parent_map[id(clade)]
                node_x[id(clade)] = node_x[id(parent)] + (
                    clade.branch_length or 0.0
                )
        for clade in reversed(preorder_clades):
            if clade.clades:
                node_y[id(clade)] = sum(
                    node_y[id(child)] for child in clade.clades
                ) / len(clade.clades)
        node_to_result = {
            id(clade): {"hybrid_score": 0.75, "significant": False}
            for clade in preorder_clades
            if clade is not species_tree.root
        }
        cmap = CountingCmap()
        config = Namespace(
            ylabel_fontsize=0,
            legend_position="none",
            show_title=False,
            axis_fontsize=None,
        )

        svc._plot_rectangular(
            FakeAxes(),
            None,
            species_tree,
            species_tree.root,
            tips,
            parent_map,
            node_x,
            node_y,
            node_to_result,
            cmap,
            lambda score: score,
            config,
            preorder_clades,
        )

        assert cmap.calls == 1


class TestCircularPlot:
    def test_circular_plot(self, tmp_path):
        output = str(tmp_path / "test_hybrid_circular.png")
        args = _make_args(plot_output=output)
        args.circular = True
        from phykit.services.tree.hybridization import Hybridization
        svc = Hybridization(args)
        svc.run()
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    def test_circular_plot_batches_hybridization_branches(self, tmp_path, monkeypatch):
        from Bio import Phylo
        from io import StringIO
        from matplotlib.axes import Axes
        from matplotlib.collections import LineCollection

        output = str(tmp_path / "test_hybrid_batched_circular.png")
        args = _make_args(plot_output=output)
        args.circular = True
        args.legend_position = "none"
        args.ylabel_fontsize = 0
        args.no_title = True
        from phykit.services.tree.hybridization import Hybridization
        svc = Hybridization(args)
        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        branch_results = [
            dict(
                taxa=["A", "B"],
                n_concordant=3,
                n_alt1=1,
                n_alt2=0,
                gcf=0.75,
                asymmetry_ratio=1.0,
                p_value=1.0,
                fdr_p=1.0,
                significant=False,
                hybrid_score=0.75,
                favored_alt=None,
            ),
            dict(
                taxa=["C", "D"],
                n_concordant=2,
                n_alt1=2,
                n_alt2=0,
                gcf=0.5,
                asymmetry_ratio=1.0,
                p_value=1.0,
                fdr_p=1.0,
                significant=False,
                hybrid_score=0.0,
                favored_alt=None,
            ),
        ]

        original_add_collection = Axes.add_collection
        line_collections = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("branch rendering should use LineCollection")

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(Axes, "plot", fail_plot)
        monkeypatch.setattr(Axes, "add_collection", capture_collection)

        svc._plot(
            species_tree,
            branch_results,
            output,
            shared_taxa=frozenset({"A", "B", "C", "D"}),
        )

        assert len(line_collections) >= 2
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_batches_significant_star_markers(
        self, tmp_path, monkeypatch, circular
    ):
        from Bio import Phylo
        from io import StringIO
        from matplotlib.axes import Axes

        output = str(tmp_path / f"test_hybrid_stars_{circular}.png")
        args = _make_args(plot_output=output)
        args.circular = circular
        args.legend_position = "none"
        args.ylabel_fontsize = 0
        args.no_title = True
        from phykit.services.tree.hybridization import Hybridization
        svc = Hybridization(args)
        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        branch_results = [
            dict(
                taxa=["A", "B"],
                n_concordant=1,
                n_alt1=4,
                n_alt2=0,
                gcf=0.2,
                asymmetry_ratio=1.0,
                p_value=0.01,
                fdr_p=0.01,
                significant=True,
                hybrid_score=0.8,
                favored_alt="alt1",
            ),
            dict(
                taxa=["C", "D"],
                n_concordant=1,
                n_alt1=0,
                n_alt2=4,
                gcf=0.2,
                asymmetry_ratio=1.0,
                p_value=0.01,
                fdr_p=0.01,
                significant=True,
                hybrid_score=0.7,
                favored_alt="alt2",
            ),
        ]

        original_scatter = Axes.scatter
        star_counts = []

        def capture_scatter(self, x, y, *args, **kwargs):
            if kwargs.get("marker") == "*" and kwargs.get("zorder") == 5:
                try:
                    star_counts.append(len(x))
                except TypeError:
                    star_counts.append(1)
            return original_scatter(self, x, y, *args, **kwargs)

        monkeypatch.setattr(Axes, "scatter", capture_scatter)

        svc._plot(
            species_tree,
            branch_results,
            output,
            shared_taxa=frozenset({"A", "B", "C", "D"}),
        )

        assert star_counts == [len(branch_results)]
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0


class TestNoSignificantBranches:
    def test_no_significant_branches(self, capsys):
        # With a very low alpha, nothing should be significant
        svc = _make_svc(alpha=1e-100)
        svc.run()
        captured = capsys.readouterr()
        assert "Estimated reticulation events: 0" in captured.out

    def test_no_significant_json(self, capsys):
        svc = _make_svc(alpha=1e-100, json=True)
        svc.run()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert data["n_reticulations"] == 0
        for branch in data["branches"]:
            assert branch["significant"] is False
            assert branch["hybrid_score"] == 0.0


class TestProcessArgs:
    def test_default_args(self):
        svc = _make_svc()
        assert svc.gene_trees_path == GENE_TREES
        assert svc.support_threshold is None
        assert svc.alpha == 0.05
        assert svc.json_output is False
        assert svc.plot_output is None

    def test_custom_args(self):
        svc = _make_svc(support=70.0, alpha=0.01)
        assert svc.support_threshold == 70.0
        assert svc.alpha == 0.01


class TestFDR:
    def test_empty_list(self):
        from phykit.services.tree.hybridization import Hybridization
        assert Hybridization._fdr([]) == []

    def test_single_pvalue(self):
        from phykit.services.tree.hybridization import Hybridization
        result = Hybridization._fdr([0.03])
        assert result == [0.03]

    def test_known_correction(self):
        from phykit.services.tree.hybridization import Hybridization
        pvals = [0.01, 0.04, 0.03]
        result = Hybridization._fdr(pvals)
        assert result[0] == pytest.approx(0.03)
        assert result[1] == pytest.approx(0.04)
        assert result[2] == pytest.approx(0.04)

    def test_known_correction_with_ties(self):
        from phykit.services.tree.hybridization import Hybridization
        pvals = [0.20, 0.01, 0.01, 0.50, 0.03, 0.80, 0.03]
        indexed = sorted(enumerate(pvals), key=lambda x: x[1])
        expected = [0.0] * len(pvals)
        previous = 1.0
        for rank_minus_1 in range(len(pvals) - 1, -1, -1):
            original_idx, p_value = indexed[rank_minus_1]
            rank = rank_minus_1 + 1
            adjusted = min(p_value * len(pvals) / rank, previous)
            adjusted = min(adjusted, 1.0)
            expected[original_idx] = adjusted
            previous = adjusted

        assert Hybridization._fdr(pvals) == pytest.approx(expected)

    def test_small_fdr_does_not_import_numpy(self):
        code = """
import sys
from phykit.services.tree.hybridization import Hybridization
corrected = Hybridization._fdr([0.20, 0.01, 0.01, 0.50, 0.03, 0.80, 0.03])
assert [round(value, 10) for value in corrected] == [0.28, 0.035, 0.035, 0.5833333333, 0.0525, 0.8, 0.0525]
assert "numpy" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)


class TestCLI:
    def test_hybridization_alias(self, capsys):
        from phykit.phykit import Phykit
        sys.argv = [
            "phykit", "hybrid",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
        ]
        Phykit()
        captured = capsys.readouterr()
        assert "Hybridization Analysis" in captured.out

    def test_reticulation_alias(self, capsys):
        from phykit.phykit import Phykit
        sys.argv = [
            "phykit", "reticulation",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
        ]
        Phykit()
        captured = capsys.readouterr()
        assert "Hybridization Analysis" in captured.out
