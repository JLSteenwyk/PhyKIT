import os
import subprocess
import sys
import tempfile

import pytest
from mock import patch
from argparse import Namespace

from phykit.errors import PhykitUserError
import phykit.services.tree.discordance_asymmetry as discordance_asymmetry_module

TREE_SIMPLE = "tests/sample_files/tree_simple.tre"
GENE_TREES = "tests/sample_files/gene_trees_simple.nwk"


def test_module_import_does_not_import_numpy_biophylo_or_matplotlib():
    code = """
import sys
import phykit.services.tree.discordance_asymmetry as module

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


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = discordance_asymmetry_module._LazyNumpy()

    first_asarray = lazy_np.asarray
    second_asarray = lazy_np.asarray

    assert first_asarray is second_asarray
    assert lazy_np.__dict__["asarray"] is first_asarray
    assert lazy_np._module is not None


def test_binomial_two_sided_p_value_matches_expected_values():
    assert discordance_asymmetry_module._binomial_two_sided_p_value(5, 10) == pytest.approx(1.0)
    assert discordance_asymmetry_module._binomial_two_sided_p_value(9, 10) == pytest.approx(0.021484375)
    assert discordance_asymmetry_module._binomial_two_sided_p_value(1, 1) == pytest.approx(1.0)


def test_large_binomial_two_sided_p_value_matches_scipy_fallback():
    from scipy.special import bdtr

    expected = min(1.0, 2.0 * float(bdtr(2, 65, 0.5)))

    assert discordance_asymmetry_module._binomial_two_sided_p_value(63, 65) == pytest.approx(expected)


def test_small_asymmetry_test_does_not_import_scipy(monkeypatch):
    original_import = __import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.stats" or name.startswith(("scipy.stats.", "scipy.special")):
            raise AssertionError("small discordance asymmetry binomial p-values should not import SciPy")
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr("builtins.__import__", fake_import)

    args = Namespace(
        tree=TREE_SIMPLE,
        gene_trees=GENE_TREES,
        verbose=False,
        json=False,
        plot_output=None,
    )
    result = discordance_asymmetry_module.DiscordanceAsymmetry(args)._test_asymmetry(9, 1)
    assert result["p_value"] == pytest.approx(0.021484375)
    assert result["favored_alt"] == "alt1"


def test_preorder_clades_direct_binary_children_avoid_reversed_iterator():
    from Bio.Phylo.BaseTree import Clade, Tree
    from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

    class NoReversedList(list):
        def __reversed__(self):
            raise AssertionError("binary preorder should not call reversed")

    left = Clade(
        name="left",
        clades=NoReversedList([Clade(name="A"), Clade(name="B")]),
    )
    right = Clade(
        name="right",
        clades=NoReversedList([Clade(name="C"), Clade(name="D")]),
    )
    root = Clade(name="root", clades=NoReversedList([left, right]))
    tree = Tree(root=root)

    assert DiscordanceAsymmetry._preorder_clades_direct(tree) == [
        root,
        left,
        left.clades[0],
        left.clades[1],
        right,
        right.clades[0],
        right.clades[1],
    ]


def test_count_split_matches_scans_gene_trees_once():
    from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

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

    assert DiscordanceAsymmetry._count_split_matches(
        gene_tree_splits, concordant, alt1, alt2
    ) == (2, 2, 2)
    assert gene_tree_splits.iterations == 1


class TestProcessArgs:
    def test_default_args(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)
        assert svc.gene_trees_path == GENE_TREES
        assert svc.verbose is False
        assert svc.json_output is False
        assert svc.plot_output is None

    def test_verbose_arg(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=True,
            json=False,
            plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)
        assert svc.verbose is True

    def test_json_arg(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=True,
            plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)
        assert svc.json_output is True

    def test_plot_output_arg(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output="/tmp/test.png",
        )
        svc = DiscordanceAsymmetry(args)
        assert svc.plot_output == "/tmp/test.png"


class TestCanonicalSplit:
    def test_equal_size_returns_lexicographically_smaller_side(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        all_taxa = frozenset({"A", "B", "C", "D"})

        assert (
            DiscordanceAsymmetry._canonical_split(
                frozenset({"C", "D"}),
                all_taxa,
            )
            == frozenset({"A", "B"})
        )
        assert (
            DiscordanceAsymmetry._canonical_split(
                frozenset({"A", "B"}),
                all_taxa,
            )
            == frozenset({"A", "B"})
        )

    def test_empty_split_still_canonicalizes_to_empty_set(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        assert (
            DiscordanceAsymmetry._canonical_split(frozenset(), frozenset())
            == frozenset()
        )


class TestParseGeneTrees:
    def test_parses_multi_newick(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)
        trees = svc._parse_gene_trees(GENE_TREES)
        assert len(trees) == 10
        for gt in trees:
            tips = [t.name for t in gt.get_terminals()]
            assert len(tips) == 8

    def test_parse_gene_trees_skips_comments_blanks_and_whitespace(self, tmp_path):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        gene_trees = tmp_path / "gene_trees.nwk"
        gene_trees.write_text(
            "   # ignored\n\n  (A:1.0,B:2.0,C:3.0);  \n\t# also ignored\n"
        )
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=str(gene_trees),
            verbose=False, json=False, plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)

        trees = svc._parse_gene_trees(str(gene_trees))

        assert len(trees) == 1

    def test_parse_gene_tree_path_list_avoids_per_row_path_objects(
        self, tmp_path, monkeypatch
    ):
        from pathlib import Path
        import phykit.services.tree.discordance_asymmetry as module
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

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
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=str(gene_trees),
            verbose=False,
            json=False,
            plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)

        trees = svc._parse_gene_trees(str(gene_trees))

        assert len(trees) == 2
        assert parent_joins == 0

    def test_missing_file_raises_error(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees="nonexistent.nwk",
            verbose=False, json=False, plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)
        with pytest.raises(PhykitUserError):
            svc._parse_gene_trees("nonexistent.nwk")


class TestCountTopologies:
    @pytest.fixture
    def svc(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return DiscordanceAsymmetry(args)

    def test_count_topologies_uses_fast_tip_name_helper_for_species_taxa(
        self, svc, mocker
    ):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        spy = mocker.spy(svc, "get_tip_names_from_tree")

        result, _ = svc._count_topologies(species_tree, gene_trees)

        assert result
        assert spy.call_count == 1

    def test_counts_with_sample_data(self, svc):
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

    def test_all_branches_present(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result, _ = svc._count_topologies(species_tree, gene_trees)
        # Expected branches for the 8-taxon tree:
        # ((raccoon,bear),((sea_lion,seal),((monkey,cat),weasel)),dog)
        # Internal branches produce these canonical splits:
        expected_branches = {
            "bear,raccoon",
            "cat,monkey",
            "cat,monkey,weasel",
            "sea_lion,seal",
            "bear,dog,raccoon",
        }
        assert set(result.keys()) == expected_branches

    def test_concordance_proportions_cross_validate(self, svc):
        """For each branch, gCF = n_concordant / (n_concordant + n_alt1 + n_alt2)
        should be > 0.5 for all branches (majority are concordant in sample data)."""
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result, _ = svc._count_topologies(species_tree, gene_trees)
        for branch_key, data in result.items():
            total = data["n_concordant"] + data["n_alt1"] + data["n_alt2"]
            if total == 0:
                continue
            gcf = data["n_concordant"] / total
            assert gcf > 0.5, (
                f"Expected gCF > 0.5 for branch {branch_key}, "
                f"got {gcf:.3f} (conc={data['n_concordant']}, "
                f"alt1={data['n_alt1']}, alt2={data['n_alt2']})"
            )

    def test_cached_clade_taxa_paths_do_not_call_get_terminals(self, svc, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

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

    def test_get_four_groups_merges_multifurcation_extras(self, svc, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        tree = Phylo.read(
            StringIO("((A:1,B:1,C:1,D:1):1,(E:1,F:1,G:1,H:1):1);"),
            "newick",
        )
        all_taxa = frozenset({"A", "B", "C", "D", "E", "F", "G", "H"})
        parent_map = svc._build_parent_map(tree)
        clade_taxa = svc._collect_clade_taxa(tree)
        node = tree.root.clades[0]

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("cached clade taxa should be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        groups = svc._get_four_groups(tree, node, parent_map, all_taxa, clade_taxa)
        assert groups == (
            frozenset({"A"}),
            frozenset({"B", "C", "D"}),
            frozenset({"E"}),
            frozenset({"F", "G", "H"}),
        )

    def test_build_parent_map_handles_mixed_child_counts(self, svc, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

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

    def test_collect_clade_taxa_handles_mixed_child_counts(self, svc, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

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
        self, svc, monkeypatch
    ):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("standard terminal helper should build directly")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        tips = svc._get_terminal_clades(tree)

        assert [tip.name for tip in tips] == ["A", "B", "C", "D", "E", "F"]


class TestTestAsymmetry:
    def _make_svc(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees="tests/sample_files/gene_trees_simple.nwk",
            verbose=False, json=False, plot_output=None,
        )
        return DiscordanceAsymmetry(args)

    def test_symmetric_case(self):
        svc = self._make_svc()
        result = svc._test_asymmetry(5, 5)
        assert result["asymmetry_ratio"] == 0.5
        assert result["p_value"] == pytest.approx(1.0)
        assert result["favored_alt"] is None

    def test_asymmetric_case(self):
        svc = self._make_svc()
        result = svc._test_asymmetry(9, 1)
        assert result["asymmetry_ratio"] == 0.9
        assert result["p_value"] < 0.05
        assert result["favored_alt"] == "alt1"

    def test_zero_discordance(self):
        svc = self._make_svc()
        result = svc._test_asymmetry(0, 0)
        assert result["asymmetry_ratio"] is None
        assert result["p_value"] is None
        assert result["favored_alt"] is None

    def test_alt2_favored(self):
        svc = self._make_svc()
        result = svc._test_asymmetry(1, 9)
        assert result["asymmetry_ratio"] == 0.9
        assert result["favored_alt"] == "alt2"

    def test_single_discordant(self):
        svc = self._make_svc()
        result = svc._test_asymmetry(1, 0)
        assert result["asymmetry_ratio"] == 1.0
        assert result["p_value"] == pytest.approx(1.0)
        assert result["favored_alt"] == "alt1"


class TestFDR:
    def test_empty_list(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        assert DiscordanceAsymmetry._fdr([]) == []

    def test_single_pvalue(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        result = DiscordanceAsymmetry._fdr([0.03])
        assert result == [0.03]

    def test_known_correction(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        # Known input: 3 p-values
        pvals = [0.01, 0.04, 0.03]
        result = DiscordanceAsymmetry._fdr(pvals)
        # FDR correction: rank p-values, adjust
        # sorted: (0, 0.01), (2, 0.03), (1, 0.04)
        # rank3 (idx=1, p=0.04): 0.04*3/3 = 0.04, prev=0.04
        # rank2 (idx=2, p=0.03): 0.03*3/2 = 0.045, min(0.045, 0.04)=0.04, prev=0.04
        # rank1 (idx=0, p=0.01): 0.01*3/1 = 0.03, min(0.03, 0.04)=0.03, prev=0.03
        # Result: [0.03, 0.04, 0.04]
        assert result[0] == pytest.approx(0.03)
        assert result[1] == pytest.approx(0.04)
        assert result[2] == pytest.approx(0.04)

    def test_all_ones(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        result = DiscordanceAsymmetry._fdr([1.0, 1.0, 1.0])
        assert all(p == 1.0 for p in result)

    def test_known_correction_with_ties(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
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

        assert DiscordanceAsymmetry._fdr(pvals) == pytest.approx(expected)

    def test_medium_fdr_matches_scalar_reference(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        pvals = [((idx * 37) % 101) / 1000 for idx in range(32)]
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

        assert DiscordanceAsymmetry._fdr(pvals) == pytest.approx(expected)

    def test_vector_fdr_caps_adjusted_values(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        pvals = [1.0 if idx % 5 else 0.001 * (idx + 1) for idx in range(64)]
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

        assert DiscordanceAsymmetry._fdr(pvals) == pytest.approx(expected)
        assert max(expected) == 1.0

    def test_small_fdr_does_not_import_numpy(self):
        code = """
import sys
from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
corrected = DiscordanceAsymmetry._fdr([0.20, 0.01, 0.01, 0.50, 0.03, 0.80, 0.03])
assert [round(value, 10) for value in corrected] == [0.28, 0.035, 0.035, 0.5833333333, 0.0525, 0.8, 0.0525]
assert "numpy" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)


class TestRun:
    def _make_svc(self, verbose=False, json_output=False, plot_output=None):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees="tests/sample_files/gene_trees_simple.nwk",
            verbose=verbose, json=json_output, plot_output=plot_output,
        )
        return DiscordanceAsymmetry(args)

    def test_run_reuses_unmodified_species_tree(self):
        svc = self._make_svc()
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
        count_topologies.assert_called_once_with(species_tree, gene_trees)

    def test_ladderized_plot_uses_copied_species_tree(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees="tests/sample_files/gene_trees_simple.nwk",
            verbose=False,
            json=False,
            plot_output="asym.png",
            ladderize=True,
        )
        svc = DiscordanceAsymmetry(args)
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

    def test_text_output_has_header(self, capsys):
        svc = self._make_svc()
        svc.run()
        captured = capsys.readouterr()
        assert "branch" in captured.out
        assert "n_conc" in captured.out
        assert "n_alt1" in captured.out
        assert "n_alt2" in captured.out
        assert "asym_ratio" in captured.out
        assert "binom_p" in captured.out
        assert "fdr_p" in captured.out
        assert "gene_flow" in captured.out

    def test_text_output_has_branches(self, capsys):
        svc = self._make_svc()
        svc.run()
        captured = capsys.readouterr()
        assert "bear,raccoon" in captured.out
        assert "cat,monkey" in captured.out

    def test_text_output_has_summary(self, capsys):
        svc = self._make_svc()
        svc.run()
        captured = capsys.readouterr()
        assert "Summary:" in captured.out
        assert "branches tested" in captured.out

    def test_json_output_structure(self, capsys):
        import json
        svc = self._make_svc(json_output=True)
        svc.run()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert "branches" in data
        assert "summary" in data
        assert isinstance(data["branches"], list)
        assert len(data["branches"]) > 0
        # Check branch entry has expected keys
        branch = data["branches"][0]
        for key in ["split", "n_concordant", "n_alt1", "n_alt2",
                     "asymmetry_ratio", "p_value", "fdr_p"]:
            assert key in branch
        # Check summary
        assert "n_gene_trees" in data["summary"]
        assert "n_branches_tested" in data["summary"]
        assert "n_significant_fdr05" in data["summary"]

    def test_json_output_n_gene_trees(self, capsys):
        import json
        svc = self._make_svc(json_output=True)
        svc.run()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert data["summary"]["n_gene_trees"] == 10

    def test_verbose_output(self, capsys):
        svc = self._make_svc(verbose=True)
        svc.run()
        captured = capsys.readouterr()
        assert "Branch:" in captured.out
        assert "gCF=" in captured.out
        assert "gDF1=" in captured.out
        assert "gDF2=" in captured.out

    def test_output_text_batches_branch_rows(self):
        svc = self._make_svc(verbose=True)
        branch_results = [
            {
                "split": ["bear", "raccoon"],
                "n_concordant": 3,
                "n_alt1": 2,
                "n_alt2": 1,
                "asymmetry_ratio": 2.0,
                "p_value": 0.125,
                "fdr_p": 0.04,
                "favored_alt": "alt1",
            }
        ]
        summary = {
            "n_branches_tested": 1,
            "n_significant_fdr05": 1,
        }

        with patch("builtins.print") as mocked_print:
            svc._output_text(branch_results, summary)

        mocked_print.assert_called_once()
        output = mocked_print.call_args.args[0]
        header = (
            f"{'branch':<30}"
            f"{'n_conc':>8}"
            f"{'n_alt1':>8}"
            f"{'n_alt2':>8}"
            f"{'asym_ratio':>12}"
            f"{'binom_p':>12}"
            f"{'fdr_p':>12}"
            f"{'gene_flow':>12}"
        )
        expected = "\n".join([
            header,
            "-" * len(header),
            (
                f"{'bear,raccoon':<30}"
                f"{3:>8}"
                f"{2:>8}"
                f"{1:>8}"
                f"{'2.000':>12}"
                f"{'0.1250':>12}"
                f"{'0.0400':>12}"
                f"{'alt1':>12}"
            ),
            "---",
            "Summary: 1 branches tested, 1 significant (FDR<0.05)",
            "",
            "Branch: bear,raccoon",
            "  gCF=0.500  gDF1=2/6  gDF2=1/6",
        ])
        assert output == expected


class TestPlot:
    def _make_svc(self, plot_output):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees="tests/sample_files/gene_trees_simple.nwk",
            verbose=False, json=False, plot_output=plot_output,
        )
        return DiscordanceAsymmetry(args)

    def test_plot_creates_file(self, tmp_path):
        output = str(tmp_path / "test_asym.png")
        svc = self._make_svc(output)
        svc.run()
        assert os.path.exists(output)

    def test_plot_file_nonempty(self, tmp_path):
        output = str(tmp_path / "test_asym.png")
        svc = self._make_svc(output)
        svc.run()
        assert os.path.getsize(output) > 0

    def test_plot_no_error_with_all_concordant(self, tmp_path):
        # Even if all gene trees are concordant (no discordance),
        # plotting should not error
        output = str(tmp_path / "test_asym_conc.png")
        svc = self._make_svc(output)
        svc.run()
        # If we get here without error, the test passes
        assert os.path.exists(output)

    def test_plot_reuses_direct_traversal_lists(self, tmp_path, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO
        import numpy as numpy_module
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        output = str(tmp_path / "test_asym_direct_traversal.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output=output,
            legend_position="none",
            ylabel_fontsize=0,
            no_title=True,
        )
        svc = DiscordanceAsymmetry(args)
        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        branch_results = [
            dict(
                split=["A", "B"],
                n_concordant=3,
                n_alt1=1,
                n_alt2=0,
                asymmetry_ratio=1.0,
                p_value=1.0,
                fdr_p=1.0,
                favored_alt=None,
            )
        ]

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("standard plotting should reuse direct traversals")

        def fail_np_mean(*args, **kwargs):
            raise AssertionError(
                "plot coordinate setup should average child y values without NumPy"
            )

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)
        monkeypatch.setattr(numpy_module, "mean", fail_np_mean)

        svc._plot(
            species_tree,
            branch_results,
            output,
            shared_taxa=frozenset({"A", "B", "C", "D"}),
        )

        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    def test_rectangular_plot_batches_asymmetry_branches(self, tmp_path, monkeypatch):
        from Bio import Phylo
        from io import StringIO
        from matplotlib.axes import Axes
        from matplotlib.collections import LineCollection
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        output = str(tmp_path / "test_asym_batched_rect.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output=output,
            legend_position="none",
            ylabel_fontsize=0,
            no_title=True,
        )
        svc = DiscordanceAsymmetry(args)
        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        branch_results = [
            dict(
                split=["A", "B"],
                n_concordant=3,
                n_alt1=1,
                n_alt2=0,
                asymmetry_ratio=1.0,
                p_value=1.0,
                fdr_p=1.0,
                favored_alt=None,
            ),
            dict(
                split=["C", "D"],
                n_concordant=2,
                n_alt1=2,
                n_alt2=0,
                asymmetry_ratio=None,
                p_value=None,
                fdr_p=None,
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
        scalar_arrays = [
            collection.get_array()
            for collection in line_collections
            if collection.get_array() is not None
        ]
        assert scalar_arrays
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    def test_rectangular_plot_uses_scalar_asymmetry_ratio_collection(
        self, tmp_path, monkeypatch
    ):
        from Bio import Phylo
        from io import StringIO
        from matplotlib.axes import Axes
        from matplotlib.collections import LineCollection
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        output = str(tmp_path / "test_asym_repeated_colors.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output=output,
            legend_position="none",
            ylabel_fontsize=0,
            no_title=True,
        )
        svc = DiscordanceAsymmetry(args)
        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        branch_results = [
            dict(
                split=["A", "B"],
                n_concordant=3,
                n_alt1=1,
                n_alt2=0,
                asymmetry_ratio=0.75,
                p_value=1.0,
                fdr_p=1.0,
                favored_alt=None,
            ),
            dict(
                split=["C", "D"],
                n_concordant=3,
                n_alt1=1,
                n_alt2=0,
                asymmetry_ratio=0.75,
                p_value=1.0,
                fdr_p=1.0,
                favored_alt=None,
            ),
        ]
        original_add_collection = Axes.add_collection
        line_collections = []

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(Axes, "add_collection", capture_collection)

        svc._plot(
            species_tree,
            branch_results,
            output,
            shared_taxa=frozenset({"A", "B", "C", "D"}),
        )

        scalar_arrays = [
            collection.get_array()
            for collection in line_collections
            if collection.get_array() is not None
        ]
        assert len(scalar_arrays) == 1
        assert list(scalar_arrays[0]) == [0.75, 0.75]
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    def test_rectangular_plot_skips_redundant_tight_layout(
        self, tmp_path, monkeypatch
    ):
        from Bio import Phylo
        from io import StringIO
        from matplotlib.figure import Figure
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        output = str(tmp_path / "test_asym_no_tight_layout.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output=output,
            legend_position="none",
            ylabel_fontsize=0,
            no_title=True,
        )
        svc = DiscordanceAsymmetry(args)
        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        branch_results = [
            dict(
                split=["A", "B"],
                n_concordant=3,
                n_alt1=1,
                n_alt2=0,
                asymmetry_ratio=1.0,
                p_value=1.0,
                fdr_p=1.0,
                favored_alt=None,
            )
        ]

        def fail_tight_layout(*args, **kwargs):
            raise AssertionError("bbox_inches='tight' handles saved bounds")

        monkeypatch.setattr(Figure, "tight_layout", fail_tight_layout)

        svc._plot(
            species_tree,
            branch_results,
            output,
            shared_taxa=frozenset({"A", "B", "C", "D"}),
        )

        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    def test_circular_plot_skips_redundant_tight_layout(
        self, tmp_path, monkeypatch
    ):
        from Bio import Phylo
        from io import StringIO
        from matplotlib.figure import Figure
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        output = str(tmp_path / "test_asym_circular_no_tight_layout.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output=output,
            circular=True,
            legend_position="none",
            ylabel_fontsize=0,
            no_title=True,
        )
        svc = DiscordanceAsymmetry(args)
        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        branch_results = [
            dict(
                split=["A", "B"],
                n_concordant=3,
                n_alt1=1,
                n_alt2=0,
                asymmetry_ratio=1.0,
                p_value=1.0,
                fdr_p=1.0,
                favored_alt=None,
            )
        ]

        def fail_tight_layout(*args, **kwargs):
            raise AssertionError("bbox_inches='tight' handles saved bounds")

        monkeypatch.setattr(Figure, "tight_layout", fail_tight_layout)

        svc._plot(
            species_tree,
            branch_results,
            output,
            shared_taxa=frozenset({"A", "B", "C", "D"}),
        )

        assert os.path.exists(output)
        assert os.path.getsize(output) > 0


class TestPlotCircular:
    def _make_svc(self, plot_output):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        args = Namespace(
            tree="tests/sample_files/tree_simple.tre",
            gene_trees="tests/sample_files/gene_trees_simple.nwk",
            verbose=False, json=False, plot_output=plot_output,
            circular=True,
        )
        return DiscordanceAsymmetry(args)

    def test_circular_plot_creates_file(self, tmp_path):
        output = str(tmp_path / "test_asym_circular.png")
        svc = self._make_svc(output)
        svc.run()
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    def test_circular_plot_batches_asymmetry_branches(self, tmp_path, monkeypatch):
        from Bio import Phylo
        from io import StringIO
        from matplotlib.axes import Axes
        from matplotlib.collections import LineCollection
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        output = str(tmp_path / "test_asym_batched_circular.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output=output,
            circular=True,
            legend_position="none",
            ylabel_fontsize=0,
            no_title=True,
        )
        svc = DiscordanceAsymmetry(args)
        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        branch_results = [
            dict(
                split=["A", "B"],
                n_concordant=3,
                n_alt1=1,
                n_alt2=0,
                asymmetry_ratio=1.0,
                p_value=1.0,
                fdr_p=1.0,
                favored_alt=None,
            ),
            dict(
                split=["C", "D"],
                n_concordant=2,
                n_alt1=2,
                n_alt2=0,
                asymmetry_ratio=None,
                p_value=None,
                fdr_p=None,
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
        scalar_arrays = [
            collection.get_array()
            for collection in line_collections
            if collection.get_array() is not None
        ]
        assert scalar_arrays
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_batches_significant_star_markers(
        self, tmp_path, monkeypatch, circular
    ):
        from Bio import Phylo
        from io import StringIO
        from matplotlib.axes import Axes
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        output = str(tmp_path / f"test_asym_stars_{circular}.png")
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output=output,
            circular=circular,
            legend_position="none",
            ylabel_fontsize=0,
            no_title=True,
        )
        svc = DiscordanceAsymmetry(args)
        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        branch_results = [
            dict(
                split=["A", "B"],
                n_concordant=1,
                n_alt1=4,
                n_alt2=0,
                asymmetry_ratio=1.0,
                p_value=0.01,
                fdr_p=0.01,
                favored_alt="alt1",
            ),
            dict(
                split=["C", "D"],
                n_concordant=1,
                n_alt1=0,
                n_alt2=4,
                asymmetry_ratio=1.0,
                p_value=0.01,
                fdr_p=0.01,
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


class TestMissingTaxa:
    """Test that taxa present in the species tree but absent from gene trees
    (or vice versa) are handled correctly."""

    @pytest.fixture
    def svc(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return DiscordanceAsymmetry(args)

    def test_species_tree_has_extra_taxon(self, svc, tmp_path):
        """Species tree has taxon 'e' absent from all gene trees.

        Species tree: ((a,e),(b,(c,d)));
        Gene trees use only {a,b,c,d}.

        The (a,e) branch should be skipped (C1 or C2 empty after filtering).
        The (c,d) and (b,(c,d)) branches should still produce correct counts.
        """
        from Bio import Phylo
        from io import StringIO

        species_tree = Phylo.read(StringIO("((a,e),(b,(c,d)));"), "newick")

        # Write gene trees: all concordant with ((a),(b,(c,d)))
        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = [
            "(a,(b,(c,d)));",
            "(a,(b,(c,d)));",
            "(a,(b,(c,d)));",
        ]
        gt_file.write_text("\n".join(gt_lines))
        gene_trees = svc._parse_gene_trees(str(gt_file))

        result, _ = svc._count_topologies(species_tree, gene_trees)

        # The branch grouping (a,e) should be skipped because 'e' is not in
        # gene trees, making C2 empty after filtering
        # We should still see branches for (c,d) and (b,(c,d)) -> complement is (a)
        assert len(result) > 0

        # (c,d) branch should be present and concordant
        assert "c,d" in result
        assert result["c,d"]["n_concordant"] == 3

        # The branch for (a) vs (b,c,d) should NOT appear because the
        # species-tree node grouping (a,e) becomes degenerate after filtering
        # (e is removed, making C2 empty)

    def test_gene_tree_has_extra_taxon(self, svc, tmp_path):
        """Gene trees have taxon 'x' absent from species tree.

        Species tree: (a,(b,(c,d)));
        Gene trees include {a,b,c,d,x}.

        The extra taxon 'x' should be ignored; results should match as if
        gene trees only had {a,b,c,d}.
        """
        from Bio import Phylo
        from io import StringIO

        species_tree = Phylo.read(StringIO("(a,(b,(c,d)));"), "newick")

        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = [
            "((a,x),(b,(c,d)));",
            "((a,x),(b,(c,d)));",
            "((a,x),(b,(c,d)));",
        ]
        gt_file.write_text("\n".join(gt_lines))
        gene_trees = svc._parse_gene_trees(str(gt_file))

        result, _ = svc._count_topologies(species_tree, gene_trees)

        # Extra taxon 'x' in gene trees should be ignored
        assert "c,d" in result
        assert result["c,d"]["n_concordant"] == 3

    def test_degenerate_branch_filtered_out(self, svc, tmp_path):
        """When filtering leaves a group empty, that branch is skipped entirely."""
        from Bio import Phylo
        from io import StringIO

        # Species tree where one child of root has only taxa missing from gene trees
        species_tree = Phylo.read(StringIO("((e,f),(a,(b,(c,d))));"), "newick")

        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = [
            "(a,(b,(c,d)));",
            "(a,(b,(c,d)));",
        ]
        gt_file.write_text("\n".join(gt_lines))
        gene_trees = svc._parse_gene_trees(str(gt_file))

        result, _ = svc._count_topologies(species_tree, gene_trees)

        # (e,f) branch should be entirely skipped (both taxa missing)
        # (a,(b,(c,d))) subtree branches should still work
        assert "c,d" in result
        assert result["c,d"]["n_concordant"] == 2

        # No branch key should reference e or f
        for key in result:
            assert "e" not in key
            assert "f" not in key


class TestBifurcatingRoot:
    """Test that bifurcating root branches get correct NNI decomposition."""

    @pytest.fixture
    def svc(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return DiscordanceAsymmetry(args)

    def test_bifurcating_root_nni_alternatives(self, svc, tmp_path):
        """Bifurcating root ((a,b),(c,(d,e))) with known NNI alternatives.

        The root branch separates {a,b} from {c,d,e}.
        Four subtrees: {a}, {b}, {c}, {d,e}.
        Concordant: {a,b} | {c,d,e}
        NNI alt 1: {a,c} | {b,d,e}   (swap b <-> c)
        NNI alt 2: {a,d,e} | {b,c}   (swap b <-> d,e)
        """
        from Bio import Phylo
        from io import StringIO

        species_tree = Phylo.read(StringIO("((a,b),(c,(d,e)));"), "newick")

        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = [
            "((a,b),(c,(d,e)));",   # concordant
            "((a,b),(c,(d,e)));",   # concordant
            "((a,c),(b,(d,e)));",   # NNI alt: {a,c} | {b,d,e}
            "((a,(d,e)),(b,c));",   # NNI alt: {b,c} | {a,d,e}
        ]
        gt_file.write_text("\n".join(gt_lines))
        gene_trees = svc._parse_gene_trees(str(gt_file))

        result, _ = svc._count_topologies(species_tree, gene_trees)

        assert "a,b" in result
        r = result["a,b"]
        assert r["n_concordant"] == 2
        # Both NNI alternatives should be detected (1 each)
        assert r["n_alt1"] + r["n_alt2"] == 2
        assert r["n_alt1"] == 1
        assert r["n_alt2"] == 1

    def test_bifurcating_root_leaf_sibling_skipped(self, svc, tmp_path):
        """When root is bifurcating and sibling is a leaf, the root branch
        should be skipped (only 3 subtrees, NNI not possible)."""
        from Bio import Phylo
        from io import StringIO

        species_tree = Phylo.read(StringIO("(a,(b,(c,d)));"), "newick")

        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = ["(a,(b,(c,d)));", "(a,(b,(c,d)));"]
        gt_file.write_text("\n".join(gt_lines))
        gene_trees = svc._parse_gene_trees(str(gt_file))

        result, _ = svc._count_topologies(species_tree, gene_trees)

        # The root branch (a | b,c,d) should be skipped because
        # sibling 'a' is a leaf — can't form 4 subtrees for NNI.
        # Only the (c,d) branch should remain.
        assert "c,d" in result
        # The branch_key "a" (root branch from bifurcating root with leaf)
        # should not appear
        assert "a" not in result


class TestSiblingSelection:
    """Test that sibling selection skips siblings with all-filtered taxa."""

    @pytest.fixture
    def svc(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return DiscordanceAsymmetry(args)

    def test_first_sibling_filtered_uses_second(self, svc, tmp_path):
        """When siblings[0] taxa are all absent from gene trees, the code
        should use the next sibling with valid taxa."""
        from Bio import Phylo
        from io import StringIO

        # Trifurcating root: (e,f) has no shared taxa, (a,b) and (c,d) do
        species_tree = Phylo.read(StringIO("((e,f),(a,b),(c,d));"), "newick")

        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = [
            "((a,b),(c,d));",
            "((a,c),(b,d));",
        ]
        gt_file.write_text("\n".join(gt_lines))
        gene_trees = svc._parse_gene_trees(str(gt_file))

        result, _ = svc._count_topologies(species_tree, gene_trees)

        # Branches for (a,b) and (c,d) should NOT be skipped — the code
        # should find a valid sibling past the filtered-out (e,f)
        assert len(result) > 0
        assert "a,b" in result
        assert result["a,b"]["n_concordant"] == 1


class TestPlotWithExtraTaxa:
    """Test that _plot correctly matches branch results when species tree
    has extra taxa not in gene trees."""

    def test_plot_with_extra_species_taxa(self, tmp_path):
        """Plot should color branches correctly even when species tree has
        taxa absent from gene trees."""
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        from Bio import Phylo
        from io import StringIO

        # Create species tree with extra taxon 'e'
        sp_file = tmp_path / "species.tre"
        sp_file.write_text("((a,e),(b,(c,d)));")

        gt_file = tmp_path / "gene_trees.nwk"
        gt_lines = [
            "(a,(b,(c,d)));",
            "(a,(b,(c,d)));",
            "(a,(b,(c,d)));",
        ]
        gt_file.write_text("\n".join(gt_lines))

        args = Namespace(
            tree=str(sp_file), gene_trees=str(gt_file),
            verbose=False, json=False, plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(str(gt_file))
        result, shared_taxa = svc._count_topologies(species_tree, gene_trees)

        # Build branch_results the same way run() does
        branch_results = []
        for branch_key in sorted(result.keys()):
            data = result[branch_key]
            test_result = svc._test_asymmetry(data["n_alt1"], data["n_alt2"])
            entry = dict(
                split=data["split"],
                n_concordant=data["n_concordant"],
                n_alt1=data["n_alt1"],
                n_alt2=data["n_alt2"],
            )
            entry.update(test_result)
            entry["fdr_p"] = None
            branch_results.append(entry)

        # Call _plot with shared_taxa — should not error and should create file
        output = str(tmp_path / "test_plot.png")
        svc._plot(species_tree, branch_results, output,
                  shared_taxa=shared_taxa)
        assert os.path.exists(output)
        assert os.path.getsize(output) > 0

    def test_plot_uses_cached_clade_taxa(self, tmp_path, monkeypatch):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        sp_file = tmp_path / "species.tre"
        sp_file.write_text("((A:1,B:1):1,(C:1,D:1):1);")
        gt_file = tmp_path / "gene_trees.nwk"
        gt_file.write_text("((A:1,B:1):1,(C:1,D:1):1);")
        args = Namespace(
            tree=str(sp_file), gene_trees=str(gt_file),
            verbose=False, annotate=False, json=False, plot_output=None,
        )
        svc = DiscordanceAsymmetry(args)
        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        branch_results = [
            dict(
                split=["A", "B"],
                n_concordant=3,
                n_alt1=1,
                n_alt2=0,
                gcf=0.75,
                asymmetry_ratio=1.0,
                p_value=1.0,
                fdr_p=1.0,
                favored_alt=None,
            )
        ]

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("plot setup should use cached clade taxa")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        output = str(tmp_path / "test_plot_cached.png")
        svc._plot(
            species_tree,
            branch_results,
            output,
            shared_taxa=frozenset({"A", "B", "C", "D"}),
        )

        assert os.path.exists(output)
        assert os.path.getsize(output) > 0


class TestTraversalEquivalence:
    """Verify that _collect_taxa and _extract_splits produce identical
    results to the Bio.Phylo get_terminals / get_nonterminals equivalents."""

    @pytest.fixture
    def svc(self):
        from phykit.services.tree.discordance_asymmetry import DiscordanceAsymmetry

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return DiscordanceAsymmetry(args)

    def test_collect_taxa_matches_get_terminals(self, svc):
        """_collect_taxa must return the same taxon names as
        Bio.Phylo get_terminals() across all gene trees."""
        gene_trees = svc._parse_gene_trees(GENE_TREES)

        # Reference: Bio.Phylo get_terminals()
        reference = set()
        for gt in gene_trees:
            for t in gt.get_terminals():
                reference.add(t.name)

        result = svc._collect_taxa(gene_trees)
        assert result == reference

    def test_extract_splits_matches_biopython(self, svc):
        """_extract_splits must produce the same canonical bipartitions as
        the Bio.Phylo get_nonterminals + get_terminals approach."""
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        species_tree = svc.read_tree_file()

        all_taxa_fs = frozenset(
            set(t.name for t in species_tree.get_terminals())
            & set().union(*(
                set(t.name for t in gt.get_terminals()) for gt in gene_trees
            ))
        )

        for gt in gene_trees:
            # Reference: original Bio.Phylo approach
            ref_splits = set()
            for clade in gt.get_nonterminals():
                tips = frozenset(
                    t.name for t in clade.get_terminals()
                    if t.name in all_taxa_fs
                )
                if len(tips) <= 1 or tips == all_taxa_fs:
                    continue
                ref_splits.add(svc._canonical_split(tips, all_taxa_fs))

            # Optimized version
            opt_splits = svc._extract_splits(gt, all_taxa_fs)

            assert opt_splits == ref_splits

    def test_collect_taxa_with_extra_taxa(self, svc, tmp_path):
        """_collect_taxa works correctly when gene trees have taxa
        not present in the species tree."""
        gt_file = tmp_path / "gene_trees.nwk"
        gt_file.write_text("((a,x),(b,(c,d)));\n(a,(b,(c,(d,y))));\n")
        gene_trees = svc._parse_gene_trees(str(gt_file))

        result = svc._collect_taxa(gene_trees)
        assert result == {"a", "b", "c", "d", "x", "y"}

    def test_extract_splits_with_filtered_taxa(self, svc, tmp_path):
        """_extract_splits correctly filters to shared taxa only."""
        from Bio import Phylo
        from io import StringIO

        gt = Phylo.read(StringIO("((a,x),(b,(c,d)));"), "newick")
        all_taxa_fs = frozenset(["a", "b", "c", "d"])

        splits = svc._extract_splits(gt, all_taxa_fs)

        # x is filtered out, so (a,x) clade becomes just {a} (size 1, skipped).
        # Remaining non-trivial clades:
        #   (c,d): tips {c,d}, complement {a,b}, both size 2 -> tiebreak -> canonical {a,b}
        #   (b,(c,d)): tips {b,c,d}, complement {a}, size 1 < 3 -> canonical {a}
        assert frozenset(["a", "b"]) in splits
        assert frozenset(["a"]) in splits
        assert len(splits) == 2


class TestCLI:
    def test_disc_asym_alias(self, capsys):
        from phykit.phykit import Phykit
        sys.argv = [
            "phykit", "disc_asym",
            "-t", "tests/sample_files/tree_simple.tre",
            "-g", "tests/sample_files/gene_trees_simple.nwk",
        ]
        Phykit()
        captured = capsys.readouterr()
        assert "branch" in captured.out
        assert "Summary:" in captured.out

    def test_da_alias(self, capsys):
        from phykit.phykit import Phykit
        sys.argv = [
            "phykit", "da",
            "-t", "tests/sample_files/tree_simple.tre",
            "-g", "tests/sample_files/gene_trees_simple.nwk",
        ]
        Phykit()
        captured = capsys.readouterr()
        assert "bear,raccoon" in captured.out
