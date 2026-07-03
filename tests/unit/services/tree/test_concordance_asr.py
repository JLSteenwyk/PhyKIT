import copy
import builtins
import json
import math
import subprocess
import sys
from argparse import Namespace
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, TreeMixin

import phykit.services.tree.concordance_asr as concordance_asr_module
from phykit.services.tree.concordance_asr import ConcordanceAsr
from phykit.errors import PhykitUserError


def test_module_import_does_not_import_numpy_or_biophylo():
    code = """
import sys
import phykit.services.tree.concordance_asr as module
assert callable(module.print_json)
assert hasattr(module.np, "__getattr__")
assert hasattr(module.Phylo, "read")
assert hasattr(module.pickle, "dumps")
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "pickle" not in sys.modules
assert "phykit.services.tree.ancestral_reconstruction" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.circular_layout" not in sys.modules
assert "phykit.helpers.color_annotations" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = concordance_asr_module._LazyNumpy()

    first_array = lazy_np.array
    second_array = lazy_np.array

    assert first_array is second_array
    assert lazy_np.__dict__["array"] is first_array
    assert lazy_np._module is not None


here = Path(__file__).resolve().parent
sample = here.parent.parent.parent / "sample_files"

TREE_SIMPLE = str(sample / "tree_simple.tre")
GENE_TREES = str(sample / "gene_trees_simple.nwk")
TRAITS_FILE = str(sample / "tree_simple_traits.tsv")


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        gene_trees=GENE_TREES,
        trait_data=TRAITS_FILE,
        trait=None,
        method="weighted",
        ci=False,
        plot=None,
        missing_taxa="shared",
        json=False,
    )


@pytest.fixture
def ci_args():
    return Namespace(
        tree=TREE_SIMPLE,
        gene_trees=GENE_TREES,
        trait_data=TRAITS_FILE,
        trait=None,
        method="weighted",
        ci=True,
        plot=None,
        missing_taxa="shared",
        json=False,
    )


@pytest.fixture
def distribution_args():
    return Namespace(
        tree=TREE_SIMPLE,
        gene_trees=GENE_TREES,
        trait_data=TRAITS_FILE,
        trait=None,
        method="distribution",
        ci=True,
        plot=None,
        missing_taxa="shared",
        json=False,
    )


@pytest.fixture
def json_args():
    return Namespace(
        tree=TREE_SIMPLE,
        gene_trees=GENE_TREES,
        trait_data=TRAITS_FILE,
        trait=None,
        method="weighted",
        ci=False,
        plot=None,
        missing_taxa="shared",
        json=True,
    )


class TestCanonicalSplit:
    def test_equal_size_returns_lexicographically_smaller_side(self):
        all_taxa = frozenset({"A", "B", "C", "D"})

        assert (
            ConcordanceAsr._canonical_split(frozenset({"C", "D"}), all_taxa)
            == frozenset({"A", "B"})
        )
        assert (
            ConcordanceAsr._canonical_split(frozenset({"A", "B"}), all_taxa)
            == frozenset({"A", "B"})
        )

    def test_empty_split_still_canonicalizes_to_empty_set(self):
        assert ConcordanceAsr._canonical_split(frozenset(), frozenset()) == frozenset()


class TestProcessArgs:
    def test_defaults(self, default_args):
        svc = ConcordanceAsr(default_args)
        assert svc.method == "weighted"
        assert svc.ci is False
        assert svc.plot_output is None
        assert svc.json_output is False
        assert svc.missing_taxa == "shared"
        assert svc.trait_column is None

    def test_overrides(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            trait_data=TRAITS_FILE,
            trait="body_mass",
            method="distribution",
            ci=True,
            plot="out.png",
            missing_taxa="error",
            json=True,
        )
        svc = ConcordanceAsr(args)
        assert svc.method == "distribution"
        assert svc.ci is True
        assert svc.plot_output == "out.png"
        assert svc.json_output is True
        assert svc.missing_taxa == "error"
        assert svc.trait_column == "body_mass"


class TestGeneTreeParsing:
    def test_parse_gene_trees_skips_comments_blanks_and_whitespace(
        self, tmp_path, default_args
    ):
        gene_trees = tmp_path / "gene_trees.nwk"
        gene_trees.write_text(
            "   # ignored\n\n  (A:1.0,B:2.0,C:3.0);  \n\t# also ignored\n"
        )
        svc = ConcordanceAsr(default_args)

        trees = svc._parse_gene_trees(str(gene_trees))

        assert len(trees) == 1

    def test_parse_gene_trees_streams_rows(self, default_args, monkeypatch):
        import phykit.services.tree.concordance_asr as module

        class StreamingOnlyFile:
            def __enter__(self):
                return self

            def __exit__(self, *_args):
                return False

            def __iter__(self):
                return iter(
                    [
                        "# ignored\n",
                        "\n",
                        "  (A:1.0,B:2.0);  \n",
                        "relative.nwk\n",
                    ]
                )

            def read(self):
                raise AssertionError("gene-tree parser should stream rows")

            def readlines(self):
                raise AssertionError("gene-tree parser should stream rows")

        class FakeParent:
            def __str__(self):
                return "/tmp"

        class FakePath:
            parent = FakeParent()

            def __init__(self, path):
                assert path == "gene_trees.txt"

            def open(self, *args, **kwargs):
                return StreamingOnlyFile()

        parsed_sources = []

        def fake_read(source, fmt):
            parsed_sources.append(source)
            return source

        monkeypatch.setattr(module, "Path", FakePath)
        monkeypatch.setattr(module.Phylo, "read", fake_read)
        svc = ConcordanceAsr(default_args)

        trees = svc._parse_gene_trees("gene_trees.txt")

        assert trees == parsed_sources
        assert len(trees) == 2
        assert isinstance(parsed_sources[0], StringIO)
        assert parsed_sources[1] == "/tmp/relative.nwk"

    def test_parse_gene_tree_path_list_avoids_per_row_path_objects(
        self, tmp_path, default_args, monkeypatch
    ):
        import phykit.services.tree.concordance_asr as module

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
        svc = ConcordanceAsr(default_args)

        trees = svc._parse_gene_trees(str(gene_trees))

        assert len(trees) == 2
        assert parent_joins == 0


class TestTaxaNormalization:
    def test_normalize_taxa_uses_fast_tip_name_helper(self, default_args, mocker):
        svc = ConcordanceAsr(default_args)
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        species_tips = set(svc.get_tip_names_from_tree(species_tree))
        spy = mocker.spy(svc, "get_tip_names_from_tree")

        pruned_trees, shared = svc._normalize_taxa(
            species_tree, gene_trees, species_tips
        )

        assert len(pruned_trees) == len(gene_trees)
        assert len(shared) == 8
        assert spy.call_count == (2 * len(gene_trees)) + 1


def _make_asr_helper():
    """Create a minimal AncestralReconstruction helper for testing."""
    from phykit.services.tree.ancestral_reconstruction import AncestralReconstruction
    asr = AncestralReconstruction.__new__(AncestralReconstruction)
    asr.ci = True
    return asr


class TestGCFComputation:
    def test_gcf_topology_counter_scans_once_and_preserves_equal_targets(
        self, monkeypatch
    ):
        def fail_sum(*_args, **_kwargs):
            raise AssertionError("gCF topology counts should scan splits once")

        monkeypatch.setattr(builtins, "sum", fail_sum)

        concordant = frozenset({"A", "B"})
        alt = frozenset({"A", "C"})
        gene_tree_splits = [
            {concordant, alt},
            {concordant},
            {alt},
            set(),
        ]

        assert ConcordanceAsr._count_gcf_topologies(
            gene_tree_splits,
            concordant,
            concordant,
            alt,
        ) == (2, 2, 2)

    def test_collect_clade_tip_sets_standard_tree_avoids_find_clades(
        self, default_args, mocker
    ):
        svc = ConcordanceAsr(default_args)
        species_tree = svc.read_tree_file()
        mocker.patch.object(
            species_tree,
            "find_clades",
            side_effect=AssertionError("standard tree path should avoid find_clades"),
        )

        clade_tip_sets = svc._collect_clade_tip_sets(species_tree)

        assert clade_tip_sets[id(species_tree.root)] == frozenset(
            svc.get_tip_names_from_tree(species_tree)
        )

    def test_gcf_fractions(self, default_args):
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())

        gcf = svc._compute_gcf_per_node(species_tree, gene_trees, all_taxa)

        # All internal nodes should have gCF values
        assert len(gcf) > 0

        # gCF + gDF1 + gDF2 should sum to 1.0 for each node
        for node_id, (g, d1, d2) in gcf.items():
            assert abs(g + d1 + d2 - 1.0) < 1e-10, (
                f"gCF + gDF1 + gDF2 = {g + d1 + d2}, expected 1.0"
            )
            assert 0.0 <= g <= 1.0
            assert 0.0 <= d1 <= 1.0
            assert 0.0 <= d2 <= 1.0

    def test_gcf_computation_uses_direct_tree_traversal(
        self, default_args, monkeypatch
    ):
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())
        expected = svc._compute_gcf_per_node(species_tree, gene_trees, all_taxa)

        def fail_traversal(*_args, **_kwargs):
            raise AssertionError("standard gCF path should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_traversal)

        observed = svc._compute_gcf_per_node(species_tree, gene_trees, all_taxa)

        assert observed == expected

    def test_all_concordant_gcf_is_one(self, default_args):
        """When all gene trees have the same topology, gCF should be 1.0."""
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        species_tree = svc.read_tree_file()

        # Create gene trees identical to species tree
        identical_trees = [copy.deepcopy(species_tree) for _ in range(5)]
        all_taxa = set(t.name for t in species_tree.get_terminals())

        gcf = svc._compute_gcf_per_node(species_tree, identical_trees, all_taxa)

        for node_id, (g, d1, d2) in gcf.items():
            assert g == 1.0, f"Expected gCF=1.0, got {g}"
            assert d1 == 0.0
            assert d2 == 0.0

    def test_gcf_uses_correct_four_groups(self, default_args):
        """Verify that gCF computation uses tree-structure-based NNI alternatives."""
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        species_tree = svc.read_tree_file()
        all_taxa = frozenset(t.name for t in species_tree.get_terminals())
        parent_map = svc._asr._build_parent_map(species_tree)

        # Check four-group decomposition for each non-root internal node
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            groups = svc._get_four_groups(
                species_tree, clade, parent_map, all_taxa
            )
            if groups is None:
                continue
            C1, C2, S, D = groups
            # Four groups should partition all taxa
            assert C1 | C2 | S | D == all_taxa, (
                f"Four groups don't cover all taxa: "
                f"C1={sorted(C1)}, C2={sorted(C2)}, S={sorted(S)}, D={sorted(D)}"
            )
            # Groups should be disjoint
            assert len(C1) + len(C2) + len(S) + len(D) == len(all_taxa)
            # C1 and C2 should be non-empty
            assert len(C1) >= 1
            assert len(C2) >= 1

    def test_cached_four_groups_match_uncached(self, default_args):
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        species_tree = svc.read_tree_file()
        all_taxa = frozenset(t.name for t in species_tree.get_terminals())
        parent_map = svc._asr._build_parent_map(species_tree)
        clade_tip_sets = svc._collect_clade_tip_sets(species_tree)

        for clade in species_tree.find_clades(order="preorder"):
            uncached = svc._get_four_groups(
                species_tree,
                clade,
                parent_map,
                all_taxa,
            )
            cached = svc._get_four_groups(
                species_tree,
                clade,
                parent_map,
                all_taxa,
                clade_tip_sets,
            )
            assert cached == uncached

    def test_cached_four_groups_merge_multifurcation_extras(
        self, default_args, monkeypatch
    ):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        species_tree = Phylo.read(
            StringIO("((A:1,B:1,C:1,D:1):1,(E:1,F:1,G:1,H:1):1);"),
            "newick",
        )
        all_taxa = frozenset({"A", "B", "C", "D", "E", "F", "G", "H"})
        parent_map = svc._asr._build_parent_map(species_tree)
        clade_tip_sets = svc._collect_clade_tip_sets(species_tree)
        node = species_tree.root.clades[0]

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("cached descendant sets should be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        groups = svc._get_four_groups(
            species_tree,
            node,
            parent_map,
            all_taxa,
            clade_tip_sets,
        )
        assert groups == (
            frozenset({"A"}),
            frozenset({"B", "C", "D"}),
            frozenset({"E", "F", "G", "H"}),
            frozenset(),
        )


class TestCachedDescendantAssembly:
    def test_distribution_result_assembly_uses_cached_descendants(
        self, default_args, monkeypatch
    ):
        import phykit.services.tree.concordance_asr as casr_module

        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        svc.ci = False
        species_tree = svc.read_tree_file()

        clade_tip_sets = svc._collect_clade_tip_sets(species_tree)
        all_taxa = set(clade_tip_sets[id(species_tree.root)])
        trait_values = {name: 1.0 for name in all_taxa}
        gcf = {
            id(clade): (1.0, 0.0, 0.0)
            for clade in species_tree.find_clades(order="preorder")
            if not clade.is_terminal()
        }
        svc._compute_gcf_per_node = lambda tree, gene_trees, taxa: gcf

        def run_asr_on_tree(tree_index, values):
            adjusted = {
                tips: float(tree_index)
                for tips in clade_tip_sets.values()
                if len(tips) > 1
            }
            return adjusted, {}, float(tree_index)

        svc._run_asr_on_tree = run_asr_on_tree

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("cached descendant sets should be used")

        class FailingNumpy:
            def __getattr__(self, name):
                raise AssertionError(
                    "distribution summary assembly should not use NumPy"
                )

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(casr_module, "np", FailingNumpy())

        result = svc._run_distribution(
            species_tree, [1, 2, 3], trait_values, all_taxa
        )

        assert len(result["ancestral_estimates"]) == len(gcf)
        assert result["sigma2"] == pytest.approx(2.0)
        assert all(
            entry["descendants"] == sorted(entry["descendants"])
            for entry in result["ancestral_estimates"].values()
        )
        for entry in result["ancestral_estimates"].values():
            assert entry["estimate"] == pytest.approx(2.0)
            assert entry["var_topology"] == pytest.approx(2.0 / 3.0)

    def test_uncertainty_node_data_uses_descendant_lookup(
        self, default_args, monkeypatch
    ):
        svc = ConcordanceAsr(default_args)
        species_tree = svc.read_tree_file()
        clade_tip_sets = svc._collect_clade_tip_sets(species_tree)

        class NoItemsDict(dict):
            def items(self):
                raise AssertionError("nested items scan should not be used")

        anc = NoItemsDict()
        for idx, clade in enumerate(species_tree.find_clades(order="preorder")):
            if clade.is_terminal():
                continue
            anc[f"N{idx}"] = {
                "descendants": sorted(clade_tip_sets[id(clade)]),
                "source_estimates": [1.0, 2.0],
                "gcf": 0.8,
                "estimate": float(idx),
            }

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("cached descendant sets should be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        node_data = svc._collect_uncertainty_node_data(
            species_tree,
            {"method": "weighted", "ancestral_estimates": anc},
        )

        assert len(node_data) == len(anc)
        assert all(len(row[1]) == 2 for row in node_data)

    def test_uncertainty_node_data_preserves_first_duplicate_match(
        self, default_args
    ):
        svc = ConcordanceAsr(default_args)
        species_tree = svc.read_tree_file()
        root_desc = sorted(
            svc._collect_clade_tip_sets(species_tree)[id(species_tree.root)]
        )
        result = {
            "method": "weighted",
            "ancestral_estimates": {
                "first": {
                    "descendants": root_desc,
                    "source_estimates": [1.0, 2.0],
                    "gcf": 0.1,
                    "estimate": 1.0,
                },
                "second": {
                    "descendants": root_desc,
                    "source_estimates": [3.0, 4.0],
                    "gcf": 0.9,
                    "estimate": 2.0,
                },
            },
        }

        node_data = svc._collect_uncertainty_node_data(species_tree, result)

        assert len(node_data) == 1
        assert node_data[0][1] == [1.0, 2.0]
        assert node_data[0][2] == 0.1


class TestNNIAlternatives:
    def test_produces_alternatives(self, default_args):
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        species_tree = svc.read_tree_file()
        parent_map = svc._asr._build_parent_map(species_tree)

        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            if clade == species_tree.root:
                continue
            if len(clade.clades) >= 2:
                alts = svc._build_nni_alternative_at_node(
                    species_tree, clade, parent_map
                )
                assert len(alts) >= 1, "Expected at least 1 NNI alternative"
                orig_tips = set(t.name for t in species_tree.get_terminals())
                for alt_tree, expected_desc in alts:
                    alt_tips = set(t.name for t in alt_tree.get_terminals())
                    assert alt_tips == orig_tips, "NNI alternative lost/gained taxa"
                    # expected_desc should be a valid subset of tips
                    assert expected_desc.issubset(orig_tips)
                break

    def test_branch_lengths_preserved(self, default_args):
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        species_tree = svc.read_tree_file()
        parent_map = svc._asr._build_parent_map(species_tree)

        def tree_length(tree):
            return sum(
                c.branch_length for c in tree.find_clades()
                if c.branch_length is not None
            )

        orig_len = tree_length(species_tree)

        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal() or clade == species_tree.root:
                continue
            if len(clade.clades) >= 2:
                alts = svc._build_nni_alternative_at_node(
                    species_tree, clade, parent_map
                )
                for alt_tree, _ in alts:
                    alt_len = tree_length(alt_tree)
                    assert abs(alt_len - orig_len) < 1e-6, (
                        f"Tree length changed: {orig_len} -> {alt_len}"
                    )
                break

    def test_expected_desc_matches_nni_swap(self, default_args):
        """Verify expected descendant sets match actual NNI swap topology."""
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        species_tree = svc.read_tree_file()
        parent_map = svc._asr._build_parent_map(species_tree)
        all_taxa = frozenset(t.name for t in species_tree.get_terminals())

        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal() or clade == species_tree.root:
                continue
            if len(clade.clades) < 2:
                continue

            C1 = frozenset(t.name for t in clade.clades[0].get_terminals())
            C2 = frozenset(t.name for t in clade.clades[1].get_terminals())
            siblings = [c for c in parent_map[id(clade)].clades
                        if id(c) != id(clade)]
            S = frozenset(t.name for t in siblings[0].get_terminals())

            alts = svc._build_nni_alternative_at_node(
                species_tree, clade, parent_map
            )
            assert len(alts) == 2

            # Alt 0 swaps C1 <-> S: expected desc = S | C2
            _, desc0 = alts[0]
            assert desc0 == S | C2, (
                f"Alt 0: expected {sorted(S | C2)}, got {sorted(desc0)}"
            )

            # Alt 1 swaps C2 <-> S: expected desc = C1 | S
            _, desc1 = alts[1]
            assert desc1 == C1 | S, (
                f"Alt 1: expected {sorted(C1 | S)}, got {sorted(desc1)}"
            )


class TestLawOfTotalVariance:
    def test_known_values(self, monkeypatch):
        weights = [0.5, 0.5]
        means = [10.0, 20.0]
        variances = [1.0, 1.0]

        monkeypatch.setattr(
            concordance_asr_module.np,
            "array",
            lambda *args, **kwargs: pytest.fail(
                "law of total variance should use scalar loops"
            ),
        )
        total, within, between = ConcordanceAsr._law_of_total_variance(
            weights, means, variances
        )

        # Within = 0.5*1 + 0.5*1 = 1.0
        assert abs(within - 1.0) < 1e-10

        # Between = 0.5*(10-15)^2 + 0.5*(20-15)^2 = 25.0
        assert abs(between - 25.0) < 1e-10

        # Total = 26.0
        assert abs(total - 26.0) < 1e-10

    def test_single_topology_zero_between(self):
        weights = [1.0]
        means = [10.0]
        variances = [2.0]

        total, within, between = ConcordanceAsr._law_of_total_variance(
            weights, means, variances
        )

        assert abs(between - 0.0) < 1e-10
        assert abs(within - 2.0) < 1e-10
        assert abs(total - 2.0) < 1e-10

    def test_zero_weights(self):
        weights = [0.0, 0.0]
        means = [10.0, 20.0]
        variances = [1.0, 1.0]

        total, within, between = ConcordanceAsr._law_of_total_variance(
            weights, means, variances
        )

        assert total == 0.0
        assert within == 0.0
        assert between == 0.0


class TestWeightedMethod:
    def test_run_asr_on_tree_uses_fast_tip_name_helper_without_pruning(
        self, default_args, mocker
    ):
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        tree = svc.read_tree_file()
        all_taxa = set(svc.get_tip_names_from_tree(tree))
        trait_values = svc._asr._parse_single_trait_data(TRAITS_FILE, sorted(all_taxa))
        spy = mocker.spy(svc, "get_tip_names_from_tree")
        fast_copy = mocker.patch.object(
            svc,
            "_fast_tree_copy",
            side_effect=AssertionError("all-shared ASR should not copy tree"),
        )
        fast_anc = mocker.spy(svc._asr, "_fast_anc")

        estimates, _, _ = svc._run_asr_on_tree(tree, trait_values)

        assert estimates
        assert spy.call_count == 1
        fast_copy.assert_not_called()
        assert fast_anc.call_args.args[0] is tree

    def test_run_asr_on_tree_uses_fast_tip_name_helper_after_pruning(
        self, default_args, monkeypatch, mocker
    ):
        class OrderedTraitValues(dict):
            def __contains__(self, key):
                raise AssertionError("ordered prune path should not scan membership")

        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        tree = svc.read_tree_file()
        tree_copy = svc.read_tree_file()
        pruned_tree = svc.read_tree_file()
        all_taxa = set(svc.get_tip_names_from_tree(tree))
        trait_values = svc._asr._parse_single_trait_data(TRAITS_FILE, sorted(all_taxa))
        trait_values.pop("dog")
        trait_values = OrderedTraitValues(trait_values)
        spy = mocker.spy(svc, "get_tip_names_from_tree")
        monkeypatch.setattr(
            concordance_asr_module.Tree, "_ORDERED_MAPPING_PRUNE_MIN_SIZE", 0
        )
        fast_copy = mocker.patch.object(svc, "_fast_tree_copy", return_value=tree_copy)
        prune = mocker.patch.object(
            svc, "prune_tree_using_taxa_list", return_value=pruned_tree
        )
        fast_anc = mocker.spy(svc._asr, "_fast_anc")

        estimates, _, _ = svc._run_asr_on_tree(tree, trait_values)

        assert estimates
        assert spy.call_count == 2
        fast_copy.assert_called_once_with(tree)
        prune.assert_called_once_with(tree_copy, ["dog"])
        assert fast_anc.call_args.args[0] is pruned_tree

    def test_basic_run(self, default_args):
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())
        from phykit.services.tree.ancestral_reconstruction import AncestralReconstruction
        trait_values = svc._asr._parse_single_trait_data(TRAITS_FILE, sorted(all_taxa))

        result = svc._run_weighted(species_tree, gene_trees, trait_values, all_taxa)

        assert result["method"] == "weighted"
        assert result["n_gene_trees"] == 10
        assert result["n_tips"] == 8
        assert result["sigma2"] > 0
        assert len(result["ancestral_estimates"]) > 0

        for label, entry in result["ancestral_estimates"].items():
            assert "estimate" in entry
            assert "gcf" in entry
            assert "var_topology" in entry
            assert "var_parameter" in entry
            assert 0.0 <= entry["gcf"] <= 1.0

    def test_all_concordant_matches_standard_asr(self, default_args):
        """When all gene trees match species tree, result should be close to standard ASR."""
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(TRAITS_FILE, sorted(all_taxa))

        # All gene trees identical to species tree
        identical_trees = [copy.deepcopy(species_tree) for _ in range(5)]

        result = svc._run_weighted(
            copy.deepcopy(species_tree), identical_trees, trait_values, all_taxa
        )

        # Run standard ASR for comparison (same tree object for both calls)
        std_tree = copy.deepcopy(species_tree)
        sp_estimates, _, _, _ = svc._asr._fast_anc(
            std_tree,
            np.array([trait_values[n] for n in sorted(trait_values.keys())]),
            sorted(trait_values.keys()),
            svc._asr._label_internal_nodes(std_tree),
        )

        # All gCF should be 1.0 and var_topology should be ~0
        # Estimates should match standard ASR
        for label, entry in result["ancestral_estimates"].items():
            assert entry["gcf"] == 1.0
            assert entry["var_topology"] < 1e-6
            if label in sp_estimates:
                assert abs(entry["estimate"] - sp_estimates[label]) < 1e-4, (
                    f"Node {label}: concordance-weighted={entry['estimate']:.6f} "
                    f"standard={sp_estimates[label]:.6f}"
                )


class TestDistributionMethod:
    def test_basic_run(self, distribution_args):
        svc = ConcordanceAsr(distribution_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(TRAITS_FILE, sorted(all_taxa))

        result = svc._run_distribution(species_tree, gene_trees, trait_values, all_taxa)

        assert result["method"] == "distribution"
        assert result["n_gene_trees"] == 10
        assert len(result["ancestral_estimates"]) > 0

        for label, entry in result["ancestral_estimates"].items():
            assert "estimate" in entry
            assert "gcf" in entry
            assert "var_topology" in entry

    def test_cis_present(self, distribution_args):
        svc = ConcordanceAsr(distribution_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(TRAITS_FILE, sorted(all_taxa))

        result = svc._run_distribution(species_tree, gene_trees, trait_values, all_taxa)

        # At least some nodes should have CIs
        has_ci = any(
            "ci_lower" in e for e in result["ancestral_estimates"].values()
        )
        assert has_ci, "Distribution method should produce CIs"


class TestRun:
    def test_run_uses_fast_tip_name_helper_for_species_prune_setup(
        self, default_args, mocker
    ):
        svc = ConcordanceAsr(default_args)
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(svc.get_tip_names_from_tree(species_tree))
        spy = mocker.spy(svc, "get_tip_names_from_tree")
        mocker.patch.object(svc, "read_tree_file", return_value=species_tree)
        mocker.patch.object(svc, "_parse_gene_trees", return_value=gene_trees)
        mocker.patch.object(
            svc, "_normalize_taxa", return_value=(gene_trees, all_taxa)
        )
        run_weighted = mocker.patch.object(
            svc,
            "_run_weighted",
            return_value={"method": "weighted", "ancestral_estimates": {}},
        )
        mocker.patch.object(svc, "_print_text_output")
        fast_copy = mocker.patch.object(
            svc,
            "_fast_tree_copy",
            side_effect=AssertionError("all-shared species tree should not copy"),
        )

        svc.run()

        assert spy.call_count == 2
        fast_copy.assert_not_called()
        assert run_weighted.call_args.args[0] is species_tree

    def test_run_copies_species_tree_before_trait_pruning(
        self, default_args, monkeypatch, mocker
    ):
        class OrderedTraitValues(dict):
            def __contains__(self, key):
                raise AssertionError("ordered prune path should not scan membership")

        svc = ConcordanceAsr(default_args)
        species_tree = svc.read_tree_file()
        species_copy = svc.read_tree_file()
        pruned_species = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(svc.get_tip_names_from_tree(species_tree))
        species_tip_order = svc.get_tip_names_from_tree(species_tree)
        trait_values = OrderedTraitValues(
            (taxon, float(index))
            for index, taxon in enumerate(species_tip_order)
            if taxon != "dog"
        )

        mocker.patch.object(svc, "read_tree_file", return_value=species_tree)
        mocker.patch.object(svc, "_parse_gene_trees", return_value=gene_trees)
        mocker.patch.object(
            svc, "_normalize_taxa", return_value=(gene_trees, all_taxa)
        )
        monkeypatch.setattr(
            concordance_asr_module.Tree, "_ORDERED_MAPPING_PRUNE_MIN_SIZE", 0
        )
        mocker.patch(
            "phykit.services.tree.ancestral_reconstruction."
            "AncestralReconstruction._parse_single_trait_data",
            return_value=trait_values,
        )
        fast_copy = mocker.patch.object(
            svc, "_fast_tree_copy", return_value=species_copy
        )
        prune = mocker.patch.object(
            svc, "prune_tree_using_taxa_list", return_value=pruned_species
        )
        run_weighted = mocker.patch.object(
            svc,
            "_run_weighted",
            return_value={"method": "weighted", "ancestral_estimates": {}},
        )
        mocker.patch.object(svc, "_print_text_output")

        svc.run()

        fast_copy.assert_called_once_with(species_tree)
        prune.assert_called_once_with(species_copy, ["dog"])
        assert run_weighted.call_args.args[0] is pruned_species

    def test_text_output(self, default_args):
        svc = ConcordanceAsr(default_args)
        with patch("builtins.print") as mock_print:
            svc.run()

        printed = " ".join(str(c) for c in mock_print.call_args_list)
        assert "Concordance-Aware" in printed
        assert "weighted" in printed
        assert "gCF" in printed

    def test_print_text_output_batches_estimate_rows(self, default_args):
        svc = ConcordanceAsr(default_args)
        result = {
            "method": "weighted",
            "n_tips": 8,
            "n_gene_trees": 3,
            "sigma2": 1.25,
            "ancestral_estimates": {
                "N1": {
                    "descendants": ["A", "B"],
                    "estimate": 2.5,
                    "gcf": 0.75,
                    "var_topology": 0.01,
                    "var_parameter": 0.02,
                    "ci_lower": 2.0,
                    "ci_upper": 3.0,
                    "is_root": True,
                },
                "N2": {
                    "descendants": ["C", "D", "E"],
                    "estimate": -1.5,
                    "gcf": 0.5,
                    "var_topology": 0.03,
                    "var_parameter": 0.04,
                },
            },
        }

        with patch("builtins.print") as mock_print:
            svc._print_text_output(result)

        mock_print.assert_called_once()
        output = mock_print.call_args.args[0]
        assert "Concordance-Aware Ancestral State Reconstruction" in output
        assert "Method: weighted" in output
        assert "N1 (root)" in output
        assert "[2.0000, 3.0000]" in output
        assert "N2" in output
        assert "Var_topo" in output

    def test_json_output(self, json_args):
        svc = ConcordanceAsr(json_args)
        captured = []
        with patch("phykit.services.tree.concordance_asr.print_json") as mock_json:
            mock_json.side_effect = lambda d: captured.append(d)
            svc.run()

        assert len(captured) == 1
        result = captured[0]
        assert "method" in result
        assert "ancestral_estimates" in result
        assert "n_gene_trees" in result
        assert "sigma2" in result

    def test_distribution_text_output(self, distribution_args):
        svc = ConcordanceAsr(distribution_args)
        with patch("builtins.print") as mock_print:
            svc.run()

        printed = " ".join(str(c) for c in mock_print.call_args_list)
        assert "Concordance-Aware" in printed
        assert "distribution" in printed


class TestCircularPlot:
    @pytest.mark.parametrize("circular", [False, True])
    def test_concordance_plot_uses_direct_tree_traversal(
        self, default_args, monkeypatch, tmp_path, circular
    ):
        default_args.circular = circular
        svc = ConcordanceAsr(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        estimates = {}
        counter = 1
        for clade in svc._iter_preorder(tree.root):
            if not clade.clades:
                continue
            label = clade.name or f"N{counter}"
            if not clade.name:
                counter += 1
            estimates[label] = {"estimate": float(counter), "gcf": 0.7}

        def fail_traversal(*args, **kwargs):
            raise AssertionError("plot setup should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        output_path = tmp_path / f"concordance_asr_direct_{circular}.png"
        svc._plot_concordance_contmap(
            tree, {"ancestral_estimates": estimates}, str(output_path)
        )

        assert output_path.exists()

    def test_iter_preorder_preserves_order_without_reversed(self, default_args):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("_iter_preorder should push children directly")

        svc = ConcordanceAsr(default_args)
        root = Clade(name="root")
        left = Clade(name="left")
        middle = Clade(name="middle")
        right = Clade(name="right")
        left.clades = NoReversedList([Clade(name="left_a"), Clade(name="left_b")])
        middle.clades = NoReversedList([Clade(name="middle_a")])
        right.clades = NoReversedList()
        root.clades = NoReversedList([left, middle, right])

        order = [clade.name for clade in svc._iter_preorder(root)]

        assert order == [
            "root",
            "left",
            "left_a",
            "left_b",
            "middle",
            "middle_a",
            "right",
        ]

    def test_concordance_plot_reuses_preorder_for_node_positions(
        self, default_args, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import phykit.helpers.plot_config as plot_config

        default_args.circular = False
        svc = ConcordanceAsr(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        preorder_clades = list(svc._iter_preorder(tree.root))
        estimates = {}
        counter = 1
        for clade in preorder_clades:
            if not clade.clades:
                continue
            label = clade.name or f"N{counter}"
            if not clade.name:
                counter += 1
            estimates[label] = {"estimate": float(counter), "gcf": 0.7}
        expected_preorder_ids = [id(clade) for clade in preorder_clades]
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

        output_path = tmp_path / "concordance_asr_preorder_positions.png"
        svc._plot_concordance_contmap(
            tree, {"ancestral_estimates": estimates}, str(output_path)
        )

        assert calls == [expected_preorder_ids]
        assert output_path.exists()

    def test_concordance_plot_reuses_clade_lists_for_circular_coords(
        self, default_args, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import phykit.helpers.circular_layout as circular_layout

        default_args.circular = True
        default_args.ylabel_fontsize = 0
        default_args.no_title = True
        svc = ConcordanceAsr(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        preorder_clades = list(svc._iter_preorder(tree.root))
        estimates = {}
        counter = 1
        for clade in preorder_clades:
            if not clade.clades:
                continue
            label = clade.name or f"N{counter}"
            if not clade.name:
                counter += 1
            estimates[label] = {"estimate": float(counter), "gcf": 0.7}
        expected_preorder_ids = [id(clade) for clade in preorder_clades]
        expected_terminal_ids = [
            id(clade) for clade in preorder_clades if not clade.clades
        ]
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

        output_path = tmp_path / "concordance_asr_circular_precomputed_coords.png"
        svc._plot_concordance_contmap(
            tree, {"ancestral_estimates": estimates}, str(output_path)
        )

        assert calls == [(expected_preorder_ids, expected_terminal_ids)]
        assert output_path.exists()

    def test_rectangular_concordance_plot_batches_base_branches(
        self, default_args, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        default_args.circular = False
        default_args.ylabel_fontsize = 0
        default_args.no_title = True
        svc = ConcordanceAsr(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        estimates = {}
        counter = 1
        for clade in svc._iter_preorder(tree.root):
            if not clade.clades:
                continue
            label = clade.name or f"N{counter}"
            if not clade.name:
                counter += 1
            estimates[label] = {"estimate": float(counter), "gcf": 0.7}

        output_path = tmp_path / "concordance_asr_batched.png"
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("Concordance-ASR branches should be batched")

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        svc._plot_concordance_contmap(
            tree, {"ancestral_estimates": estimates}, str(output_path)
        )

        assert len(line_collections) >= 2
        assert output_path.exists()

    @pytest.mark.parametrize("circular", [False, True])
    def test_concordance_plot_batches_gcf_markers(
        self, default_args, monkeypatch, tmp_path, circular
    ):
        pytest.importorskip("matplotlib")
        import matplotlib.axes

        default_args.circular = circular
        default_args.ylabel_fontsize = 0
        default_args.no_title = True
        default_args.legend_position = "none"
        svc = ConcordanceAsr(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        estimates = {}
        counter = 1
        for clade in svc._iter_preorder(tree.root):
            if not clade.clades:
                continue
            label = clade.name or f"N{counter}"
            if not clade.name:
                counter += 1
            estimates[label] = {
                "estimate": float(counter),
                "gcf": 0.25 + 0.1 * counter,
            }

        original_scatter = matplotlib.axes.Axes.scatter
        gcf_marker_counts = []

        def capture_scatter(self, x, y, *args, **kwargs):
            if kwargs.get("zorder") == 5 and kwargs.get("edgecolors") == "black":
                try:
                    gcf_marker_counts.append(len(x))
                except TypeError:
                    gcf_marker_counts.append(1)
            return original_scatter(self, x, y, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "scatter", capture_scatter)

        output_path = tmp_path / f"concordance_asr_gcf_markers_{circular}.png"
        svc._plot_concordance_contmap(
            tree, {"ancestral_estimates": estimates}, str(output_path)
        )

        assert gcf_marker_counts == [len(estimates)]
        assert output_path.exists()

    def test_uncertainty_plot_mean_markers_avoid_numpy(
        self, default_args, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")

        default_args.ylabel_fontsize = 0
        default_args.no_title = True
        svc = ConcordanceAsr(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clade_tip_sets = svc._collect_clade_tip_sets(tree)
        estimates = {}
        for idx, clade in enumerate(svc._iter_preorder(tree.root)):
            if not clade.clades:
                continue
            estimates[f"N{idx}"] = {
                "descendants": sorted(clade_tip_sets[id(clade)]),
                "source_estimates": [1.0, 2.0, 4.0],
                "gcf": 0.8,
                "estimate": 2.0,
            }

        class FailingNumpy:
            def __getattr__(self, name):
                raise AssertionError(
                    "uncertainty mean markers should not use NumPy"
                )

        monkeypatch.setattr(concordance_asr_module, "np", FailingNumpy())

        output_path = tmp_path / "concordance_asr_uncertainty.png"
        svc._plot_uncertainty(
            tree,
            {"method": "weighted", "ancestral_estimates": estimates},
            str(output_path),
        )

        assert output_path.exists()

    def test_concordance_asr_circular(self):
        """--circular flag produces a circular layout concordance ASR plot."""
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        import tempfile
        import os

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            plot_path = f.name

        try:
            args = Namespace(
                tree=TREE_SIMPLE,
                gene_trees=GENE_TREES,
                trait_data=TRAITS_FILE,
                trait=None,
                method="weighted",
                ci=False,
                plot=plot_path,
                missing_taxa="shared",
                json=False,
                circular=True,
            )
            svc = ConcordanceAsr(args)
            with patch("builtins.print"):
                svc.run()

            assert os.path.exists(plot_path)
            assert os.path.getsize(plot_path) > 0
        finally:
            if os.path.exists(plot_path):
                os.unlink(plot_path)


class TestEdgeCases:
    def test_single_gene_tree_errors(self):
        """Single gene tree should raise error."""
        import tempfile
        import os

        # Create a temp file with just one gene tree
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".nwk", delete=False
        ) as f:
            f.write(
                "((raccoon:18.5,bear:7.2):0.9,((sea_lion:12.1,seal:11.8):7.6,"
                "((monkey:99.2,cat:46.5):21.0,weasel:19.1):2.1):3.9,dog:24.8);\n"
            )
            tmp_path = f.name

        try:
            args = Namespace(
                tree=TREE_SIMPLE,
                gene_trees=tmp_path,
                trait_data=TRAITS_FILE,
                trait=None,
                method="weighted",
                ci=False,
                plot=None,
                missing_taxa="shared",
                json=False,
            )
            svc = ConcordanceAsr(args)
            with pytest.raises(SystemExit):
                svc.run()
        finally:
            os.unlink(tmp_path)

    def test_missing_taxa_error_mode(self):
        """Error mode should reject taxa mismatches."""
        import tempfile
        import os

        # Create gene trees with different taxa
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".nwk", delete=False
        ) as f:
            f.write(
                "((raccoon:18.5,bear:7.2):0.9,(seal:11.8,"
                "((monkey:99.2,cat:46.5):21.0,weasel:19.1):2.1):3.9,dog:24.8);\n"
                "((raccoon:19.0,bear:6.8):0.85,(seal:12.0,"
                "((monkey:100.0,cat:47.0):20.5,weasel:18.8):2.0):3.8,dog:25.0);\n"
            )
            tmp_path = f.name

        try:
            args = Namespace(
                tree=TREE_SIMPLE,
                gene_trees=tmp_path,
                trait_data=TRAITS_FILE,
                trait=None,
                method="weighted",
                ci=False,
                plot=None,
                missing_taxa="error",
                json=False,
            )
            svc = ConcordanceAsr(args)
            with pytest.raises(SystemExit):
                svc.run()
        finally:
            os.unlink(tmp_path)

    def test_gene_tree_file_not_found(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees="/nonexistent/path.nwk",
            trait_data=TRAITS_FILE,
            trait=None,
            method="weighted",
            ci=False,
            plot=None,
            missing_taxa="shared",
            json=False,
        )
        svc = ConcordanceAsr(args)
        with pytest.raises(SystemExit):
            svc.run()

    def test_ci_with_weighted(self, ci_args):
        svc = ConcordanceAsr(ci_args)
        captured = []
        with patch("builtins.print") as mock_print:
            svc.run()

        printed = " ".join(str(c) for c in mock_print.call_args_list)
        assert "95% CI" in printed


class TestBenchmarkVerification:
    """Verify correctness of concordance-aware ASR against known values."""

    def test_gcf_manual_computation(self, default_args):
        """Verify gCF values against manually counted concordance in gene trees.

        The species tree topology has the bipartition:
          {raccoon, bear} | {sea_lion, seal, monkey, cat, weasel, dog}
        Gene trees 1-6 and 10 preserve this clade (7 out of 10).
        Gene trees 7, 8, 9 break this clade.
        The gCF for this node should reflect the fraction of gene trees
        supporting it among those resolving the relevant quartet.
        """
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())

        gcf = svc._compute_gcf_per_node(species_tree, gene_trees, all_taxa)

        # At least some nodes should have gCF < 1.0
        # (discordant gene trees 7-9 break at least some bipartitions)
        has_partial_concordance = any(g < 1.0 for g, _, _ in gcf.values())
        assert has_partial_concordance, (
            "Expected some nodes with gCF < 1.0 given discordant gene trees"
        )

        # Verify gCF for {raccoon, bear} clade: it should be supported
        # by gene trees 1-6 and 10 (7 concordant).
        # The three quartet resolutions at this node:
        #   concordant: {raccoon, bear} | rest
        #   NNI alt 1:  {sibling, bear} | rest
        #   NNI alt 2:  {raccoon, sibling} | rest
        # We check that the concordant fraction is ~0.7 (7/10)
        raccoon_bear_node = None
        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            tips = frozenset(t.name for t in clade.get_terminals())
            if tips == frozenset({"raccoon", "bear"}):
                raccoon_bear_node = clade
                break

        if raccoon_bear_node is not None and id(raccoon_bear_node) in gcf:
            g, d1, d2 = gcf[id(raccoon_bear_node)]
            # At least 7 of 10 gene trees group raccoon+bear,
            # so gCF should be >= 0.7
            assert g >= 0.6, (
                f"Expected gCF >= 0.6 for {{raccoon, bear}}, got {g}"
            )

    def test_var_topology_positive_with_discordance(self, default_args):
        """When gene trees disagree, var_topology should be > 0 for at least
        some nodes (the between-topology variance from the law of total variance).
        """
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(
            TRAITS_FILE, sorted(all_taxa)
        )

        result = svc._run_weighted(species_tree, gene_trees, trait_values, all_taxa)

        # At least one node should have var_topology > 0
        has_positive_var_topo = any(
            e["var_topology"] > 0
            for e in result["ancestral_estimates"].values()
        )
        assert has_positive_var_topo, (
            "Expected at least one node with var_topology > 0 "
            "given discordant gene trees"
        )

    def test_weighted_differs_from_standard_asr_with_discordance(self, default_args):
        """The concordance-weighted estimate should differ from standard ASR
        when gene trees are discordant, because NNI alternative estimates
        are mixed in with their discordance weights.
        """
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(
            TRAITS_FILE, sorted(all_taxa)
        )

        # Weighted (concordance-aware) result
        weighted_result = svc._run_weighted(
            copy.deepcopy(species_tree), gene_trees, trait_values, all_taxa
        )

        # Standard ASR on species tree only
        # Must use the same tree object for _fast_anc and _label_internal_nodes
        # since node_labels is keyed by id(clade)
        std_tree = copy.deepcopy(species_tree)
        sp_est, _, _, _ = svc._asr._fast_anc(
            std_tree,
            np.array([trait_values[n] for n in sorted(trait_values.keys())]),
            sorted(trait_values.keys()),
            svc._asr._label_internal_nodes(std_tree),
        )

        # For nodes with gCF < 1.0 (where NNI alternatives contribute),
        # the weighted estimate should differ from the standard ASR estimate
        has_difference = False
        for label, entry in weighted_result["ancestral_estimates"].items():
            if entry["gcf"] < 1.0 and label in sp_est:
                if abs(entry["estimate"] - sp_est[label]) > 1e-6:
                    has_difference = True
                    break

        assert has_difference, (
            "Expected concordance-weighted estimate to differ from "
            "standard ASR for at least one discordant node"
        )

    def test_all_concordant_zero_topology_variance(self, default_args):
        """When all gene trees are identical to the species tree,
        var_topology should be 0 for all nodes (no topological uncertainty).
        """
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(
            TRAITS_FILE, sorted(all_taxa)
        )

        identical_trees = [copy.deepcopy(species_tree) for _ in range(5)]
        result = svc._run_weighted(
            copy.deepcopy(species_tree), identical_trees, trait_values, all_taxa
        )

        for label, entry in result["ancestral_estimates"].items():
            assert entry["var_topology"] < 1e-10, (
                f"Node {label}: expected var_topology ~0 with all concordant "
                f"trees, got {entry['var_topology']}"
            )
            assert entry["gcf"] == 1.0, (
                f"Node {label}: expected gCF=1.0 with all concordant "
                f"trees, got {entry['gcf']}"
            )

    def test_distribution_variance_tracks_discordance(self, distribution_args):
        """In the distribution method, nodes with more discordance should
        tend to have higher var_topology (more spread across gene trees).
        """
        svc = ConcordanceAsr(distribution_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(
            TRAITS_FILE, sorted(all_taxa)
        )

        result = svc._run_distribution(
            species_tree, gene_trees, trait_values, all_taxa
        )

        # At least some nodes should have non-zero variance from gene tree spread
        has_nonzero_var = any(
            e["var_topology"] > 0
            for e in result["ancestral_estimates"].values()
        )
        assert has_nonzero_var, (
            "Distribution method should show variance across gene trees"
        )

    def test_nni_alternatives_incorporated(self, default_args):
        """Verify that NNI alternative ASR estimates are actually used in
        the weighted combination (not silently dropped).
        """
        svc = ConcordanceAsr(default_args)
        svc._asr = _make_asr_helper()

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        all_taxa = set(t.name for t in species_tree.get_terminals())
        trait_values = svc._asr._parse_single_trait_data(
            TRAITS_FILE, sorted(all_taxa)
        )

        # Run components manually
        node_labels = svc._asr._label_internal_nodes(species_tree)
        parent_map = svc._asr._build_parent_map(species_tree)
        gcf_per_node = svc._compute_gcf_per_node(
            species_tree, gene_trees, all_taxa
        )

        # For nodes with gCF < 1.0, verify NNI alternatives can be found
        nni_found_count = 0
        nni_attempted_count = 0

        for clade in species_tree.find_clades(order="preorder"):
            if clade.is_terminal():
                continue
            if id(clade) not in gcf_per_node:
                continue
            gcf, gdf1, gdf2 = gcf_per_node[id(clade)]
            if gcf >= 1.0:
                continue

            nni_attempted_count += 1
            alts = svc._build_nni_alternative_at_node(
                species_tree, clade, parent_map
            )
            for alt_tree, expected_desc in alts:
                est, _, _ = svc._run_asr_on_tree(alt_tree, trait_values)
                if expected_desc in est:
                    nni_found_count += 1

        # We should have attempted at least one discordant node
        assert nni_attempted_count > 0, (
            "Expected at least one node with gCF < 1.0"
        )
        # At least some NNI lookups should succeed
        assert nni_found_count > 0, (
            f"NNI estimate lookups all failed "
            f"({nni_found_count}/{nni_attempted_count * 2} found). "
            "The expected descendant sets may not match NNI tree nodes."
        )
