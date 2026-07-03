import builtins
import importlib
import os
import subprocess
import sys
import tempfile

import numpy as np
import pytest
from mock import patch
from argparse import Namespace

from phykit.errors import PhykitUserError

TREE_SIMPLE = "tests/sample_files/tree_simple.tre"
GENE_TREES = "tests/sample_files/gene_trees_simple.nwk"


def test_module_import_does_not_import_numpy_biophylo_or_scipy():
    code = """
import builtins
import sys
module_name = "phykit.services.tree.evo_tempo_map"
sys.modules.pop(module_name, None)
sys.modules.pop("pathlib", None)
original_import = builtins.__import__
def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
    if name == "pathlib" or name.startswith("pathlib."):
        raise AssertionError("evo_tempo_map import should not import pathlib")
    return original_import(name, globals, locals, fromlist, level)
builtins.__import__ = guarded_import
try:
    import phykit.services.tree.evo_tempo_map as module
finally:
    builtins.__import__ = original_import
assert callable(module.print_json)
assert callable(module.Path)
assert hasattr(module.np, "__getattr__")
assert hasattr(module.Phylo, "read")
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "pathlib" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "scipy.stats" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_module_import_does_not_import_scipy_stats(monkeypatch):
    module_name = "phykit.services.tree.evo_tempo_map"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.stats" or name.startswith("scipy.stats."):
            raise AssertionError(
                "evo_tempo_map module import should not import scipy.stats"
            )
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", guarded_import)
    try:
        importlib.import_module(module_name)
    finally:
        imported = sys.modules.pop(module_name, None)
        if previous is not None:
            sys.modules[module_name] = previous
        parent_name, _, child_name = module_name.rpartition(".")
        parent = sys.modules.get(parent_name)
        if parent is not None:
            if previous is not None:
                setattr(parent, child_name, previous)
            elif getattr(parent, child_name, None) is imported:
                delattr(parent, child_name)


class TestProcessArgs:
    def test_default_args(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output=None,
        )
        svc = EvoTempoMap(args)
        assert svc.gene_trees_path == GENE_TREES
        assert svc.verbose is False
        assert svc.json_output is False
        assert svc.plot_output is None

    def test_verbose_arg(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=True,
            json=False,
            plot_output=None,
        )
        svc = EvoTempoMap(args)
        assert svc.verbose is True

    def test_json_arg(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=True,
            plot_output=None,
        )
        svc = EvoTempoMap(args)
        assert svc.json_output is True

    def test_plot_output_arg(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output="/tmp/test.png",
        )
        svc = EvoTempoMap(args)
        assert svc.plot_output == "/tmp/test.png"


class TestCanonicalSplit:
    def test_equal_size_returns_lexicographically_smaller_side(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        all_taxa = frozenset({"A", "B", "C", "D"})

        assert (
            EvoTempoMap._canonical_split(frozenset({"C", "D"}), all_taxa)
            == frozenset({"A", "B"})
        )
        assert (
            EvoTempoMap._canonical_split(frozenset({"A", "B"}), all_taxa)
            == frozenset({"A", "B"})
        )

    def test_empty_split_still_canonicalizes_to_empty_set(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        assert EvoTempoMap._canonical_split(frozenset(), frozenset()) == frozenset()


class TestParseGeneTrees:
    def test_parses_multi_newick(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        trees = svc._parse_gene_trees(GENE_TREES)
        assert len(trees) == 10
        for gt in trees:
            tips = [t.name for t in gt.get_terminals()]
            assert len(tips) == 8

    def test_parse_gene_trees_skips_comments_blanks_and_whitespace(self, tmp_path):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        gene_trees = tmp_path / "gene_trees.nwk"
        gene_trees.write_text(
            "   # ignored\n\n  (A:1.0,B:2.0,C:3.0);  \n\t# also ignored\n"
        )
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=str(gene_trees),
            verbose=False, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)

        trees = svc._parse_gene_trees(str(gene_trees))

        assert len(trees) == 1

    def test_parse_gene_tree_path_list_avoids_per_row_path_objects(
        self, tmp_path, monkeypatch
    ):
        import phykit.services.tree.evo_tempo_map as module
        from pathlib import Path
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

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
        svc = EvoTempoMap(args)

        trees = svc._parse_gene_trees(str(gene_trees))

        assert len(trees) == 2
        assert parent_joins == 0

    def test_file_not_found_error(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees="nonexistent.nwk",
            verbose=False, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        with pytest.raises(PhykitUserError):
            svc._parse_gene_trees("nonexistent.nwk")

    def test_validates_branch_lengths(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        no_bl = "((A,B),(C,D));\n((E,F),(G,H));\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".nwk", delete=False
        ) as f:
            f.write(no_bl)
            tmp_path = f.name
        try:
            args = Namespace(
                tree=TREE_SIMPLE, gene_trees=tmp_path,
                verbose=False, json=False, plot_output=None,
            )
            svc = EvoTempoMap(args)
            with pytest.raises(PhykitUserError) as excinfo:
                svc._parse_and_validate_gene_trees()
            assert "branch lengths" in excinfo.value.messages[1]
        finally:
            os.unlink(tmp_path)

    def test_validates_branch_lengths_with_direct_preorder(self, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees="unused.nwk",
            verbose=False,
            json=False,
            plot_output=None,
        )
        svc = EvoTempoMap(args)
        trees = [
            Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick"),
            Phylo.read(StringIO("((A:2,C:2):2,(B:2,D:2):2);"), "newick"),
        ]

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError(
                "standard gene-tree validation should use direct preorder"
            )

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)
        monkeypatch.setattr(svc, "_parse_gene_trees", lambda _path: trees)

        assert svc._parse_and_validate_gene_trees() == trees

    def test_requires_at_least_2_gene_trees(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        single = "((raccoon:19.2,bear:6.8):0.85,((sea_lion:12.0,seal:12.0):7.5,((monkey:100.9,cat:47.1):20.6,weasel:18.9):2.1):3.9,dog:25.5);\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".nwk", delete=False
        ) as f:
            f.write(single)
            tmp_path = f.name
        try:
            args = Namespace(
                tree=TREE_SIMPLE, gene_trees=tmp_path,
                verbose=False, json=False, plot_output=None,
            )
            svc = EvoTempoMap(args)
            with pytest.raises(PhykitUserError) as excinfo:
                svc._parse_and_validate_gene_trees()
            assert "At least 2" in excinfo.value.messages[0]
        finally:
            os.unlink(tmp_path)


class TestClassifyGeneTrees:
    @pytest.fixture
    def svc(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return EvoTempoMap(args)

    def test_classify_returns_dict_per_branch(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._classify_gene_trees(species_tree, gene_trees)
        assert len(result) > 0
        for branch_label, data in result.items():
            assert "concordant_lengths" in data
            assert "discordant_lengths" in data
            assert "n_concordant" in data
            assert "n_discordant" in data
            assert "split" in data

    def test_concordant_plus_discordant_leq_total(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._classify_gene_trees(species_tree, gene_trees)
        n_gene_trees = len(gene_trees)
        for branch_label, data in result.items():
            total = data["n_concordant"] + data["n_discordant"]
            assert total <= n_gene_trees

    def test_all_concordant_when_identical_trees(self, svc):
        species_tree = svc.read_tree_file()
        identical_trees = [svc.read_tree_file() for _ in range(3)]
        result = svc._classify_gene_trees(species_tree, identical_trees)
        for branch_label, data in result.items():
            assert data["n_concordant"] == 3
            assert data["n_discordant"] == 0

    def test_branch_lengths_are_positive(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._classify_gene_trees(species_tree, gene_trees)
        for branch_label, data in result.items():
            for bl in data["concordant_lengths"]:
                assert bl > 0
            for bl in data["discordant_lengths"]:
                assert bl > 0

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
            frozenset({"C", "D"}),
            frozenset(),
        )

        splits = svc._extract_bipartitions_with_lengths(tree, all_taxa, clade_taxa)
        assert frozenset({"A", "B"}) in splits

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
            frozenset({"E", "F", "G", "H"}),
            frozenset(),
        )

    def test_build_parent_map_handles_mixed_child_counts(self, svc, monkeypatch):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_iter_preorder(*args, **kwargs):
            raise AssertionError("parent map should build directly")

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("standard parent map should not use find_clades")

        monkeypatch.setattr(
            type(svc), "_iter_preorder_clades", staticmethod(fail_iter_preorder)
        )
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

    def test_iter_postorder_clades_matches_biopython_order(self, svc):
        from Bio import Phylo
        from io import StringIO

        tree = Phylo.read(
            StringIO("((A:1,B:2,C:3):4,(D:5,E:6):7,F:8):0;"),
            "newick",
        )

        direct = svc._iter_postorder_clades(tree)
        reference = list(tree.find_clades(order="postorder"))

        assert [id(clade) for clade in direct] == [
            id(clade) for clade in reference
        ]

    def test_iter_preorder_clades_preserves_order_without_reversed(self, svc):
        from Bio.Phylo.BaseTree import Clade, Tree

        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("preorder helper should not call reversed()")

        left = Clade(name="left")
        left.clades = NoReversedList([Clade(name="left_a"), Clade(name="left_b")])
        middle = Clade(name="middle")
        middle.clades = NoReversedList([Clade(name="middle_a")])
        right = Clade(name="right")
        right.clades = NoReversedList()
        root = Clade(name="root")
        root.clades = NoReversedList([left, middle, right])

        direct = svc._iter_preorder_clades(Tree(root=root))

        assert [clade.name for clade in direct] == [
            "root",
            "left",
            "left_a",
            "left_b",
            "middle",
            "middle_a",
            "right",
        ]

    def test_extract_bipartitions_without_cache_uses_direct_postorder_pass(
        self, svc, monkeypatch
    ):
        from Bio import Phylo
        from io import StringIO

        tree = Phylo.read(StringIO("((A:1,B:1):2,(C:1,D:1):3);"), "newick")
        all_taxa = frozenset({"A", "B", "C", "D"})

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("standard trees should use direct traversal")

        def fail_get_nonterminals(*args, **kwargs):
            raise AssertionError("no-cache path should not do a second traversal")

        monkeypatch.setattr(tree, "find_clades", fail_find_clades)
        monkeypatch.setattr(tree, "get_nonterminals", fail_get_nonterminals)

        splits = svc._extract_bipartitions_with_lengths(tree, all_taxa)

        assert splits[frozenset({"A", "B"})] == 3

    def test_classify_shared_taxa_setup_uses_fast_tip_names(
        self, svc, monkeypatch
    ):
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin
        from io import StringIO

        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick"
        )
        gene_trees = [
            Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick"),
            Phylo.read(StringIO("((A:1,C:1):1,(B:1,D:1):1);"), "newick"),
        ]

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("no-prune setup should use fast terminal names")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        result = svc._classify_gene_trees(species_tree, gene_trees)

        assert result

    def test_cross_check_with_concordance_asr_gcf(self, svc):
        """The gCF (fraction concordant) at each branch should match
        concordance_asr's computation."""
        from phykit.services.tree.concordance_asr import ConcordanceAsr
        from phykit.services.tree.ancestral_reconstruction import AncestralReconstruction

        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        etm_result = svc._classify_gene_trees(species_tree, gene_trees)

        # Get gCF from concordance_asr
        casr_args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            trait_data="tests/sample_files/tree_simple_traits.tsv",
            trait=None, method="weighted", ci=False, plot=None,
            missing_taxa="shared", json=False,
        )
        casr = ConcordanceAsr(casr_args)
        # Initialize the _asr helper that _compute_gcf_per_node depends on
        casr._asr = AncestralReconstruction.__new__(AncestralReconstruction)
        casr._asr.ci = True
        casr_species_tree = casr.read_tree_file()
        casr_gene_trees = casr._parse_gene_trees(GENE_TREES)
        all_taxa = sorted(set(
            t.name for t in casr_species_tree.get_terminals()
        ))
        gcf_result = casr._compute_gcf_per_node(
            casr_species_tree, casr_gene_trees, all_taxa
        )

        # For each species tree branch, compute gCF from ETM and compare
        for branch_key, etm_data in etm_result.items():
            n_conc = etm_data["n_concordant"]
            n_disc = etm_data["n_discordant"]
            total = n_conc + n_disc
            if total == 0:
                continue
            etm_gcf = n_conc / total

            # Find matching node in concordance_asr by descendant set
            split_set = frozenset(etm_data["split"])
            matched = False
            for node_id, (gcf, gdf1, gdf2) in gcf_result.items():
                for clade in casr_species_tree.find_clades():
                    if id(clade) == node_id:
                        node_tips = frozenset(
                            t.name for t in clade.get_terminals()
                        )
                        # canonical split might be the complement
                        complement = frozenset(all_taxa) - node_tips
                        canon = min(node_tips, complement, key=len) if len(node_tips) != len(complement) else min(node_tips, complement, key=lambda s: sorted(s))
                        if canon == split_set:
                            assert abs(etm_gcf - gcf) < 0.01, (
                                f"gCF mismatch at {split_set}: "
                                f"ETM={etm_gcf:.3f}, concordance_asr={gcf:.3f}"
                            )
                            matched = True
                            break
                if matched:
                    break


class TestStatisticalTests:
    @pytest.fixture
    def svc(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return EvoTempoMap(args)

    def test_mann_whitney_returns_U_and_pvalue(self, svc):
        conc = [10.0, 12.0, 11.0, 13.0, 10.5]
        disc = [5.0, 6.0, 4.5, 5.5]
        result = svc._test_branch(conc, disc)
        assert "mann_whitney_U" in result
        assert "mann_whitney_p" in result
        assert "permutation_p" in result
        assert result["mann_whitney_p"] < 0.05

    def test_permutation_p_between_0_and_1(self, svc):
        conc = [10.0, 12.0, 11.0, 13.0, 10.5]
        disc = [5.0, 6.0, 4.5, 5.5]
        result = svc._test_branch(conc, disc)
        assert 0.0 <= result["permutation_p"] <= 1.0

    def test_no_difference_gives_high_p(self, svc):
        rng = np.random.default_rng(42)
        vals = rng.normal(10, 1, 20).tolist()
        conc = vals[:10]
        disc = vals[10:]
        result = svc._test_branch(conc, disc)
        assert result["mann_whitney_p"] > 0.05

    def test_insufficient_data_returns_none(self, svc):
        conc = [10.0, 12.0, 11.0]
        disc = [5.0]
        result = svc._test_branch(conc, disc)
        assert result["mann_whitney_p"] is None
        assert result["permutation_p"] is None

    def test_empty_discordant_returns_none(self, svc):
        conc = [10.0, 12.0, 11.0]
        disc = []
        result = svc._test_branch(conc, disc)
        assert result["mann_whitney_p"] is None
        assert result["permutation_p"] is None

    def test_insufficient_data_does_not_import_numpy_or_scipy(self):
        code = """
import sys
from phykit.services.tree.evo_tempo_map import EvoTempoMap

svc = EvoTempoMap.__new__(EvoTempoMap)
result = svc._test_branch([10.0], [5.0])
assert result["concordant_mean"] == 10.0
assert result["discordant_mean"] == 5.0
assert result["mann_whitney_p"] is None
assert result["permutation_p"] is None
assert "numpy" not in sys.modules
assert "scipy" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_summary_stats_correct(self, svc):
        conc = [10.0, 20.0, 30.0]
        disc = [5.0, 15.0, 25.0]
        result = svc._test_branch(conc, disc)
        assert abs(result["concordant_mean"] - 20.0) < 1e-10
        assert abs(result["concordant_median"] - 20.0) < 1e-10
        assert abs(result["discordant_mean"] - 15.0) < 1e-10
        assert abs(result["discordant_median"] - 15.0) < 1e-10


class TestFDRCorrection:
    def test_fdr_correction(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        p_values = [0.01, 0.04, 0.03, 0.20]
        corrected = EvoTempoMap._fdr(p_values)
        assert len(corrected) == 4
        for orig, corr in zip(p_values, corrected):
            assert corr >= orig
        for c in corrected:
            assert c <= 1.0

    def test_fdr_empty_list(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        assert EvoTempoMap._fdr([]) == []

    def test_fdr_single_value(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        assert EvoTempoMap._fdr([0.05]) == [0.05]

    def test_fdr_matches_scalar_reference_with_ties(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        p_values = [0.20, 0.01, 0.01, 0.50, 0.03, 0.80, 0.03]
        indexed = sorted(enumerate(p_values), key=lambda x: x[1])
        expected = [0.0] * len(p_values)
        previous = 1.0
        for rank_minus_1 in range(len(p_values) - 1, -1, -1):
            original_idx, p_value = indexed[rank_minus_1]
            rank = rank_minus_1 + 1
            adjusted = min(p_value * len(p_values) / rank, previous)
            adjusted = min(adjusted, 1.0)
            expected[original_idx] = adjusted
            previous = adjusted

        assert EvoTempoMap._fdr(p_values) == pytest.approx(expected)

    def test_medium_fdr_matches_scalar_reference(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        p_values = [((idx * 37) % 101) / 1000 for idx in range(32)]
        indexed = sorted(enumerate(p_values), key=lambda x: x[1])
        expected = [0.0] * len(p_values)
        previous = 1.0
        for rank_minus_1 in range(len(p_values) - 1, -1, -1):
            original_idx, p_value = indexed[rank_minus_1]
            rank = rank_minus_1 + 1
            adjusted = min(p_value * len(p_values) / rank, previous)
            adjusted = min(adjusted, 1.0)
            expected[original_idx] = adjusted
            previous = adjusted

        assert EvoTempoMap._fdr(p_values) == pytest.approx(expected)

    def test_vector_fdr_caps_adjusted_values(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        p_values = [1.0 if idx % 5 else 0.001 * (idx + 1) for idx in range(64)]
        indexed = sorted(enumerate(p_values), key=lambda x: x[1])
        expected = [0.0] * len(p_values)
        previous = 1.0
        for rank_minus_1 in range(len(p_values) - 1, -1, -1):
            original_idx, p_value = indexed[rank_minus_1]
            rank = rank_minus_1 + 1
            adjusted = min(p_value * len(p_values) / rank, previous)
            adjusted = min(adjusted, 1.0)
            expected[original_idx] = adjusted
            previous = adjusted

        assert EvoTempoMap._fdr(p_values) == pytest.approx(expected)
        assert max(expected) == 1.0

    def test_small_fdr_does_not_import_numpy(self):
        code = """
import sys
from phykit.services.tree.evo_tempo_map import EvoTempoMap
corrected = EvoTempoMap._fdr([0.20, 0.01, 0.01, 0.50, 0.03, 0.80, 0.03])
assert [round(value, 10) for value in corrected] == [0.28, 0.035, 0.035, 0.5833333333, 0.0525, 0.8, 0.0525]
assert "numpy" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_fdr_matches_relative_rate_test(self):
        """Cross-check: our FDR should match relative_rate_test._fdr."""
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        from phykit.services.tree.relative_rate_test import RelativeRateTest
        p_values = [0.001, 0.01, 0.05, 0.10, 0.50]
        etm_corrected = EvoTempoMap._fdr(p_values)
        rrt_corrected = RelativeRateTest._fdr(p_values)
        for a, b in zip(etm_corrected, rrt_corrected):
            assert abs(a - b) < 1e-10


class TestGlobalTreeness:
    @pytest.fixture
    def svc(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        return EvoTempoMap(args)

    def test_treeness_matches_base_class(self, svc):
        """Cross-check: treeness computed here should match
        Tree.calculate_treeness() for each gene tree."""
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        for gt in gene_trees:
            expected = svc.calculate_treeness(gt)
            actual = svc._compute_treeness(gt)
            assert abs(actual - expected) < 1e-10, (
                f"Treeness mismatch: expected {expected}, got {actual}"
            )

    def test_treeness_zero_length_returns_zero(self, svc):
        from Bio.Phylo.Newick import Clade, Tree

        tree = Tree(
            root=Clade(
                branch_length=0.0,
                clades=[
                    Clade(branch_length=0.0, name="a"),
                    Clade(branch_length=0.0, name="b"),
                ],
            )
        )

        assert svc._compute_treeness(tree) == 0.0

    def test_global_treeness_returns_structure(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._compute_global_treeness(species_tree, gene_trees)
        assert "treeness_concordant" in result
        assert "treeness_discordant" in result
        assert "treeness_U_p" in result
        assert "n" in result["treeness_concordant"]
        assert "mean" in result["treeness_concordant"]

    def test_all_concordant_treeness(self, svc):
        species_tree = svc.read_tree_file()
        identical_trees = [svc.read_tree_file() for _ in range(3)]
        result = svc._compute_global_treeness(species_tree, identical_trees)
        assert result["treeness_concordant"]["n"] == 3
        assert result["treeness_discordant"]["n"] == 0
        assert result["treeness_U_p"] is None  # Can't test with 0 discordant

    def test_all_concordant_treeness_avoids_numpy_and_mann_whitney(
        self, svc, monkeypatch
    ):
        import phykit.services.tree.evo_tempo_map as etm_module

        class FailingNumpy:
            def __getattr__(self, name):
                raise AssertionError("treeness summaries should not use NumPy")

        def fail_mann_whitney(*args, **kwargs):
            raise AssertionError("single-group treeness should not run Mann-Whitney")

        monkeypatch.setattr(etm_module, "np", FailingNumpy())
        monkeypatch.setattr(etm_module, "_mannwhitneyu", fail_mann_whitney)

        species_tree = svc.read_tree_file()
        identical_trees = [svc.read_tree_file() for _ in range(3)]
        result = svc._compute_global_treeness(species_tree, identical_trees)

        assert result["treeness_concordant"]["n"] == 3
        assert result["treeness_discordant"]["n"] == 0
        assert result["treeness_U_p"] is None

    def test_global_treeness_uses_cached_clade_taxa(self, svc, monkeypatch):
        from io import StringIO

        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin

        species_tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick"
        )
        gene_trees = [
            Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick"),
            Phylo.read(StringIO("((A:1,C:1):1,(B:1,D:1):1);"), "newick"),
        ]

        def fail_get_terminals(*args, **kwargs):
            raise AssertionError("cached clade taxa should be used")

        def fail_find_clades(*args, **kwargs):
            raise AssertionError("standard trees should use direct traversal")

        def fail_get_nonterminals(*args, **kwargs):
            raise AssertionError("cached split sets should use direct traversal")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)
        monkeypatch.setattr(TreeMixin, "get_nonterminals", fail_get_nonterminals)

        result = svc._compute_global_treeness(species_tree, gene_trees)

        assert result["treeness_concordant"]["n"] == 1
        assert result["treeness_discordant"]["n"] == 1

    def test_concordant_plus_discordant_equals_total(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._compute_global_treeness(species_tree, gene_trees)
        total = result["treeness_concordant"]["n"] + result["treeness_discordant"]["n"]
        assert total == 10  # 10 gene trees in the sample file

    def test_treeness_values_between_0_and_1(self, svc):
        species_tree = svc.read_tree_file()
        gene_trees = svc._parse_gene_trees(GENE_TREES)
        result = svc._compute_global_treeness(species_tree, gene_trees)
        for group in ["treeness_concordant", "treeness_discordant"]:
            if result[group]["mean"] is not None:
                assert 0 < result[group]["mean"] < 1
            if result[group]["median"] is not None:
                assert 0 < result[group]["median"] < 1


class TestPlot:
    def test_plot_creates_file(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        with tempfile.NamedTemporaryFile(
            suffix=".png", delete=False
        ) as f:
            tmp_path = f.name
        try:
            args = Namespace(
                tree=TREE_SIMPLE, gene_trees=GENE_TREES,
                verbose=False, json=False, plot_output=tmp_path,
            )
            svc = EvoTempoMap(args)
            svc.run()
            assert os.path.exists(tmp_path)
            assert os.path.getsize(tmp_path) > 0
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)

    def test_plot_batches_strip_points(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes

        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            verbose=False,
            json=False,
            plot_output=str(tmp_path / "tempo.png"),
            no_title=True,
        )
        svc = EvoTempoMap(args)
        branch_results = [
            {
                "split": [f"T{i}", f"U{i}"],
                "_concordant_lengths": [1.0 + i, 1.1 + i],
                "_discordant_lengths": [1.3 + i],
                "fdr_p": 0.01 if i == 0 else 0.5,
            }
            for i in range(4)
        ]

        original_text = matplotlib.axes.Axes.text
        original_scatter = matplotlib.axes.Axes.scatter
        scatter_calls = []

        def fail_star_text(self, *args, **kwargs):
            if len(args) >= 3 and args[2] == "*":
                raise AssertionError("significance stars should be batched")
            return original_text(self, *args, **kwargs)

        def capture_scatter(self, x, y, *args, **kwargs):
            scatter_calls.append((x, y, kwargs))
            return original_scatter(self, x, y, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "text", fail_star_text)
        monkeypatch.setattr(matplotlib.axes.Axes, "scatter", capture_scatter)

        svc._plot(branch_results, str(tmp_path / "tempo.png"))

        strip_calls = [
            call for call in scatter_calls
            if call[2].get("color") in {"#4C72B0", "#DD8452"}
        ]
        star_calls = [
            call for call in scatter_calls
            if call[2].get("marker") == "$*$"
        ]
        assert len(strip_calls) == 2
        conc_call, disc_call = strip_calls
        assert len(conc_call[0]) == len(conc_call[1]) == 8
        assert len(disc_call[0]) == len(disc_call[1]) == 4
        assert conc_call[2]["color"] == "#4C72B0"
        assert disc_call[2]["color"] == "#DD8452"
        assert len(star_calls) == 1
        assert len(star_calls[0][0]) == len(star_calls[0][1]) == 1
        assert (tmp_path / "tempo.png").exists()

    def test_plot_not_created_without_flag(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        svc.run()  # Should not crash even though no plot output


class TestRun:
    def test_run_reuses_unmodified_species_tree(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap

        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        species_tree = object()
        gene_trees = [object()] * 5
        classification = {
            "A,B": {
                "split": ["A", "B"],
                "n_concordant": 3,
                "n_discordant": 2,
                "concordant_lengths": [1.0, 1.1],
                "discordant_lengths": [1.3, 1.4],
            }
        }
        global_stats = {
            "treeness_concordant": {"n": 3, "mean": 0.5, "median": 0.5},
            "treeness_discordant": {"n": 2, "mean": 0.4, "median": 0.4},
            "mann_whitney_p": None,
        }
        branch_stats = {
            "concordant_mean": None,
            "concordant_median": None,
            "concordant_std": None,
            "discordant_mean": None,
            "discordant_median": None,
            "discordant_std": None,
            "mann_whitney_U": None,
            "mann_whitney_p": None,
            "permutation_p": None,
        }

        with patch.object(
            svc,
            "read_tree_file",
            side_effect=AssertionError("run should not copy the cached species tree"),
        ), patch.object(
            svc, "read_tree_file_unmodified", return_value=species_tree
        ) as read_unmodified, patch.object(
            svc, "_parse_and_validate_gene_trees", return_value=gene_trees
        ), patch.object(
            svc, "_classify_gene_trees", return_value=classification
        ) as classify, patch.object(
            svc, "_compute_global_treeness", return_value=global_stats
        ) as global_treeness, patch.object(
            svc, "_test_branch", return_value=branch_stats
        ), patch.object(
            svc, "_output_text"
        ):
            svc.run()

        read_unmodified.assert_called_once_with()
        classify.assert_called_once_with(species_tree, gene_trees)
        global_treeness.assert_called_once_with(species_tree, gene_trees)

    def test_run_text_output(self, capsys):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        svc.run()
        captured = capsys.readouterr()
        assert "n_conc" in captured.out
        assert "n_disc" in captured.out
        assert "Global treeness" in captured.out
        assert "Branches tested" in captured.out

    def test_run_json_output(self, capsys):
        import json
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=False, json=True, plot_output=None,
        )
        svc = EvoTempoMap(args)
        svc.run()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert "branches" in data
        assert "global_" in data
        assert "n_gene_trees" in data["global_"]
        for branch in data["branches"]:
            assert "split" in branch
            assert "n_concordant" in branch
            assert "mann_whitney_p" in branch
            assert "fdr_p" in branch

    def test_run_verbose_output(self, capsys):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=True, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        svc.run()
        captured = capsys.readouterr()
        # Verbose should print per-branch raw lengths
        assert "Concordant lengths" in captured.out or "concordant" in captured.out.lower()

    def test_output_text_batches_branch_rows(self):
        from phykit.services.tree.evo_tempo_map import EvoTempoMap
        args = Namespace(
            tree=TREE_SIMPLE, gene_trees=GENE_TREES,
            verbose=True, json=False, plot_output=None,
        )
        svc = EvoTempoMap(args)
        branch_results = [
            {
                "split": ["bear", "raccoon"],
                "n_concordant": 3,
                "n_discordant": 2,
                "concordant_median": 0.4,
                "discordant_median": 0.7,
                "mann_whitney_p": 0.125,
                "permutation_p": 0.25,
                "fdr_p": 0.5,
                "_concordant_lengths": [1.0, 2.0],
                "_discordant_lengths": [3.0],
            }
        ]
        global_stats = {
            "treeness_concordant": {"median": 0.4, "n": 2},
            "treeness_discordant": {"median": 0.7, "n": 1},
            "n_branches_tested": 1,
            "n_significant_fdr05": 0,
        }

        with patch("builtins.print") as mocked_print:
            svc._output_text(branch_results, global_stats)

        mocked_print.assert_called_once()
        output = mocked_print.call_args.args[0]
        header = (
            f"{'branch':<30}"
            f"{'n_conc':>8}"
            f"{'n_disc':>8}"
            f"{'med_conc':>12}"
            f"{'med_disc':>12}"
            f"{'U_pval':>12}"
            f"{'perm_pval':>12}"
            f"{'fdr_p':>12}"
        )
        expected = "\n".join([
            header,
            "-" * len(header),
            (
                f"{'bear,raccoon':<30}"
                f"{3:>8}"
                f"{2:>8}"
                f"{'0.400000':>12}"
                f"{'0.700000':>12}"
                f"{'0.125000':>12}"
                f"{'0.250000':>12}"
                f"{'0.500000':>12}"
            ),
            "---",
            "Global treeness: concordant=0.400000 (n=2), discordant=0.700000 (n=1)",
            "Branches tested: 1, significant (FDR<0.05): 0",
            "",
            "Branch: bear,raccoon",
            "  Concordant lengths: [1.0, 2.0]",
            "  Discordant lengths: [3.0]",
        ])
        assert output == expected


class TestCLIIntegration:
    def test_cli_runs_without_error(self, capsys):
        from phykit.phykit import Phykit

        testargs = [
            "phykit",
            "evo_tempo_map",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        captured = capsys.readouterr()
        assert "n_conc" in captured.out

    def test_cli_alias_etm(self, capsys):
        from phykit.phykit import Phykit

        testargs = [
            "phykit",
            "etm",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        captured = capsys.readouterr()
        assert "n_conc" in captured.out

    def test_cli_json_flag(self, capsys):
        import json
        from phykit.phykit import Phykit

        testargs = [
            "phykit",
            "evo_tempo_map",
            "-t", TREE_SIMPLE,
            "-g", GENE_TREES,
            "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        captured = capsys.readouterr()
        data = json.loads(captured.out)
        assert "branches" in data
