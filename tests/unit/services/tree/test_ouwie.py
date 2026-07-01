import builtins
import importlib
import json
import subprocess
import sys

import numpy as np
import pytest
from argparse import Namespace
from Bio.Phylo.BaseTree import TreeMixin
from Bio.Phylo.Newick import Clade, Tree as NewickTree
from pathlib import Path
from unittest.mock import patch

import phykit.services.tree.ouwie as ouwie_module
from phykit.services.tree.ouwie import OUwie, ALL_MODELS
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")
REGIMES_FILE = str(SAMPLE_FILES / "tree_simple_regimes.tsv")


def test_module_import_does_not_import_scipy_linalg_or_optimize(monkeypatch):
    module_name = "phykit.services.tree.ouwie"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if (
            name == "scipy.linalg"
            or name.startswith("scipy.linalg.")
            or name == "scipy.optimize"
            or name.startswith("scipy.optimize.")
        ):
            raise AssertionError(
                "ouwie module import should not import SciPy linalg/optimize"
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


def test_module_import_does_not_import_numpy_or_scipy():
    code = """
import sys
import phykit.services.tree.ouwie as module
assert callable(module.print_json)
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "pickle" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy.linalg" not in sys.modules
assert "scipy.optimize" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_merge_nonempty_child_state_sets_does_not_slice_children():
    class NoSliceList(list):
        def __getitem__(self, key):
            if isinstance(key, slice):
                raise AssertionError("child state-set merge should not slice")
            return super().__getitem__(key)

    left = object()
    middle = object()
    right = object()
    empty = object()
    state_sets = {
        id(left): {"A", "B"},
        id(middle): {"A", "C"},
        id(right): {"B", "C"},
        id(empty): set(),
    }

    merged = ouwie_module._merge_nonempty_child_state_sets(
        state_sets,
        NoSliceList([left, middle, right, empty]),
    )

    assert merged == {"A", "B", "C"}
    assert (
        ouwie_module._merge_nonempty_child_state_sets(
            state_sets,
            NoSliceList([left, middle]),
        )
        == {"A"}
    )
    assert (
        ouwie_module._merge_nonempty_child_state_sets(
            state_sets,
            NoSliceList([empty]),
        )
        == set()
    )


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        regime_data=REGIMES_FILE,
        models=None,
        json=False,
    )


@pytest.fixture
def svc(default_args):
    return OUwie(default_args)


@pytest.fixture
def precomputed(svc):
    """Pre-compute all common data structures."""
    import copy
    tree = svc.read_tree_file()
    tree_tips = svc.get_tip_names_from_tree(tree)
    traits = svc._parse_trait_file(TRAITS_FILE, tree_tips)
    regime_assignments = svc._parse_regime_file(REGIMES_FILE, tree_tips)

    shared = set(traits.keys()) & set(regime_assignments.keys())
    traits = {k: traits[k] for k in shared}
    regime_assignments = {k: regime_assignments[k] for k in shared}

    tree_copy = copy.deepcopy(tree)
    tip_names_in_tree = [t.name for t in tree_copy.get_terminals()]
    tips_to_prune = [t for t in tip_names_in_tree if t not in shared]
    if tips_to_prune:
        tree_copy = svc.prune_tree_using_taxa_list(tree_copy, tips_to_prune)

    ordered_names = sorted(traits.keys())
    x = np.array([traits[name] for name in ordered_names])

    regimes = sorted(set(regime_assignments.values()))
    parent_map = svc._build_parent_map(tree_copy)
    branch_regimes = svc._assign_branch_regimes(
        tree_copy, regime_assignments, parent_map
    )
    root_regime = svc._get_root_regime(tree_copy, regime_assignments, parent_map)

    lineage_info = svc._build_lineage_info(
        tree_copy, ordered_names, parent_map, branch_regimes
    )

    per_regime_vcv = svc._build_per_regime_vcv(
        tree_copy, ordered_names, regimes, branch_regimes, parent_map
    )
    vcv_total = sum(per_regime_vcv.values())
    tree_height = float(np.max(np.diag(vcv_total)))

    return dict(
        svc=svc, tree=tree_copy, x=x, ordered_names=ordered_names,
        regimes=regimes, lineage_info=lineage_info,
        per_regime_vcv=per_regime_vcv, vcv_total=vcv_total,
        root_regime=root_regime, tree_height=tree_height,
        n_regimes=len(regimes),
    )


# ── TestProcessArgs ──────────────────────────────────────────────────

class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE, models=None, json=False,
        )
        svc = OUwie(args)
        assert svc.selected_models == ALL_MODELS

    def test_subset_models(self):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE, models="BM1,OUM,OUMVA", json=False,
        )
        svc = OUwie(args)
        assert svc.selected_models == ["BM1", "OUM", "OUMVA"]

    def test_invalid_model(self):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE, models="BM1,BOGUS", json=False,
        )
        with pytest.raises(PhykitUserError):
            OUwie(args)

    def test_shared_trait_regime_taxa_returns_overlap(self):
        traits = {"C": 3.0, "A": 1.0, "B": 2.0}
        regimes = {"B": "r1", "D": "r2", "A": "r1"}

        assert OUwie._shared_trait_regime_taxa(traits, regimes) == {"A", "B"}


class TestTraitAndRegimeParsing:
    def test_trait_file_skips_comments_and_blanks(self, tmp_path, default_args):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "# comment\n\nraccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\n"
            "seal\t4.0\nmonkey\t5.0\ncat\t6.0\nweasel\t7.0\ndog\t8.0\n"
        )
        svc = OUwie(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        traits = svc._parse_trait_file(str(trait_file), tree_tips)

        assert len(traits) == 8
        assert traits["raccoon"] == pytest.approx(1.0)

    def test_all_shared_trait_file_emits_no_warnings(
        self, tmp_path, default_args, capsys
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "raccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\nseal\t4.0\n"
            "monkey\t5.0\ncat\t6.0\nweasel\t7.0\ndog\t8.0\n"
        )
        svc = OUwie(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        traits = svc._parse_trait_file(str(trait_file), tree_tips)

        stderr = capsys.readouterr().err
        assert set(traits) == set(tree_tips)
        assert stderr == ""

    def test_trait_file_extra_columns_error(self, tmp_path, default_args):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "raccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\ndog\t4.0\nextra\t1.0\t2.0\n"
        )
        svc = OUwie(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_trait_file(str(trait_file), tree_tips)

        assert (
            "Line 5 in trait file has 3 columns; expected 2."
            in exc_info.value.messages
        )

    def test_trait_file_non_numeric_error(self, tmp_path, default_args):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("raccoon\t1.0\nbear\tbad\nsea_lion\t3.0\ndog\t4.0\n")
        svc = OUwie(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_trait_file(str(trait_file), tree_tips)

        assert (
            "Non-numeric trait value 'bad' for taxon 'bear' on line 2."
            in exc_info.value.messages
        )

    def test_all_shared_regime_file_emits_no_warnings(
        self, tmp_path, default_args, capsys
    ):
        regime_file = tmp_path / "regimes.tsv"
        regime_file.write_text(
            "raccoon\tA\nbear\tA\nsea_lion\tB\nseal\tB\n"
            "monkey\tA\ncat\tA\nweasel\tB\ndog\tB\n"
        )
        svc = OUwie(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        regimes = svc._parse_regime_file(str(regime_file), tree_tips)

        stderr = capsys.readouterr().err
        assert set(regimes) == set(tree_tips)
        assert stderr == ""

    def test_regime_file_extra_columns_error(self, tmp_path, default_args):
        regime_file = tmp_path / "regimes.tsv"
        regime_file.write_text(
            "raccoon\tA\nbear\tA\nsea_lion\tB\ndog\tA\nextra\tA\tB\n"
        )
        svc = OUwie(default_args)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_regime_file(str(regime_file), tree_tips)

        assert (
            "Line 5 in regime file has 3 columns; expected 2."
            in exc_info.value.messages
        )


# ── TestWeightMatrix ─────────────────────────────────────────────────

class TestWeightMatrix:
    def test_rows_sum_to_one(self, precomputed):
        d = precomputed
        alpha = 0.5
        W = d["svc"]._build_weight_matrix_single_alpha(
            d["ordered_names"], d["lineage_info"], d["regimes"],
            alpha, d["root_regime"],
        )
        row_sums = W.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)

    def test_single_regime_reduces_to_ones(self, precomputed):
        d = precomputed
        alpha = 1.0
        # If all taxa were same regime, W would be n x 1 of ones
        # With 2 regimes, rows should still sum to 1
        W = d["svc"]._build_weight_matrix_single_alpha(
            d["ordered_names"], d["lineage_info"], d["regimes"],
            alpha, d["root_regime"],
        )
        assert W.shape == (len(d["ordered_names"]), len(d["regimes"]))
        assert np.all(W >= -1e-10)  # Non-negative weights

    def test_multi_alpha_rows_sum_to_one(self, precomputed):
        d = precomputed
        alphas_dict = {r: 0.5 for r in d["regimes"]}
        W = d["svc"]._build_weight_matrix_multi_alpha(
            d["ordered_names"], d["lineage_info"], d["regimes"],
            alphas_dict, d["root_regime"],
        )
        row_sums = W.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)

    @pytest.mark.parametrize("alpha", [1e-12, 0.5])
    def test_single_alpha_weight_cache_matches_direct_builder(
        self, precomputed, alpha
    ):
        d = precomputed
        cache = d["svc"]._prepare_single_alpha_weight_cache(
            d["ordered_names"],
            d["lineage_info"],
            d["regimes"],
            d["root_regime"],
        )

        expected = d["svc"]._build_weight_matrix_single_alpha(
            d["ordered_names"],
            d["lineage_info"],
            d["regimes"],
            alpha,
            d["root_regime"],
        )
        observed = d["svc"]._build_weight_matrix_single_alpha_from_cache(
            cache,
            alpha,
        )

        np.testing.assert_allclose(observed, expected)


# ── TestVCV ──────────────────────────────────────────────────────────

class TestRegimeAssignment:
    def test_parent_map_uses_direct_traversal(self, svc, monkeypatch):
        tree = svc.read_tree_file()

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("parent-map helper should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        parent_map = svc._build_parent_map(tree)

        seen_children = []
        stack = [tree.root]
        while stack:
            clade = stack.pop()
            for child in clade.clades:
                assert parent_map[id(child)] is clade
                seen_children.append(child)
                stack.append(child)

        assert seen_children
        assert id(tree.root) not in parent_map

    def test_parent_map_handles_mixed_child_counts(self, svc, monkeypatch):
        tree = NewickTree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="A"),
                    Clade(
                        branch_length=1.0,
                        clades=[
                            Clade(branch_length=1.0, name="B"),
                            Clade(branch_length=1.0, name="C"),
                        ],
                    ),
                    Clade(
                        branch_length=1.0,
                        clades=[
                            Clade(branch_length=1.0, name="D"),
                            Clade(branch_length=1.0, name="E"),
                            Clade(branch_length=1.0, name="F"),
                        ],
                    ),
                ],
            )
        )

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("parent-map helper should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        parent_map = svc._build_parent_map(tree)
        binary = tree.root.clades[1]
        trifurcation = tree.root.clades[2]

        assert parent_map[id(tree.root.clades[0])] is tree.root
        assert parent_map[id(binary)] is tree.root
        assert parent_map[id(trifurcation)] is tree.root
        assert parent_map[id(binary.clades[0])] is binary
        assert parent_map[id(binary.clades[1])] is binary
        assert parent_map[id(trifurcation.clades[0])] is trifurcation
        assert parent_map[id(trifurcation.clades[1])] is trifurcation
        assert parent_map[id(trifurcation.clades[2])] is trifurcation
        assert id(tree.root) not in parent_map

    def test_combined_branch_and_root_regime_uses_direct_traversal(
        self, svc, monkeypatch
    ):
        tree = svc.read_tree_file()
        tip_names = svc.get_tip_names_from_tree(tree)
        tip_regimes = {
            name: ("r1" if idx % 2 == 0 else "r2")
            for idx, name in enumerate(tip_names)
        }
        parent_map = svc._build_parent_map(tree)

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("combined helper should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        branch_regimes, root_regime = svc._assign_branch_regimes_and_root(
            tree,
            tip_regimes,
            parent_map,
        )

        assert root_regime in {"r1", "r2"}
        assert branch_regimes
        assert all(regime in {"r1", "r2"} for regime in branch_regimes.values())

    def test_combined_branch_and_root_regime_handles_mixed_child_counts(self, svc):
        binary = Clade(
            clades=[
                Clade(name="B", branch_length=1.0),
                Clade(name="C", branch_length=1.0),
            ],
            branch_length=1.0,
        )
        trifurcation = Clade(
            clades=[
                Clade(name="D", branch_length=1.0),
                Clade(name="E", branch_length=1.0),
                Clade(name="F", branch_length=1.0),
            ],
            branch_length=1.0,
        )
        tree = NewickTree(
            root=Clade(
                clades=[
                    Clade(name="A", branch_length=1.0),
                    binary,
                    trifurcation,
                ]
            )
        )
        tip_regimes = {
            "A": "r1",
            "B": "r2",
            "C": "r2",
            "D": "r3",
            "E": "r3",
            "F": "r1",
        }
        parent_map = svc._build_parent_map(tree)

        branch_regimes, root_regime = svc._assign_branch_regimes_and_root(
            tree,
            tip_regimes,
            parent_map,
        )

        assert root_regime == "r1"
        assert branch_regimes[id(tree.root.clades[0])] == "r1"
        assert branch_regimes[id(binary)] == "r2"
        assert branch_regimes[id(binary.clades[0])] == "r2"
        assert branch_regimes[id(binary.clades[1])] == "r2"
        assert branch_regimes[id(trifurcation)] == "r1"
        assert branch_regimes[id(trifurcation.clades[0])] == "r3"
        assert branch_regimes[id(trifurcation.clades[1])] == "r3"
        assert branch_regimes[id(trifurcation.clades[2])] == "r1"

    def test_combined_branch_and_root_regime_uses_lexicographic_tiebreak(self, svc):
        left = Clade(
            clades=[
                Clade(name="A", branch_length=1.0),
                Clade(name="B", branch_length=1.0),
            ],
            branch_length=1.0,
        )
        right = Clade(
            clades=[
                Clade(name="C", branch_length=1.0),
                Clade(name="D", branch_length=1.0),
            ],
            branch_length=1.0,
        )
        tree = NewickTree(root=Clade(clades=[left, right]))
        tip_regimes = {
            "A": "zeta",
            "B": "beta",
            "C": "alpha",
            "D": "alpha",
        }
        parent_map = svc._build_parent_map(tree)

        branch_regimes, root_regime = svc._assign_branch_regimes_and_root(
            tree,
            tip_regimes,
            parent_map,
        )

        assert root_regime == "alpha"
        assert branch_regimes[id(left)] == "beta"
        assert branch_regimes[id(right)] == "alpha"


class TestVCV:
    def test_build_root_to_tip_paths_uses_set_for_tip_filtering(self, precomputed):
        class NoContainsList(list):
            def __contains__(self, _item):
                raise AssertionError("tip filtering should not scan ordered_names")

        d = precomputed
        ordered_names = NoContainsList(d["ordered_names"])
        paths = d["svc"]._build_root_to_tip_paths(
            d["tree"],
            ordered_names,
            d["svc"]._build_parent_map(d["tree"]),
        )

        assert set(paths) == set(d["ordered_names"])
        assert all(paths[name] for name in d["ordered_names"])

    def test_build_lineage_info_uses_set_for_tip_filtering(self, precomputed):
        class NoContainsList(list):
            def __contains__(self, _item):
                raise AssertionError("tip filtering should not scan ordered_names")

        d = precomputed
        ordered_names = NoContainsList(d["ordered_names"])
        lineage_info = d["svc"]._build_lineage_info(
            d["tree"],
            ordered_names,
            d["svc"]._build_parent_map(d["tree"]),
            {},
        )

        assert set(lineage_info) == set(d["ordered_names"])
        assert all(lineage_info[name] for name in d["ordered_names"])

    def test_per_regime_vcv_matches_explicit_pairwise_prefixes(self, precomputed):
        d = precomputed
        ordered_names = d["ordered_names"]
        regimes = d["regimes"]
        paths = d["svc"]._build_root_to_tip_paths(
            d["tree"], ordered_names, d["svc"]._build_parent_map(d["tree"])
        )
        branch_ids = []
        for name in ordered_names:
            for clade_id, _ in paths[name]:
                if clade_id not in branch_ids:
                    branch_ids.append(clade_id)
        branch_regimes = {
            clade_id: regimes[idx % len(regimes)]
            for idx, clade_id in enumerate(branch_ids)
        }

        expected = {r: np.zeros((len(ordered_names), len(ordered_names))) for r in regimes}
        for i, name_i in enumerate(ordered_names):
            path_i = paths[name_i]
            for clade_id, branch_length in path_i:
                regime = branch_regimes.get(clade_id, regimes[0])
                expected[regime][i, i] += branch_length

            for j in range(i + 1, len(ordered_names)):
                path_j = paths[ordered_names[j]]
                for branch_i, branch_j in zip(path_i, path_j):
                    if branch_i[0] != branch_j[0]:
                        break
                    clade_id, branch_length = branch_i
                    regime = branch_regimes.get(clade_id, regimes[0])
                    expected[regime][i, j] += branch_length
                    expected[regime][j, i] += branch_length

        observed = d["svc"]._build_per_regime_vcv(
            d["tree"], ordered_names, regimes, branch_regimes,
            d["svc"]._build_parent_map(d["tree"]),
        )

        for regime in regimes:
            np.testing.assert_allclose(observed[regime], expected[regime])

    def test_per_regime_vcv_single_tip_branches_skip_block_indexing(
        self, svc, monkeypatch
    ):
        ordered_names = ["A", "B", "C"]
        regimes = ["r1", "r2"]
        paths = {
            "A": [("A_tip", 1.0)],
            "B": [("B_tip", 2.0)],
            "C": [("C_tip", 3.0)],
        }
        branch_regimes = {
            "A_tip": "r1",
            "B_tip": "r2",
            "C_tip": "r1",
        }

        def fail_ix(*_args, **_kwargs):
            raise AssertionError("single-tip branches should update diagonals directly")

        monkeypatch.setattr(
            svc, "_build_root_to_tip_paths",
            lambda _tree, _ordered_names, _parent_map: paths,
        )
        monkeypatch.setattr(ouwie_module.np, "ix_", fail_ix)

        observed = svc._build_per_regime_vcv(
            None, ordered_names, regimes, branch_regimes, {}
        )

        np.testing.assert_allclose(observed["r1"], np.diag([1.0, 0.0, 3.0]))
        np.testing.assert_allclose(observed["r2"], np.diag([0.0, 2.0, 0.0]))

    def test_per_regime_vcv_reuses_lineage_info(self, svc, monkeypatch):
        ordered_names = ["A", "B"]
        regimes = ["r1", "r2"]
        lineage_info = {
            "A": [
                ("shared", 2.0, "r1", 0.0, 2.0),
                ("A_tip", 1.0, "r2", 2.0, 3.0),
            ],
            "B": [
                ("shared", 2.0, "r1", 0.0, 2.0),
                ("B_tip", 3.0, "r2", 2.0, 5.0),
            ],
        }

        monkeypatch.setattr(
            svc,
            "_build_root_to_tip_paths",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("lineage-info cache should avoid rebuilding paths")
            ),
        )

        observed = svc._build_per_regime_vcv(
            None,
            ordered_names,
            regimes,
            {},
            {},
            lineage_info=lineage_info,
        )

        np.testing.assert_allclose(observed["r1"], [[2.0, 2.0], [2.0, 2.0]])
        np.testing.assert_allclose(observed["r2"], [[1.0, 0.0], [0.0, 3.0]])

    def test_single_alpha_cache_uses_shared_branch_lengths(self, svc):
        ordered_names = ["A", "B", "C"]
        lineage_info = {
            "A": [
                ("root_left", 2.0, "r1", 0.0, 2.0),
                ("A_tip", 1.0, "r1", 2.0, 3.0),
            ],
            "B": [
                ("root_left", 2.0, "r1", 0.0, 2.0),
                ("B_tip", 3.0, "r2", 2.0, 5.0),
            ],
            "C": [
                ("C_tip", 4.0, "r2", 0.0, 4.0),
            ],
        }

        cache = svc._prepare_single_alpha_ou_cache(ordered_names, lineage_info)

        expected_shared = np.array(
            [
                [3.0, 2.0, 0.0],
                [2.0, 5.0, 0.0],
                [0.0, 0.0, 4.0],
            ]
        )
        expected_unique = np.array(
            [
                [0.0, 4.0, 7.0],
                [4.0, 0.0, 9.0],
                [7.0, 9.0, 0.0],
            ]
        )

        np.testing.assert_allclose(cache["tip_heights"], [3.0, 5.0, 4.0])
        np.testing.assert_allclose(cache["shared_paths"], expected_shared)
        np.testing.assert_allclose(cache["unique_paths"], expected_unique)

    def test_single_alpha_cache_single_tip_branches_skip_block_indexing(
        self, svc, monkeypatch
    ):
        ordered_names = ["A", "B", "C"]
        lineage_info = {
            "A": [("A_tip", 1.0, "r1", 0.0, 1.0)],
            "B": [("B_tip", 2.0, "r2", 0.0, 2.0)],
            "C": [("C_tip", 3.0, "r1", 0.0, 3.0)],
        }

        def fail_ix(*_args, **_kwargs):
            raise AssertionError("single-tip branches should update diagonals directly")

        monkeypatch.setattr(ouwie_module.np, "ix_", fail_ix)

        cache = svc._prepare_single_alpha_ou_cache(ordered_names, lineage_info)

        np.testing.assert_allclose(cache["tip_heights"], [1.0, 2.0, 3.0])
        np.testing.assert_allclose(cache["shared_paths"], np.diag([1.0, 2.0, 3.0]))
        expected_unique = np.array([
            [0.0, 3.0, 4.0],
            [3.0, 0.0, 5.0],
            [4.0, 5.0, 0.0],
        ])
        np.testing.assert_allclose(cache["unique_paths"], expected_unique)

    def test_single_alpha_cache_matches_explicit_path_builder(self, precomputed):
        d = precomputed
        alpha = 0.7
        sigma2 = 1.3

        def reference():
            ordered_names = d["ordered_names"]
            lineage_info = d["lineage_info"]
            n = len(ordered_names)
            vcv = np.zeros((n, n))
            tip_heights = {
                name: sum(bl for _, bl, _, _, _ in lineage_info[name])
                for name in ordered_names
            }
            for i in range(n):
                for j in range(i, n):
                    t_i = tip_heights[ordered_names[i]]
                    t_j = tip_heights[ordered_names[j]]
                    if i == j:
                        shared = t_i
                    else:
                        shared = 0.0
                        path_i = lineage_info[ordered_names[i]]
                        path_j = lineage_info[ordered_names[j]]
                        for branch_i, branch_j in zip(path_i, path_j):
                            if branch_i[0] != branch_j[0]:
                                break
                            shared += branch_i[1]
                    if alpha < 1e-10:
                        value = sigma2 * shared
                    else:
                        value = (
                            sigma2 / (2.0 * alpha)
                        ) * np.exp(
                            -alpha * (t_i + t_j - 2.0 * shared)
                        ) * (
                            1.0 - np.exp(-2.0 * alpha * shared)
                        )
                    vcv[i, j] = value
                    vcv[j, i] = value
            return vcv

        cache = d["svc"]._prepare_single_alpha_ou_cache(
            d["ordered_names"], d["lineage_info"]
        )
        observed = d["svc"]._build_ou_vcv_single_alpha_from_cache(
            cache, alpha, sigma2
        )

        np.testing.assert_allclose(observed, reference())

    def test_positive_definite(self, precomputed):
        d = precomputed
        alpha = 0.5
        V = d["svc"]._build_ou_vcv_single_alpha(
            d["ordered_names"], d["lineage_info"], alpha, 1.0,
        )
        eigenvalues = np.linalg.eigvalsh(V)
        assert np.all(eigenvalues > 0)

    def test_H_matrices_sum_correctly(self, precomputed):
        d = precomputed
        alpha = 0.5
        H = d["svc"]._build_ou_H_matrices(
            d["ordered_names"], d["lineage_info"], d["regimes"], alpha,
        )
        # Sum of H_r * sigma2=1 should equal single-alpha VCV with sigma2=1
        V_from_H = sum(H[r] for r in d["regimes"])
        V_direct = d["svc"]._build_ou_vcv_single_alpha(
            d["ordered_names"], d["lineage_info"], alpha, 1.0,
        )
        np.testing.assert_allclose(V_from_H, V_direct, atol=1e-10)

    def test_H_matrices_match_explicit_pairwise_prefixes(self, precomputed):
        d = precomputed
        alpha = 0.7
        ordered_names = d["ordered_names"]
        lineage_info = d["lineage_info"]
        regimes = d["regimes"]
        tip_heights = {
            name: sum(bl for _, bl, _, _, _ in lineage_info[name])
            for name in ordered_names
        }
        expected = {
            regime: np.zeros((len(ordered_names), len(ordered_names)))
            for regime in regimes
        }

        for i, name_i in enumerate(ordered_names):
            path_i = lineage_info[name_i]
            t_i = tip_heights[name_i]
            for _, branch_length, regime, _d_start, d_end in path_i:
                if regime not in expected:
                    continue
                contribution = np.exp(-2.0 * alpha * (t_i - d_end)) * (
                    1.0 - np.exp(-2.0 * alpha * branch_length)
                ) / (2.0 * alpha)
                expected[regime][i, i] += contribution

            for j in range(i + 1, len(ordered_names)):
                path_j = lineage_info[ordered_names[j]]
                t_j = tip_heights[ordered_names[j]]
                for branch_i, branch_j in zip(path_i, path_j):
                    if branch_i[0] != branch_j[0]:
                        break
                    _, branch_length, regime, _d_start, d_end = branch_i
                    if regime not in expected:
                        continue
                    contribution = np.exp(-alpha * (t_i + t_j - 2.0 * d_end)) * (
                        1.0 - np.exp(-2.0 * alpha * branch_length)
                    ) / (2.0 * alpha)
                    expected[regime][i, j] += contribution
                    expected[regime][j, i] += contribution

        observed = d["svc"]._build_ou_H_matrices(
            ordered_names, lineage_info, regimes, alpha,
        )

        for regime in regimes:
            np.testing.assert_allclose(observed[regime], expected[regime])

    def test_H_matrices_single_tip_branches_skip_block_indexing(
        self, svc, monkeypatch
    ):
        ordered_names = ["A", "B", "C"]
        lineage_info = {
            "A": [("A_tip", 1.0, "r1", 0.0, 1.0)],
            "B": [("B_tip", 2.0, "r2", 0.0, 2.0)],
            "C": [("C_tip", 3.0, "r1", 0.0, 3.0)],
        }

        def fail_ix(*_args, **_kwargs):
            raise AssertionError("single-tip branches should update diagonals directly")

        monkeypatch.setattr(ouwie_module.np, "ix_", fail_ix)

        observed = svc._build_ou_H_matrices(
            ordered_names, lineage_info, ["r1", "r2"], alpha=0.7,
        )

        def coeff(bl):
            return (1.0 - np.exp(-2.0 * 0.7 * bl)) / (2.0 * 0.7)

        np.testing.assert_allclose(
            observed["r1"], np.diag([coeff(1.0), 0.0, coeff(3.0)])
        )
        np.testing.assert_allclose(observed["r2"], np.diag([0.0, coeff(2.0), 0.0]))

    def test_multi_alpha_positive_definite(self, precomputed):
        d = precomputed
        alphas_dict = {r: 0.5 for r in d["regimes"]}
        sigmas_dict = {r: 1.0 for r in d["regimes"]}
        V = d["svc"]._build_ou_vcv_multi_alpha(
            d["ordered_names"], d["lineage_info"], d["regimes"],
            alphas_dict, sigmas_dict,
        )
        eigenvalues = np.linalg.eigvalsh(V)
        assert np.all(eigenvalues > 0)

    def test_multi_alpha_vcv_matches_explicit_pairwise_prefixes(self, precomputed):
        d = precomputed
        ordered_names = d["ordered_names"]
        lineage_info = d["lineage_info"]
        regimes = d["regimes"]
        alphas_dict = {
            regime: 0.2 + idx * 0.4
            for idx, regime in enumerate(regimes)
        }
        sigmas_dict = {
            regime: 0.5 + idx * 0.7
            for idx, regime in enumerate(regimes)
        }

        def decay_to_tip(path):
            decays = np.ones(len(path))
            for b in range(len(path) - 2, -1, -1):
                _, next_branch_length, next_regime, _, _ = path[b + 1]
                decays[b] = decays[b + 1] * np.exp(
                    -alphas_dict.get(next_regime, 0.01) * next_branch_length
                )
            return decays

        path_decays = {
            name: decay_to_tip(lineage_info[name])
            for name in ordered_names
        }
        expected = np.zeros((len(ordered_names), len(ordered_names)))
        for i, name_i in enumerate(ordered_names):
            path_i = lineage_info[name_i]
            decay_i = path_decays[name_i]
            for branch_idx, (
                _clade_id, branch_length, regime, _d_start, _d_end,
            ) in enumerate(path_i):
                alpha = alphas_dict.get(regime, 0.01)
                sigma2 = sigmas_dict.get(regime, 0.01)
                if alpha < 1e-10:
                    contribution = sigma2 * branch_length * decay_i[branch_idx] ** 2
                else:
                    contribution = sigma2 / (2.0 * alpha) * (
                        1.0 - np.exp(-2.0 * alpha * branch_length)
                    ) * decay_i[branch_idx] ** 2
                expected[i, i] += contribution

            for j in range(i + 1, len(ordered_names)):
                path_j = lineage_info[ordered_names[j]]
                decay_j = path_decays[ordered_names[j]]
                covariance = 0.0
                for branch_idx, (branch_i, branch_j) in enumerate(zip(path_i, path_j)):
                    if branch_i[0] != branch_j[0]:
                        break
                    _clade_id, branch_length, regime, _d_start, _d_end = branch_i
                    alpha = alphas_dict.get(regime, 0.01)
                    sigma2 = sigmas_dict.get(regime, 0.01)
                    if alpha < 1e-10:
                        contribution = (
                            sigma2 * branch_length
                            * decay_i[branch_idx] * decay_j[branch_idx]
                        )
                    else:
                        contribution = sigma2 / (2.0 * alpha) * (
                            1.0 - np.exp(-2.0 * alpha * branch_length)
                        ) * decay_i[branch_idx] * decay_j[branch_idx]
                    covariance += contribution
                expected[i, j] = covariance
                expected[j, i] = covariance

        observed = d["svc"]._build_ou_vcv_multi_alpha(
            ordered_names, lineage_info, regimes, alphas_dict, sigmas_dict,
        )

        np.testing.assert_allclose(observed, expected)

    def test_multi_alpha_single_tip_branches_skip_block_indexing(
        self, svc, monkeypatch
    ):
        ordered_names = ["A", "B", "C"]
        regimes = ["r1", "r2"]
        lineage_info = {
            "A": [("A_tip", 1.0, "r1", 0.0, 1.0)],
            "B": [("B_tip", 2.0, "r2", 0.0, 2.0)],
            "C": [("C_tip", 3.0, "r1", 0.0, 3.0)],
        }
        alphas = {"r1": 0.7, "r2": 0.4}
        sigmas = {"r1": 1.5, "r2": 0.8}

        def fail_ix(*_args, **_kwargs):
            raise AssertionError("single-tip branches should update diagonals directly")

        monkeypatch.setattr(ouwie_module.np, "ix_", fail_ix)

        observed = svc._build_ou_vcv_multi_alpha(
            ordered_names, lineage_info, regimes, alphas, sigmas,
        )

        def coeff(regime, bl):
            alpha = alphas[regime]
            return sigmas[regime] / (2.0 * alpha) * (
                1.0 - np.exp(-2.0 * alpha * bl)
            )

        expected = np.diag([
            coeff("r1", 1.0),
            coeff("r2", 2.0),
            coeff("r1", 3.0),
        ])
        np.testing.assert_allclose(observed, expected)


# ── TestBM1 ──────────────────────────────────────────────────────────

class TestBM1:
    def test_k_params_equals_2(self, precomputed):
        d = precomputed
        res = d["svc"]._fit_bm1(d["x"], d["vcv_total"])
        assert res["k_params"] == 2

    def test_ll_is_finite(self, precomputed):
        d = precomputed
        res = d["svc"]._fit_bm1(d["x"], d["vcv_total"])
        assert np.isfinite(res["log_likelihood"])

    def test_cholesky_bm_likelihood_matches_inverse(self, precomputed):
        d = precomputed
        fast = d["svc"]._concentrated_ll_bm_cholesky(d["x"], d["vcv_total"])
        inverse = d["svc"]._concentrated_ll_bm_inverse(d["x"], d["vcv_total"])
        np.testing.assert_allclose(fast, inverse)

    def test_cholesky_bm_likelihood_uses_single_solve(
        self, precomputed, monkeypatch
    ):
        d = precomputed
        original = ouwie_module.cho_solve
        calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original(*args, **kwargs)

        monkeypatch.setattr(ouwie_module, "cho_solve", counting_cho_solve)

        fast = d["svc"]._concentrated_ll_bm_cholesky(d["x"], d["vcv_total"])
        inverse = d["svc"]._concentrated_ll_bm_inverse(d["x"], d["vcv_total"])

        assert calls == 1
        np.testing.assert_allclose(fast, inverse)


# ── TestBMS ──────────────────────────────────────────────────────────

class TestBMS:
    def test_ll_ge_bm1(self, precomputed):
        d = precomputed
        bm1 = d["svc"]._fit_bm1(d["x"], d["vcv_total"])
        bms = d["svc"]._fit_bms(d["x"], d["per_regime_vcv"], d["regimes"])
        # BMS has more parameters, should fit at least as well
        assert bms["log_likelihood"] >= bm1["log_likelihood"] - 1e-6

    def test_per_regime_sigma2_positive(self, precomputed):
        d = precomputed
        bms = d["svc"]._fit_bms(d["x"], d["per_regime_vcv"], d["regimes"])
        for r in d["regimes"]:
            assert bms["params"]["sigma2"][r] > 0

    def test_cholesky_vcv_likelihood_matches_inverse(self, precomputed):
        d = precomputed
        ones = np.ones(len(d["x"]))
        V = d["per_regime_vcv"][d["regimes"][0]] * 0.8
        for regime in d["regimes"][1:]:
            V += d["per_regime_vcv"][regime] * 1.2

        fast = d["svc"]._bm_log_likelihood_from_vcv_cholesky(d["x"], V, ones)
        inverse = d["svc"]._bm_log_likelihood_from_vcv_inverse(d["x"], V, ones)
        np.testing.assert_allclose(fast, inverse)

    def test_cholesky_vcv_likelihood_uses_single_solve(
        self, precomputed, monkeypatch
    ):
        d = precomputed
        ones = np.ones(len(d["x"]))
        V = d["per_regime_vcv"][d["regimes"][0]] * 0.8
        for regime in d["regimes"][1:]:
            V += d["per_regime_vcv"][regime] * 1.2
        original = ouwie_module.cho_solve
        calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original(*args, **kwargs)

        monkeypatch.setattr(ouwie_module, "cho_solve", counting_cho_solve)

        fast = d["svc"]._bm_log_likelihood_from_vcv_cholesky(d["x"], V, ones)
        inverse = d["svc"]._bm_log_likelihood_from_vcv_inverse(d["x"], V, ones)

        assert calls == 1
        np.testing.assert_allclose(fast, inverse)


# ── TestOU1 ──────────────────────────────────────────────────────────

class TestOU1:
    def test_fit_ou1_prepares_single_alpha_cache_once(
        self, precomputed, monkeypatch
    ):
        d = precomputed
        original = d["svc"]._prepare_single_alpha_ou_cache
        calls = 0

        def counting_prepare(*args, **kwargs):
            nonlocal calls
            calls += 1
            if calls > 1:
                raise AssertionError("single-alpha OU cache should be reused")
            return original(*args, **kwargs)

        monkeypatch.setattr(
            d["svc"], "_prepare_single_alpha_ou_cache", counting_prepare
        )

        result = d["svc"]._fit_ou1(
            d["x"], d["ordered_names"], d["lineage_info"], d["tree_height"],
        )

        assert calls == 1
        assert np.isfinite(result["log_likelihood"])

    def test_alpha_positive(self, precomputed):
        d = precomputed
        res = d["svc"]._fit_ou1(
            d["x"], d["ordered_names"], d["lineage_info"], d["tree_height"],
        )
        assert res["params"]["alpha"] > 0

    def test_ll_is_finite(self, precomputed):
        d = precomputed
        res = d["svc"]._fit_ou1(
            d["x"], d["ordered_names"], d["lineage_info"], d["tree_height"],
        )
        assert np.isfinite(res["log_likelihood"])

    def test_cholesky_profile_likelihood_matches_inverse(self, precomputed):
        d = precomputed
        V = d["svc"]._build_ou_vcv_single_alpha(
            d["ordered_names"], d["lineage_info"], 0.7, 1.0,
        )
        W = d["svc"]._build_weight_matrix_single_alpha(
            d["ordered_names"], d["lineage_info"], d["regimes"],
            0.7, d["root_regime"],
        )

        theta_fast, sig2_fast, ll_fast = d["svc"]._ou_profile_likelihood_cholesky(
            d["x"], W, V
        )
        theta_inv, sig2_inv, ll_inv = d["svc"]._ou_profile_likelihood_inverse(
            d["x"], W, V
        )

        np.testing.assert_allclose(theta_fast, theta_inv)
        np.testing.assert_allclose(sig2_fast, sig2_inv)
        np.testing.assert_allclose(ll_fast, ll_inv)

    def test_cholesky_profile_likelihood_reuses_solved_columns(
        self, precomputed, monkeypatch
    ):
        d = precomputed
        V = d["svc"]._build_ou_vcv_single_alpha(
            d["ordered_names"], d["lineage_info"], 0.7, 1.0,
        )
        W = d["svc"]._build_weight_matrix_single_alpha(
            d["ordered_names"], d["lineage_info"], d["regimes"],
            0.7, d["root_regime"],
        )
        original = ouwie_module.cho_solve
        calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original(*args, **kwargs)

        monkeypatch.setattr(ouwie_module, "cho_solve", counting_cho_solve)

        theta, sig2, ll = d["svc"]._ou_profile_likelihood_cholesky(
            d["x"], W, V
        )

        assert calls == 1
        assert np.all(np.isfinite(theta))
        assert np.isfinite(sig2)
        assert np.isfinite(ll)

    def test_cholesky_gls_log_likelihood_matches_inverse(self, precomputed):
        d = precomputed
        alpha = 0.7
        sigma2s = [0.5 + i for i in range(len(d["regimes"]))]
        H = d["svc"]._build_ou_H_matrices(
            d["ordered_names"], d["lineage_info"], d["regimes"], alpha,
        )
        V = np.zeros((len(d["x"]), len(d["x"])))
        for r_idx, regime in enumerate(d["regimes"]):
            V += sigma2s[r_idx] * H[regime]
        W = d["svc"]._build_weight_matrix_single_alpha(
            d["ordered_names"], d["lineage_info"], d["regimes"],
            alpha, d["root_regime"],
        )

        theta_fast, ll_fast = d["svc"]._ou_gls_log_likelihood_cholesky(
            d["x"], W, V
        )
        theta_inv, ll_inv = d["svc"]._ou_gls_log_likelihood_inverse(
            d["x"], W, V
        )

        np.testing.assert_allclose(theta_fast, theta_inv)
        np.testing.assert_allclose(ll_fast, ll_inv)

    def test_cholesky_gls_log_likelihood_reuses_solved_columns(
        self, precomputed, monkeypatch
    ):
        d = precomputed
        alpha = 0.7
        sigma2s = [0.5 + i for i in range(len(d["regimes"]))]
        H = d["svc"]._build_ou_H_matrices(
            d["ordered_names"], d["lineage_info"], d["regimes"], alpha,
        )
        V = np.zeros((len(d["x"]), len(d["x"])))
        for r_idx, regime in enumerate(d["regimes"]):
            V += sigma2s[r_idx] * H[regime]
        W = d["svc"]._build_weight_matrix_single_alpha(
            d["ordered_names"], d["lineage_info"], d["regimes"],
            alpha, d["root_regime"],
        )
        original = ouwie_module.cho_solve
        calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original(*args, **kwargs)

        monkeypatch.setattr(ouwie_module, "cho_solve", counting_cho_solve)

        theta, ll = d["svc"]._ou_gls_log_likelihood_cholesky(d["x"], W, V)

        assert calls == 1
        assert np.all(np.isfinite(theta))
        assert np.isfinite(ll)


# ── TestOUM ──────────────────────────────────────────────────────────

class TestOUM:
    def test_fit_oum_reuses_single_alpha_weight_cache(
        self, precomputed, monkeypatch
    ):
        d = precomputed

        def fail_direct_builder(*args, **kwargs):
            raise AssertionError("OUM should reuse the single-alpha weight cache")

        monkeypatch.setattr(
            d["svc"], "_build_weight_matrix_single_alpha", fail_direct_builder
        )

        oum = d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )

        assert np.isfinite(oum["log_likelihood"])

    def test_ll_ge_ou1(self, precomputed):
        d = precomputed
        ou1 = d["svc"]._fit_ou1(
            d["x"], d["ordered_names"], d["lineage_info"], d["tree_height"],
        )
        oum = d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oum["log_likelihood"] >= ou1["log_likelihood"] - 1e-6

    def test_per_regime_theta(self, precomputed):
        d = precomputed
        oum = d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert isinstance(oum["params"]["theta"], dict)
        assert len(oum["params"]["theta"]) == len(d["regimes"])


# ── TestOUMV ─────────────────────────────────────────────────────────

class TestOUMV:
    def test_fit_oumv_reuses_single_alpha_weight_cache(
        self, precomputed, monkeypatch
    ):
        d = precomputed

        def fail_direct_builder(*args, **kwargs):
            raise AssertionError("OUMV should reuse the single-alpha weight cache")

        monkeypatch.setattr(
            d["svc"], "_build_weight_matrix_single_alpha", fail_direct_builder
        )

        oumv = d["svc"]._fit_oumv(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )

        assert np.isfinite(oumv["log_likelihood"])

    def test_ll_ge_oum(self, precomputed):
        d = precomputed
        oum = d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        oumv = d["svc"]._fit_oumv(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oumv["log_likelihood"] >= oum["log_likelihood"] - 1e-6

    def test_per_regime_sigma2(self, precomputed):
        d = precomputed
        oumv = d["svc"]._fit_oumv(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert isinstance(oumv["params"]["sigma2"], dict)
        for r in d["regimes"]:
            assert oumv["params"]["sigma2"][r] > 0


# ── TestOUMA ─────────────────────────────────────────────────────────

class TestOUMA:
    def test_ll_ge_oum(self, precomputed):
        d = precomputed
        oum = d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        ouma = d["svc"]._fit_ouma(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert ouma["log_likelihood"] >= oum["log_likelihood"] - 1e-6

    def test_per_regime_alpha(self, precomputed):
        d = precomputed
        ouma = d["svc"]._fit_ouma(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert isinstance(ouma["params"]["alpha"], dict)
        for r in d["regimes"]:
            assert ouma["params"]["alpha"][r] > 0


# ── TestOUMVA ────────────────────────────────────────────────────────

class TestOUMVA:
    def test_ll_ge_oumv_and_ouma(self, precomputed):
        d = precomputed
        oumv = d["svc"]._fit_oumv(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        ouma = d["svc"]._fit_ouma(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        oumva = d["svc"]._fit_oumva(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oumva["log_likelihood"] >= oumv["log_likelihood"] - 1e-6
        assert oumva["log_likelihood"] >= ouma["log_likelihood"] - 1e-6

    def test_k_params_equals_3R(self, precomputed):
        d = precomputed
        R = len(d["regimes"])
        oumva = d["svc"]._fit_oumva(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oumva["k_params"] == 3 * R


# ── TestModelNesting ─────────────────────────────────────────────────

class TestModelNesting:
    """Verify that nested models have LL(simpler) <= LL(complex)."""

    def test_bm1_le_bms(self, precomputed):
        d = precomputed
        bm1 = d["svc"]._fit_bm1(d["x"], d["vcv_total"])
        bms = d["svc"]._fit_bms(d["x"], d["per_regime_vcv"], d["regimes"])
        assert bms["log_likelihood"] >= bm1["log_likelihood"] - 1e-6

    def test_ou1_le_oum(self, precomputed):
        d = precomputed
        ou1 = d["svc"]._fit_ou1(
            d["x"], d["ordered_names"], d["lineage_info"], d["tree_height"],
        )
        oum = d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oum["log_likelihood"] >= ou1["log_likelihood"] - 1e-6

    def test_oum_le_oumv(self, precomputed):
        d = precomputed
        oum = d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        oumv = d["svc"]._fit_oumv(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oumv["log_likelihood"] >= oum["log_likelihood"] - 1e-6

    def test_oum_le_ouma(self, precomputed):
        d = precomputed
        oum = d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        ouma = d["svc"]._fit_ouma(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert ouma["log_likelihood"] >= oum["log_likelihood"] - 1e-6

    def test_oumv_le_oumva(self, precomputed):
        d = precomputed
        oumv = d["svc"]._fit_oumv(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        oumva = d["svc"]._fit_oumva(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oumva["log_likelihood"] >= oumv["log_likelihood"] - 1e-6

    def test_ouma_le_oumva(self, precomputed):
        d = precomputed
        ouma = d["svc"]._fit_ouma(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        oumva = d["svc"]._fit_oumva(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        )
        assert oumva["log_likelihood"] >= ouma["log_likelihood"] - 1e-6


# ── TestModelComparison ──────────────────────────────────────────────

class TestModelComparison:
    def test_aicc_weights_sum_to_one(self, precomputed, monkeypatch):
        d = precomputed
        results = []
        results.append(d["svc"]._fit_bm1(d["x"], d["vcv_total"]))
        results.append(d["svc"]._fit_ou1(
            d["x"], d["ordered_names"], d["lineage_info"], d["tree_height"],
        ))
        results.append(d["svc"]._fit_oum(
            d["x"], d["ordered_names"], d["lineage_info"],
            d["regimes"], d["root_regime"], d["tree_height"],
        ))
        n = len(d["x"])
        monkeypatch.setattr(
            ouwie_module.np,
            "exp",
            lambda *args, **kwargs: pytest.fail("AICc weights should use math.exp"),
        )
        results = d["svc"]._compute_model_comparison(results, n)
        weights = [r["aicc_weight"] for r in results]
        np.testing.assert_allclose(sum(weights), 1.0, atol=1e-10)

    def test_sorted_by_aicc(self, precomputed):
        d = precomputed
        results = []
        results.append(d["svc"]._fit_bm1(d["x"], d["vcv_total"]))
        results.append(d["svc"]._fit_ou1(
            d["x"], d["ordered_names"], d["lineage_info"], d["tree_height"],
        ))
        n = len(d["x"])
        results = d["svc"]._compute_model_comparison(results, n)
        aiccs = [r["aicc"] for r in results]
        assert aiccs == sorted(aiccs)

    def test_multiregime_r_squared_uses_unweighted_sigma2_average(self, precomputed):
        d = precomputed
        results = [
            {
                "model": "BM1",
                "k_params": 2,
                "log_likelihood": -10.0,
                "params": {"sigma2": 4.0},
            },
            {
                "model": "OUMV",
                "k_params": 4,
                "log_likelihood": -9.0,
                "params": {"sigma2": {"r1": 1.0, "r2": 3.0}},
            },
        ]

        compared = d["svc"]._compute_model_comparison(results, n=20)
        oumv = next(r for r in compared if r["model"] == "OUMV")

        assert oumv["r_squared"] == 0.5

    def test_all_lls_finite(self, precomputed):
        d = precomputed
        results = []
        for model in ALL_MODELS:
            res = d["svc"]._fit_model(
                model, d["x"], d["vcv_total"], d["per_regime_vcv"],
                d["ordered_names"], d["regimes"], d["lineage_info"],
                d["root_regime"], d["tree_height"], d["n_regimes"],
            )
            results.append(res)
        for r in results:
            assert np.isfinite(r["log_likelihood"]), f"{r['model']} LL not finite"


# ── TestRun ──────────────────────────────────────────────────────────

class TestRun:
    @patch("builtins.print")
    def test_text_output(self, mocked_print, svc):
        svc.run()
        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "OUwie Model Comparison" in all_output
        assert "BM1" in all_output
        assert "OUM" in all_output

    def test_print_text_output_batches_model_and_parameter_rows(self, svc, mocker):
        printed = mocker.patch("builtins.print")
        results = [
            {
                "model": "OUM",
                "k_params": 8,
                "log_likelihood": -12.3456,
                "aic": 30.1,
                "aicc": 31.2,
                "delta_aicc": 0.0,
                "aicc_weight": 0.75,
                "bic": 35.4,
                "delta_bic": 0.0,
                "r_squared": 0.8,
                "params": {
                    "theta": {"aquatic": 1.25, "terrestrial": 2.5},
                    "z0": 1.1,
                    "alpha": {"aquatic": 0.3, "terrestrial": 0.4},
                    "sigma2": {"aquatic": 0.05, "terrestrial": 0.06},
                },
            },
            {
                "model": "BM1",
                "k_params": 2,
                "log_likelihood": -18.0,
                "aic": 40.0,
                "aicc": 41.0,
                "delta_aicc": 9.8,
                "aicc_weight": 0.01,
                "bic": 42.0,
                "delta_bic": 6.6,
                "r_squared": 0.1,
                "params": {},
            },
        ]

        svc._print_text_output(results, 8, ["aquatic", "terrestrial"])

        printed.assert_called_once_with(
            "OUwie Model Comparison (Multi-Regime OU)\n"
            "\n"
            "Number of tips: 8\n"
            "Regimes: 2 (aquatic, terrestrial)\n"
            "\n"
            "Model   k    LL          AIC       AICc      dAICc    AICcW    "
            "BIC       dBIC     R2     \n"
            "OUM     8    -12.346     30.10     31.20     0.00     0.750    "
            "35.40     0.00     0.800  \n"
            "BM1     2    -18.000     40.00     41.00     9.80     0.010    "
            "42.00     6.60     0.100  \n"
            "\n"
            "Best model (AICc): OUM\n"
            "Best model (BIC):  OUM\n"
            "\n"
            "Parameter estimates (OUM):\n"
            "  Theta (trait optima):\n"
            "    aquatic: 1.2500\n"
            "    terrestrial: 2.5000\n"
            "  z0: 1.1000\n"
            "  Alpha (selection strength):\n"
            "    aquatic: 0.3000\n"
            "    terrestrial: 0.4000\n"
            "  Sigma-squared (diffusion rate):\n"
            "    aquatic: 0.0500\n"
            "    terrestrial: 0.0600"
        )

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE, models=None, json=True,
        )
        svc = OUwie(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_tips"] == 8
        assert "models" in payload
        assert "BM1" in payload["models"]
        assert "best_model_aicc" in payload

    def test_run_uses_fast_tip_name_helper_for_prune_setup(self, default_args, mocker):
        svc = OUwie(default_args)
        spy = mocker.spy(svc, "get_tip_names_from_tree")
        mocker.patch("builtins.print")

        svc.run()

        assert spy.call_count == 1

    def test_run_uses_unmodified_tree_read(self, default_args, mocker, monkeypatch):
        svc = OUwie(default_args)
        tree = OUwie(default_args).read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        traits = {name: float(index) for index, name in enumerate(tree_tips)}
        regimes_by_tip = {
            name: ("x" if index % 2 == 0 else "y")
            for index, name in enumerate(tree_tips)
        }
        per_regime_vcv = {
            "x": np.array([[1.0, 0.1], [0.1, 3.0]]),
            "y": np.array([[2.0, 0.2], [0.2, 1.0]]),
        }
        fit_result = {
            "model": "BM1",
            "log_likelihood": -1.0,
            "aic": 2.0,
            "aicc": 2.5,
            "bic": 3.0,
            "n_params": 2,
        }

        read_unmodified = mocker.patch.object(
            svc, "read_tree_file_unmodified", return_value=tree
        )
        mocker.patch.object(
            svc,
            "read_tree_file",
            side_effect=AssertionError("run should use read_tree_file_unmodified"),
        )
        validate = mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(svc, "get_tip_names_from_tree", return_value=tree_tips)
        mocker.patch.object(svc, "_parse_trait_file", return_value=traits)
        mocker.patch.object(svc, "_parse_regime_file", return_value=regimes_by_tip)
        mocker.patch.object(svc, "_build_parent_map", return_value={})
        mocker.patch.object(
            svc, "_assign_branch_regimes_and_root", return_value=({}, "x")
        )
        lineage_info = mocker.patch.object(svc, "_build_lineage_info", return_value={})
        mocker.patch.object(svc, "_build_per_regime_vcv", return_value=per_regime_vcv)
        fit_model = mocker.patch.object(svc, "_fit_model", return_value=fit_result)
        mocker.patch.object(
            svc, "_compute_model_comparison", side_effect=lambda results, *_args: results
        )
        mocker.patch.object(svc, "_print_text_output")
        fast_copy = mocker.patch.object(
            svc,
            "_fast_copy",
            side_effect=AssertionError("all-shared analysis should not copy tree"),
        )
        monkeypatch.setattr(
            ouwie_module.np,
            "diag",
            lambda _matrix: (_ for _ in ()).throw(
                AssertionError("tree height should read the VCV diagonal view")
            ),
        )

        svc.run()

        read_unmodified.assert_called_once_with()
        validate.assert_called_once()
        fast_copy.assert_not_called()
        assert lineage_info.call_args.args[0] is tree
        assert fit_model.call_args.args[8] == pytest.approx(4.0)

    def test_run_copies_before_pruning_missing_tree_tips(self, default_args, mocker):
        svc = OUwie(default_args)
        tree = OUwie(default_args).read_tree_file()
        tree_copy = OUwie(default_args).read_tree_file()
        pruned_tree = OUwie(default_args).read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        shared_tips = tree_tips[:-1]
        traits = {name: float(index) for index, name in enumerate(shared_tips)}
        regimes_by_tip = {
            name: ("x" if index % 2 == 0 else "y")
            for index, name in enumerate(shared_tips)
        }
        fit_result = {
            "model": "BM1",
            "log_likelihood": -1.0,
            "aic": 2.0,
            "aicc": 2.5,
            "bic": 3.0,
            "n_params": 2,
        }

        mocker.patch.object(svc, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(svc, "get_tip_names_from_tree", return_value=tree_tips)
        mocker.patch.object(svc, "_parse_trait_file", return_value=traits)
        mocker.patch.object(svc, "_parse_regime_file", return_value=regimes_by_tip)
        fast_copy = mocker.patch.object(svc, "_fast_copy", return_value=tree_copy)
        prune = mocker.patch.object(
            svc, "prune_tree_using_taxa_list", return_value=pruned_tree
        )
        mocker.patch.object(svc, "_build_parent_map", return_value={})
        mocker.patch.object(
            svc, "_assign_branch_regimes_and_root", return_value=({}, "x")
        )
        lineage_info = mocker.patch.object(svc, "_build_lineage_info", return_value={})
        mocker.patch.object(
            svc, "_build_per_regime_vcv", return_value={"x": np.eye(2), "y": np.eye(2)}
        )
        mocker.patch.object(svc, "_fit_model", return_value=fit_result)
        mocker.patch.object(
            svc, "_compute_model_comparison", side_effect=lambda results, *_args: results
        )
        mocker.patch.object(svc, "_print_text_output")

        svc.run()

        fast_copy.assert_called_once_with(tree)
        prune.assert_called_once_with(tree_copy, [tree_tips[-1]])
        assert lineage_info.call_args.args[0] is pruned_tree

    @patch("builtins.print")
    def test_subset_models(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE, models="BM1,OUM", json=True,
        )
        svc = OUwie(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert set(payload["models"].keys()) == {"BM1", "OUM"}


# ── TestEffectSize ──────────────────────────────────────────────────

class TestEffectSize:
    @patch("builtins.print")
    def test_r2_in_json(self, mocked_print):
        """All models have r_squared in JSON output."""
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE, models=None, json=True,
        )
        svc = OUwie(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        for model_name, model_data in payload["models"].items():
            assert "r_squared" in model_data, (
                f"r_squared missing from {model_name}"
            )

    @patch("builtins.print")
    def test_bm1_r2_zero(self, mocked_print):
        """BM1 is baseline so its R² should be exactly 0."""
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE, models=None, json=True,
        )
        svc = OUwie(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["models"]["BM1"]["r_squared"] == 0.0

    @patch("builtins.print")
    def test_r2_finite(self, mocked_print):
        """R² values are finite for BM1 and OU1."""
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE, models="BM1,OU1", json=True,
        )
        svc = OUwie(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        for model_name in ["BM1", "OU1"]:
            r2 = payload["models"][model_name]["r_squared"]
            assert np.isfinite(r2), (
                f"R² for {model_name} is not finite: {r2}"
            )

    @patch("builtins.print")
    def test_r2_without_bm1_in_models(self, mocked_print):
        """When BM1 is not in selected models, R² is still computed via silent fit.

        The design doc states: 'BM1 not in selected models (ouwie): Fit
        silently as baseline.' This verifies that excluding BM1 still
        produces finite R² values, not NaN.
        """
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE, models="OU1,OUM", json=True,
        )
        svc = OUwie(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert "BM1" not in payload["models"]
        for model_name, model_data in payload["models"].items():
            assert np.isfinite(model_data["r_squared"]), (
                f"R² for {model_name} is not finite when BM1 excluded"
            )
