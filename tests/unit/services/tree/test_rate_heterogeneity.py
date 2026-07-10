import builtins
import importlib
import json
import os
import subprocess
import sys
import tempfile
from io import StringIO

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, TreeMixin

import phykit.services.tree.rate_heterogeneity as rate_heterogeneity_module
from phykit.services.tree.rate_heterogeneity import RateHeterogeneity
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")
REGIMES_FILE = str(SAMPLE_FILES / "tree_simple_regimes.tsv")


def test_module_import_does_not_import_numpy_or_scipy():
    code = """
import sys
import phykit.services.tree.rate_heterogeneity as module

assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "pickle" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "phykit.helpers.circular_layout" not in sys.modules
assert "phykit.helpers.color_annotations" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy.linalg" not in sys.modules
assert "scipy.optimize" not in sys.modules
assert "scipy.stats" not in sys.modules
assert module._CHO_FACTOR is None
assert module._CHO_SOLVE is None
assert module._MINIMIZE is None
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = rate_heterogeneity_module._LazyNumpy()

    array_attr = lazy_np.array

    assert lazy_np.__dict__["array"] is array_attr
    assert lazy_np.array is array_attr
    assert lazy_np._module is not None


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

    merged = rate_heterogeneity_module._merge_nonempty_child_state_sets(
        state_sets,
        NoSliceList([left, middle, right, empty]),
    )

    assert merged == {"A", "B", "C"}
    assert (
        rate_heterogeneity_module._merge_nonempty_child_state_sets(
            state_sets,
            NoSliceList([left, middle]),
        )
        == {"A"}
    )
    assert (
        rate_heterogeneity_module._merge_nonempty_child_state_sets(
            state_sets,
            NoSliceList([empty]),
        )
        == set()
    )


def test_module_import_does_not_import_scipy_linalg_or_optimize(monkeypatch):
    module_name = "phykit.services.tree.rate_heterogeneity"
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
                "rate_heterogeneity module import should not import SciPy linalg/optimize"
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


def test_repeated_cholesky_calls_cache_scipy_linalg_imports(monkeypatch):
    previous_cho_factor = rate_heterogeneity_module._CHO_FACTOR
    previous_cho_solve = rate_heterogeneity_module._CHO_SOLVE
    rate_heterogeneity_module._CHO_FACTOR = None
    rate_heterogeneity_module._CHO_SOLVE = None
    original_import = builtins.__import__
    scipy_linalg_imports = 0

    def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
        nonlocal scipy_linalg_imports
        if name == "scipy.linalg":
            scipy_linalg_imports += 1
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", counting_import)
    try:
        rng = np.random.default_rng(20260628)
        svc = RateHeterogeneity.__new__(RateHeterogeneity)
        n = 12
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        y = rng.normal(size=n)

        svc._fit_single_rate_cholesky(y, C)
        first_call_imports = scipy_linalg_imports
        svc._fit_single_rate_cholesky(y, C)
    finally:
        rate_heterogeneity_module._CHO_FACTOR = previous_cho_factor
        rate_heterogeneity_module._CHO_SOLVE = previous_cho_solve

    assert first_call_imports > 0
    assert scipy_linalg_imports == first_call_imports


def test_repeated_minimize_calls_cache_scipy_optimize_imports(monkeypatch):
    previous_minimize = rate_heterogeneity_module._MINIMIZE
    rate_heterogeneity_module._MINIMIZE = None
    original_import = builtins.__import__
    scipy_optimize_imports = 0

    def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
        nonlocal scipy_optimize_imports
        if name == "scipy.optimize":
            scipy_optimize_imports += 1
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", counting_import)
    try:
        first = rate_heterogeneity_module.minimize(
            lambda values: (values[0] - 0.25) ** 2,
            x0=np.array([0.0]),
            method="Nelder-Mead",
        )
        first_call_imports = scipy_optimize_imports
        second = rate_heterogeneity_module.minimize(
            lambda values: (values[0] - 0.75) ** 2,
            x0=np.array([1.0]),
            method="Nelder-Mead",
        )
    finally:
        rate_heterogeneity_module._MINIMIZE = previous_minimize

    assert first.x[0] == pytest.approx(0.25, abs=1e-4)
    assert second.x[0] == pytest.approx(0.75, abs=1e-4)
    assert first_call_imports > 0
    assert scipy_optimize_imports == first_call_imports


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        regime_data=REGIMES_FILE,
        nsim=0,
        seed=None,
        plot=None,
        json=False,
    )


@pytest.fixture
def json_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        regime_data=REGIMES_FILE,
        nsim=0,
        seed=None,
        plot=None,
        json=True,
    )


@pytest.fixture
def sim_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        regime_data=REGIMES_FILE,
        nsim=10,
        seed=42,
        plot=None,
        json=False,
    )


@pytest.fixture
def service(default_args):
    return RateHeterogeneity(default_args)


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            regime_data="r.tsv",
            nsim=0,
            seed=None,
            plot=None,
            json=False,
        )
        svc = RateHeterogeneity.__new__(RateHeterogeneity)
        parsed = svc.process_args(args)
        assert parsed["tree_file_path"] == "t.tre"
        assert parsed["trait_data_path"] == "d.tsv"
        assert parsed["regime_data_path"] == "r.tsv"
        assert parsed["n_sim"] == 0
        assert parsed["seed"] is None
        assert parsed["json_output"] is False

    def test_overrides(self):
        args = Namespace(
            tree="t.tre",
            trait_data="d.tsv",
            regime_data="r.tsv",
            nsim=100,
            seed=42,
            plot="out.png",
            json=True,
        )
        svc = RateHeterogeneity.__new__(RateHeterogeneity)
        parsed = svc.process_args(args)
        assert parsed["n_sim"] == 100
        assert parsed["seed"] == 42
        assert parsed["plot_output"] == "out.png"
        assert parsed["json_output"] is True


class TestRegimeParsing:
    def test_valid_file(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        regimes = service._parse_regime_file(REGIMES_FILE, tree_tips)
        assert len(regimes) == 8
        assert regimes["sea_lion"] == "aquatic"
        assert regimes["raccoon"] == "terrestrial"

    def test_comments_and_blanks(self, service, tmp_path):
        regime_file = tmp_path / "regimes.tsv"
        regime_file.write_text(
            "# comment\n\nraccoon\tA\nbear\tA\nsea_lion\tB\n"
            "seal\tB\nmonkey\tA\ncat\tA\nweasel\tA\ndog\tA\n"
        )
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        regimes = service._parse_regime_file(str(regime_file), tree_tips)
        assert len(regimes) == 8
        assert regimes["raccoon"] == "A"

    def test_ordered_all_shared_regime_file_skips_sets(
        self, service, tmp_path, monkeypatch
    ):
        regime_file = tmp_path / "regimes.tsv"
        regime_file.write_text("A\tr1\nB\tr2\nC\tr1\n")

        def fail_set(*args, **kwargs):
            raise AssertionError("ordered exact regime path should not build sets")

        monkeypatch.setattr(builtins, "set", fail_set)
        regimes = service._parse_regime_file(str(regime_file), ["A", "B", "C"])

        assert regimes == {"A": "r1", "B": "r2", "C": "r1"}
        assert builtins.set is fail_set

    def test_missing_file(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        with pytest.raises(PhykitUserError):
            service._parse_regime_file("/nonexistent/file.tsv", tree_tips)

    def test_extra_columns_error(self, service, tmp_path):
        regime_file = tmp_path / "regimes.tsv"
        regime_file.write_text(
            "raccoon\tA\nbear\tA\nsea_lion\tB\ndog\tA\nextra\tA\tB\n"
        )
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            service._parse_regime_file(str(regime_file), tree_tips)

        assert (
            "Line 5 in regime file has 3 columns; expected 2."
            in exc_info.value.messages
        )

    def test_missing_tab_regime_error_reports_one_column(self, service, tmp_path):
        regime_file = tmp_path / "regimes.tsv"
        regime_file.write_text("raccoon\tA\nbear\nsea_lion\tB\ndog\tA\n")
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            service._parse_regime_file(str(regime_file), tree_tips)

        assert (
            "Line 2 in regime file has 1 columns; expected 2."
            in exc_info.value.messages
        )

    def test_equal_sized_taxon_mismatch_still_warns(self, service, tmp_path, capsys):
        regime_file = tmp_path / "regimes.tsv"
        regime_file.write_text(
            "raccoon\tA\nbear\tA\nsea_lion\tB\nseal\tB\nmonkey\tA\ncat\tA\n"
            "weasel\tA\noff_tree\tA\n"
        )
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)

        regimes = service._parse_regime_file(str(regime_file), tree_tips)

        stderr = capsys.readouterr().err
        assert set(regimes) == set(tree_tips) - {"dog"}
        assert "1 taxa in tree but not in regime file: dog" in stderr
        assert "1 taxa in regime file but not in tree: off_tree" in stderr


class TestTraitParsing:
    def test_comments_and_blanks(self, service, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "   # comment\n\nraccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\n"
            "seal\t4.0\nmonkey\t5.0\ncat\t6.0\nweasel\t7.0\ndog\t8.0\n"
        )
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        traits = service._parse_trait_file(str(trait_file), tree_tips)
        assert len(traits) == 8
        assert traits["raccoon"] == pytest.approx(1.0)

    def test_all_shared_trait_file_emits_no_warnings(
        self, service, tmp_path, capsys
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "raccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\nseal\t4.0\n"
            "monkey\t5.0\ncat\t6.0\nweasel\t7.0\ndog\t8.0\n"
        )
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)

        traits = service._parse_trait_file(str(trait_file), tree_tips)

        stderr = capsys.readouterr().err
        assert set(traits) == set(tree_tips)
        assert stderr == ""

    def test_ordered_all_shared_trait_file_skips_sets(
        self, service, tmp_path, monkeypatch
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("A\t1.0\nB\t2.0\nC\t3.0\n")

        def fail_set(*args, **kwargs):
            raise AssertionError("ordered exact trait path should not build sets")

        monkeypatch.setattr(builtins, "set", fail_set)
        traits = service._parse_trait_file(str(trait_file), ["A", "B", "C"])

        assert traits == {"A": 1.0, "B": 2.0, "C": 3.0}
        assert builtins.set is fail_set

    def test_extra_columns_error(self, service, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "raccoon\t1.0\nbear\t2.0\nsea_lion\t3.0\ndog\t4.0\nextra\t1.0\t2.0\n"
        )
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            service._parse_trait_file(str(trait_file), tree_tips)

        assert (
            "Line 5 in trait file has 3 columns; expected 2."
            in exc_info.value.messages
        )

    def test_missing_tab_error_reports_one_column(self, service, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("raccoon\t1.0\nbear\nsea_lion\t3.0\ndog\t4.0\n")
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            service._parse_trait_file(str(trait_file), tree_tips)

        assert (
            "Line 2 in trait file has 1 columns; expected 2."
            in exc_info.value.messages
        )

    def test_non_numeric_trait_error(self, service, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text("raccoon\t1.0\nbear\tbad\nsea_lion\t3.0\ndog\t4.0\n")
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)

        with pytest.raises(PhykitUserError) as exc_info:
            service._parse_trait_file(str(trait_file), tree_tips)

        assert (
            "Non-numeric trait value 'bad' for taxon 'bear' on line 2."
            in exc_info.value.messages
        )


class TestFitchParsimony:
    def test_iter_postorder_matches_biopython_order(self, service):
        tree = Phylo.read(
            StringIO("((A:1,B:2,C:3):4,(D:5,E:6):7,F:8):0;"),
            "newick",
        )

        direct = list(service._iter_postorder(tree.root))
        reference = list(tree.find_clades(order="postorder"))

        assert [id(clade) for clade in direct] == [
            id(clade) for clade in reference
        ]

    def test_setup_helpers_use_direct_tree_traversal(self, service, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        regimes = ["aquatic", "terrestrial"]
        regime_assignments = {
            "A": "aquatic",
            "B": "aquatic",
            "C": "terrestrial",
            "D": "terrestrial",
        }

        def fail_traversal(*args, **kwargs):
            raise AssertionError("optimized path should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        parent_map = service._build_parent_map(tree)
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        paths = service._build_root_to_tip_paths(tree, ordered_names, parent_map)
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )

        assert len(parent_map) == 6
        assert len(branch_regimes) == 6
        assert sum(bl for _, bl in paths["A"]) == pytest.approx(2.0)
        assert all(matrix.shape == (4, 4) for matrix in per_regime_vcv.values())

    def test_parent_map_handles_mixed_child_counts_without_preorder(
        self, service, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_iter_preorder(*args, **kwargs):
            raise AssertionError("parent map should build directly")

        monkeypatch.setattr(service, "_iter_preorder", fail_iter_preorder)

        parent_map = service._build_parent_map(tree)

        terminal, binary, trifurcating = tree.root.clades
        assert id(tree.root) not in parent_map
        assert parent_map[id(terminal)] is tree.root
        assert parent_map[id(binary)] is tree.root
        assert parent_map[id(trifurcating)] is tree.root
        assert all(parent_map[id(child)] is binary for child in binary.clades)
        assert all(
            parent_map[id(child)] is trifurcating for child in trifurcating.clades
        )

    def test_iter_preorder_preserves_order_without_reversed(self, service):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("_iter_preorder should push children directly")

        root = Clade(name="root")
        left = Clade(name="left")
        middle = Clade(name="middle")
        right = Clade(name="right")
        left.clades = NoReversedList([Clade(name="left_a"), Clade(name="left_b")])
        middle.clades = NoReversedList([Clade(name="middle_a")])
        right.clades = NoReversedList()
        root.clades = NoReversedList([left, middle, right])

        order = [clade.name for clade in service._iter_preorder(root)]

        assert order == [
            "root",
            "left",
            "left_a",
            "left_b",
            "middle",
            "middle_a",
            "right",
        ]

    def test_branch_regimes_handle_mixed_child_counts(self, service):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )
        regime_assignments = {
            "A": "r1",
            "B": "r2",
            "C": "r2",
            "D": "r3",
            "E": "r3",
            "F": "r1",
        }
        parent_map = service._build_parent_map(tree)

        branch_regimes = service._assign_branch_regimes(
            tree,
            regime_assignments,
            parent_map,
        )

        terminal, binary, trifurcating = tree.root.clades
        assert branch_regimes[id(terminal)] == "r1"
        assert branch_regimes[id(binary)] == "r2"
        assert branch_regimes[id(binary.clades[0])] == "r2"
        assert branch_regimes[id(binary.clades[1])] == "r2"
        assert branch_regimes[id(trifurcating)] == "r1"
        assert branch_regimes[id(trifurcating.clades[0])] == "r3"
        assert branch_regimes[id(trifurcating.clades[1])] == "r3"
        assert branch_regimes[id(trifurcating.clades[2])] == "r1"

    def test_branch_regime_ambiguity_uses_lexicographic_tiebreak(self, service):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        regime_assignments = {
            "A": "zeta",
            "B": "beta",
            "C": "alpha",
            "D": "alpha",
        }
        parent_map = service._build_parent_map(tree)

        branch_regimes = service._assign_branch_regimes(
            tree,
            regime_assignments,
            parent_map,
        )

        ambiguous_child = tree.root.clades[0]
        alpha_child = tree.root.clades[1]
        assert branch_regimes[id(ambiguous_child)] == "beta"
        assert branch_regimes[id(alpha_child)] == "alpha"

    def test_tip_states_preserved(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        parent_map = service._build_parent_map(tree)
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        # Tips should have their correct regime
        for tip in tree.get_terminals():
            if tip.name in regime_assignments:
                assert branch_regimes[id(tip)] == regime_assignments[tip.name]

    def test_correct_branch_assignments(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        parent_map = service._build_parent_map(tree)
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        # All non-root nodes should have a regime
        for clade in tree.find_clades(order="preorder"):
            if clade != tree.root:
                assert id(clade) in branch_regimes
                assert branch_regimes[id(clade)] in {"aquatic", "terrestrial"}


class TestPerRegimeVCV:
    def test_build_per_regime_vcv_matches_pairwise_reference(
        self, service, monkeypatch
    ):
        ordered_names = ["A", "B", "C"]
        regimes = ["aquatic", "terrestrial"]
        paths = {
            "A": [(1, 1.0), (2, 0.5)],
            "B": [(1, 1.0), (3, 0.75)],
            "C": [(4, 1.25)],
        }
        branch_regimes = {
            1: "aquatic",
            2: "terrestrial",
            3: "aquatic",
            4: "terrestrial",
        }

        monkeypatch.setattr(
            service,
            "_build_root_to_tip_paths",
            lambda *args, **kwargs: paths,
        )

        expected = {regime: np.zeros((3, 3)) for regime in regimes}
        for i, name_i in enumerate(ordered_names):
            path_i = paths[name_i]
            for clade_id, branch_length in path_i:
                regime = branch_regimes.get(clade_id, regimes[0])
                expected[regime][i, i] += branch_length
            for j in range(i + 1, len(ordered_names)):
                path_j = paths[ordered_names[j]]
                for (clade_i, branch_i), (clade_j, _) in zip(path_i, path_j):
                    if clade_i != clade_j:
                        break
                    regime = branch_regimes.get(clade_i, regimes[0])
                    expected[regime][i, j] += branch_i
                    expected[regime][j, i] += branch_i

        observed = service._build_per_regime_vcv(
            tree=None,
            ordered_names=ordered_names,
            regimes=regimes,
            branch_regimes=branch_regimes,
            parent_map={},
        )

        for regime in regimes:
            np.testing.assert_allclose(observed[regime], expected[regime])

    def test_single_tip_branches_skip_block_indexing(self, service, monkeypatch):
        ordered_names = ["A", "B", "C"]
        regimes = ["aquatic", "terrestrial"]
        paths = {
            "A": [("A_tip", 1.0)],
            "B": [("B_tip", 2.0)],
            "C": [("C_tip", 3.0)],
        }
        branch_regimes = {
            "A_tip": "aquatic",
            "B_tip": "terrestrial",
            "C_tip": "aquatic",
        }

        def fail_ix(*_args, **_kwargs):
            raise AssertionError("single-tip branches should update diagonals directly")

        monkeypatch.setattr(
            service,
            "_build_root_to_tip_paths",
            lambda *args, **kwargs: paths,
        )
        monkeypatch.setattr(rate_heterogeneity_module.np, "ix_", fail_ix)

        observed = service._build_per_regime_vcv(
            tree=None,
            ordered_names=ordered_names,
            regimes=regimes,
            branch_regimes=branch_regimes,
            parent_map={},
        )

        np.testing.assert_allclose(
            observed["aquatic"], np.diag([1.0, 0.0, 3.0])
        )
        np.testing.assert_allclose(
            observed["terrestrial"], np.diag([0.0, 2.0, 0.0])
        )

    def test_sum_equals_total(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )

        # Sum of per-regime VCVs should approximate total VCV
        C_sum = sum(per_regime_vcv.values())

        # Build total VCV using root-to-tip distances
        n = len(ordered_names)
        C_total = np.zeros((n, n))
        root_to_tip = {}
        for name in ordered_names:
            root_to_tip[name] = tree.distance(tree.root, name)
        for i in range(n):
            for j in range(i, n):
                if i == j:
                    C_total[i, j] = root_to_tip[ordered_names[i]]
                else:
                    d_ij = tree.distance(ordered_names[i], ordered_names[j])
                    shared_path = (
                        root_to_tip[ordered_names[i]]
                        + root_to_tip[ordered_names[j]]
                        - d_ij
                    ) / 2.0
                    C_total[i, j] = shared_path
                    C_total[j, i] = shared_path

        np.testing.assert_allclose(C_sum, C_total, rtol=1e-6)

    def test_individual_values_nonnegative(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )

        for r in regimes:
            assert np.all(per_regime_vcv[r] >= -1e-10)


class TestSingleRate:
    def test_sigma2_positive(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        C_total = sum(per_regime_vcv.values())

        sigma2, anc, ll = service._fit_single_rate(y, C_total)
        assert sigma2 > 0

    def test_log_likelihood_finite(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        C_total = sum(per_regime_vcv.values())

        sigma2, anc, ll = service._fit_single_rate(y, C_total)
        assert np.isfinite(ll)

    def test_cholesky_matches_inverse(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        C_total = sum(per_regime_vcv.values())

        fast = service._fit_single_rate_cholesky(y, C_total)
        inverse = service._fit_single_rate_inverse(y, C_total)

        np.testing.assert_allclose(fast, inverse)

    def test_cholesky_uses_single_solve(self, service, monkeypatch):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        C_total = sum(per_regime_vcv.values())
        original = rate_heterogeneity_module.cho_solve
        calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original(*args, **kwargs)

        monkeypatch.setattr(
            rate_heterogeneity_module, "cho_solve", counting_cho_solve
        )

        fast = service._fit_single_rate_cholesky(y, C_total)
        inverse = service._fit_single_rate_inverse(y, C_total)

        assert calls == 1
        np.testing.assert_allclose(fast, inverse)


class TestMultiRate:
    def test_log_likelihood_geq_single(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        C_total = sum(per_regime_vcv.values())

        _, _, ll_single = service._fit_single_rate(y, C_total)
        _, _, ll_multi = service._fit_multi_rate(y, per_regime_vcv, regimes)
        # Multi-rate model should fit at least as well
        assert ll_multi >= ll_single - 1e-6

    def test_per_regime_sigma2_positive(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )

        sigma2s, _, _ = service._fit_multi_rate(y, per_regime_vcv, regimes)
        assert np.all(sigma2s > 0)

    def test_log_likelihood_cholesky_matches_inverse(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        V = np.zeros((len(y), len(y)))
        for sigma2, regime in zip((0.7, 1.4), regimes):
            V += sigma2 * per_regime_vcv[regime]
        ones = np.ones(len(y))

        fast = service._multi_rate_log_likelihood_cholesky(y, V, ones)
        inverse = service._multi_rate_log_likelihood_inverse(y, V, ones)

        np.testing.assert_allclose(fast, inverse)

    def test_log_likelihood_cholesky_uses_single_solve(self, service, monkeypatch):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        V = np.zeros((len(y), len(y)))
        for sigma2, regime in zip((0.7, 1.4), regimes):
            V += sigma2 * per_regime_vcv[regime]
        ones = np.ones(len(y))
        original = rate_heterogeneity_module.cho_solve
        calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original(*args, **kwargs)

        monkeypatch.setattr(
            rate_heterogeneity_module, "cho_solve", counting_cho_solve
        )

        fast = service._multi_rate_log_likelihood_cholesky(y, V, ones)
        inverse = service._multi_rate_log_likelihood_inverse(y, V, ones)

        assert calls == 1
        np.testing.assert_allclose(fast, inverse)


class TestLRT:
    def test_chi2_p_in_range(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        C_total = sum(per_regime_vcv.values())

        _, _, ll_single = service._fit_single_rate(y, C_total)
        _, _, ll_multi = service._fit_multi_rate(y, per_regime_vcv, regimes)

        lrt_stat = max(2.0 * (ll_multi - ll_single), 0.0)
        from scipy.stats import chi2
        p = float(chi2.sf(lrt_stat, df=len(regimes) - 1))
        assert 0.0 <= p <= 1.0

    def test_lrt_stat_nonnegative(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        C_total = sum(per_regime_vcv.values())

        _, _, ll_single = service._fit_single_rate(y, C_total)
        _, _, ll_multi = service._fit_multi_rate(y, per_regime_vcv, regimes)

        lrt_stat = max(2.0 * (ll_multi - ll_single), 0.0)
        assert lrt_stat >= 0


class TestRValidation:
    """Validate against R 4.4.0 phytools::brownie.lite with stem=TRUE.

    R code:
        library(phytools)
        tree <- read.tree("tree_simple.tre")
        x <- setNames(c(1.04,2.39,2.30,1.88,0.60,0.56,-0.30,1.18),
                       c("raccoon","bear","sea_lion","seal","monkey","cat","weasel","dog"))
        aquatic_mrca <- findMRCA(tree, c("sea_lion", "seal"))
        tree2 <- paintSubTree(tree, node=aquatic_mrca, state="aquatic",
                              anc.state="terrestrial", stem=TRUE)
        brownie.lite(tree2, x)

    PhyKIT uses Fitch parsimony for regime assignment, which matches
    R's paintSubTree(stem=TRUE) for this tree topology.
    """

    def _get_results(self, service):
        tree = service.read_tree_file()
        tree_tips = service.get_tip_names_from_tree(tree)
        trait_values = service._parse_trait_file(TRAITS_FILE, tree_tips)
        regime_assignments = service._parse_regime_file(REGIMES_FILE, tree_tips)
        ordered_names = sorted(trait_values.keys())
        y = np.array([trait_values[name] for name in ordered_names])

        parent_map = service._build_parent_map(tree)
        regimes = sorted(set(regime_assignments.values()))
        branch_regimes = service._assign_branch_regimes(
            tree, regime_assignments, parent_map
        )
        per_regime_vcv = service._build_per_regime_vcv(
            tree, ordered_names, regimes, branch_regimes, parent_map
        )
        C_total = sum(per_regime_vcv.values())

        sigma2_s, anc_s, ll_s = service._fit_single_rate(y, C_total)
        sigma2_m, anc_m, ll_m = service._fit_multi_rate(y, per_regime_vcv, regimes)

        return dict(
            sigma2_single=sigma2_s,
            anc_single=anc_s,
            ll_single=ll_s,
            sigma2_multi=dict(zip(regimes, sigma2_m)),
            anc_multi=anc_m,
            ll_multi=ll_m,
        )

    def test_single_rate_sigma2(self, service):
        r = self._get_results(service)
        # R: 0.0384065703382948
        np.testing.assert_allclose(r["sigma2_single"], 0.0384065703, rtol=1e-6)

    def test_single_rate_log_likelihood(self, service):
        r = self._get_results(service)
        # R: -11.5696774876938
        np.testing.assert_allclose(r["ll_single"], -11.5696774877, rtol=1e-6)

    def test_multi_rate_log_likelihood(self, service):
        r = self._get_results(service)
        # R: -11.204608463063
        np.testing.assert_allclose(r["ll_multi"], -11.2046, rtol=1e-3)


class TestRun:
    @patch("builtins.print")
    def test_text_output(self, mocked_print, default_args):
        svc = RateHeterogeneity(default_args)
        svc.run()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Rate Heterogeneity Test" in all_output
        assert "Single-rate model" in all_output
        assert "Multi-rate model" in all_output
        assert "Likelihood ratio test" in all_output

    def test_print_text_output_batches_regime_rows(self, default_args, mocker):
        svc = RateHeterogeneity(default_args)
        printed = mocker.patch("builtins.print")

        svc._print_text_output(
            n_tips=8,
            regimes=["aquatic", "terrestrial"],
            sigma2_single=0.0384065703,
            anc_single=1.23456,
            ll_single=-11.569677,
            aic_single=27.139354,
            sigma2_multi_dict={"aquatic": 0.008811, "terrestrial": 0.050019},
            anc_multi=1.5,
            ll_multi=-11.204608,
            aic_multi=30.409216,
            lrt_stat=0.730138,
            df=1,
            chi2_p=0.3929,
            sim_p=0.25,
            n_sim_done=10,
            r2_regime=-0.0341227314,
        )

        printed.assert_called_once_with(
            "Rate Heterogeneity Test (Multi-rate Brownian Motion)\n"
            "\n"
            "Regimes: 2 (aquatic, terrestrial)\n"
            "Number of tips: 8\n"
            "\n"
            "Single-rate model (H0):\n"
            "  Sigma-squared:      0.0384\n"
            "  Ancestral state:    1.2346\n"
            "  Log-likelihood:     -11.57\n"
            "  AIC:                27.14\n"
            "\n"
            "Multi-rate model (H1):\n"
            "  Regime               Sigma-squared\n"
            "  aquatic                     0.0088\n"
            "  terrestrial                 0.0500\n"
            "  Ancestral state:    1.5000\n"
            "  Log-likelihood:     -11.20\n"
            "  AIC:                30.41\n"
            "\n"
            "Likelihood ratio test:\n"
            "  LRT statistic:      0.7301\n"
            "  Degrees of freedom: 1\n"
            "  Chi-squared p-value: 0.3929\n"
            "  Simulated p-value:  0.2500 (10 simulations)\n"
            "\n"
            "Effect size:\n"
            "  R2_regime: -0.0341"
        )

    @patch("builtins.print")
    def test_json_output(self, mocked_print, json_args):
        svc = RateHeterogeneity(json_args)
        svc.run()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "n_tips" in payload
        assert payload["n_tips"] == 8
        assert "regimes" in payload
        assert "single_rate" in payload
        assert "multi_rate" in payload
        assert "lrt" in payload
        assert "sigma2" in payload["single_rate"]
        assert "sigma2" in payload["multi_rate"]
        assert "chi2_p_value" in payload["lrt"]

    def test_count_regime_tips_scans_assignments_once(self):
        class CountingValues:
            def __init__(self, values):
                self.values = values
                self.iterations = 0

            def __iter__(self):
                self.iterations += 1
                return iter(self.values)

        class Assignments:
            def __init__(self):
                self.iterable = CountingValues(["b", "a", "b", "c", "a"])

            def values(self):
                return self.iterable

        assignments = Assignments()

        counts = RateHeterogeneity._count_regime_tips(
            assignments,
            ["a", "b", "c", "missing"],
        )

        assert counts == {"a": 2, "b": 2, "c": 1, "missing": 0}
        assert assignments.iterable.iterations == 1

    def test_prepare_shared_trait_regime_data_reuses_all_shared_dicts(self):
        tree_tips = ["A", "B", "C"]
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0}
        regime_assignments = {"A": "r1", "B": "r2", "C": "r1"}

        (
            shared_traits,
            shared_regimes,
            tips_to_prune,
            ordered_names,
            regimes,
        ) = RateHeterogeneity._prepare_shared_trait_regime_data(
            tree_tips,
            trait_values,
            regime_assignments,
        )

        assert shared_traits is trait_values
        assert shared_regimes is regime_assignments
        assert tips_to_prune == []
        assert ordered_names == ["A", "B", "C"]
        assert regimes == ["r1", "r2"]

    def test_shared_trait_regime_taxa_scans_smaller_mapping(self):
        class NoIterTraitValues(dict):
            def __iter__(self):
                raise AssertionError("large trait mapping should not be scanned")

        trait_values = NoIterTraitValues(
            {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0, "E": 5.0}
        )
        regime_assignments = {"B": "r1", "D": "r2"}

        assert RateHeterogeneity._shared_trait_regime_taxa(
            trait_values,
            regime_assignments,
        ) == {"B", "D"}

    def test_shared_trait_regime_taxa_scans_traits_when_smaller(self):
        class NoIterRegimeAssignments(dict):
            def __iter__(self):
                raise AssertionError("large regime mapping should not be scanned")

        trait_values = {"B": 2.0, "D": 4.0}
        regime_assignments = NoIterRegimeAssignments(
            {"A": "r1", "B": "r1", "C": "r2", "D": "r2", "E": "r3"}
        )

        assert RateHeterogeneity._shared_trait_regime_taxa(
            trait_values,
            regime_assignments,
        ) == {"B", "D"}

    def test_prepare_shared_trait_regime_data_exact_keys_prunes_tree_only_tips(
        self, monkeypatch
    ):
        class OrderedTraitValues(dict):
            def __contains__(self, key):
                raise AssertionError("ordered prune path should not scan membership")

        tree_tips = ["A", "B", "C", "D"]
        trait_values = OrderedTraitValues({"A": 1.0, "B": 2.0, "C": 3.0})
        regime_assignments = {"A": "r1", "B": "r2", "C": "r1"}

        monkeypatch.setattr(
            rate_heterogeneity_module.Tree, "_ORDERED_MAPPING_PRUNE_MIN_SIZE", 0
        )
        (
            shared_traits,
            shared_regimes,
            tips_to_prune,
            ordered_names,
            regimes,
        ) = RateHeterogeneity._prepare_shared_trait_regime_data(
            tree_tips,
            trait_values,
            regime_assignments,
        )

        assert shared_traits is trait_values
        assert shared_regimes is regime_assignments
        assert tips_to_prune == ["D"]
        assert ordered_names == ["A", "B", "C"]
        assert regimes == ["r1", "r2"]

    def test_prepare_shared_trait_regime_data_filters_partial_overlap(self):
        tree_tips = ["A", "B", "C", "D"]
        trait_values = {"A": 1.0, "B": 2.0, "C": 3.0, "trait_only": 4.0}
        regime_assignments = {
            "A": "r1",
            "B": "r2",
            "C": "r1",
            "regime_only": "r3",
        }

        (
            shared_traits,
            shared_regimes,
            tips_to_prune,
            ordered_names,
            regimes,
        ) = RateHeterogeneity._prepare_shared_trait_regime_data(
            tree_tips,
            trait_values,
            regime_assignments,
        )

        assert shared_traits == {"A": 1.0, "B": 2.0, "C": 3.0}
        assert shared_regimes == {"A": "r1", "B": "r2", "C": "r1"}
        assert tips_to_prune == ["D"]
        assert ordered_names == ["A", "B", "C"]
        assert regimes == ["r1", "r2"]

    def test_prepare_shared_trait_regime_data_partial_overlap_scans_smaller_mapping(self):
        class NoIterTraitValues(dict):
            def __iter__(self):
                raise AssertionError("large trait mapping should not be scanned")

        tree_tips = ["A", "B", "C", "D", "E"]
        trait_values = NoIterTraitValues(
            {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0, "E": 5.0}
        )
        regime_assignments = {"B": "r1", "C": "r1", "D": "r2", "regime_only": "r3"}

        (
            shared_traits,
            shared_regimes,
            tips_to_prune,
            ordered_names,
            regimes,
        ) = RateHeterogeneity._prepare_shared_trait_regime_data(
            tree_tips,
            trait_values,
            regime_assignments,
        )

        assert shared_traits == {"B": 2.0, "C": 3.0, "D": 4.0}
        assert shared_regimes == {"B": "r1", "C": "r1", "D": "r2"}
        assert tips_to_prune == ["A", "E"]
        assert ordered_names == ["B", "C", "D"]
        assert regimes == ["r1", "r2"]

    def test_sum_vcv_matrices_accumulates_without_mutating_inputs(self):
        x = np.array([[1.0, 2.0], [3.0, 4.0]])
        y = np.array([[0.5, 1.5], [2.5, 3.5]])
        z = np.array([[2.0, 0.0], [0.0, 2.0]])
        original_x = x.copy()
        original_y = y.copy()
        original_z = z.copy()

        observed = RateHeterogeneity._sum_vcv_matrices({"x": x, "y": y, "z": z})

        np.testing.assert_allclose(observed, x + y + z)
        np.testing.assert_allclose(x, original_x)
        np.testing.assert_allclose(y, original_y)
        np.testing.assert_allclose(z, original_z)
        assert observed is not x

    @patch("builtins.print")
    def test_two_regime_lrt_p_value_does_not_import_scipy_stats(
        self, mocked_print, json_args, monkeypatch
    ):
        svc = RateHeterogeneity(json_args)
        real_import = __import__

        def fail_scipy_stats_import(name, *args, **kwargs):
            if name == "scipy.stats" or name.startswith("scipy.stats."):
                raise AssertionError("df=1 LRT p-value should not import scipy.stats")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr("builtins.__import__", fail_scipy_stats_import)

        svc.run()

        payload = json.loads(mocked_print.call_args.args[0])
        assert 0 <= payload["lrt"]["chi2_p_value"] <= 1

    def test_three_regime_chi2_p_value_does_not_import_scipy_stats(
        self, monkeypatch
    ):
        real_import = __import__

        def fail_scipy_stats_import(name, *args, **kwargs):
            if name == "scipy.stats" or name.startswith("scipy.stats."):
                raise AssertionError("df=2 LRT p-value should not import scipy.stats")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr("builtins.__import__", fail_scipy_stats_import)

        assert rate_heterogeneity_module._chi2_sf(2.0, df=2) == pytest.approx(
            0.36787944117144233
        )

    def test_run_uses_fast_tip_name_helper_for_prune_setup(self, default_args, mocker):
        svc = RateHeterogeneity(default_args)
        spy = mocker.spy(svc, "get_tip_names_from_tree")
        mocker.patch("builtins.print")

        svc.run()

        assert spy.call_count == 1

    def test_run_uses_unmodified_tree_read(self, default_args, mocker):
        svc = RateHeterogeneity(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tip_names = ["A", "B", "C", "D"]
        traits = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        regimes = {"A": "x", "B": "x", "C": "y", "D": "y"}
        per_regime_vcv = {"x": np.eye(4), "y": np.eye(4)}

        read_unmodified = mocker.patch.object(
            svc, "read_tree_file_unmodified", return_value=tree
        )
        mocker.patch.object(
            svc,
            "read_tree_file",
            side_effect=AssertionError("run should use read_tree_file_unmodified"),
        )
        validate = mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(svc, "get_tip_names_from_tree", return_value=tip_names)
        mocker.patch.object(svc, "_parse_trait_file", return_value=traits)
        mocker.patch.object(svc, "_parse_regime_file", return_value=regimes)
        mocker.patch.object(svc, "_iter_preorder", return_value=[])
        mocker.patch.object(svc, "_build_parent_map", return_value={})
        mocker.patch.object(svc, "_assign_branch_regimes", return_value={})
        build_vcv = mocker.patch.object(
            svc, "_build_per_regime_vcv", return_value=per_regime_vcv
        )
        mocker.patch.object(svc, "_fit_single_rate", return_value=(1.0, 0.0, -10.0))
        mocker.patch.object(
            svc, "_fit_multi_rate", return_value=(np.array([1.0, 2.0]), 0.0, -9.0)
        )
        mocker.patch.object(svc, "_print_text_output")
        fast_copy = mocker.patch.object(
            svc,
            "_fast_copy",
            side_effect=AssertionError("all-shared analysis should not copy tree"),
        )

        svc.run()

        read_unmodified.assert_called_once_with()
        validate.assert_called_once()
        fast_copy.assert_not_called()
        assert build_vcv.call_args.args[0] is tree

    def test_run_copies_before_pruning_missing_tree_tips(self, default_args, mocker):
        svc = RateHeterogeneity(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree_copy = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        pruned_tree = Phylo.read(StringIO("((A:1,B:1):1,C:1);"), "newick")
        tip_names = ["A", "B", "C", "D"]
        traits = {"A": 1.0, "B": 2.0, "C": 3.0}
        regimes = {"A": "x", "B": "x", "C": "y"}

        mocker.patch.object(svc, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(svc, "get_tip_names_from_tree", return_value=tip_names)
        mocker.patch.object(svc, "_parse_trait_file", return_value=traits)
        mocker.patch.object(svc, "_parse_regime_file", return_value=regimes)
        fast_copy = mocker.patch.object(svc, "_fast_copy", return_value=tree_copy)
        prune = mocker.patch.object(
            svc, "prune_tree_using_taxa_list", return_value=pruned_tree
        )
        mocker.patch.object(svc, "_iter_preorder", return_value=[])
        mocker.patch.object(svc, "_build_parent_map", return_value={})
        mocker.patch.object(svc, "_assign_branch_regimes", return_value={})
        build_vcv = mocker.patch.object(
            svc, "_build_per_regime_vcv", return_value={"x": np.eye(3), "y": np.eye(3)}
        )
        mocker.patch.object(svc, "_fit_single_rate", return_value=(1.0, 0.0, -10.0))
        mocker.patch.object(
            svc, "_fit_multi_rate", return_value=(np.array([1.0, 2.0]), 0.0, -9.0)
        )
        mocker.patch.object(svc, "_print_text_output")

        svc.run()

        fast_copy.assert_called_once_with(tree)
        prune.assert_called_once_with(tree_copy, ["D"])
        assert build_vcv.call_args.args[0] is pruned_tree

    def test_run_copies_before_ladderized_plot(self, default_args, mocker, tmp_path):
        default_args.plot = str(tmp_path / "rate.png")
        default_args.ladderize = True
        svc = RateHeterogeneity(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tree_copy = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        tip_names = ["A", "B", "C", "D"]
        traits = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        regimes = {"A": "x", "B": "x", "C": "y", "D": "y"}

        mocker.patch.object(svc, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(svc, "get_tip_names_from_tree", return_value=tip_names)
        mocker.patch.object(svc, "_parse_trait_file", return_value=traits)
        mocker.patch.object(svc, "_parse_regime_file", return_value=regimes)
        fast_copy = mocker.patch.object(svc, "_fast_copy", return_value=tree_copy)
        prune = mocker.patch.object(
            svc,
            "prune_tree_using_taxa_list",
            side_effect=AssertionError("all-shared ladderized plot should not prune"),
        )
        ladderize = mocker.spy(tree_copy, "ladderize")
        mocker.patch.object(svc, "_iter_preorder", return_value=[])
        mocker.patch.object(svc, "_build_parent_map", return_value={})
        mocker.patch.object(svc, "_assign_branch_regimes", return_value={})
        build_vcv = mocker.patch.object(
            svc, "_build_per_regime_vcv", return_value={"x": np.eye(4), "y": np.eye(4)}
        )
        mocker.patch.object(svc, "_fit_single_rate", return_value=(1.0, 0.0, -10.0))
        mocker.patch.object(
            svc, "_fit_multi_rate", return_value=(np.array([1.0, 2.0]), 0.0, -9.0)
        )
        plot = mocker.patch.object(svc, "_plot_regime_tree")
        mocker.patch.object(svc, "_print_text_output")

        svc.run()

        fast_copy.assert_called_once_with(tree)
        prune.assert_not_called()
        ladderize.assert_called_once_with()
        assert build_vcv.call_args.args[0] is tree_copy
        assert plot.call_args.args[0] is tree_copy

    @patch("builtins.print")
    def test_with_simulations(self, mocked_print, sim_args):
        sim_args.json = True
        svc = RateHeterogeneity(sim_args)
        svc.run()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["lrt"]["n_sim"] == 10
        assert payload["lrt"]["sim_p_value"] is not None
        assert 0 <= payload["lrt"]["sim_p_value"] <= 1


class TestEffectSize:
    @patch("builtins.print")
    def test_r2_regime_in_json(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE,
            nsim=0,
            seed=None,
            plot=None,
            json=True,
        )
        svc = RateHeterogeneity(args)
        svc.run()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "r_squared_regime" in payload
        assert isinstance(payload["r_squared_regime"], float)
        assert np.isfinite(payload["r_squared_regime"])

    @patch("builtins.print")
    def test_r2_regime_in_text_output(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE,
            nsim=0,
            seed=None,
            plot=None,
            json=False,
        )
        svc = RateHeterogeneity(args)
        svc.run()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "R2_regime" in all_output

    @patch("builtins.print")
    def test_r2_regime_bounded(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE,
            nsim=0,
            seed=None,
            plot=None,
            json=True,
        )
        svc = RateHeterogeneity(args)
        svc.run()

        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["r_squared_regime"] <= 1.0

    @patch("builtins.print")
    def test_r2_regime_matches_reference(self, mocked_print):
        """R²_regime must match reference values from Fitch-parsimony regime assignment.

        Reference (PhyKIT, tip-count-weighted):
          sigma2_single = 0.0384065703  (matches R brownie.lite)
          sigma2_aquatic = 0.0088113062
          sigma2_terrestrial = 0.0500190412
          weighted avg = (2/8)*0.00881 + (6/8)*0.05002 = 0.0397170574
          R²_regime = 1 - 0.0397170574 / 0.0384065703 = -0.0341227314

        Note: R brownie.lite uses stochastic mapping (make.simmap) for
        regime assignment, so its per-regime σ² values differ from PhyKIT's
        Fitch parsimony. The single-rate σ² matches exactly.
        """
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE,
            nsim=0,
            seed=None,
            plot=None,
            json=True,
        )
        svc = RateHeterogeneity(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        # Single-rate σ² matches R brownie.lite exactly
        assert payload["single_rate"]["sigma2"] == pytest.approx(
            0.0384065703, abs=1e-6
        )
        # R²_regime is negative (multi-rate weighted avg > single rate)
        assert payload["r_squared_regime"] == pytest.approx(
            -0.0341227314, abs=1e-4
        )


class TestPlot:
    @pytest.mark.parametrize("circular", [False, True])
    def test_plot_regime_tree_uses_direct_tree_traversal(
        self, default_args, monkeypatch, tmp_path, circular
    ):
        default_args.circular = circular
        svc = RateHeterogeneity(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clades = list(svc._iter_preorder(tree.root))
        parent_map = svc._build_parent_map(tree, clades)
        regimes = ["aquatic", "terrestrial"]
        branch_regimes = {
            id(clade): regimes[index % 2]
            for index, clade in enumerate(clades)
            if clade != tree.root
        }

        def fail_traversal(*args, **kwargs):
            raise AssertionError("plot setup should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)

        output_path = tmp_path / f"rate_heterogeneity_direct_{circular}.png"
        svc._plot_regime_tree(
            tree, branch_regimes, regimes, parent_map, str(output_path)
        )

        assert output_path.exists()

    def test_plot_regime_tree_reuses_preorder_for_node_positions(
        self, default_args, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import phykit.helpers.plot_config as plot_config

        default_args.circular = False
        svc = RateHeterogeneity(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clades = list(svc._iter_preorder(tree.root))
        parent_map = svc._build_parent_map(tree, clades)
        regimes = ["aquatic", "terrestrial"]
        branch_regimes = {
            id(clade): regimes[index % 2]
            for index, clade in enumerate(clades)
            if clade != tree.root
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

        output_path = tmp_path / "rate_heterogeneity_preorder_positions.png"
        svc._plot_regime_tree(
            tree, branch_regimes, regimes, parent_map, str(output_path)
        )

        assert calls == [expected_preorder_ids]
        assert output_path.exists()

    def test_plot_regime_tree_reuses_clade_lists_for_circular_coords(
        self, default_args, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import phykit.helpers.circular_layout as circular_layout

        default_args.circular = True
        default_args.ylabel_fontsize = 0
        default_args.no_title = True
        svc = RateHeterogeneity(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clades = list(svc._iter_preorder(tree.root))
        parent_map = svc._build_parent_map(tree, clades)
        regimes = ["aquatic", "terrestrial"]
        branch_regimes = {
            id(clade): regimes[index % 2]
            for index, clade in enumerate(clades)
            if clade != tree.root
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

        output_path = tmp_path / "rate_heterogeneity_circular_precomputed_coords.png"
        svc._plot_regime_tree(
            tree, branch_regimes, regimes, parent_map, str(output_path)
        )

        assert calls == [(expected_preorder_ids, expected_terminal_ids)]
        assert output_path.exists()

    def test_rectangular_plot_batches_regime_branches(
        self, default_args, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import matplotlib.axes

        default_args.circular = False
        svc = RateHeterogeneity(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clades = list(svc._iter_preorder(tree.root))
        parent_map = svc._build_parent_map(tree, clades)
        regimes = ["aquatic", "terrestrial"]
        branch_regimes = {
            id(clade): regimes[index % 2]
            for index, clade in enumerate(clades)
            if clade != tree.root
        }

        def fail_plot(*args, **kwargs):
            raise AssertionError("rectangular regime branches should be batched")

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)

        output_path = tmp_path / "rate_heterogeneity_batched.png"
        svc._plot_regime_tree(
            tree, branch_regimes, regimes, parent_map, str(output_path)
        )

        assert output_path.exists()

    def test_rectangular_plot_skips_redundant_tight_layout(
        self, default_args, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        from matplotlib.figure import Figure

        default_args.circular = False
        default_args.legend_position = "none"
        default_args.ylabel_fontsize = 0
        default_args.no_title = True
        svc = RateHeterogeneity(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clades = list(svc._iter_preorder(tree.root))
        parent_map = svc._build_parent_map(tree, clades)
        regimes = ["aquatic", "terrestrial"]
        branch_regimes = {
            id(clade): regimes[index % 2]
            for index, clade in enumerate(clades)
            if clade != tree.root
        }

        def fail_tight_layout(*args, **kwargs):
            raise AssertionError("bbox_inches='tight' handles saved bounds")

        monkeypatch.setattr(Figure, "tight_layout", fail_tight_layout)

        output_path = tmp_path / "rate_heterogeneity_no_tight_layout.png"
        svc._plot_regime_tree(
            tree, branch_regimes, regimes, parent_map, str(output_path)
        )

        assert output_path.exists()
        assert output_path.stat().st_size > 0

    def test_circular_plot_batches_regime_branches(
        self, default_args, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        default_args.circular = True
        default_args.ylabel_fontsize = 0
        default_args.no_title = True
        svc = RateHeterogeneity(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clades = list(svc._iter_preorder(tree.root))
        parent_map = svc._build_parent_map(tree, clades)
        regimes = ["aquatic", "terrestrial"]
        branch_regimes = {
            id(clade): regimes[index % 2]
            for index, clade in enumerate(clades)
            if clade != tree.root
        }
        original_add_collection = matplotlib.axes.Axes.add_collection
        line_collections = []

        def fail_plot(*args, **kwargs):
            raise AssertionError("circular regime branches should be batched")

        def count_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(
            matplotlib.axes.Axes, "add_collection", count_collection
        )

        output_path = tmp_path / "rate_heterogeneity_circular_batched.png"
        svc._plot_regime_tree(
            tree, branch_regimes, regimes, parent_map, str(output_path)
        )

        assert len(line_collections) >= 2
        assert output_path.exists()

    def test_circular_plot_skips_redundant_tight_layout(
        self, default_args, monkeypatch, tmp_path
    ):
        pytest.importorskip("matplotlib")
        from matplotlib.figure import Figure

        default_args.circular = True
        default_args.ylabel_fontsize = 0
        default_args.no_title = True
        svc = RateHeterogeneity(default_args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        clades = list(svc._iter_preorder(tree.root))
        parent_map = svc._build_parent_map(tree, clades)
        regimes = ["aquatic", "terrestrial"]
        branch_regimes = {
            id(clade): regimes[index % 2]
            for index, clade in enumerate(clades)
            if clade != tree.root
        }

        def fail_tight_layout(*args, **kwargs):
            raise AssertionError("bbox_inches='tight' handles saved bounds")

        monkeypatch.setattr(Figure, "tight_layout", fail_tight_layout)

        output_path = tmp_path / "rate_heterogeneity_circular_no_tight_layout.png"
        svc._plot_regime_tree(
            tree, branch_regimes, regimes, parent_map, str(output_path)
        )

        assert output_path.exists()
        assert output_path.stat().st_size > 0

    def test_plot_file_created(self, default_args):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            tmppath = f.name

        try:
            default_args.plot = tmppath
            svc = RateHeterogeneity(default_args)
            svc.run()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)

    def test_plot_circular(self):
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
                regime_data=REGIMES_FILE,
                nsim=0,
                seed=None,
                plot=tmppath,
                json=False,
                circular=True,
            )
            svc = RateHeterogeneity(args)
            svc.run()
            assert os.path.exists(tmppath)
            assert os.path.getsize(tmppath) > 0
        finally:
            if os.path.exists(tmppath):
                os.unlink(tmppath)
