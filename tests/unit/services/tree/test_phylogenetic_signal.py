import builtins
import importlib
import subprocess
import sys
import pytest
import json
import numpy as np
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from phykit.services.tree.phylogenetic_signal import PhylogeneticSignal
from phykit.helpers.trait_parsing import parse_multi_trait_file
import phykit.services.tree.phylogenetic_signal as ps_module
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits.tsv")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


def test_module_import_does_not_import_numpy_pgls_or_scipy():
    code = """
import sys
import phykit.services.tree.phylogenetic_signal as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert callable(module.parse_multi_trait_file)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.trait_parsing" not in sys.modules
assert "numpy" not in sys.modules
assert "phykit.helpers.pgls_utils" not in sys.modules
assert "scipy.linalg" not in sys.modules
assert "scipy.optimize" not in sys.modules
assert module._CHO_FACTOR is None
assert module._CHO_SOLVE is None
assert module._MINIMIZE_SCALAR is None
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_module_import_does_not_import_scipy_linalg_or_optimize(monkeypatch):
    module_name = "phykit.services.tree.phylogenetic_signal"
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
                "phylogenetic_signal module import should not import SciPy linalg/optimize"
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


def test_permutation_p_value_ge_counts_extreme_permutations(monkeypatch):
    permutations = np.array([0.5, 1.5, 2.5, 3.5])
    original_count_nonzero = ps_module.np.count_nonzero
    calls = []

    def counting_count_nonzero(values):
        calls.append(values.copy())
        return original_count_nonzero(values)

    monkeypatch.setattr(ps_module.np, "count_nonzero", counting_count_nonzero)

    p_value = ps_module._permutation_p_value_ge(permutations, 2.0)

    assert p_value == pytest.approx(0.5)
    assert len(calls) == 1
    np.testing.assert_array_equal(calls[0], np.array([False, False, True, True]))
    assert np.isnan(ps_module._permutation_p_value_ge(np.array([]), 2.0))


def test_repeated_cholesky_calls_cache_scipy_linalg_imports(monkeypatch):
    previous_cho_factor = ps_module._CHO_FACTOR
    previous_cho_solve = ps_module._CHO_SOLVE
    ps_module._CHO_FACTOR = None
    ps_module._CHO_SOLVE = None
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
        svc = PhylogeneticSignal.__new__(PhylogeneticSignal)
        n = 12
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        x = rng.normal(size=n)

        svc._log_likelihood(x, C)
        first_call_imports = scipy_linalg_imports
        svc._log_likelihood(x, C)
    finally:
        ps_module._CHO_FACTOR = previous_cho_factor
        ps_module._CHO_SOLVE = previous_cho_solve

    assert first_call_imports > 0
    assert scipy_linalg_imports == first_call_imports


def test_repeated_minimize_scalar_calls_cache_scipy_optimize_imports(monkeypatch):
    previous_minimize_scalar = ps_module._MINIMIZE_SCALAR
    ps_module._MINIMIZE_SCALAR = None
    original_import = builtins.__import__
    scipy_optimize_imports = 0

    def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
        nonlocal scipy_optimize_imports
        if name == "scipy.optimize":
            scipy_optimize_imports += 1
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", counting_import)
    try:
        first = ps_module.minimize_scalar(
            lambda value: (value - 0.25) ** 2,
            bounds=(0.0, 1.0),
            method="bounded",
        )
        first_call_imports = scipy_optimize_imports
        second = ps_module.minimize_scalar(
            lambda value: (value - 0.75) ** 2,
            bounds=(0.0, 1.0),
            method="bounded",
        )
    finally:
        ps_module._MINIMIZE_SCALAR = previous_minimize_scalar

    assert first.x == pytest.approx(0.25, abs=1e-5)
    assert second.x == pytest.approx(0.75, abs=1e-5)
    assert first_call_imports > 0
    assert scipy_optimize_imports == first_call_imports


def test_inverse_from_spd_matrix_cholesky_matches_inverse():
    rng = np.random.default_rng(20260628)
    A = rng.normal(size=(8, 8))
    matrix = A @ A.T + np.eye(8) * 0.25

    observed = ps_module._inverse_from_spd_matrix(matrix)

    np.testing.assert_allclose(observed, np.linalg.inv(matrix))


def test_inverse_from_spd_matrix_avoids_inverse_for_spd_matrix(monkeypatch):
    rng = np.random.default_rng(20260628)
    A = rng.normal(size=(6, 6))
    matrix = A @ A.T + np.eye(6)

    def fail_inverse(_matrix):
        raise AssertionError("SPD covariance inversion should use Cholesky")

    monkeypatch.setattr(ps_module.np.linalg, "inv", fail_inverse)

    observed = ps_module._inverse_from_spd_matrix(matrix)

    assert observed.shape == matrix.shape
    assert np.all(np.isfinite(observed))


def test_inverse_from_spd_matrix_keeps_inverse_fallback():
    matrix = np.array([[1.0, 2.0], [2.0, 1.0]])

    observed = ps_module._inverse_from_spd_matrix(matrix)

    np.testing.assert_allclose(observed, np.linalg.inv(matrix))


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        method="blombergs_k",
        permutations=1000,
        json=False,
    )


@pytest.fixture
def lambda_args():
    return Namespace(
        tree=TREE_SIMPLE,
        trait_data=TRAITS_FILE,
        method="lambda",
        permutations=1000,
        json=False,
    )


class TestProcessArgs:
    def test_defaults(self):
        args = Namespace(tree="t.tre", trait_data="d.tsv")
        svc = PhylogeneticSignal(args)
        assert svc.tree_file_path == "t.tre"
        assert svc.trait_data_path == "d.tsv"
        assert svc.method == "blombergs_k"
        assert svc.permutations == 1000
        assert svc.json_output is False

    def test_json_true(self):
        args = Namespace(tree="t.tre", trait_data="d.tsv", json=True, method="lambda", permutations=500)
        svc = PhylogeneticSignal(args)
        assert svc.json_output is True
        assert svc.method == "lambda"
        assert svc.permutations == 500


class TestParseTraitFile:
    def test_valid_file(self, default_args):
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        traits = svc._parse_trait_file(TRAITS_FILE, tips)
        assert len(traits) == 8
        assert traits["raccoon"] == pytest.approx(1.04)
        assert traits["weasel"] == pytest.approx(-0.30)

    def test_missing_file(self, default_args):
        svc = PhylogeneticSignal(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_trait_file("/nonexistent/path.tsv", ["a", "b", "c"])

    def test_non_numeric_value(self, default_args, tmp_path):
        trait_file = tmp_path / "bad.tsv"
        trait_file.write_text("taxon1\tabc\n")
        svc = PhylogeneticSignal(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_trait_file(str(trait_file), ["taxon1"])

    def test_wrong_columns(self, default_args, tmp_path):
        trait_file = tmp_path / "bad.tsv"
        trait_file.write_text("taxon1\t1.0\textra\n")
        svc = PhylogeneticSignal(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_trait_file(str(trait_file), ["taxon1"])

    def test_comments_and_blanks(self, default_args, tmp_path):
        trait_file = tmp_path / "good.tsv"
        trait_file.write_text("# comment\n\ntaxon1\t1.0\ntaxon2\t2.0\ntaxon3\t3.0\n")
        svc = PhylogeneticSignal(default_args)
        traits = svc._parse_trait_file(str(trait_file), ["taxon1", "taxon2", "taxon3"])
        assert len(traits) == 3

    def test_all_shared_trait_file_emits_no_warnings(
        self, default_args, tmp_path, capsys
    ):
        trait_file = tmp_path / "good.tsv"
        trait_file.write_text("taxon1\t1.0\ntaxon2\t2.0\ntaxon3\t3.0\n")
        svc = PhylogeneticSignal(default_args)

        traits = svc._parse_trait_file(
            str(trait_file), ["taxon1", "taxon2", "taxon3"]
        )

        stderr = capsys.readouterr().err
        assert traits == {"taxon1": 1.0, "taxon2": 2.0, "taxon3": 3.0}
        assert stderr == ""

    def test_taxon_mismatch_warns(self, default_args, tmp_path, capsys):
        trait_file = tmp_path / "partial.tsv"
        trait_file.write_text("taxon1\t1.0\ntaxon2\t2.0\ntaxon3\t3.0\nextra\t4.0\n")
        svc = PhylogeneticSignal(default_args)
        traits = svc._parse_trait_file(
            str(trait_file), ["taxon1", "taxon2", "taxon3", "taxon4"]
        )
        _, err = capsys.readouterr()
        assert "Warning" in err
        assert len(traits) == 3

    def test_too_few_shared_taxa(self, default_args, tmp_path):
        trait_file = tmp_path / "few.tsv"
        trait_file.write_text("taxon1\t1.0\ntaxon2\t2.0\n")
        svc = PhylogeneticSignal(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_trait_file(str(trait_file), ["taxon1", "taxon2"])


class TestValidateTree:
    def test_too_few_tips(self, default_args, mocker):
        svc = PhylogeneticSignal(default_args)
        from Bio import Phylo
        from io import StringIO
        tree = Phylo.read(StringIO("(A:1.0,B:2.0);"), "newick")
        with pytest.raises(PhykitUserError):
            svc.validate_tree(tree, min_tips=3, require_branch_lengths=True)

    def test_missing_branch_lengths(self, default_args):
        svc = PhylogeneticSignal(default_args)
        from Bio import Phylo
        from io import StringIO
        tree = Phylo.read(StringIO("(A,B,C);"), "newick")
        with pytest.raises(PhykitUserError):
            svc.validate_tree(tree, min_tips=3, require_branch_lengths=True)

    def test_valid_tree_passes(self, default_args):
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        svc.validate_tree(tree, min_tips=3, require_branch_lengths=True)


class TestBuildVCVMatrix:
    def test_symmetric(self, default_args):
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        vcv = svc._build_vcv_matrix(tree, tips)
        np.testing.assert_array_almost_equal(vcv, vcv.T)

    def test_diagonal_is_root_to_tip(self, default_args):
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        vcv = svc._build_vcv_matrix(tree, tips)
        for i, name in enumerate(tips):
            expected = tree.distance(tree.root, name)
            assert vcv[i, i] == pytest.approx(expected, rel=1e-6)

    def test_correct_shape(self, default_args):
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        vcv = svc._build_vcv_matrix(tree, tips)
        assert vcv.shape == (8, 8)


class TestBlombergsK:
    def test_returns_expected_keys(self, default_args):
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        traits = svc._parse_trait_file(TRAITS_FILE, tips)
        ordered = sorted(traits.keys())
        x = np.array([traits[n] for n in ordered])
        vcv = svc._build_vcv_matrix(tree, ordered)
        result = svc._blombergs_k(x, vcv, 100)
        assert "K" in result
        assert "p_value" in result
        assert "permutations" in result
        assert result["permutations"] == 100

    def test_k_matches_phytools(self, default_args):
        """K must match R phytools::phylosig(method='K') reference value."""
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        traits = svc._parse_trait_file(TRAITS_FILE, tips)
        ordered = sorted(traits.keys())
        x = np.array([traits[n] for n in ordered])
        vcv = svc._build_vcv_matrix(tree, ordered)
        result = svc._blombergs_k(x, vcv, 1000)
        # R reference: phylosig(tree, x, method="K") = 0.584216
        assert result["K"] == pytest.approx(0.584216, abs=1e-4)
        assert 0 <= result["p_value"] <= 1

    def test_batched_permutations_match_scalar_loop(self):
        rng = np.random.default_rng(20260622)
        svc = PhylogeneticSignal.__new__(PhylogeneticSignal)
        n_perm = 12
        n = 8
        A = rng.normal(size=(n, n))
        vcv = A @ A.T + np.eye(n)
        C_inv = np.linalg.inv(vcv)
        x = rng.normal(size=n)
        ones = np.ones(n)
        sum_C_inv = float(ones @ C_inv @ ones)
        expected_ratio = (np.trace(vcv) - n / sum_C_inv) / (n - 1)

        scalar = np.empty(n_perm)
        perm_rng = np.random.default_rng(seed=42)
        for perm_i in range(n_perm):
            x_perm = perm_rng.permutation(x)
            a_p = float((ones @ C_inv @ x_perm) / sum_C_inv)
            e_p = x_perm - a_p
            obs_p = float(e_p @ e_p) / float(e_p @ C_inv @ e_p)
            scalar[perm_i] = obs_p / expected_ratio

        batched = svc._blombergs_k_permutations(
            x,
            C_inv,
            sum_C_inv,
            expected_ratio,
            n_perm,
            batch_size=5,
        )

        assert batched == pytest.approx(scalar)


class TestPagelsLambda:
    def test_returns_expected_keys(self, default_args):
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        traits = svc._parse_trait_file(TRAITS_FILE, tips)
        ordered = sorted(traits.keys())
        x = np.array([traits[n] for n in ordered])
        vcv = svc._build_vcv_matrix(tree, ordered)
        result = svc._pagels_lambda(x, vcv)
        assert "lambda" in result
        assert "log_likelihood" in result
        assert "p_value" in result

    def test_lambda_matches_phytools(self, default_args):
        """Lambda and logL must match R phytools::phylosig(method='lambda') reference."""
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        traits = svc._parse_trait_file(TRAITS_FILE, tips)
        ordered = sorted(traits.keys())
        x = np.array([traits[n] for n in ordered])
        vcv = svc._build_vcv_matrix(tree, ordered)
        result = svc._pagels_lambda(x, vcv)
        # R reference: phylosig(tree, x, method="lambda", test=TRUE)
        #   lambda = 0.999927, logL = -11.5697, P = 0.716569
        assert result["lambda"] == pytest.approx(0.9999, abs=1e-3)
        assert result["log_likelihood"] == pytest.approx(-11.5697, abs=1e-3)
        assert result["p_value"] == pytest.approx(0.7166, abs=1e-2)
        assert 0 <= result["lambda"] <= 1

    def test_lambda_p_value_does_not_import_scipy_stats(
        self, default_args, monkeypatch
    ):
        svc = PhylogeneticSignal(default_args)
        tree = svc.read_tree_file()
        tips = sorted(svc.get_tip_names_from_tree(tree))
        traits = svc._parse_trait_file(TRAITS_FILE, tips)
        ordered = sorted(traits.keys())
        x = np.array([traits[n] for n in ordered])
        vcv = svc._build_vcv_matrix(tree, ordered)
        real_import = __import__

        def fail_scipy_stats_import(name, *args, **kwargs):
            if name == "scipy.stats" or name.startswith("scipy.stats."):
                raise AssertionError("lambda p-value should not import scipy.stats")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr("builtins.__import__", fail_scipy_stats_import)

        result = svc._pagels_lambda(x, vcv)

        assert 0 <= result["p_value"] <= 1

    def test_log_likelihood_cholesky_matches_inverse(self):
        rng = np.random.default_rng(20260623)
        svc = PhylogeneticSignal.__new__(PhylogeneticSignal)
        n = 12
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        x = rng.normal(size=n)

        fast_ll, fast_a = svc._log_likelihood(x, C)
        inverse_ll, inverse_a = svc._log_likelihood_inverse(x, C)

        assert fast_ll == pytest.approx(inverse_ll)
        assert fast_a == pytest.approx(inverse_a)

    def test_log_likelihood_uses_single_cholesky_solve(self, monkeypatch):
        rng = np.random.default_rng(20260628)
        svc = PhylogeneticSignal.__new__(PhylogeneticSignal)
        n = 12
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        x = rng.normal(size=n)
        original_cho_solve = ps_module.cho_solve
        solve_calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal solve_calls
            solve_calls += 1
            return original_cho_solve(*args, **kwargs)

        monkeypatch.setattr(ps_module, "cho_solve", counting_cho_solve)

        fast_ll, fast_a = svc._log_likelihood(x, C)
        inverse_ll, inverse_a = svc._log_likelihood_inverse(x, C)

        assert solve_calls == 1
        assert fast_ll == pytest.approx(inverse_ll)
        assert fast_a == pytest.approx(inverse_a)

    def test_log_likelihood_singular_matrix_keeps_inverse_error(self):
        svc = PhylogeneticSignal.__new__(PhylogeneticSignal)
        C = np.ones((4, 4))
        x = np.arange(4, dtype=float)

        with pytest.raises(np.linalg.LinAlgError):
            svc._log_likelihood(x, C)


class TestRun:
    def test_run_uses_unmodified_tree_read(self, mocker):
        args = Namespace(
            tree="dummy.tre",
            trait_data="traits.tsv",
            method="blombergs_k",
            permutations=100,
            json=False,
            gene_trees=None,
            multivariate=False,
        )
        tree = object()
        svc = PhylogeneticSignal(args)
        read_tree = mocker.patch.object(
            svc, "read_tree_file_unmodified", return_value=tree
        )
        mocker.patch.object(
            svc,
            "read_tree_file",
            side_effect=AssertionError("copying tree reader should not be used"),
        )
        mocker.patch.object(svc, "validate_tree")
        mocker.patch.object(svc, "get_tip_names_from_tree", return_value=["a", "b", "c"])
        mocker.patch.object(
            svc,
            "_parse_trait_file",
            return_value={"a": 1.0, "b": 2.0, "c": 3.0},
        )
        mocker.patch(
            "phykit.services.tree.vcv_utils.build_vcv_matrix",
            return_value=np.eye(3),
        )
        mocker.patch.object(svc, "_compute_r2_phylo", return_value=0.5)
        mocker.patch.object(
            svc,
            "_blombergs_k",
            return_value={"K": 1.0, "p_value": 0.25, "permutations": 100},
        )

        svc.run()

        read_tree.assert_called_once_with()
        svc.validate_tree.assert_called_once_with(
            tree,
            min_tips=3,
            require_branch_lengths=True,
            context="phylogenetic signal analysis",
        )

    def test_blombergs_k_text_output(self, default_args, capsys):
        svc = PhylogeneticSignal(default_args)
        svc.run()
        out, _ = capsys.readouterr()
        parts = out.strip().split("\t")
        assert len(parts) == 3
        float(parts[0])  # K
        float(parts[1])  # p_value
        float(parts[2])  # r_squared_phylo

    def test_lambda_text_output(self, lambda_args, capsys):
        svc = PhylogeneticSignal(lambda_args)
        svc.run()
        out, _ = capsys.readouterr()
        parts = out.strip().split("\t")
        assert len(parts) == 4
        float(parts[0])  # lambda
        float(parts[1])  # log_likelihood
        float(parts[2])  # p_value
        float(parts[3])  # r_squared_phylo

    def test_json_output_blombergs_k(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "K" in data
        assert "p_value" in data
        assert "permutations" in data

    def test_json_output_lambda(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="lambda",
            permutations=1000,
            json=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "lambda" in data
        assert "log_likelihood" in data
        assert "p_value" in data


GENE_TREES_FILE = str(SAMPLE_FILES / "gene_trees_simple.nwk")


class TestDiscordanceVCV:
    def test_run_with_gene_trees(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "K" in data
        assert "p_value" in data
        assert "vcv_metadata" in data
        assert data["vcv_metadata"]["n_gene_trees"] == 10

    def test_all_concordant_gene_trees(self, capsys):
        """When all gene trees are identical to species tree, K should be
        very close to species-tree-only K."""
        # Run without gene trees
        args_no_gt = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
        )
        svc = PhylogeneticSignal(args_no_gt)
        svc.run()
        out, _ = capsys.readouterr()
        k_no_gt = json.loads(out)["K"]

        # Run with species tree as gene trees (all concordant)
        import tempfile, shutil
        from Bio import Phylo
        tree = Phylo.read(TREE_SIMPLE, "newick")
        with tempfile.NamedTemporaryFile(mode="w", suffix=".nwk", delete=False) as f:
            for _ in range(3):
                Phylo.write(tree, f, "newick")
            concordant_path = f.name

        try:
            args_gt = Namespace(
                tree=TREE_SIMPLE,
                trait_data=TRAITS_FILE,
                method="blombergs_k",
                permutations=100,
                json=True,
                gene_trees=concordant_path,
            )
            svc = PhylogeneticSignal(args_gt)
            svc.run()
            out, _ = capsys.readouterr()
            k_gt = json.loads(out)["K"]
            assert k_gt == pytest.approx(k_no_gt, rel=0.01)
        finally:
            import os
            os.unlink(concordant_path)

    def test_discordance_vcv_differs_from_species(self, capsys):
        """With discordant gene trees, K should differ from species-tree-only K."""
        args_no_gt = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
        )
        svc = PhylogeneticSignal(args_no_gt)
        svc.run()
        out, _ = capsys.readouterr()
        k_no_gt = json.loads(out)["K"]

        args_gt = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticSignal(args_gt)
        svc.run()
        out, _ = capsys.readouterr()
        k_gt = json.loads(out)["K"]
        # They should differ (gene trees 7-9 are discordant)
        assert k_gt != pytest.approx(k_no_gt, abs=1e-6)

    def test_lambda_with_gene_trees(self, capsys):
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="lambda",
            permutations=1000,
            json=True,
            gene_trees=GENE_TREES_FILE,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "lambda" in data
        assert "vcv_metadata" in data


class TestEffectSize:
    def test_r2_phylo_cholesky_matches_inverse(self):
        rng = np.random.default_rng(20260623)
        svc = PhylogeneticSignal.__new__(PhylogeneticSignal)
        n = 12
        A = rng.normal(size=(n, n))
        vcv = A @ A.T + np.eye(n)
        x = rng.normal(size=n)

        assert svc._compute_r2_phylo(x, vcv) == pytest.approx(
            svc._compute_r2_phylo_inverse(x, vcv)
        )

    def test_r2_phylo_uses_single_cholesky_solve(self, monkeypatch):
        rng = np.random.default_rng(20260628)
        svc = PhylogeneticSignal.__new__(PhylogeneticSignal)
        n = 12
        A = rng.normal(size=(n, n))
        vcv = A @ A.T + np.eye(n)
        x = rng.normal(size=n)
        original_cho_solve = ps_module.cho_solve
        solve_calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal solve_calls
            solve_calls += 1
            return original_cho_solve(*args, **kwargs)

        monkeypatch.setattr(ps_module, "cho_solve", counting_cho_solve)

        assert svc._compute_r2_phylo(x, vcv) == pytest.approx(
            svc._compute_r2_phylo_inverse(x, vcv)
        )
        assert solve_calls == 1

    def test_r2_phylo_constant_trait_returns_nan(self):
        svc = PhylogeneticSignal.__new__(PhylogeneticSignal)
        vcv = np.eye(4)
        x = np.ones(4)

        assert np.isnan(svc._compute_r2_phylo(x, vcv))

    def test_r2_phylo_blombergs_k(self, capsys):
        """JSON output for blombergs_k contains r_squared_phylo as a finite float."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "r_squared_phylo" in data
        assert isinstance(data["r_squared_phylo"], float)
        assert np.isfinite(data["r_squared_phylo"])

    def test_r2_phylo_lambda(self, capsys):
        """JSON output for lambda contains r_squared_phylo as a finite float."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="lambda",
            permutations=1000,
            json=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "r_squared_phylo" in data
        assert isinstance(data["r_squared_phylo"], float)
        assert np.isfinite(data["r_squared_phylo"])

    def test_r2_phylo_positive_for_phylogenetic_trait(self, capsys):
        """Body mass trait has strong phylogenetic signal, so R2_phylo > 0."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert data["r_squared_phylo"] > 0

    def test_r2_phylo_in_text_output(self, capsys):
        """Text output includes R2_phylo value."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=False,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        parts = out.strip().split("\t")
        # Third column is r_squared_phylo
        assert len(parts) == 3
        r2 = float(parts[2])
        assert np.isfinite(r2)

    def test_r2_phylo_matches_r_reference(self, capsys):
        """R²_phylo must match R phytools/manual GLS reference value.

        R reference (tests/r_validation/validate_signal_r2.R):
          sigma2_bm_manual  = 0.0384065703
          sigma2_wn_manual  = 0.7667234375
          r2_phylo_manual   = 0.9499081827
        """
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert data["r_squared_phylo"] == pytest.approx(0.9499081827, abs=1e-6)

    def test_r2_phylo_negative_for_anti_phylogenetic_trait(self, tmp_path, capsys):
        """R²_phylo < 0 when closely related species have very different values.

        The design doc states: 'R²_phylo < 0: possible (phylogeny actively
        misleads), valid. Report as-is.' This test verifies that negative R²
        values are correctly computed and not clipped to 0 or NaN.

        We use a balanced tree with long internal branches (strong expected
        phylogenetic signal) but assign opposite trait values to sister taxa,
        maximally violating the BM assumption.
        """
        # Tree with strong phylogenetic structure (long shared paths)
        tree_file = tmp_path / "balanced.tre"
        tree_file.write_text("((A:0.1,B:0.1):10,(C:0.1,D:0.1):10);")

        # Anti-phylogenetic: sisters have opposite values
        trait_file = tmp_path / "anti_phylo_traits.tsv"
        trait_file.write_text("A\t10.0\nB\t-10.0\nC\t5.0\nD\t-5.0\n")

        args = Namespace(
            tree=str(tree_file),
            trait_data=str(trait_file),
            method="blombergs_k",
            permutations=100,
            json=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        # R²_phylo should be strongly negative: phylogeny misleads
        assert data["r_squared_phylo"] < 0, (
            f"Expected negative R²_phylo for anti-phylogenetic trait, "
            f"got {data['r_squared_phylo']}"
        )
        # Should be a real number, not NaN
        assert np.isfinite(data["r_squared_phylo"])
        # Reference: R²_phylo = -9.0 for this exact setup
        assert data["r_squared_phylo"] == pytest.approx(-9.0, abs=0.1)


class TestKmult:
    def test_kmult_returns_positive(self):
        """K_mult should be > 0 for real trait data."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            multivariate=True,
        )
        svc = PhylogeneticSignal(args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits_multi = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered = sorted(traits_multi.keys())
        from phykit.services.tree.vcv_utils import build_vcv_matrix
        vcv = build_vcv_matrix(tree, ordered)
        p = len(trait_names)
        Y = np.array([[traits_multi[name][j] for j in range(p)] for name in ordered])
        result = svc._kmult(Y, vcv, 100)
        assert result["K_mult"] > 0

    def test_kmult_permutation_pvalue(self):
        """p-value from K_mult permutation test must be between 0 and 1."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            multivariate=True,
        )
        svc = PhylogeneticSignal(args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits_multi = parse_multi_trait_file(MULTI_TRAITS_FILE, tips)
        ordered = sorted(traits_multi.keys())
        from phykit.services.tree.vcv_utils import build_vcv_matrix
        vcv = build_vcv_matrix(tree, ordered)
        p = len(trait_names)
        Y = np.array([[traits_multi[name][j] for j in range(p)] for name in ordered])
        result = svc._kmult(Y, vcv, 100)
        assert 0 <= result["p_value"] <= 1

    def test_kmult_batched_permutations_match_scalar_loop(self):
        rng = np.random.default_rng(20260623)
        n_perm = 11
        n = 7
        p = 3
        A = rng.normal(size=(n, n))
        vcv = A @ A.T + np.eye(n)
        C_inv = np.linalg.inv(vcv)
        Y = rng.normal(size=(n, p))
        ones = np.ones(n)
        sum_C_inv = float(ones @ C_inv @ ones)
        expected_ratio = (np.trace(vcv) - n / sum_C_inv) / (n - 1)

        scalar = np.empty(n_perm)
        perm_rng = np.random.default_rng(seed=42)
        for perm_i in range(n_perm):
            Y_perm = Y[perm_rng.permutation(n)]
            a_p = (ones @ C_inv @ Y_perm) / sum_C_inv
            E_p = Y_perm - a_p
            mse_obs_p = np.trace(E_p.T @ E_p)
            mse_phylo_p = np.trace(E_p.T @ C_inv @ E_p)
            obs_ratio_p = 0.0 if mse_phylo_p == 0 else mse_obs_p / mse_phylo_p
            scalar[perm_i] = obs_ratio_p / expected_ratio

        batched = PhylogeneticSignal._kmult_permutations(
            Y,
            C_inv,
            sum_C_inv,
            expected_ratio,
            n_perm,
            batch_size=4,
        )

        assert batched == pytest.approx(scalar)

    def test_kmult_single_trait_matches_k(self):
        """K_mult with 1 trait should approximately equal univariate K."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            multivariate=False,
        )
        svc = PhylogeneticSignal(args)
        tree = svc.read_tree_file()
        tips = svc.get_tip_names_from_tree(tree)
        traits = svc._parse_trait_file(TRAITS_FILE, tips)
        ordered = sorted(traits.keys())
        x = np.array([traits[n] for n in ordered])
        from phykit.services.tree.vcv_utils import build_vcv_matrix
        vcv = build_vcv_matrix(tree, ordered)

        # Univariate K
        result_k = svc._blombergs_k(x, vcv, 100)

        # K_mult with single trait (n x 1 matrix)
        Y = x.reshape(-1, 1)
        result_kmult = svc._kmult(Y, vcv, 100)

        assert result_kmult["K_mult"] == pytest.approx(result_k["K"], rel=1e-6)

    def test_multivariate_flag_uses_multi_trait_file(self, capsys):
        """--multivariate flag reads multi-column TSV and runs K_mult."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            multivariate=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "K_mult" in data
        assert "p_value" in data
        assert "n_traits" in data
        assert data["n_traits"] == 3
        assert "permutations" in data
        assert data["permutations"] == 100

    def test_multivariate_text_output(self, capsys):
        """Text output for K_mult should have 4 tab-separated fields."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=False,
            multivariate=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        parts = out.strip().split("\t")
        assert len(parts) == 4
        float(parts[0])  # K_mult
        float(parts[1])  # p_value
        int(parts[2])    # n_traits
        int(parts[3])    # permutations

    def test_multivariate_with_lambda_raises_error(self):
        """Using --multivariate with --method lambda should raise an error."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="lambda",
            permutations=100,
            json=False,
            multivariate=True,
        )
        svc = PhylogeneticSignal(args)
        with pytest.raises(PhykitUserError):
            svc.run()

    def test_multivariate_json_output(self, capsys):
        """JSON output for K_mult should contain all expected fields."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            multivariate=True,
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "K_mult" in data
        assert isinstance(data["K_mult"], float)
        assert "p_value" in data
        assert isinstance(data["p_value"], float)
        assert "n_traits" in data
        assert data["n_traits"] == 3
        assert "permutations" in data
        assert data["permutations"] == 100

    def test_multivariate_with_gene_trees(self, capsys):
        """K_mult should work with discordance-aware VCV (--gene-trees)."""
        args = Namespace(
            tree=TREE_SIMPLE,
            trait_data=MULTI_TRAITS_FILE,
            method="blombergs_k",
            permutations=100,
            json=True,
            multivariate=True,
            gene_trees=str(SAMPLE_FILES / "gene_trees_simple.nwk"),
        )
        svc = PhylogeneticSignal(args)
        svc.run()
        out, _ = capsys.readouterr()
        data = json.loads(out)
        assert "K_mult" in data
        assert "vcv_metadata" in data
