import builtins
import importlib
import subprocess
import sys
from argparse import Namespace
from io import StringIO
import json

import pytest
import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

from phykit.errors import PhykitUserError
from phykit.services.tree.phylo_anova import PhyloAnova
import phykit.services.tree.phylo_anova as phylo_anova_module


TREE = "tests/sample_files/tree_simple.tre"
TRAITS_SINGLE = "tests/sample_files/tree_simple_anova_single.tsv"
TRAITS_MULTI = "tests/sample_files/tree_simple_anova_traits.tsv"


def test_module_import_does_not_import_scipy_linalg(monkeypatch):
    module_name = "phykit.services.tree.phylo_anova"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.linalg" or name.startswith("scipy.linalg."):
            raise AssertionError(
                "phylo_anova module import should not import SciPy linalg"
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


def test_module_import_does_not_import_numpy_or_scipy_linalg():
    code = """
import sys
import phykit.services.tree.phylo_anova as module
assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy.linalg" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_permutation_p_value_and_z_matches_legacy_reductions():
    permutations = np.array([0.5, 1.5, 2.5, 3.5])
    observed = 2.0

    p_value, z_score = phylo_anova_module._permutation_p_value_and_z(
        observed, permutations
    )

    assert p_value == pytest.approx(float(np.mean(permutations >= observed)))
    assert z_score == pytest.approx(
        (observed - np.mean(permutations)) / np.std(permutations)
    )


def test_permutation_p_value_and_z_computes_std_once(monkeypatch):
    permutations = np.array([1.0, 2.0, 3.0, 4.0])
    original_std = phylo_anova_module.np.std
    std_calls = 0

    def counting_std(*args, **kwargs):
        nonlocal std_calls
        std_calls += 1
        return original_std(*args, **kwargs)

    monkeypatch.setattr(phylo_anova_module.np, "std", counting_std)

    p_value, z_score = phylo_anova_module._permutation_p_value_and_z(
        2.5, permutations
    )

    assert std_calls == 1
    assert p_value == pytest.approx(0.5)
    assert np.isfinite(z_score)


def test_permutation_p_value_and_z_counts_extreme_permutations(monkeypatch):
    permutations = np.array([0.5, 1.5, 2.5, 3.5])
    original_count_nonzero = phylo_anova_module.np.count_nonzero
    calls = []

    def counting_count_nonzero(values):
        calls.append(values.copy())
        return original_count_nonzero(values)

    monkeypatch.setattr(
        phylo_anova_module.np,
        "count_nonzero",
        counting_count_nonzero,
    )

    p_value, z_score = phylo_anova_module._permutation_p_value_and_z(
        2.0,
        permutations,
    )

    assert p_value == pytest.approx(0.5)
    assert np.isfinite(z_score)
    assert len(calls) == 1
    np.testing.assert_array_equal(calls[0], np.array([False, False, True, True]))


def test_permutation_p_value_and_z_zero_variance_returns_zero_z():
    p_value, z_score = phylo_anova_module._permutation_p_value_and_z(
        1.0, np.ones(4)
    )

    assert p_value == pytest.approx(1.0)
    assert z_score == 0.0


def _make_args(**overrides):
    defaults = dict(
        tree=TREE,
        traits=TRAITS_SINGLE,
        group_column=None,
        method="auto",
        permutations=99,
        pairwise=False,
        plot_output=None,
        plot_type="auto",
        seed=42,
        json=False,
        fig_width=None,
        fig_height=None,
        dpi=150,
        no_title=False,
        title=None,
        legend_position=None,
        ylabel_fontsize=None,
        xlabel_fontsize=None,
        title_fontsize=None,
        axis_fontsize=None,
        colors=None,
        ladderize=False,
        cladogram=False,
        circular=False,
        color_file=None,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


class TestPhyloAnovaParsing:
    def test_parse_trait_file_streams_without_readlines(self, monkeypatch):
        class StreamingOnlyFile(StringIO):
            def __enter__(self):
                return self

            def __exit__(self, *_args):
                self.close()

            def readlines(self):
                raise AssertionError("PhyloANOVA parser should stream rows")

        data = (
            "# ignored before header\n"
            "\n"
            "taxon\tgroup\tbody_mass\tlength\n"
            "A\tg1\t1.0\t10.0\n"
            "  # ignored between rows\n"
            "B\tg2\t2.0\t20.0\n"
            "C\tg1\t3.0\t30.0\n"
        )

        def fake_open(path):
            assert path == "traits.tsv"
            return StreamingOnlyFile(data)

        monkeypatch.setattr("builtins.open", fake_open)
        svc = PhyloAnova.__new__(PhyloAnova)
        svc.group_column = None

        header, traits = svc._parse_trait_file("traits.tsv", ["A", "B", "C"])

        assert header == ["group", "body_mass", "length"]
        assert traits == {
            "A": ["g1", 1.0, 10.0],
            "B": ["g2", 2.0, 20.0],
            "C": ["g1", 3.0, 30.0],
        }

    def test_parse_trait_file_preserves_logical_line_numbers(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "\n".join(
                [
                    "# ignored before header",
                    "taxon\tgroup\tbody_mass",
                    "A\tg1\t1.0",
                    "  # ignored between rows",
                    "B\tg2\tbad",
                ]
            )
            + "\n"
        )
        svc = PhyloAnova.__new__(PhyloAnova)
        svc.group_column = None

        with pytest.raises(PhykitUserError) as excinfo:
            svc._parse_trait_file(str(trait_file), ["A", "B", "C"])

        assert (
            "Non-numeric trait value 'bad' for taxon 'B' "
            "(trait 'body_mass') on line 3."
        ) in excinfo.value.messages

    def test_parse_trait_file_handles_nondefault_group_column(
        self, tmp_path, capsys
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "taxon\tbody_mass\tgroup\tlength\n"
            "A\t1.0\tg1\t10.0\n"
            "B\t2.0\tg2\t20.0\n"
            "C\t3.0\tg1\t30.0\n"
        )
        svc = PhyloAnova.__new__(PhyloAnova)
        svc.group_column = "group"

        header, traits = svc._parse_trait_file(str(trait_file), ["A", "B", "C"])

        assert header == ["body_mass", "group", "length"]
        assert traits == {
            "A": [1.0, "g1", 10.0],
            "B": [2.0, "g2", 20.0],
            "C": [3.0, "g1", 30.0],
        }
        assert capsys.readouterr().err == ""


class TestPhyloAnova:
    def test_groups_and_response_matrix_preserves_nondefault_group_column(self):
        traits = {
            "A": [1.0, "g1", 10.0],
            "B": [2.0, "g2", 20.0],
            "C": [3.0, "g1", 30.0],
        }

        groups, Y = PhyloAnova._groups_and_response_matrix(
            traits,
            ["A", "B", "C"],
            group_idx=1,
            n_cols=3,
        )

        assert groups == ["g1", "g2", "g1"]
        np.testing.assert_array_equal(
            Y,
            np.array(
                [
                    [1.0, 10.0],
                    [2.0, 20.0],
                    [3.0, 30.0],
                ]
            ),
        )

    def test_groups_and_response_matrix_keeps_single_response_2d(self):
        traits = {
            "A": ["g1", 1.0],
            "B": ["g2", 2.0],
        }

        groups, Y = PhyloAnova._groups_and_response_matrix(
            traits,
            ["A", "B"],
            group_idx=0,
            n_cols=2,
        )

        assert groups == ["g1", "g2"]
        assert Y.shape == (2, 1)
        np.testing.assert_array_equal(Y, np.array([[1.0], [2.0]]))

    def test_groups_and_response_matrix_precomputes_response_getter(
        self, monkeypatch
    ):
        calls = []
        original_itemgetter = phylo_anova_module.itemgetter

        def tracking_itemgetter(*items):
            calls.append(items)
            return original_itemgetter(*items)

        monkeypatch.setattr(phylo_anova_module, "itemgetter", tracking_itemgetter)
        traits = {
            f"t{idx}": [float(idx), f"g{idx % 2}", float(idx + 1)]
            for idx in range(10)
        }

        groups, Y = PhyloAnova._groups_and_response_matrix(
            traits,
            sorted(traits),
            group_idx=1,
            n_cols=3,
        )

        assert calls == [(0, 2)]
        assert len(groups) == 10
        assert Y.shape == (10, 2)

    def test_build_design_matrix_uses_cached_group_lookup(self):
        class NoIndexList(list):
            def index(self, *_args, **_kwargs):
                raise AssertionError("design matrix should cache group indices")

        groups = ["g1", "g2", "g3", "g1", "g3"]
        unique_groups = NoIndexList(["g1", "g2", "g3"])

        design = PhyloAnova._build_design_matrix(
            groups, unique_groups, len(groups)
        )

        np.testing.assert_array_equal(
            design,
            np.array(
                [
                    [1.0, 0.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [1.0, 0.0, 1.0],
                    [1.0, 0.0, 0.0],
                    [1.0, 0.0, 1.0],
                ]
            ),
        )

    def test_phylomorphospace_overlay_uses_direct_traversal(self, monkeypatch):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        pc = np.array([
            [0.0, 0.0],
            [2.0, 0.0],
            [0.0, 2.0],
            [2.0, 2.0],
        ])

        expected_parent_map = {}
        for clade in tree.find_clades(order="preorder"):
            for child in clade.clades:
                expected_parent_map[id(child)] = clade
        expected_coords = {
            id(tip): pc[ordered_names.index(tip.name)]
            for tip in tree.get_terminals()
        }
        for clade in tree.find_clades(order="postorder"):
            if id(clade) in expected_coords:
                continue
            expected_coords[id(clade)] = np.mean(
                [expected_coords[id(child)] for child in clade.clades],
                axis=0,
            )

        def fail_generic_traversal(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_generic_traversal)
        monkeypatch.setattr(TreeMixin, "get_terminals", fail_generic_traversal)

        parent_map, node_coords, preorder = (
            PhyloAnova._prepare_phylomorphospace_overlay(
                tree, ordered_names, pc
            )
        )

        assert set(parent_map) == set(expected_parent_map)
        assert [clade.name for clade in preorder if not clade.clades] == ordered_names
        for cid, expected in expected_coords.items():
            np.testing.assert_allclose(node_coords[cid], expected)

    def test_phylomorphospace_overlay_averages_available_child_coords(self):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "C", "D"]
        pc = np.array([
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
        ])

        parent_map, node_coords, preorder = (
            PhyloAnova._prepare_phylomorphospace_overlay(
                tree, ordered_names, pc
            )
        )

        tips = {tip.name: tip for tip in tree.get_terminals()}
        left_internal = parent_map[id(tips["A"])]
        right_internal = parent_map[id(tips["C"])]

        np.testing.assert_allclose(node_coords[id(left_internal)], [0.0, 0.0])
        np.testing.assert_allclose(node_coords[id(right_internal)], [2.0, 1.0])
        np.testing.assert_allclose(node_coords[id(tree.root)], [1.0, 0.5])
        assert preorder[0] is tree.root

    def test_cholesky_transform_matches_explicit_inverse(self):
        """Triangular-solve whitening matches the previous L inverse path."""
        rng = np.random.default_rng(123)
        A = rng.normal(size=(8, 8))
        vcv = A @ A.T + np.eye(8)
        Y = rng.normal(size=(8, 3))
        X_full = np.column_stack(
            [np.ones(8), [0, 1, 0, 1, 0, 1, 0, 1]]
        )
        X_red = np.ones((8, 1))

        fast = PhyloAnova._cholesky_transform(vcv, Y, X_full, X_red)
        L = np.linalg.cholesky(vcv)
        L_inv = np.linalg.inv(L)
        expected = (L_inv @ Y, L_inv @ X_full, L_inv @ X_red)

        for observed, expected_matrix in zip(fast, expected):
            np.testing.assert_allclose(observed, expected_matrix)

    def test_anova_deterministic_values(self, capsys):
        """Deterministic SS/MS/F match R validation values."""
        args = _make_args(permutations=99, seed=42)
        svc = PhyloAnova(args)
        svc.run()
        captured = capsys.readouterr()
        # SS_model = 0.0382, SS_resid = 0.2691, F = 0.3548
        assert "0.0382" in captured.out
        assert "0.2691" in captured.out
        assert "0.3548" in captured.out

    def test_anova_json_output(self, capsys):
        """JSON output has correct structure and deterministic values."""
        args = _make_args(json=True, permutations=99, seed=42)
        svc = PhyloAnova(args)
        svc.run()
        captured = capsys.readouterr()
        # Find the JSON line
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        assert payload["method"] == "anova"
        assert payload["n_taxa"] == 8
        assert payload["n_groups"] == 3
        tab = payload["anova_table"]
        assert tab["ss_model"] == pytest.approx(0.038184, abs=1e-4)
        assert tab["ss_resid"] == pytest.approx(0.269069, abs=1e-4)
        assert tab["f_stat"] == pytest.approx(0.354776, abs=1e-3)
        assert tab["df_model"] == 2
        assert tab["df_resid"] == 5

    def test_anova_rrpp_projection_ss_matches_materialized_residuals(self):
        rng = np.random.default_rng(20260628)
        n = 32
        n_groups = 4
        groups = np.array([idx % n_groups for idx in range(n)])
        Y = rng.normal(size=n)
        X_full = np.ones((n, n_groups))
        for group_idx in range(1, n_groups):
            X_full[:, group_idx] = groups == group_idx
        X_red = np.ones((n, 1))
        args = _make_args(permutations=37, seed=17)
        svc = PhyloAnova(args)

        observed = svc._run_anova(Y, X_full, X_red, n, n_groups)

        Q_red = svc._projection_basis(X_red)
        Q_full = svc._projection_basis(X_full)
        hat_red = svc._project(Q_red, Y)
        resid_red = Y - hat_red
        hat_full = svc._project(Q_full, Y)
        resid_full = Y - hat_full
        ss_total = float(resid_red @ resid_red)
        ss_resid = float(resid_full @ resid_full)
        ss_model = ss_total - ss_resid
        df_model = n_groups - 1
        df_resid = n - n_groups
        ms_model = ss_model / df_model
        ms_resid = ss_resid / df_resid
        f_obs = ms_model / ms_resid

        legacy_rng = np.random.RandomState(args.seed)
        f_perms = np.zeros(args.permutations)
        for perm_idx in range(args.permutations):
            permuted = legacy_rng.permutation(n)
            y_perm = hat_red + resid_red[permuted]
            hat_full_p = svc._project(Q_full, y_perm)
            resid_full_p = y_perm - hat_full_p
            ss_resid_p = float(resid_full_p @ resid_full_p)
            ss_model_p = ss_total - ss_resid_p
            ms_model_p = ss_model_p / df_model
            ms_resid_p = ss_resid_p / df_resid
            f_perms[perm_idx] = (
                ms_model_p / ms_resid_p if ms_resid_p > 0 else 0.0
            )

        assert observed["test_statistic_name"] == "F"
        assert observed["ss_model"] == pytest.approx(ss_model)
        assert observed["ss_resid"] == pytest.approx(ss_resid)
        assert observed["ss_total"] == pytest.approx(ss_total)
        assert observed["f_stat"] == pytest.approx(f_obs)
        assert observed["z_score"] == pytest.approx(
            (f_obs - np.mean(f_perms)) / np.std(f_perms)
        )
        assert observed["p_value"] == pytest.approx(
            float(np.mean(f_perms >= f_obs))
        )

    def test_manova_auto_detection(self, capsys):
        """Multi-trait file auto-detects MANOVA."""
        args = _make_args(traits=TRAITS_MULTI, permutations=99, seed=42)
        svc = PhyloAnova(args)
        svc.run()
        captured = capsys.readouterr()
        assert "MANOVA" in captured.out
        assert "Pillai" in captured.out

    def test_manova_deterministic_values(self, capsys):
        """MANOVA Pillai's trace matches R validation."""
        args = _make_args(
            traits=TRAITS_MULTI, permutations=99, seed=42, json=True,
        )
        svc = PhyloAnova(args)
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        tab = payload["anova_table"]
        assert tab["ss_model"] == pytest.approx(0.049484, abs=1e-4)
        assert tab["ss_resid"] == pytest.approx(0.415777, abs=1e-4)
        assert tab["pillai_trace"] == pytest.approx(0.794354, abs=1e-3)

    def test_pillai_trace_uses_array_sum(self, monkeypatch):
        eigenvalues = np.array([0.2, 0.5, 1.0, 1.5])
        expected = float(np.sum(eigenvalues / (1.0 + eigenvalues)))

        monkeypatch.setattr(
            phylo_anova_module.np,
            "sum",
            lambda *args, **kwargs: pytest.fail(
                "Pillai trace should use ndarray.sum"
            ),
        )

        assert phylo_anova_module._pillai_trace_from_eigenvalues(
            eigenvalues
        ) == pytest.approx(expected)

    def test_singular_value_total_variance_uses_dot_product(self, monkeypatch):
        singular_values = np.array([1.0, 2.0, 4.0, 8.0])
        expected = float(np.sum(singular_values ** 2))

        monkeypatch.setattr(
            phylo_anova_module.np,
            "sum",
            lambda *args, **kwargs: pytest.fail(
                "phylomorphospace PCA variance should use np.dot"
            ),
        )

        assert phylo_anova_module._singular_value_total_variance(
            singular_values
        ) == pytest.approx(expected)

    def test_manova_rrpp_projection_sscp_matches_materialized_residuals(self):
        rng = np.random.default_rng(20260628)
        n = 36
        n_groups = 4
        n_traits = 5
        groups = np.array([idx % n_groups for idx in range(n)])
        Y = rng.normal(size=(n, n_traits))
        X_full = np.ones((n, n_groups))
        for group_idx in range(1, n_groups):
            X_full[:, group_idx] = groups == group_idx
        X_red = np.ones((n, 1))
        args = _make_args(permutations=31, seed=19)
        svc = PhyloAnova(args)

        observed = svc._run_manova(Y, X_full, X_red, n, n_groups)

        Q_red = svc._projection_basis(X_red)
        Q_full = svc._projection_basis(X_full)
        hat_red = svc._project(Q_red, Y)
        resid_red = Y - hat_red
        hat_full = svc._project(Q_full, Y)
        resid_full = Y - hat_full
        SS_total = resid_red.T @ resid_red
        SS_resid = resid_full.T @ resid_full
        SS_model = SS_total - SS_resid
        df_model = n_groups - 1
        df_resid = n - n_groups

        SS_resid_inv = np.linalg.inv(SS_resid)
        H_E_inv = SS_model @ SS_resid_inv
        eigenvalues = np.real(np.linalg.eigvals(H_E_inv))
        pillai = float(np.sum(eigenvalues / (1.0 + eigenvalues)))

        s = min(df_model, n_traits)
        m = (abs(df_model - n_traits) - 1) / 2.0
        nn = (df_resid - n_traits - 1) / 2.0
        f_approx = (pillai / s) * ((2 * nn + s + 1) / (2 * m + s + 1))

        legacy_rng = np.random.RandomState(args.seed)
        pillai_perms = np.zeros(args.permutations)
        for perm_idx in range(args.permutations):
            permuted = legacy_rng.permutation(n)
            Y_perm = hat_red + resid_red[permuted]
            hat_full_p = svc._project(Q_full, Y_perm)
            resid_full_p = Y_perm - hat_full_p
            SS_resid_p = resid_full_p.T @ resid_full_p
            SS_model_p = SS_total - SS_resid_p
            SS_resid_inv_p = np.linalg.inv(SS_resid_p)
            H_E_inv_p = SS_model_p @ SS_resid_inv_p
            eig_p = np.real(np.linalg.eigvals(H_E_inv_p))
            pillai_perms[perm_idx] = float(
                np.sum(eig_p / (1.0 + eig_p))
            )

        assert observed["test_statistic_name"] == "Pillai's trace"
        assert observed["ss_model"] == pytest.approx(float(np.trace(SS_model)))
        assert observed["ss_resid"] == pytest.approx(float(np.trace(SS_resid)))
        assert observed["ss_total"] == pytest.approx(float(np.trace(SS_total)))
        assert observed["f_stat"] == pytest.approx(f_approx)
        assert observed["pillai_trace"] == pytest.approx(pillai)
        assert observed["z_score"] == pytest.approx(
            (pillai - np.mean(pillai_perms)) / np.std(pillai_perms)
        )
        assert observed["p_value"] == pytest.approx(
            float(np.mean(pillai_perms >= pillai))
        )

    def test_pairwise_output(self, capsys):
        """--pairwise produces pairwise comparisons."""
        args = _make_args(pairwise=True, permutations=99, seed=42)
        svc = PhyloAnova(args)
        svc.run()
        captured = capsys.readouterr()
        assert "Pairwise comparisons" in captured.out
        assert "carnivore vs herbivore" in captured.out
        assert "carnivore vs omnivore" in captured.out
        assert "herbivore vs omnivore" in captured.out

    def _stub_run_tail(self, monkeypatch, svc, captured):
        def build_vcv_matrix(tree, ordered_names):
            captured["tree"] = tree
            return np.eye(len(ordered_names))

        monkeypatch.setattr(
            "phykit.services.tree.vcv_utils.build_vcv_matrix",
            build_vcv_matrix,
        )
        monkeypatch.setattr(
            svc,
            "_parse_trait_file",
            lambda *_args, **_kwargs: (
                ["group", "trait"],
                {
                    "A": ["g1", 1.0],
                    "B": ["g1", 2.0],
                    "C": ["g2", 3.0],
                    "D": ["g2", 4.0],
                },
            ),
        )
        monkeypatch.setattr(
            svc,
            "_cholesky_transform",
            lambda _vcv, y, x_full, x_red: (y, x_full, x_red),
        )
        monkeypatch.setattr(
            svc,
            "_run_anova",
            lambda *_args, **_kwargs: {
                "ss_model": 0.0,
                "ss_resid": 0.0,
                "ss_total": 0.0,
                "df_model": 1,
                "df_resid": 2,
                "df_total": 3,
                "ms_model": 0.0,
                "ms_resid": 0.0,
                "f_stat": 0.0,
                "z_score": 0.0,
                "p_value": 1.0,
                "test_statistic_name": "F",
            },
        )
        monkeypatch.setattr(svc, "_print_results", lambda *_args, **_kwargs: None)

    def test_run_uses_read_only_tree_without_copy_when_branch_lengths_present(
        self, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        svc = PhyloAnova(_make_args(permutations=0))
        captured = {}

        monkeypatch.setattr(svc, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(
            svc,
            "read_tree_file",
            lambda: (_ for _ in ()).throw(
                AssertionError("mutable tree read should not be used")
            ),
        )
        monkeypatch.setattr(
            svc,
            "_fast_copy",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("tree should not be copied")
            ),
        )
        self._stub_run_tail(monkeypatch, svc, captured)

        svc.run()

        assert captured["tree"] is tree
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}

    def test_default_branch_length_scan_handles_mixed_child_counts(self, monkeypatch):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree scan should not use find_clades")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        assert PhyloAnova._needs_default_branch_lengths(tree) is False
        tree.root.clades[2].clades[1].branch_length = None
        assert PhyloAnova._needs_default_branch_lengths(tree) is True

    def test_run_missing_branch_lengths_copies_before_validation(
        self, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("((A,B):1,(C:1,D:1):1);"),
            "newick",
        )
        svc = PhyloAnova(_make_args(permutations=0))
        original_fast_copy = svc._fast_copy
        copied_trees = []
        captured = {}

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        monkeypatch.setattr(svc, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(svc, "_fast_copy", copy_spy)
        self._stub_run_tail(monkeypatch, svc, captured)

        svc.run()

        assert len(copied_trees) == 1
        assert captured["tree"] is copied_trees[0]
        assert tree.root.clades[0].clades[0].branch_length is None
        assert copied_trees[0].root.clades[0].clades[0].branch_length == 1e-8

    def test_print_results_batches_text_rows(self, mocker):
        args = _make_args(pairwise=True, permutations=99, seed=42)
        svc = PhyloAnova(args)
        printed = mocker.patch("builtins.print")
        result = {
            "df_model": 2,
            "ss_model": 0.1234,
            "ms_model": 0.0617,
            "f_stat": 4.25,
            "z_score": 2.1,
            "p_value": 0.03,
            "df_resid": 7,
            "ss_resid": 1.25,
            "ms_resid": 0.1786,
            "df_total": 9,
            "ss_total": 1.3734,
            "pillai_trace": 0.4567,
        }
        pairwise_results = [
            {
                "group1": "carnivore",
                "group2": "herbivore",
                "distance": 1.25,
                "z_score": 0.5,
                "p_value": 0.04,
            },
            {
                "group1": "carnivore",
                "group2": "omnivore",
                "distance": 2.5,
                "z_score": 1.5,
                "p_value": 0.20,
            },
        ]

        svc._print_results(
            "manova",
            result,
            ["carnivore", "herbivore", "omnivore"],
            ["mass", "length"],
            pairwise_results,
            n=10,
            n_groups=3,
        )

        printed.assert_called_once_with(
            "Phylogenetic MANOVA (RRPP, 99 permutations)\n"
            "Response traits: mass, length\n"
            "Groups (3): carnivore, herbivore, omnivore\n"
            "N taxa: 10\n"
            "\n"
            "Source             Df           SS           MS          F          Z    p-value\n"
            "----------------------------------------------------------------------------\n"
            "group               2       0.1234       0.0617     4.2500     2.1000     0.0300\n"
            "residuals           7       1.2500       0.1786\n"
            "total               9       1.3734\n"
            "\n"
            "Pillai's trace: 0.4567\n"
            "\n"
            "Pairwise comparisons:\n"
            "  Comparison                              d          Z    p-value\n"
            "  --------------------------------------------------------------\n"
            "  carnivore vs herbivore             1.2500     0.5000     0.0400\n"
            "  carnivore vs omnivore              2.5000     1.5000     0.2000"
        )

    def test_pairwise_json(self, capsys):
        """Pairwise results in JSON output."""
        args = _make_args(pairwise=True, permutations=99, seed=42, json=True)
        svc = PhyloAnova(args)
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        assert "pairwise" in payload
        assert len(payload["pairwise"]) == 3
        groups = {(pw["group1"], pw["group2"]) for pw in payload["pairwise"]}
        assert ("carnivore", "herbivore") in groups

    def test_pairwise_rrpp_matches_materialized_permutations(self):
        rng = np.random.default_rng(123)
        n = 30
        Y = rng.normal(size=(n, 3))
        groups = np.array(["a", "b", "c"] * 10)
        unique_groups = ["a", "b", "c"]
        X_red = np.ones((n, 1))
        args = _make_args(permutations=25, seed=7)
        svc = PhyloAnova(args)

        optimized = svc._run_pairwise(Y, list(groups), unique_groups, X_red)

        legacy_rng = np.random.RandomState(args.seed)
        beta_red = np.linalg.lstsq(X_red, Y, rcond=None)[0]
        hat_red = X_red @ beta_red
        resid_red = Y - hat_red
        expected = []
        for i in range(len(unique_groups)):
            for j in range(i + 1, len(unique_groups)):
                g1, g2 = unique_groups[i], unique_groups[j]
                mask = (groups == g1) | (groups == g2)
                idx = np.where(mask)[0]
                Y_sub = Y[idx]
                groups_sub = groups[idx]
                mean1 = Y_sub[groups_sub == g1].mean(axis=0)
                mean2 = Y_sub[groups_sub == g2].mean(axis=0)
                d_obs = float(np.linalg.norm(mean1 - mean2))
                resid_sub = resid_red[idx]
                hat_sub = hat_red[idx]
                d_perms = np.zeros(args.permutations)
                for p in range(args.permutations):
                    perm_idx = legacy_rng.permutation(len(idx))
                    Y_perm = hat_sub + resid_sub[perm_idx]
                    mean1_p = Y_perm[groups_sub == g1].mean(axis=0)
                    mean2_p = Y_perm[groups_sub == g2].mean(axis=0)
                    d_perms[p] = float(np.linalg.norm(mean1_p - mean2_p))
                p_val = float(np.mean(d_perms >= d_obs))
                z_val = (
                    (d_obs - np.mean(d_perms)) / np.std(d_perms)
                    if np.std(d_perms) > 0
                    else 0.0
                )
                expected.append(
                    dict(group1=g1, group2=g2, distance=d_obs, z_score=z_val, p_value=p_val)
                )

        assert len(optimized) == len(expected)
        for observed, reference in zip(optimized, expected):
            assert observed["group1"] == reference["group1"]
            assert observed["group2"] == reference["group2"]
            assert observed["distance"] == pytest.approx(reference["distance"])
            assert observed["z_score"] == pytest.approx(reference["z_score"])
            assert observed["p_value"] == pytest.approx(reference["p_value"])

    def test_pairwise_rrpp_matches_materialized_unbalanced_groups(self):
        rng = np.random.default_rng(321)
        groups = np.array(["a"] * 7 + ["b"] * 11 + ["c"] * 5)
        unique_groups = ["a", "b", "c"]
        n = len(groups)
        Y = rng.normal(size=(n, 4))
        X_red = np.ones((n, 1))
        args = _make_args(permutations=17, seed=11)
        svc = PhyloAnova(args)

        optimized = svc._run_pairwise(Y, list(groups), unique_groups, X_red)

        legacy_rng = np.random.RandomState(args.seed)
        beta_red = np.linalg.lstsq(X_red, Y, rcond=None)[0]
        hat_red = X_red @ beta_red
        resid_red = Y - hat_red
        expected = []
        for i in range(len(unique_groups)):
            for j in range(i + 1, len(unique_groups)):
                g1, g2 = unique_groups[i], unique_groups[j]
                mask = (groups == g1) | (groups == g2)
                idx = np.where(mask)[0]
                Y_sub = Y[idx]
                groups_sub = groups[idx]
                mean1 = Y_sub[groups_sub == g1].mean(axis=0)
                mean2 = Y_sub[groups_sub == g2].mean(axis=0)
                d_obs = float(np.linalg.norm(mean1 - mean2))
                resid_sub = resid_red[idx]
                hat_sub = hat_red[idx]
                d_perms = np.zeros(args.permutations)
                for p in range(args.permutations):
                    perm_idx = legacy_rng.permutation(len(idx))
                    Y_perm = hat_sub + resid_sub[perm_idx]
                    mean1_p = Y_perm[groups_sub == g1].mean(axis=0)
                    mean2_p = Y_perm[groups_sub == g2].mean(axis=0)
                    d_perms[p] = float(np.linalg.norm(mean1_p - mean2_p))
                p_val = float(np.mean(d_perms >= d_obs))
                z_val = (
                    (d_obs - np.mean(d_perms)) / np.std(d_perms)
                    if np.std(d_perms) > 0
                    else 0.0
                )
                expected.append(
                    dict(group1=g1, group2=g2, distance=d_obs, z_score=z_val, p_value=p_val)
                )

        for observed, reference in zip(optimized, expected):
            assert observed["group1"] == reference["group1"]
            assert observed["group2"] == reference["group2"]
            assert observed["distance"] == pytest.approx(reference["distance"])
            assert observed["z_score"] == pytest.approx(reference["z_score"])
            assert observed["p_value"] == pytest.approx(reference["p_value"])

    def test_pairwise_rrpp_generates_legacy_permutation_sequence(self, monkeypatch):
        rng = np.random.default_rng(123)
        n = 24
        Y = rng.normal(size=(n, 2))
        groups = np.array(["a", "b", "c"] * 8)
        unique_groups = ["a", "b", "c"]
        X_red = np.ones((n, 1))
        args = _make_args(permutations=5, seed=7)
        svc = PhyloAnova(args)

        expected_rng = np.random.RandomState(args.seed)
        expected_perms = []
        for i in range(len(unique_groups)):
            for j in range(i + 1, len(unique_groups)):
                pair_size = np.count_nonzero(
                    (groups == unique_groups[i]) | (groups == unique_groups[j])
                )
                for _ in range(args.permutations):
                    expected_perms.append(expected_rng.permutation(pair_size))

        actual_perms = []
        original_random_state = np.random.RandomState

        class TrackingRandomState:
            def __init__(self, *state_args, **state_kwargs):
                self._rng = original_random_state(*state_args, **state_kwargs)

            def permutation(self, *perm_args, **perm_kwargs):
                perm = self._rng.permutation(*perm_args, **perm_kwargs)
                actual_perms.append(perm.copy())
                return perm

        monkeypatch.setattr(np.random, "RandomState", TrackingRandomState)

        def fail_where(*_args, **_kwargs):
            raise AssertionError("pairwise RRPP should use flat mask indices")

        monkeypatch.setattr(phylo_anova_module.np, "where", fail_where)

        svc._run_pairwise(Y, list(groups), unique_groups, X_red)

        assert len(actual_perms) == len(expected_perms)
        for observed, expected in zip(actual_perms, expected_perms):
            np.testing.assert_array_equal(observed, expected)

    def test_method_anova_with_multi_trait_raises(self):
        """--method anova with multiple traits raises error."""
        args = _make_args(traits=TRAITS_MULTI, method="anova")
        svc = PhyloAnova(args)
        with pytest.raises(SystemExit) as exc_info:
            svc.run()
        assert exc_info.value.code == 2

    def test_missing_trait_file_raises(self):
        """Nonexistent trait file raises error."""
        args = _make_args(traits="nonexistent.tsv")
        svc = PhyloAnova(args)
        with pytest.raises(SystemExit) as exc_info:
            svc.run()
        assert exc_info.value.code == 2

    def test_bad_group_column_raises(self):
        """Nonexistent group column raises error."""
        args = _make_args(group_column="nonexistent")
        svc = PhyloAnova(args)
        with pytest.raises(SystemExit) as exc_info:
            svc.run()
        assert exc_info.value.code == 2

    def test_plot_boxplot_creates_file(self, tmp_path):
        """Boxplot is created."""
        out = str(tmp_path / "anova_boxplot.png")
        args = _make_args(
            plot_output=out, plot_type="boxplot",
            permutations=99, seed=42,
        )
        svc = PhyloAnova(args)
        svc.run()
        assert (tmp_path / "anova_boxplot.png").exists()

    def test_plot_boxplot_uses_vectorized_group_masks(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        groups = ["a", "b", "a", "c", "b", "c"]
        converted_groups = 0
        original_asarray = phylo_anova_module.np.asarray

        def tracking_asarray(value, *args, **kwargs):
            nonlocal converted_groups
            if value is groups:
                converted_groups += 1
            return original_asarray(value, *args, **kwargs)

        monkeypatch.setattr(phylo_anova_module.np, "asarray", tracking_asarray)
        svc = PhyloAnova(
            _make_args(plot_output=str(tmp_path / "anova_boxplot.png"))
        )

        svc._plot_boxplot(
            None,
            ["t1", "t2", "t3", "t4", "t5", "t6"],
            np.array([[0.1], [1.1], [0.2], [2.1], [1.2], [2.2]]),
            groups,
            ["a", "b", "c"],
            ["trait"],
            str(tmp_path / "anova_boxplot.png"),
        )

        assert converted_groups == 1
        assert (tmp_path / "anova_boxplot.png").exists()

    def test_plot_phylomorphospace_creates_file(self, tmp_path):
        """Phylomorphospace is created for MANOVA."""
        out = str(tmp_path / "manova_morpho.png")
        args = _make_args(
            traits=TRAITS_MULTI,
            plot_output=out, plot_type="phylomorphospace",
            permutations=99, seed=42,
        )
        svc = PhyloAnova(args)
        svc.run()
        assert (tmp_path / "manova_morpho.png").exists()

    def test_plot_phylomorphospace_batches_branch_segments(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import LineCollection

        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ordered_names = ["A", "B", "C", "D"]
        Y = np.array([
            [0.0, 0.0],
            [2.0, 0.0],
            [0.0, 2.0],
            [2.0, 2.0],
        ])
        groups = ["g1", "g1", "g2", "g2"]
        args = _make_args(
            traits=TRAITS_MULTI,
            plot_output=str(tmp_path / "manova_morpho.png"),
            no_title=True,
        )
        svc = PhyloAnova(args)

        original_add_collection = matplotlib.axes.Axes.add_collection
        original_asarray = phylo_anova_module.np.asarray
        line_collections = []
        converted_groups = 0

        def fail_plot(*args, **kwargs):
            raise AssertionError("phylomorphospace branches should use LineCollection")

        def capture_collection(self, collection, *args, **kwargs):
            if isinstance(collection, LineCollection):
                line_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        def tracking_asarray(value, *args, **kwargs):
            nonlocal converted_groups
            if value is groups:
                converted_groups += 1
            return original_asarray(value, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "plot", fail_plot)
        monkeypatch.setattr(matplotlib.axes.Axes, "add_collection", capture_collection)
        monkeypatch.setattr(phylo_anova_module.np, "asarray", tracking_asarray)

        svc._plot_phylomorphospace(
            tree,
            ordered_names,
            Y,
            groups,
            ["g1", "g2"],
            ["x", "y"],
            str(tmp_path / "manova_morpho.png"),
        )

        assert any(
            len(collection.get_segments()) == 6
            for collection in line_collections
        )
        assert converted_groups == 1
        assert (tmp_path / "manova_morpho.png").exists()

    def test_reproducible_with_seed(self, capsys):
        """Same seed gives same results."""
        args1 = _make_args(permutations=99, seed=123, json=True)
        svc1 = PhyloAnova(args1)
        svc1.run()
        out1 = capsys.readouterr().out

        args2 = _make_args(permutations=99, seed=123, json=True)
        svc2 = PhyloAnova(args2)
        svc2.run()
        out2 = capsys.readouterr().out

        # Extract JSON lines
        j1 = j2 = None
        for line in out1.strip().split("\n"):
            if line.startswith("{"):
                j1 = json.loads(line)
        for line in out2.strip().split("\n"):
            if line.startswith("{"):
                j2 = json.loads(line)
        assert j1["anova_table"]["p_value"] == j2["anova_table"]["p_value"]
        assert j1["anova_table"]["z_score"] == j2["anova_table"]["z_score"]
