import builtins
import importlib
import math
import subprocess
import sys
from argparse import Namespace
from io import StringIO
import json

import numpy as np
import pytest
from Bio import Phylo

from phykit.services.tree.phylo_path import PhyloPath
import phykit.services.tree.phylo_path as phylo_path_module
from phykit.errors import PhykitUserError


TREE = "tests/sample_files/tree_simple.tre"
TRAITS = "tests/sample_files/tree_simple_multi_traits.tsv"
MODELS = "tests/sample_files/phylo_path_models.txt"


def test_module_import_does_not_import_scipy_linalg(monkeypatch):
    module_name = "phykit.services.tree.phylo_path"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.linalg" or name.startswith("scipy.linalg."):
            raise AssertionError(
                "phylo_path module import should not import SciPy linalg"
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
import phykit.services.tree.phylo_path as module
assert callable(module.print_json)
assert callable(module.trait_matrix_from_rows)
assert hasattr(module.np, "__getattr__")
assert "phykit.helpers.trait_parsing" not in sys.modules
assert "json" not in sys.modules
assert "typing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy.linalg" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def _make_args(**overrides):
    defaults = dict(
        tree=TREE, traits=TRAITS, models=MODELS,
        best_only=False, plot_output=None, csv=None, json=False,
        fig_width=None, fig_height=None, dpi=150, no_title=False,
        title=None, legend_position=None, ylabel_fontsize=None,
        xlabel_fontsize=None, title_fontsize=None, axis_fontsize=None,
        colors=None, ladderize=False, cladogram=False, circular=False,
        color_file=None,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


class TestPhyloPath:
    def test_probability_helpers_match_expected_values(self):
        assert phylo_path_module._chi2_sf(5.991464547107979, 2) == pytest.approx(0.05)
        assert phylo_path_module._chi2_sf(6.0, 6) == pytest.approx(
            np.exp(-3.0) * (1.0 + 3.0 + 4.5)
        )
        assert phylo_path_module._t_two_tailed_p_value(2.2281388519649385, 10) == pytest.approx(0.05)

    def test_fishers_c_uses_scalar_math_log(self, monkeypatch):
        def fail_np_log(*_args, **_kwargs):
            raise AssertionError("Fisher's C scalar p-values should use math.log")

        monkeypatch.setattr(phylo_path_module.np, "log", fail_np_log, raising=False)

        p_values = [0.5, 0.25, 0.125, 0.75]
        expected = -2.0 * sum(math.log(p_value) for p_value in p_values)

        assert phylo_path_module._fishers_c_from_p_values(p_values) == pytest.approx(
            expected
        )

    def test_relative_likelihood_uses_scalar_math_exp(self, monkeypatch):
        def fail_np_exp(*_args, **_kwargs):
            raise AssertionError("model weights should use scalar math.exp")

        monkeypatch.setattr(phylo_path_module.np, "exp", fail_np_exp, raising=False)

        assert phylo_path_module._relative_likelihood_from_delta(
            2.5
        ) == pytest.approx(math.exp(-1.25))

    def test_probability_helpers_do_not_import_scipy_stats(self, monkeypatch):
        original_import = __import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "scipy.stats" or name.startswith("scipy.stats."):
                raise AssertionError("phylo_path probability helpers should not import scipy.stats")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr("builtins.__import__", fake_import)

        assert phylo_path_module._chi2_sf(5.991464547107979, 2) == pytest.approx(0.05)
        assert phylo_path_module._t_two_tailed_p_value(2.2281388519649385, 10) == pytest.approx(0.05)

    def test_even_df_chi2_helper_does_not_import_scipy_special(self, monkeypatch):
        previous_chdtrc = phylo_path_module._CHDTRC
        phylo_path_module._CHDTRC = None
        original_import = builtins.__import__

        def fail_special_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "scipy.special":
                raise AssertionError("even-df chi-square p-values should not import scipy.special")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fail_special_import)
        try:
            assert phylo_path_module._chi2_sf(5.991464547107979, 2) == pytest.approx(0.05)
            assert phylo_path_module._chi2_sf(6.0, 6) == pytest.approx(
                np.exp(-3.0) * (1.0 + 3.0 + 4.5)
            )
        finally:
            phylo_path_module._CHDTRC = previous_chdtrc

    def test_probability_helpers_cache_scipy_special_imports(self, monkeypatch):
        previous_chdtrc = phylo_path_module._CHDTRC
        previous_stdtr = phylo_path_module._STDTR
        phylo_path_module._CHDTRC = None
        phylo_path_module._STDTR = None
        original_import = builtins.__import__
        special_imports = 0

        def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
            nonlocal special_imports
            if name == "scipy.special":
                special_imports += 1
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", counting_import)
        try:
            phylo_path_module._chi2_sf(7.814727903251179, 3)
            phylo_path_module._chi2_sf(11.344866730144373, 3)
            phylo_path_module._t_two_tailed_p_value(2.2281388519649385, 10)
            phylo_path_module._t_two_tailed_p_value(1.8124611228107335, 10)
        finally:
            phylo_path_module._CHDTRC = previous_chdtrc
            phylo_path_module._STDTR = previous_stdtr

        assert special_imports == 2

    def test_parse_trait_file_streams_comments_and_blanks(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "# comment\n\n"
            "taxon\tbody_mass\tbrain_size\n"
            "A\t1.0\t10.0\n"
            "# ignored between rows\n"
            "B\t2.0\t20.0\n"
            "C\t3.0\t30.0\n"
            "D\t4.0\t40.0\n"
            "off_tree\t5.0\t50.0\n"
        )
        svc = PhyloPath(_make_args())

        trait_names, traits = svc._parse_trait_file(
            str(trait_file), ["A", "B", "C", "D"]
        )

        assert trait_names == ["body_mass", "brain_size"]
        assert traits == {
            "A": [1.0, 10.0],
            "B": [2.0, 20.0],
            "C": [3.0, 30.0],
            "D": [4.0, 40.0],
        }

    def test_parse_trait_file_all_shared_taxa(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "taxon\tbody_mass\tbrain_size\n"
            "A\t1.0\t10.0\n"
            "B\t2.0\t20.0\n"
            "C\t3.0\t30.0\n"
            "D\t4.0\t40.0\n"
        )
        svc = PhyloPath(_make_args())

        trait_names, traits = svc._parse_trait_file(
            str(trait_file), ["D", "B", "A", "C"]
        )

        assert trait_names == ["body_mass", "brain_size"]
        assert traits == {
            "A": [1.0, 10.0],
            "B": [2.0, 20.0],
            "C": [3.0, 30.0],
            "D": [4.0, 40.0],
        }

    def test_parse_trait_file_non_numeric_error(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "taxon\tbody_mass\tbrain_size\n"
            "A\t1.0\t10.0\n"
            "B\tbad\t20.0\n"
            "C\t3.0\t30.0\n"
            "D\t4.0\t40.0\n"
        )
        svc = PhyloPath(_make_args())

        with pytest.raises(PhykitUserError) as exc_info:
            svc._parse_trait_file(str(trait_file), ["A", "B", "C", "D"])

        assert (
            "Non-numeric value 'bad' for taxon 'B' on line 3."
            in exc_info.value.messages
        )

    def test_parse_models_file_streams_without_readlines(self, monkeypatch):
        class StreamingOnlyFile:
            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, traceback):
                return False

            def __iter__(self):
                return iter(
                    [
                        "# ignored\n",
                        "\n",
                        "direct: body_mass->brain_size\n",
                        "indirect: body_mass->body_length, body_length->brain_size\n",
                    ]
                )

            def readlines(self):
                raise AssertionError("model parser should stream rows")

        monkeypatch.setattr("builtins.open", lambda path: StreamingOnlyFile())
        svc = PhyloPath(_make_args())

        models = svc._parse_models_file("models.txt")

        assert models == {
            "direct": [("body_mass", "brain_size")],
            "indirect": [
                ("body_mass", "body_length"),
                ("body_length", "brain_size"),
            ],
        }

    def test_is_dag_handles_wide_acyclic_graph(self):
        variables = [f"v{i}" for i in range(200)]
        adj = {v: set() for v in variables}
        for i in range(100):
            adj[variables[i]].add(variables[100 + i])

        assert PhyloPath._is_dag(adj, variables) is True

    def test_is_dag_rejects_cycle(self):
        variables = ["A", "B", "C"]
        adj = {"A": {"B"}, "B": {"C"}, "C": {"A"}}

        assert PhyloPath._is_dag(adj, variables) is False

    def test_topological_order_preserves_existing_queue_order(self):
        variables = [f"v{i}" for i in range(10)]
        adj = {v: set() for v in variables}
        for i in range(5):
            adj[variables[i]].add(variables[5 + i])

        assert PhyloPath._topological_order(adj, variables) == [
            "v0",
            "v1",
            "v2",
            "v3",
            "v4",
            "v5",
            "v6",
            "v7",
            "v8",
            "v9",
        ]

    def test_basis_set_reuses_precomputed_parent_sets(self, monkeypatch):
        svc = PhyloPath(_make_args())
        variables = ["A", "B", "C", "D"]
        adj = {
            "A": {"B", "C"},
            "B": set(),
            "C": set(),
            "D": set(),
        }

        monkeypatch.setattr(
            PhyloPath,
            "_get_parents",
            staticmethod(
                lambda *_args, **_kwargs: (_ for _ in ()).throw(
                    AssertionError("basis set should not scan parents per pair")
                )
            ),
        )

        assert svc._basis_set(adj, variables) == [
            ("A", "D", set()),
            ("D", "B", {"A"}),
            ("D", "C", {"A"}),
            ("B", "C", {"A"}),
        ]

    def test_model_average_indexes_coefficients_by_edge(self):
        class NoContainsDict(dict):
            def __contains__(self, key):
                raise AssertionError("model averaging should not scan every model per edge")

        svc = PhyloPath(_make_args())
        model_results = {
            "m1": {
                "weight": 1.0,
                "coefficients": NoContainsDict(
                    {
                        "b -> c": {"coef": 4.0, "se": 0.5},
                        "a -> b": {"coef": 1.0, "se": 0.2},
                    }
                ),
            },
            "m2": {
                "weight": 3.0,
                "coefficients": NoContainsDict(
                    {"a -> b": {"coef": 5.0, "se": 0.4}}
                ),
            },
            "m3": {
                "weight": 0.0,
                "coefficients": NoContainsDict(
                    {"zero -> edge": {"coef": 9.0, "se": 0.1}}
                ),
            },
        }

        averaged = svc._model_average(model_results, ["a", "b", "c"])

        assert list(averaged) == ["a -> b", "b -> c"]
        assert averaged["a -> b"]["coef"] == pytest.approx(4.0)
        assert averaged["a -> b"]["se"] == pytest.approx(3.13 ** 0.5)
        assert averaged["a -> b"]["lower"] == pytest.approx(
            4.0 - 1.96 * (3.13 ** 0.5)
        )
        assert averaged["a -> b"]["upper"] == pytest.approx(
            4.0 + 1.96 * (3.13 ** 0.5)
        )
        assert averaged["b -> c"]["coef"] == pytest.approx(4.0)
        assert averaged["b -> c"]["se"] == pytest.approx(0.5)

    def test_basic_run(self, capsys):
        """Runs without error and prints model comparison."""
        svc = PhyloPath(_make_args())
        svc.run()
        captured = capsys.readouterr()
        assert "Phylogenetic Path Analysis" in captured.out
        assert "Model comparison" in captured.out
        assert "Path coefficients" in captured.out

    def test_three_models_ranked(self, capsys):
        """All three models appear in output."""
        svc = PhyloPath(_make_args())
        svc.run()
        captured = capsys.readouterr()
        assert "direct" in captured.out
        assert "mediated" in captured.out
        assert "full" in captured.out

    def test_best_model_identified(self, capsys):
        """Best model is identified."""
        svc = PhyloPath(_make_args())
        svc.run()
        captured = capsys.readouterr()
        assert "Best model:" in captured.out

    def test_model_averaging_default(self, capsys):
        """Default output uses model averaging."""
        svc = PhyloPath(_make_args())
        svc.run()
        captured = capsys.readouterr()
        assert "model-averaged" in captured.out

    def test_run_resolves_model_variables_with_first_index_map(self, monkeypatch):
        svc = PhyloPath(_make_args())
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        captured = {}

        monkeypatch.setattr(svc, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(svc, "validate_tree", lambda *_args, **_kwargs: None)
        monkeypatch.setattr(svc, "get_tip_names_from_tree", lambda _tree: ["A", "B", "C", "D"])
        monkeypatch.setattr(
            svc,
            "_parse_trait_file",
            lambda *_args, **_kwargs: (
                ["x", "y", "x", "y"],
                {
                    "A": [1.0, 2.0, 99.0, 99.0],
                    "B": [2.0, 3.0, 99.0, 99.0],
                    "C": [3.0, 4.0, 99.0, 99.0],
                    "D": [4.0, 5.0, 99.0, 99.0],
                },
            ),
        )
        monkeypatch.setattr(
            svc,
            "_parse_models_file",
            lambda *_args, **_kwargs: {"direct": [("x", "y")]},
        )

        def trait_matrix(traits, ordered_names):
            matrix = np.array([traits[name] for name in ordered_names], dtype=float)
            captured["matrix"] = matrix
            return matrix

        monkeypatch.setattr(
            "phykit.services.tree.phylo_path.trait_matrix_from_rows",
            trait_matrix,
        )
        monkeypatch.setattr(
            "phykit.services.tree.vcv_utils.build_vcv_matrix",
            lambda _tree, ordered_names: np.eye(len(ordered_names)),
        )
        monkeypatch.setattr(svc, "_is_dag", lambda *_args, **_kwargs: True)
        monkeypatch.setattr(svc, "_basis_set", lambda *_args, **_kwargs: [])

        def estimate_path_coefficients(edges, variables, z_data, *_args, **_kwargs):
            captured["z_data"] = z_data
            return {"x -> y": {"coef": 1.0, "se": 0.1, "p": 0.5}}

        monkeypatch.setattr(
            svc,
            "_estimate_path_coefficients",
            estimate_path_coefficients,
        )
        monkeypatch.setattr(svc, "_model_average", lambda results, variables: {})
        monkeypatch.setattr(svc, "_print_text", lambda *_args, **_kwargs: None)

        svc.run()

        expected_x = np.array([1.0, 2.0, 3.0, 4.0])
        expected_y = np.array([2.0, 3.0, 4.0, 5.0])
        expected_x = (expected_x - expected_x.mean()) / expected_x.std()
        expected_y = (expected_y - expected_y.mean()) / expected_y.std()

        # The duplicate x/y columns contain 99s; first-index selection keeps them out.
        np.testing.assert_allclose(captured["z_data"]["x"], expected_x)
        np.testing.assert_allclose(captured["z_data"]["y"], expected_y)

    def test_run_standardizes_constant_variable_without_division(self, monkeypatch):
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        svc = PhyloPath(_make_args())
        captured = {}

        monkeypatch.setattr(svc, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(svc, "validate_tree", lambda *_args, **_kwargs: None)
        monkeypatch.setattr(svc, "get_tip_names_from_tree", lambda _tree: ["A", "B", "C", "D"])
        monkeypatch.setattr(
            svc,
            "_parse_trait_file",
            lambda *_args, **_kwargs: (
                ["x", "y"],
                {
                    "A": [1.0, 2.0],
                    "B": [1.0, 3.0],
                    "C": [1.0, 4.0],
                    "D": [1.0, 5.0],
                },
            ),
        )
        monkeypatch.setattr(
            svc,
            "_parse_models_file",
            lambda *_args, **_kwargs: {"direct": [("x", "y")]},
        )
        monkeypatch.setattr(
            "phykit.services.tree.vcv_utils.build_vcv_matrix",
            lambda _tree, ordered_names: np.eye(len(ordered_names)),
        )
        monkeypatch.setattr(svc, "_is_dag", lambda *_args, **_kwargs: True)
        monkeypatch.setattr(svc, "_basis_set", lambda *_args, **_kwargs: [])

        def estimate_path_coefficients(edges, variables, z_data, *_args, **_kwargs):
            captured["z_data"] = z_data
            return {}

        monkeypatch.setattr(
            svc,
            "_estimate_path_coefficients",
            estimate_path_coefficients,
        )
        monkeypatch.setattr(svc, "_model_average", lambda *_args, **_kwargs: {})
        monkeypatch.setattr(svc, "_print_text", lambda *_args, **_kwargs: None)

        svc.run()

        expected_y = np.array([2.0, 3.0, 4.0, 5.0])
        expected_y = (expected_y - expected_y.mean()) / expected_y.std()
        np.testing.assert_allclose(captured["z_data"]["x"], np.zeros(4))
        np.testing.assert_allclose(captured["z_data"]["y"], expected_y)

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
                ["x", "y", "z"],
                {
                    "A": [1.0, 2.0, 3.0],
                    "B": [2.0, 3.0, 4.0],
                    "C": [3.0, 4.0, 5.0],
                    "D": [4.0, 5.0, 6.0],
                },
            ),
        )
        monkeypatch.setattr(
            svc,
            "_parse_models_file",
            lambda *_args, **_kwargs: {
                "direct": [("x", "y"), ("y", "z")],
                "alt": [("x", "z"), ("z", "y")],
            },
        )
        def dsep_test(*args, **_kwargs):
            captured.setdefault("z_data", args[3])
            return 0.5

        monkeypatch.setattr(svc, "_dsep_test", dsep_test)
        monkeypatch.setattr(svc, "_estimate_path_coefficients", lambda *_args, **_kwargs: {})
        monkeypatch.setattr(svc, "_model_average", lambda *_args, **_kwargs: {})
        monkeypatch.setattr(svc, "_print_text", lambda *_args, **_kwargs: None)

    def test_run_uses_read_only_tree_without_copy_when_branch_lengths_present(
        self, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )
        svc = PhyloPath(_make_args())
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
        expected = [-1.3416407865, -0.4472135955, 0.4472135955, 1.3416407865]
        assert captured["z_data"]["x"].tolist() == pytest.approx(expected)
        assert captured["z_data"]["y"].tolist() == pytest.approx(expected)
        assert captured["z_data"]["z"].tolist() == pytest.approx(expected)

    def test_default_branch_length_scan_handles_mixed_child_counts(
        self, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree scan should not use find_clades")

        monkeypatch.setattr(tree, "find_clades", fail_find_clades)

        assert PhyloPath._needs_default_branch_lengths(tree) is False
        tree.root.clades[2].clades[1].branch_length = None
        assert PhyloPath._needs_default_branch_lengths(tree) is True

    def test_run_missing_branch_lengths_copies_before_validation(
        self, monkeypatch
    ):
        tree = Phylo.read(
            StringIO("((A,B):1,(C:1,D:1):1);"),
            "newick",
        )
        svc = PhyloPath(_make_args())
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

    def test_print_text_batches_report_rows(self, mocker):
        svc = PhyloPath(_make_args())
        printed = mocker.patch("builtins.print")
        ranked = [
            (
                "direct",
                {
                    "k": 2,
                    "q": 1,
                    "C": 3.25,
                    "CICc": 4.5,
                    "delta": 0.0,
                    "weight": 0.75,
                    "p": 0.02,
                },
            ),
            (
                "mediated",
                {
                    "k": 3,
                    "q": 2,
                    "C": 6.25,
                    "CICc": 8.5,
                    "delta": 4.0,
                    "weight": 0.25,
                    "p": 0.20,
                },
            ),
        ]
        avg_coefs = {
            "body_mass -> brain_size": {
                "coef": 1.25,
                "se": 0.25,
                "lower": 0.75,
                "upper": 1.75,
                "p": 0.004,
            },
            "longevity -> brain_size": {
                "coef": -0.5,
                "se": 0.5,
                "lower": -1.5,
                "upper": 0.5,
                "p": 0.20,
            },
        }

        svc._print_text(ranked, avg_coefs, "model-averaged", 8)

        printed.assert_called_once_with(
            "Phylogenetic Path Analysis (n = 8 taxa)\n"
            "\n"
            "Model comparison (ranked by CICc):\n"
            "  Model                   k    q          C       CICc    delta   weight        p\n"
            "  --------------------------------------------------------------------------\n"
            "  direct                  2    1     3.2500     4.5000   0.0000   0.7500   0.0200\n"
            "  mediated                3    2     6.2500     8.5000   4.0000   0.2500   0.2000\n"
            "\n"
            "Best model: direct (CICc weight = 0.7500)\n"
            "\n"
            "Path coefficients (model-averaged):\n"
            "  Path                            coef         SE      lower      upper\n"
            "  -------------------------------------------------------------------\n"
            "  body_mass -> brain_size       1.2500     0.2500     0.7500     1.7500 **\n"
            "  longevity -> brain_size      -0.5000     0.5000    -1.5000     0.5000 "
        )

    def test_best_only_flag(self, capsys):
        """--best-only shows only best model coefficients."""
        svc = PhyloPath(_make_args(best_only=True))
        svc.run()
        captured = capsys.readouterr()
        assert "best model:" in captured.out

    def test_json_structure(self, capsys):
        """JSON output has correct structure."""
        svc = PhyloPath(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        assert payload["n_taxa"] == 8
        assert payload["n_models"] == 3
        assert "model_comparison" in payload
        assert "path_coefficients" in payload
        assert len(payload["model_comparison"]) == 3

    def test_ciccs_are_finite(self, capsys):
        """CICc values should be finite numbers."""
        svc = PhyloPath(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        for m in payload["model_comparison"]:
            assert m["CICc"] > 0

    def test_weights_sum_to_one(self, capsys):
        """Model weights should sum to ~1."""
        svc = PhyloPath(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        total_w = sum(m["weight"] for m in payload["model_comparison"])
        assert total_w == pytest.approx(1.0, abs=0.01)

    def test_path_coefficients_present(self, capsys):
        """Path coefficients are in output."""
        svc = PhyloPath(_make_args(json=True))
        svc.run()
        captured = capsys.readouterr()
        for line in captured.out.strip().split("\n"):
            if line.startswith("{"):
                payload = json.loads(line)
                break
        coefs = payload["path_coefficients"]
        assert "body_mass -> brain_size" in coefs
        # body_mass->brain_size should be positive and strong
        assert coefs["body_mass -> brain_size"]["coef"] > 0.5

    def test_design_matrix_from_z_matches_column_stack_reference(self):
        n = 5
        z_data = {
            "x1": np.arange(n, dtype=float),
            "x2": np.arange(n, dtype=float) + 10.0,
        }

        observed = PhyloPath._design_matrix_from_z(z_data, ["x1", "x2"], n)
        expected = np.column_stack(
            [np.ones(n), z_data["x1"], z_data["x2"]]
        )

        np.testing.assert_allclose(observed, expected)

    def test_estimate_path_coefficients_preallocates_design_matrix(
        self, monkeypatch
    ):
        svc = PhyloPath(_make_args())
        n = 5
        z_data = {
            "x1": np.arange(n, dtype=float),
            "x2": np.arange(n, dtype=float) + 10.0,
            "y": np.arange(n, dtype=float) * 2.0,
        }
        captured = {}

        def fail_column_stack(*_args, **_kwargs):
            raise AssertionError("path coefficient setup should preallocate X")

        def fake_pgls_fit(y, X, _vcv, _n):
            captured["X"] = X.copy()
            return (
                np.array([0.0, 0.25, 0.5]),
                np.array([1.0, 0.1, 0.2]),
                1.0,
            )

        monkeypatch.setattr(phylo_path_module.np, "column_stack", fail_column_stack)
        monkeypatch.setattr(svc, "_pgls_fit", fake_pgls_fit)

        coefficients = svc._estimate_path_coefficients(
            [("x1", "y"), ("x2", "y")],
            ["x1", "x2", "y"],
            z_data,
            np.eye(n),
            n,
        )

        expected_X = np.empty((n, 3), dtype=np.float64)
        expected_X[:, 0] = 1.0
        expected_X[:, 1] = z_data["x1"]
        expected_X[:, 2] = z_data["x2"]
        np.testing.assert_allclose(captured["X"], expected_X)
        assert coefficients["x1 -> y"]["coef"] == pytest.approx(0.25)
        assert coefficients["x2 -> y"]["coef"] == pytest.approx(0.5)

    def test_fit_gls_from_vcv_cholesky_matches_inverse(self):
        rng = np.random.default_rng(20260623)
        svc = PhyloPath.__new__(PhyloPath)
        n = 16
        A = rng.normal(size=(n, n))
        vcv = A @ A.T + np.eye(n)
        X = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
        beta = np.array([0.4, -0.2, 1.3])
        y = X @ beta + rng.normal(size=n) * 0.1

        fast = svc._fit_gls_from_vcv(y, X, vcv)
        inverse = svc._fit_gls_from_vcv_inverse(y, X, vcv)

        for fast_value, inverse_value in zip(fast, inverse):
            assert fast_value == pytest.approx(inverse_value)

    def test_fit_gls_from_vcv_uses_single_cholesky_solve(self, monkeypatch):
        rng = np.random.default_rng(20260628)
        svc = PhyloPath.__new__(PhyloPath)
        n = 16
        A = rng.normal(size=(n, n))
        vcv = A @ A.T + np.eye(n)
        X = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
        beta = np.array([0.4, -0.2, 1.3])
        y = X @ beta + rng.normal(size=n) * 0.1
        original_cho_solve = phylo_path_module.cho_solve
        solve_calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal solve_calls
            solve_calls += 1
            return original_cho_solve(*args, **kwargs)

        monkeypatch.setattr(phylo_path_module, "cho_solve", counting_cho_solve)

        fast = svc._fit_gls_from_vcv(y, X, vcv)
        inverse = svc._fit_gls_from_vcv_inverse(y, X, vcv)

        assert solve_calls == 1
        for fast_value, inverse_value in zip(fast, inverse):
            assert fast_value == pytest.approx(inverse_value)

    def test_fit_gls_from_vcv_cholesky_avoids_normal_matrix_inverse(
        self, monkeypatch
    ):
        rng = np.random.default_rng(20260628)
        svc = PhyloPath.__new__(PhyloPath)
        n = 16
        A = rng.normal(size=(n, n))
        vcv = A @ A.T + np.eye(n)
        X = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
        beta = np.array([0.4, -0.2, 1.3])
        y = X @ beta + rng.normal(size=n) * 0.1

        def fail_inv(*_args, **_kwargs):
            raise AssertionError("fast Cholesky path should use solve")

        monkeypatch.setattr(phylo_path_module.np.linalg, "inv", fail_inv)

        beta_hat, residuals, sigma2, var_beta = svc._fit_gls_from_vcv(
            y, X, vcv
        )

        assert beta_hat.shape == (3,)
        assert residuals.shape == (n,)
        assert sigma2 > 0
        assert var_beta.shape == (3, 3)

    def test_csv_output(self, tmp_path):
        """CSV file is created."""
        csv_path = str(tmp_path / "path.csv")
        svc = PhyloPath(_make_args(csv=csv_path))
        svc.run()
        with open(csv_path) as f:
            content = f.read()
        assert "model,k,q" in content
        assert "path,coef,se" in content

    def test_plot_output(self, tmp_path):
        """Plot file is created."""
        plot_path = str(tmp_path / "dag.png")
        svc = PhyloPath(_make_args(plot_output=plot_path))
        svc.run()
        assert (tmp_path / "dag.png").exists()

    def test_plot_dag_batches_node_circles(self, tmp_path, monkeypatch):
        pytest.importorskip("matplotlib")
        import matplotlib.axes
        from matplotlib.collections import PatchCollection
        from matplotlib.patches import Circle

        plot_path = str(tmp_path / "dag.png")
        svc = PhyloPath(_make_args(plot_output=plot_path))
        patch_collections = []
        original_add_collection = matplotlib.axes.Axes.add_collection
        original_add_patch = matplotlib.axes.Axes.add_patch

        def tracking_add_collection(self, collection, *args, **kwargs):
            if isinstance(collection, PatchCollection):
                patch_collections.append(collection)
            return original_add_collection(self, collection, *args, **kwargs)

        def fail_circle_add_patch(self, patch, *args, **kwargs):
            if isinstance(patch, Circle):
                raise AssertionError("DAG node circles should be batched")
            return original_add_patch(self, patch, *args, **kwargs)

        monkeypatch.setattr(
            matplotlib.axes.Axes,
            "add_collection",
            tracking_add_collection,
        )
        monkeypatch.setattr(
            matplotlib.axes.Axes,
            "add_patch",
            fail_circle_add_patch,
        )

        svc._plot_dag({}, ["a", "b", "c", "d"], plot_path)

        assert (tmp_path / "dag.png").exists()
        assert len(patch_collections) == 1
        assert len(patch_collections[0].get_paths()) == 4

    def test_missing_models_file_raises(self):
        """Nonexistent models file raises error."""
        svc = PhyloPath(_make_args(models="nonexistent.txt"))
        with pytest.raises(SystemExit) as exc_info:
            svc.run()
        assert exc_info.value.code == 2

    def test_single_model_raises(self, tmp_path):
        """Only 1 model should raise error (need >= 2)."""
        m = tmp_path / "one_model.txt"
        m.write_text("only: body_mass->brain_size\n")
        svc = PhyloPath(_make_args(models=str(m)))
        with pytest.raises(SystemExit) as exc_info:
            svc.run()
        assert exc_info.value.code == 2
