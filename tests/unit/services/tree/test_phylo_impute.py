import builtins
import importlib
import json
import os
import subprocess
import sys
import tempfile

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.phylo_impute import PhyloImpute
import phykit.services.tree.phylo_impute as phylo_impute_module
from phykit.errors import PhykitUserError


here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
MISSING_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_traits_missing.tsv")
MULTI_TRAITS_FILE = str(SAMPLE_FILES / "tree_simple_multi_traits.tsv")


def test_module_import_does_not_import_scipy_linalg(monkeypatch):
    module_name = "phykit.services.tree.phylo_impute"
    previous = sys.modules.pop(module_name, None)
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "scipy.linalg" or name.startswith("scipy.linalg."):
            raise AssertionError(
                "phylo_impute module import should not import SciPy linalg"
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
import phykit.services.tree.phylo_impute as module
assert callable(module.print_json)
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "numpy" not in sys.modules
assert "scipy.linalg" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = phylo_impute_module._LazyNumpy()

    arange_attr = lazy_np.arange

    assert lazy_np.__dict__["arange"] is arange_attr
    assert lazy_np.arange is arange_attr
    assert lazy_np._module is not None


def test_repeated_cholesky_helpers_cache_scipy_linalg_imports(monkeypatch):
    previous_cho_factor = phylo_impute_module._CHO_FACTOR
    previous_cho_solve = phylo_impute_module._CHO_SOLVE
    phylo_impute_module._CHO_FACTOR = None
    phylo_impute_module._CHO_SOLVE = None
    original_import = builtins.__import__
    scipy_linalg_imports = 0

    def counting_import(name, globals=None, locals=None, fromlist=(), level=0):
        nonlocal scipy_linalg_imports
        if name == "scipy.linalg":
            scipy_linalg_imports += 1
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", counting_import)
    try:
        matrix = np.eye(4)
        rhs = np.ones((4, 2))

        factor = phylo_impute_module.cho_factor(
            matrix,
            lower=True,
            check_finite=False,
        )
        phylo_impute_module.cho_solve(factor, rhs, check_finite=False)
        first_call_imports = scipy_linalg_imports

        factor = phylo_impute_module.cho_factor(
            matrix,
            lower=True,
            check_finite=False,
        )
        phylo_impute_module.cho_solve(factor, rhs, check_finite=False)
    finally:
        phylo_impute_module._CHO_FACTOR = previous_cho_factor
        phylo_impute_module._CHO_SOLVE = previous_cho_solve

    assert first_call_imports > 0
    assert scipy_linalg_imports == first_call_imports


def _make_args(output_path, **overrides):
    defaults = dict(
        tree=TREE_SIMPLE,
        trait_data=MISSING_TRAITS_FILE,
        output=output_path,
        gene_trees=None,
        json=False,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


def _build_service(output_path, **overrides):
    args = _make_args(output_path, **overrides)
    return PhyloImpute(args)


def test_build_data_matrix_preserves_order_and_missing_values():
    trait_data = {
        "taxon_b": [3.5, np.nan, 7.0],
        "taxon_a": [1.0, 2.0, np.nan],
        "taxon_c": [np.nan, 4.25, 5.5],
    }

    matrix = PhyloImpute._build_data_matrix(
        ["taxon_a", "taxon_b", "taxon_c"],
        trait_data,
    )

    expected = np.array([
        [1.0, 2.0, np.nan],
        [3.5, np.nan, 7.0],
        [np.nan, 4.25, 5.5],
    ])
    assert np.allclose(matrix, expected, equal_nan=True)


def test_run_uses_unmodified_tree_read(mocker, tmp_path):
    out = str(tmp_path / "imputed.tsv")
    svc = _build_service(out)
    tree = object()
    trait_names = ["x", "y"]
    trait_data = {"A": [1.0, 2.0], "B": [3.0, 4.0]}

    read_unmodified = mocker.patch.object(
        svc, "read_tree_file_unmodified", return_value=tree
    )
    mocker.patch.object(
        svc,
        "read_tree_file",
        side_effect=AssertionError("run should use read_tree_file_unmodified"),
    )
    validate = mocker.patch.object(svc, "validate_tree")
    mocker.patch.object(svc, "get_tip_names_from_tree", return_value=["A", "B"])
    mocker.patch.object(
        svc,
        "_parse_trait_file_with_na",
        return_value=(trait_names, trait_data, []),
    )
    build_vcv = mocker.patch(
        "phykit.services.tree.vcv_utils.build_vcv_matrix",
        return_value=np.eye(2),
    )
    mocker.patch.object(
        svc,
        "_estimate_complete_case_stats",
        return_value=(np.zeros(2), np.eye(2)),
    )
    mocker.patch.object(svc, "_write_output_tsv")
    mocker.patch.object(svc, "_print_text")

    svc.run()

    read_unmodified.assert_called_once_with()
    validate.assert_called_once_with(
        tree,
        min_tips=3,
        require_branch_lengths=True,
        context="phylogenetic imputation",
    )
    assert build_vcv.call_args.args[0] is tree


def test_run_reuses_missing_mask_for_imputation_indices(mocker, tmp_path):
    out = str(tmp_path / "imputed.tsv")
    svc = _build_service(out)
    tree = object()
    trait_names = ["x", "y", "z"]
    trait_data = {
        "A": [1.0, 2.0, 3.0],
        "B": [4.0, np.nan, 6.0],
        "C": [7.0, 8.0, 9.0],
        "D": [np.nan, 11.0, np.nan],
    }
    calls = []

    mocker.patch.object(svc, "read_tree_file_unmodified", return_value=tree)
    mocker.patch.object(svc, "validate_tree")
    mocker.patch.object(svc, "get_tip_names_from_tree", return_value=list(trait_data))
    mocker.patch.object(
        svc,
        "_parse_trait_file_with_na",
        return_value=(trait_names, trait_data, []),
    )
    mocker.patch(
        "phykit.services.tree.vcv_utils.build_vcv_matrix",
        return_value=np.eye(4),
    )
    mocker.patch.object(
        svc,
        "_estimate_complete_case_stats",
        return_value=(np.zeros(3), np.eye(3)),
    )
    write_output = mocker.patch.object(svc, "_write_output_tsv")
    mocker.patch.object(svc, "_print_text")

    def fake_impute(Y, vcv, taxon_idx, observed_traits, missing_traits, a_hat, sigma):
        calls.append((taxon_idx, observed_traits.copy(), missing_traits.copy()))
        values = {int(j): float(100 + taxon_idx * 10 + j) for j in missing_traits}
        ses = {int(j): 0.5 for j in missing_traits}
        return values, ses

    mocker.patch.object(svc, "_impute_taxon", side_effect=fake_impute)

    svc.run()

    assert len(calls) == 2
    assert calls[0][0] == 1
    np.testing.assert_array_equal(calls[0][1], np.array([0, 2]))
    np.testing.assert_array_equal(calls[0][2], np.array([1]))
    assert calls[1][0] == 3
    np.testing.assert_array_equal(calls[1][1], np.array([1]))
    np.testing.assert_array_equal(calls[1][2], np.array([0, 2]))

    written_matrix = write_output.call_args.args[3]
    assert written_matrix[1, 1] == pytest.approx(111.0)
    assert written_matrix[3, 0] == pytest.approx(130.0)
    assert written_matrix[3, 2] == pytest.approx(132.0)


class TestPhyloImpute:
    def test_parse_na_values(self, tmp_path):
        """NA, ?, and empty values are all parsed as NaN."""
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        tree = svc.read_tree_file()
        tree_tips = svc.get_tip_names_from_tree(tree)
        trait_names, traits, missing_info = svc._parse_trait_file_with_na(
            MISSING_TRAITS_FILE, tree_tips
        )
        # bear has NA for brain_size, sea_lion has ? for longevity,
        # monkey has NA for longevity
        missing_taxa_traits = [(m["taxon"], m["trait"]) for m in missing_info]
        assert ("bear", "brain_size") in missing_taxa_traits
        assert ("sea_lion", "longevity") in missing_taxa_traits
        assert ("monkey", "longevity") in missing_taxa_traits

        # Check that the values are NaN
        assert np.isnan(traits["bear"][1])  # brain_size index 1
        assert np.isnan(traits["sea_lion"][2])  # longevity index 2
        assert np.isnan(traits["monkey"][2])  # longevity index 2

    def test_parse_trait_file_with_na_skips_comments_and_filters_shared_taxa(
        self, tmp_path, capsys
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "\n".join(
                [
                    "# ignored before header",
                    "taxon\tbody_mass\tbrain_size",
                    "A\t1.0\tNA",
                    "# ignored between rows",
                    "B\t2.0\t?",
                    "C\t3.0\t4.0",
                    "off_tree\t5.0\tNA",
                ]
            )
            + "\n"
        )
        svc = _build_service(str(tmp_path / "imputed.tsv"))

        trait_names, traits, missing_info = svc._parse_trait_file_with_na(
            str(trait_file), ["A", "B", "C", "missing"]
        )

        assert trait_names == ["body_mass", "brain_size"]
        assert set(traits) == {"A", "B", "C"}
        assert np.isnan(traits["A"][1])
        assert np.isnan(traits["B"][1])
        assert missing_info == [
            {"taxon": "A", "trait": "brain_size", "line": 2},
            {"taxon": "B", "trait": "brain_size", "line": 3},
        ]
        err = capsys.readouterr().err
        assert "1 taxa in tree but not in trait file: missing" in err
        assert "1 taxa in trait file but not in tree: off_tree" in err

    def test_parse_trait_file_validates_but_does_not_store_off_tree_rows(
        self, tmp_path, capsys
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "\n".join(
                [
                    "taxon\tbody_mass\tbrain_size",
                    "A\t1.0\tNA",
                    "B\t2.0\t3.0",
                    "C\t4.0\t5.0",
                    "off_tree\t6.0\t?",
                ]
            )
            + "\n"
        )
        svc = _build_service(str(tmp_path / "imputed.tsv"))

        trait_names, traits, missing_info = svc._parse_trait_file_with_na(
            str(trait_file), ["A", "B", "C"]
        )

        assert trait_names == ["body_mass", "brain_size"]
        assert set(traits) == {"A", "B", "C"}
        assert "off_tree" not in traits
        assert missing_info == [
            {"taxon": "A", "trait": "brain_size", "line": 2},
        ]
        assert all(entry["taxon"] != "off_tree" for entry in missing_info)
        err = capsys.readouterr().err
        assert "1 taxa in trait file but not in tree: off_tree" in err

    def test_parse_trait_file_all_numeric_rows_have_no_missing_metadata(
        self, tmp_path
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "\n".join(
                [
                    "taxon\tbody_mass\tbrain_size",
                    "A\t1.0\t2.0",
                    "B\t3.0\t4.0",
                    "C\t5.0\t6.0",
                ]
            )
            + "\n"
        )
        svc = _build_service(str(tmp_path / "imputed.tsv"))

        trait_names, traits, missing_info = svc._parse_trait_file_with_na(
            str(trait_file), ["A", "B", "C"]
        )

        assert trait_names == ["body_mass", "brain_size"]
        assert traits == {
            "A": [1.0, 2.0],
            "B": [3.0, 4.0],
            "C": [5.0, 6.0],
        }
        assert missing_info == []

    def test_parse_trait_file_all_shared_taxa_emit_no_warnings(self, tmp_path, capsys):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "\n".join(
                [
                    "taxon\tbody_mass\tbrain_size",
                    "A\t1.0\t2.0",
                    "B\t3.0\tNA",
                    "C\t5.0\t6.0",
                ]
            )
            + "\n"
        )
        svc = _build_service(str(tmp_path / "imputed.tsv"))

        trait_names, traits, missing_info = svc._parse_trait_file_with_na(
            str(trait_file), ["C", "A", "B"]
        )

        assert trait_names == ["body_mass", "brain_size"]
        assert set(traits) == {"A", "B", "C"}
        assert np.isnan(traits["B"][1])
        assert missing_info == [
            {"taxon": "B", "trait": "brain_size", "line": 3},
        ]
        assert capsys.readouterr().err == ""

    def test_parse_trait_file_still_validates_off_tree_non_numeric_values(
        self, tmp_path
    ):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "\n".join(
                [
                    "taxon\tbody_mass\tbrain_size",
                    "A\t1.0\t2.0",
                    "B\t2.0\t3.0",
                    "C\t4.0\t5.0",
                    "off_tree\tbad\t6.0",
                ]
            )
            + "\n"
        )
        svc = _build_service(str(tmp_path / "imputed.tsv"))

        with pytest.raises(PhykitUserError) as excinfo:
            svc._parse_trait_file_with_na(str(trait_file), ["A", "B", "C"])

        assert (
            "Non-numeric trait value 'bad' for taxon 'off_tree' "
            "(trait 'body_mass') on line 5."
        ) in excinfo.value.messages

    def test_parse_trait_file_with_na_preserves_logical_row_errors(self, tmp_path):
        trait_file = tmp_path / "traits.tsv"
        trait_file.write_text(
            "\n".join(
                [
                    "taxon\tbody_mass\tbrain_size",
                    "A\t1.0\t2.0",
                    "# ignored between rows",
                    "B\tbad\t4.0",
                ]
            )
            + "\n"
        )
        svc = _build_service(str(tmp_path / "imputed.tsv"))

        with pytest.raises(PhykitUserError) as excinfo:
            svc._parse_trait_file_with_na(str(trait_file), ["A", "B", "C"])

        assert (
            "Non-numeric trait value 'bad' for taxon 'B' "
            "(trait 'body_mass') on line 3."
        ) in excinfo.value.messages

    def test_subset_to_shared_taxa_reuses_shared_set(self, tmp_path):
        class CountingShared(list):
            def __init__(self, values):
                super().__init__(values)
                self.iterations = 0

            def __iter__(self):
                self.iterations += 1
                return super().__iter__()

        ordered_names = ["A", "B", "C", "D"]
        shared = CountingShared(["A", "C"])
        Y = np.arange(8).reshape(4, 2)
        svc = _build_service(str(tmp_path / "imputed.tsv"))

        subset_names, subset_Y = svc._subset_to_shared_taxa(
            ordered_names, Y, shared
        )

        assert shared.iterations == 1
        assert subset_names == ["A", "C"]
        assert subset_Y.tolist() == [[0, 1], [4, 5]]

    def test_subset_to_shared_taxa_returns_original_for_ordered_shared_taxa(self):
        class CountingShared(list):
            def __init__(self, values):
                super().__init__(values)
                self.iterations = 0

            def __iter__(self):
                self.iterations += 1
                return super().__iter__()

        ordered_names = ["A", "B", "C", "D"]
        shared = CountingShared(["A", "B", "C", "D"])
        Y = np.arange(8).reshape(4, 2)

        subset_names, subset_Y = PhyloImpute._subset_to_shared_taxa(
            ordered_names, Y, shared
        )

        assert shared.iterations == 0
        assert subset_names is ordered_names
        assert subset_Y is Y

    def test_imputed_values_finite(self, tmp_path):
        """All imputed values should be finite numbers."""
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        svc.run()

        # Read output and verify no NaN/inf
        with open(out) as f:
            lines = f.readlines()
        for line in lines[1:]:
            parts = line.strip().split("\t")
            for val_str in parts[1:]:
                val = float(val_str)
                assert np.isfinite(val), f"Non-finite value {val_str} in output"

    def test_se_positive(self, tmp_path, capsys):
        """All standard errors should be positive."""
        out = str(tmp_path / "imputed.tsv")
        args = _make_args(out, json=True)
        svc = PhyloImpute(args)
        svc.run()

        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        for entry in payload["imputed"]:
            assert entry["se"] > 0, (
                f"SE for {entry['taxon']}:{entry['trait']} = {entry['se']} <= 0"
            )

    def test_run_counts_complete_cases_with_count_nonzero(
        self, tmp_path, monkeypatch
    ):
        svc = PhyloImpute.__new__(PhyloImpute)
        svc.trait_data_path = "traits.tsv"
        svc.output_path = str(tmp_path / "imputed.tsv")
        svc.gene_trees_path = None
        svc.json_output = True

        tree = object()
        svc.read_tree_file_unmodified = lambda: tree
        svc.validate_tree = lambda *_args, **_kwargs: None
        svc.get_tip_names_from_tree = lambda _tree: ["a", "b", "c"]
        svc._parse_trait_file_with_na = lambda _path, _tips: (
            ["trait_1", "trait_2"],
            {
                "a": [1.0, 2.0],
                "b": [3.0, 4.0],
                "c": [5.0, float("nan")],
            },
            [{"taxon": "c", "trait": "trait_2", "line": 4}],
        )
        svc._estimate_complete_case_stats = lambda _Y, _C: (
            np.array([2.0, 3.0]),
            np.eye(2),
        )
        svc._impute_taxon = lambda *_args: (
            np.array([5.0, 6.0]),
            np.array([0.1, 0.2]),
        )
        svc._write_output_tsv = lambda *_args: None
        svc._print_json = lambda *_args: None

        monkeypatch.setattr(
            "phykit.services.tree.vcv_utils.build_vcv_matrix",
            lambda _tree, names: np.eye(len(names)),
        )

        observed_masks = []
        original_count_nonzero = np.count_nonzero

        def count_nonzero_spy(values, *args, **kwargs):
            if getattr(values, "shape", None) == (3,):
                observed_masks.append(values.copy())
            return original_count_nonzero(values, *args, **kwargs)

        def fail_sum(*_args, **_kwargs):
            raise AssertionError("complete-case mask should use count_nonzero")

        monkeypatch.setattr(
            phylo_impute_module.np,
            "count_nonzero",
            count_nonzero_spy,
            raising=False,
        )
        monkeypatch.setattr(phylo_impute_module.np, "sum", fail_sum, raising=False)

        svc.run()

        assert len(observed_masks) == 1
        assert observed_masks[0].tolist() == [True, True, False]

    def test_ci_contains_estimate(self, tmp_path, capsys):
        """CI lower < estimate < CI upper for all imputed values."""
        out = str(tmp_path / "imputed.tsv")
        args = _make_args(out, json=True)
        svc = PhyloImpute(args)
        svc.run()

        captured = capsys.readouterr()
        payload = json.loads(captured.out)
        for entry in payload["imputed"]:
            assert entry["ci_lower"] < entry["value"] < entry["ci_upper"], (
                f"CI [{entry['ci_lower']}, {entry['ci_upper']}] does not "
                f"contain estimate {entry['value']} for "
                f"{entry['taxon']}:{entry['trait']}"
            )

    def test_output_file_created(self, tmp_path):
        """Output TSV file should be created."""
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        svc.run()
        assert os.path.exists(out)

    def test_output_no_missing(self, tmp_path):
        """Output TSV should have no NA/missing values."""
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        svc.run()

        with open(out) as f:
            lines = f.readlines()

        missing_markers = {"NA", "na", "Na", "?", ""}
        for line_idx, line in enumerate(lines[1:], 2):
            parts = line.strip().split("\t")
            for i, val_str in enumerate(parts[1:]):
                assert val_str.strip() not in missing_markers, (
                    f"Missing value on line {line_idx}, column {i + 1}"
                )

    def test_complete_data_unchanged(self, tmp_path):
        """Taxa with no missing values should have unchanged values."""
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        svc.run()

        # Read original
        with open(MISSING_TRAITS_FILE) as f:
            orig_lines = f.readlines()
        orig_data = {}
        for line in orig_lines[1:]:
            parts = line.strip().split("\t")
            taxon = parts[0]
            vals = parts[1:]
            # Only track taxa with no missing values
            if all(v.strip() not in {"NA", "na", "Na", "?", ""} for v in vals):
                orig_data[taxon] = [float(v) for v in vals]

        # Read output
        with open(out) as f:
            out_lines = f.readlines()
        out_data = {}
        for line in out_lines[1:]:
            parts = line.strip().split("\t")
            out_data[parts[0]] = [float(v) for v in parts[1:]]

        for taxon, orig_vals in orig_data.items():
            for j, (orig_v, out_v) in enumerate(zip(orig_vals, out_data[taxon])):
                assert abs(orig_v - out_v) < 1e-4, (
                    f"Value changed for {taxon} trait {j}: "
                    f"{orig_v} -> {out_v}"
                )

    def test_json_output(self, tmp_path, capsys):
        """JSON output should have correct structure."""
        out = str(tmp_path / "imputed.tsv")
        args = _make_args(out, json=True)
        svc = PhyloImpute(args)
        svc.run()

        captured = capsys.readouterr()
        payload = json.loads(captured.out)

        assert "n_taxa" in payload
        assert "n_traits" in payload
        assert "trait_names" in payload
        assert "n_missing" in payload
        assert "vcv_type" in payload
        assert "imputed" in payload
        assert "output_file" in payload
        assert payload["n_traits"] == 3
        assert payload["n_missing"] == 3
        assert payload["vcv_type"] == "BM"
        assert len(payload["imputed"]) == 3

        for entry in payload["imputed"]:
            assert "taxon" in entry
            assert "trait" in entry
            assert "value" in entry
            assert "se" in entry
            assert "ci_lower" in entry
            assert "ci_upper" in entry

    def test_complete_case_stats_cholesky_matches_inverse(self, tmp_path):
        rng = np.random.default_rng(20260623)
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        n = 16
        p = 4
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        Y = rng.normal(size=(n, p))

        fast_mean, fast_cov = svc._estimate_complete_case_stats_cholesky(Y, C)
        inverse_mean, inverse_cov = svc._estimate_complete_case_stats_inverse(Y, C)

        assert fast_mean == pytest.approx(inverse_mean)
        assert fast_cov == pytest.approx(inverse_cov)

    def test_complete_case_stats_cholesky_uses_single_solve(
        self, tmp_path, monkeypatch
    ):
        import phykit.services.tree.phylo_impute as phylo_impute_module

        rng = np.random.default_rng(20260628)
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        n = 16
        p = 4
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        Y = rng.normal(size=(n, p))
        original_cho_solve = phylo_impute_module.cho_solve
        solve_calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal solve_calls
            solve_calls += 1
            return original_cho_solve(*args, **kwargs)

        monkeypatch.setattr(phylo_impute_module, "cho_solve", counting_cho_solve)

        fast_mean, fast_cov = svc._estimate_complete_case_stats_cholesky(Y, C)
        inverse_mean, inverse_cov = svc._estimate_complete_case_stats_inverse(Y, C)

        assert solve_calls == 1
        assert fast_mean == pytest.approx(inverse_mean)
        assert fast_cov == pytest.approx(inverse_cov)

    def test_complete_case_stats_inverse_matches_scalar_reference(self, tmp_path):
        rng = np.random.default_rng(20260624)
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        n = 14
        p = 5
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        Y = rng.normal(size=(n, p))

        observed_mean, observed_cov = svc._estimate_complete_case_stats_inverse(Y, C)

        C_inv = np.linalg.inv(C)
        ones = np.ones(n)
        denom = ones @ C_inv @ ones
        expected_mean = np.array([(ones @ C_inv @ Y[:, j]) / denom for j in range(p)])
        residuals = Y - expected_mean
        expected_cov = (residuals.T @ C_inv @ residuals) / (n - 1)

        assert observed_mean == pytest.approx(expected_mean)
        assert observed_cov == pytest.approx(expected_cov)

    def test_impute_taxon_cholesky_matches_inverse(self, tmp_path):
        rng = np.random.default_rng(20260623)
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        n = 16
        p = 7
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        B = rng.normal(size=(p, p))
        Sigma = B @ B.T + np.eye(p)
        Y = rng.normal(size=(n, p))
        taxon_idx = 3
        missing_traits = [1, 4, 6]
        observed_traits = [j for j in range(p) if j not in missing_traits]
        Y[taxon_idx, missing_traits] = np.nan
        Y[[1, 5], 1] = np.nan
        a_hat = rng.normal(size=p)

        fast_values, fast_ses = svc._impute_taxon_cholesky(
            Y, C, taxon_idx, observed_traits, missing_traits, a_hat, Sigma
        )
        inverse_values, inverse_ses = svc._impute_taxon_inverse(
            Y, C, taxon_idx, observed_traits, missing_traits, a_hat, Sigma
        )

        assert fast_values == pytest.approx(inverse_values)
        assert fast_ses == pytest.approx(inverse_ses)

    def test_other_taxon_indices_vectorized_order(self):
        np.testing.assert_array_equal(
            PhyloImpute._other_taxon_indices(5, 0),
            np.array([1, 2, 3, 4]),
        )
        np.testing.assert_array_equal(
            PhyloImpute._other_taxon_indices(5, 2),
            np.array([0, 1, 3, 4]),
        )
        np.testing.assert_array_equal(
            PhyloImpute._other_taxon_indices(5, 4),
            np.array([0, 1, 2, 3]),
        )

    def test_impute_taxon_paths_use_other_taxon_index_helper(
        self, tmp_path, monkeypatch
    ):
        rng = np.random.default_rng(20260630)
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        n = 8
        p = 4
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        B = rng.normal(size=(p, p))
        Sigma = B @ B.T + np.eye(p)
        Y = rng.normal(size=(n, p))
        taxon_idx = 3
        missing_traits = [1]
        observed_traits = [0, 2, 3]
        Y[taxon_idx, missing_traits] = np.nan
        a_hat = rng.normal(size=p)
        original = PhyloImpute._other_taxon_indices
        calls = []

        def tracking_other_indices(n_taxa, idx):
            calls.append((n_taxa, idx))
            return original(n_taxa, idx)

        monkeypatch.setattr(
            PhyloImpute,
            "_other_taxon_indices",
            staticmethod(tracking_other_indices),
        )

        svc._impute_taxon_cholesky(
            Y, C, taxon_idx, observed_traits, missing_traits, a_hat, Sigma
        )
        svc._impute_taxon_inverse(
            Y, C, taxon_idx, observed_traits, missing_traits, a_hat, Sigma
        )

        assert calls == [(n, taxon_idx), (n, taxon_idx)]

    def test_impute_taxon_cholesky_combines_phylo_rhs(self, tmp_path, monkeypatch):
        import phykit.services.tree.phylo_impute as phylo_impute_module

        rng = np.random.default_rng(20260628)
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        n = 16
        p = 7
        A = rng.normal(size=(n, n))
        C = A @ A.T + np.eye(n)
        B = rng.normal(size=(p, p))
        Sigma = B @ B.T + np.eye(p)
        Y = rng.normal(size=(n, p))
        taxon_idx = 3
        missing_traits = [1, 4, 6]
        observed_traits = [j for j in range(p) if j not in missing_traits]
        Y[taxon_idx, missing_traits] = np.nan
        Y[[1, 5], 1] = np.nan
        Y[[2, 8], 4] = np.nan
        Y[[6, 11], 6] = np.nan
        a_hat = rng.normal(size=p)
        original_cho_solve = phylo_impute_module.cho_solve
        solve_calls = 0

        def counting_cho_solve(*args, **kwargs):
            nonlocal solve_calls
            solve_calls += 1
            return original_cho_solve(*args, **kwargs)

        monkeypatch.setattr(phylo_impute_module, "cho_solve", counting_cho_solve)

        fast_values, fast_ses = svc._impute_taxon_cholesky(
            Y, C, taxon_idx, observed_traits, missing_traits, a_hat, Sigma
        )
        inverse_values, inverse_ses = svc._impute_taxon_inverse(
            Y, C, taxon_idx, observed_traits, missing_traits, a_hat, Sigma
        )

        assert solve_calls == 7
        assert fast_values == pytest.approx(inverse_values)
        assert fast_ses == pytest.approx(inverse_ses)

    def test_single_trait_works(self, tmp_path):
        """Single-trait imputation should work."""
        # Create a single-trait file with a missing value
        single_trait = tmp_path / "single_trait.tsv"
        single_trait.write_text(
            "taxon\tbody_mass\n"
            "raccoon\t1.04\n"
            "bear\tNA\n"
            "sea_lion\t2.30\n"
            "seal\t1.88\n"
            "monkey\t0.60\n"
            "cat\t0.56\n"
            "weasel\t-0.30\n"
            "dog\t1.18\n"
        )

        out = str(tmp_path / "imputed_single.tsv")
        args = _make_args(out, trait_data=str(single_trait))
        svc = PhyloImpute(args)
        svc.run()

        assert os.path.exists(out)
        with open(out) as f:
            lines = f.readlines()
        # Should have 9 lines (header + 8 taxa)
        assert len(lines) == 9

        # All values should be finite
        for line in lines[1:]:
            parts = line.strip().split("\t")
            val = float(parts[1])
            assert np.isfinite(val)

    def test_all_missing_one_trait(self, tmp_path):
        """If all taxa are missing one trait, use phylogenetic mean."""
        all_missing = tmp_path / "all_missing.tsv"
        all_missing.write_text(
            "taxon\tbody_mass\tnew_trait\n"
            "raccoon\t1.04\tNA\n"
            "bear\t2.39\tNA\n"
            "sea_lion\t2.30\tNA\n"
            "seal\t1.88\tNA\n"
            "monkey\t0.60\tNA\n"
            "cat\t0.56\tNA\n"
            "weasel\t-0.30\tNA\n"
            "dog\t1.18\tNA\n"
        )

        out = str(tmp_path / "imputed_all_missing.tsv")
        args = _make_args(out, trait_data=str(all_missing))
        svc = PhyloImpute(args)

        # This should fail because no taxa have complete data for all traits
        with pytest.raises(PhykitUserError):
            svc.run()

    def test_text_output_format(self, tmp_path, capsys):
        """Text output should contain expected fields."""
        out = str(tmp_path / "imputed.tsv")
        args = _make_args(out, json=False)
        svc = PhyloImpute(args)
        svc.run()

        captured = capsys.readouterr()
        text = captured.out
        assert "Phylogenetic Imputation" in text
        assert "Taxa:" in text
        assert "Traits:" in text
        assert "Missing values:" in text
        assert "VCV:" in text
        assert "Imputed values:" in text
        assert "Output:" in text

    def test_write_output_tsv_preserves_format(self, tmp_path):
        out = tmp_path / "imputed.tsv"
        Y = np.array([
            [1.2345678, 2.0, 3.5],
            [4.0, 5.6789012, 6.25],
        ])

        PhyloImpute._write_output_tsv(
            str(out),
            ["height", "weight", "depth"],
            ["taxon_a", "taxon_b"],
            Y,
        )

        assert out.read_text() == (
            "taxon\theight\tweight\tdepth\n"
            "taxon_a\t1.234568\t2.000000\t3.500000\n"
            "taxon_b\t4.000000\t5.678901\t6.250000\n"
        )

    def test_write_output_tsv_uses_chunked_savetxt(self, tmp_path, monkeypatch):
        out = tmp_path / "imputed.tsv"
        Y = np.arange(1001 * 2, dtype=float).reshape(1001, 2) / 10.0
        ordered_names = [f"taxon_{i}" for i in range(Y.shape[0])]
        chunk_shapes = []
        real_savetxt = np.savetxt

        def tracking_savetxt(handle, block, *args, **kwargs):
            chunk_shapes.append(block.shape)
            return real_savetxt(handle, block, *args, **kwargs)

        monkeypatch.setattr(phylo_impute_module.np, "savetxt", tracking_savetxt)

        PhyloImpute._write_output_tsv(
            str(out),
            ["height", "weight"],
            ordered_names,
            Y,
        )

        assert chunk_shapes == [(1000, 3), (1, 3)]
        lines = out.read_text().splitlines()
        assert lines[0] == "taxon\theight\tweight"
        assert lines[1] == "taxon_0\t0.000000\t0.100000"
        assert lines[-1] == "taxon_1000\t200.000000\t200.100000"

    def test_print_text_batches_imputed_rows(self, tmp_path, mocker):
        out = str(tmp_path / "imputed.tsv")
        svc = _build_service(out)
        printed = mocker.patch("builtins.print")

        svc._print_text(
            4,
            2,
            ["height", "weight"],
            "BM",
            [
                {
                    "taxon": "A",
                    "trait": "height",
                    "value": 1.23456,
                    "se": 0.2,
                    "ci_lower": 0.8,
                    "ci_upper": 1.6,
                },
                {
                    "taxon": "LongTaxon",
                    "trait": "weight",
                    "value": 2.0,
                    "se": 0.33333,
                    "ci_lower": 1.25,
                    "ci_upper": 2.75,
                },
            ],
            "traits.tsv",
            "imputed.tsv",
        )

        printed.assert_called_once_with(
            "Phylogenetic Imputation\n"
            "Data: traits.tsv\n"
            "Taxa: 4\n"
            "Traits: 2 (height, weight)\n"
            "Missing values: 2\n"
            "VCV: BM (standard)\n"
            "\n"
            "Imputed values:\n"
            "  Taxon      Trait      Imputed          SE                95% CI\n"
            "  A            height      1.2346      0.2000      [0.8000, 1.6000]\n"
            "  LongTaxon    weight      2.0000      0.3333      [1.2500, 2.7500]\n"
            "\n"
            "Output: imputed.tsv"
        )
