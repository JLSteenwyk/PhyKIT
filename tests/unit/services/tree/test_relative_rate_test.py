import pytest
import numpy as np
import subprocess
import sys
from argparse import Namespace
from unittest.mock import call
import phykit.services.tree.relative_rate_test as relative_rate_test_module
from phykit.services.tree.relative_rate_test import RelativeRateTest


def _make_args(**overrides):
    defaults = dict(
        tree="tree.tre",
        alignment="alignment.fa",
        alignment_list=None,
        verbose=False,
        json=False,
        plot_output=None,
        fig_width=4,
        fig_height=4,
        dpi=80,
        no_title=True,
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


def test_module_import_does_not_import_numpy():
    code = """
import builtins
import sys
module_name = "phykit.services.tree.relative_rate_test"
sys.modules.pop(module_name, None)
sys.modules.pop("pathlib", None)
original_import = builtins.__import__

def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
    if name == "pathlib" or name.startswith("pathlib."):
        raise AssertionError("relative_rate_test import should not import pathlib")
    return original_import(name, globals, locals, fromlist, level)

builtins.__import__ = guarded_import
try:
    import phykit.services.tree.relative_rate_test as module
finally:
    builtins.__import__ = original_import

assert hasattr(module.np, "__getattr__")
assert callable(module.print_json)
assert callable(module.get_alignment_and_format_helper)
assert callable(module.Path)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "numpy" not in sys.modules
assert "matplotlib" not in sys.modules
assert "scipy" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.files" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = relative_rate_test_module._LazyNumpy()

    frombuffer_attr = lazy_np.frombuffer

    assert lazy_np.__dict__["frombuffer"] is frombuffer_attr
    assert lazy_np.frombuffer is frombuffer_attr
    assert lazy_np._module is not None


class TestHeatmapMatrices:
    def test_build_heatmap_matrices_are_symmetric(self):
        results = [
            {"taxon1": "B", "taxon2": "A", "p_fdr": 0.01},
            {"taxon1": "C", "taxon2": "A", "p_fdr": 0.2},
            {"taxon1": "B", "taxon2": "C", "p_fdr": 0.0},
        ]

        taxa, matrix, p_fdr_matrix = RelativeRateTest._build_heatmap_matrices(results)
        idx = {taxon: i for i, taxon in enumerate(taxa)}

        assert taxa == ["A", "B", "C"]
        assert matrix[idx["A"], idx["B"]] == pytest.approx(2.0)
        assert matrix[idx["B"], idx["A"]] == pytest.approx(2.0)
        assert matrix[idx["B"], idx["C"]] == 10.0
        assert p_fdr_matrix[idx["A"], idx["B"]] == 0.01
        assert p_fdr_matrix[idx["B"], idx["A"]] == 0.01
        assert np.isnan(matrix[idx["A"], idx["A"]])
        assert np.isnan(p_fdr_matrix[idx["A"], idx["A"]])

    def test_significant_heatmap_cells_no_hits_skips_where(self, monkeypatch):
        p_fdr_matrix = np.ones((4, 4))
        np.fill_diagonal(p_fdr_matrix, np.nan)

        def fail_where(*args, **kwargs):
            raise AssertionError("no significant cells should skip np.where")

        monkeypatch.setattr(
            relative_rate_test_module.np, "where", fail_where, raising=False
        )

        sig_rows, sig_cols = RelativeRateTest._significant_heatmap_cells(p_fdr_matrix)

        assert sig_rows.size == 0
        assert sig_cols.size == 0

    def test_significant_heatmap_cells_sparse_hits_use_flat_indices(self, monkeypatch):
        p_fdr_matrix = np.ones((4, 4))
        np.fill_diagonal(p_fdr_matrix, np.nan)
        p_fdr_matrix[0, 2] = 0.01
        p_fdr_matrix[1, 3] = 0.02
        p_fdr_matrix[3, 0] = 0.03

        def fail_where(*args, **kwargs):
            raise AssertionError("sparse significant cells should avoid np.where")

        monkeypatch.setattr(
            relative_rate_test_module.np, "where", fail_where, raising=False
        )

        sig_rows, sig_cols = RelativeRateTest._significant_heatmap_cells(p_fdr_matrix)

        np.testing.assert_array_equal(sig_rows, np.array([0, 1, 3]))
        np.testing.assert_array_equal(sig_cols, np.array([2, 3, 0]))

    def test_plot_heatmap_batches_significance_markers(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        import matplotlib.axes

        svc = RelativeRateTest(_make_args())
        output_path = tmp_path / "rrt_heatmap.png"
        results = [
            {"taxon1": "A", "taxon2": "B", "p_fdr": 0.01},
            {"taxon1": "A", "taxon2": "C", "p_fdr": 0.50},
            {"taxon1": "B", "taxon2": "C", "p_fdr": 0.02},
        ]
        original_text = matplotlib.axes.Axes.text
        original_scatter = matplotlib.axes.Axes.scatter
        scatter_sizes = []

        def fail_star_text(self, *args, **kwargs):
            if len(args) >= 3 and args[2] == "*":
                raise AssertionError("significance markers should be batched")
            return original_text(self, *args, **kwargs)

        def count_scatter(self, x, y, *args, **kwargs):
            if kwargs.get("marker") == "$*$":
                scatter_sizes.append(len(x))
            return original_scatter(self, x, y, *args, **kwargs)

        monkeypatch.setattr(matplotlib.axes.Axes, "text", fail_star_text)
        monkeypatch.setattr(matplotlib.axes.Axes, "scatter", count_scatter)

        svc._plot_heatmap(results, str(output_path))

        assert scatter_sizes == [4]
        assert output_path.exists()

    def test_plot_heatmap_skips_redundant_tight_layout(self, monkeypatch, tmp_path):
        pytest.importorskip("matplotlib")
        from matplotlib.figure import Figure

        svc = RelativeRateTest(_make_args())
        output_path = tmp_path / "rrt_no_tight_layout.png"
        results = [
            {"taxon1": "A", "taxon2": "B", "p_fdr": 0.01},
            {"taxon1": "A", "taxon2": "C", "p_fdr": 0.50},
            {"taxon1": "B", "taxon2": "C", "p_fdr": 0.02},
        ]

        def fail_tight_layout(self, *args, **kwargs):
            raise AssertionError("bbox_inches='tight' handles saved bounds")

        monkeypatch.setattr(Figure, "tight_layout", fail_tight_layout)
        svc._plot_heatmap(results, str(output_path))

        assert output_path.exists()
        assert output_path.stat().st_size > 0


class TestRun:
    def test_run_uses_unmodified_tree_read(self, mocker):
        args = Namespace(
            tree="tree.tre",
            alignment="alignment.fa",
            alignment_list=None,
            verbose=False,
            json=False,
            plot_output=None,
        )
        svc = RelativeRateTest(args)
        tree = object()
        results = [{"taxon1": "A", "taxon2": "B", "p_fdr": 1.0}]

        read_unmodified = mocker.patch.object(
            svc, "read_tree_file_unmodified", return_value=tree
        )
        mocker.patch.object(
            svc,
            "read_tree_file",
            side_effect=AssertionError("run should use read_tree_file_unmodified"),
        )
        identify = mocker.patch.object(svc, "_identify_outgroup", return_value="out")
        run_single = mocker.patch.object(svc, "_run_single", return_value=results)
        output_single = mocker.patch.object(svc, "_output_single")

        svc.run()

        read_unmodified.assert_called_once_with()
        identify.assert_called_once_with(tree)
        run_single.assert_called_once_with("alignment.fa", tree, "out")
        output_single.assert_called_once_with("out", results)

    def test_batch_list_skips_comments_blanks_and_resolves_relative_paths(
        self, tmp_path, mocker
    ):
        list_path = tmp_path / "alignments.txt"
        list_path.write_text(
            "   # ignored\n\n  first.fa  \n\t# also ignored\n  nested/second.fa\n"
        )
        args = _make_args(alignment=None, alignment_list=str(list_path))
        svc = RelativeRateTest(args)
        tree = object()
        results = [{"taxon1": "A", "taxon2": "B", "p_fdr": 1.0}]

        mocker.patch.object(svc, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(svc, "_identify_outgroup", return_value="out")
        run_single = mocker.patch.object(svc, "_run_single", return_value=results)
        output_batch = mocker.patch.object(svc, "_output_batch")

        svc.run()

        assert run_single.call_args_list == [
            call(str(tmp_path / "first.fa"), tree, "out"),
            call(str(tmp_path / "nested/second.fa"), tree, "out"),
        ]
        output_batch.assert_called_once_with(
            "out", 2, {("A", "B"): [results[0], results[0]]}
        )

    def test_batch_list_resolution_avoids_per_row_path_objects(
        self, tmp_path, mocker, monkeypatch
    ):
        import phykit.services.tree.relative_rate_test as module
        from pathlib import Path

        list_path = tmp_path / "alignments.txt"
        list_path.write_text("first.fa\nnested/second.fa\n")
        path_constructs = 0
        absolute_checks = 0
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
                nonlocal path_constructs
                path_constructs += 1
                self._path = Path(path)

            @property
            def parent(self):
                return CountingParent(self._path.parent)

            def exists(self):
                return self._path.exists()

            def open(self, *args, **kwargs):
                return self._path.open(*args, **kwargs)

            def is_absolute(self):
                nonlocal absolute_checks
                absolute_checks += 1
                return self._path.is_absolute()

        monkeypatch.setattr(module, "Path", CountingPath)
        args = _make_args(alignment=None, alignment_list=str(list_path))
        svc = RelativeRateTest(args)
        tree = object()
        results = [{"taxon1": "A", "taxon2": "B", "p_fdr": 1.0}]

        mocker.patch.object(svc, "read_tree_file_unmodified", return_value=tree)
        mocker.patch.object(svc, "_identify_outgroup", return_value="out")
        run_single = mocker.patch.object(svc, "_run_single", return_value=results)
        mocker.patch.object(svc, "_output_batch")

        svc.run()

        assert run_single.call_args_list == [
            call(str(tmp_path / "first.fa"), tree, "out"),
            call(str(tmp_path / "nested/second.fa"), tree, "out"),
        ]
        assert path_constructs == 1
        assert absolute_checks == 0
        assert parent_joins == 0


class TestPairwiseTajima:
    def test_vectorized_pairwise_matches_legacy(self):
        svc = RelativeRateTest.__new__(RelativeRateTest)
        ingroup = ["A", "B", "C", "D"]
        seq_dict = {
            "A": "ACGT-NX?ACGT",
            "B": "ACTTACGTACGT",
            "C": "TCGTACGTACGA",
            "D": "NNNN----ACGT",
        }
        seq_out = "ACGTACGTACGT"

        fast = svc._run_pairwise_tests_vectorized(seq_dict, ingroup, seq_out)
        legacy = svc._run_pairwise_tests_legacy(seq_dict, ingroup, seq_out)

        assert fast == legacy

    def test_vectorized_pairwise_preserves_pair_stat_order(self):
        svc = RelativeRateTest.__new__(RelativeRateTest)
        ingroup = ["A", "B", "C"]
        seq_dict = {
            "A": "AAAACC",
            "B": "AACCCC",
            "C": "CCCCAA",
        }
        seq_out = "AAAAAA"

        results = svc._run_pairwise_tests_vectorized(seq_dict, ingroup, seq_out)

        assert results == [
            {
                "m1": 0,
                "m2": 2,
                "chi2": 2.0,
                "p_value": pytest.approx(0.15729920705028516),
                "taxon1": "A",
                "taxon2": "B",
            },
            {
                "m1": 2,
                "m2": 4,
                "chi2": pytest.approx(0.6666666666666666),
                "p_value": pytest.approx(0.4142161782425252),
                "taxon1": "A",
                "taxon2": "C",
            },
            {
                "m1": 2,
                "m2": 2,
                "chi2": 0.0,
                "p_value": 1.0,
                "taxon1": "B",
                "taxon2": "C",
            },
        ]

    def test_large_vectorized_pairwise_uses_df1_erfc_pvalues(self, monkeypatch):
        svc = RelativeRateTest.__new__(RelativeRateTest)
        n_taxa = 46
        ingroup = [f"t{idx}" for idx in range(n_taxa)]
        seq_out = "A" * n_taxa
        seq_dict = {}
        for idx, taxon in enumerate(ingroup):
            seq = ["A"] * n_taxa
            seq[idx] = "C"
            seq_dict[taxon] = "".join(seq)

        real_import = __import__

        def fail_stats_import(name, *args, **kwargs):
            if name == "scipy.stats" or name.startswith("scipy.stats."):
                raise AssertionError("df=1 vector p-values should avoid scipy.stats")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr("builtins.__import__", fail_stats_import)
        monkeypatch.setattr(
            relative_rate_test_module.np,
            "triu_indices",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("large pairwise path should scan rows directly")
            ),
        )

        results = svc._run_pairwise_tests_vectorized(seq_dict, ingroup, seq_out)

        assert len(results) == 1035
        assert all(row["m1"] == 1 for row in results)
        assert all(row["m2"] == 1 for row in results)
        assert all(row["chi2"] == 0.0 for row in results)
        assert all(row["p_value"] == 1.0 for row in results)

    def test_large_pairwise_results_preserve_no_informative_rows(self):
        matrix = np.array(
            [
                [0, 0, 0, 0],
                [0, 0, 2, 0],
                [0, 1, 0, 3],
                [0, 0, 4, 0],
            ],
            dtype=np.int32,
        )
        rows = RelativeRateTest._pairwise_results_large(
            matrix,
            ["A", "B", "C", "D"],
        )

        assert rows[:3] == [
            {
                "m1": 0,
                "m2": 0,
                "chi2": 0.0,
                "p_value": 1.0,
                "taxon1": "A",
                "taxon2": "B",
            },
            {
                "m1": 0,
                "m2": 0,
                "chi2": 0.0,
                "p_value": 1.0,
                "taxon1": "A",
                "taxon2": "C",
            },
            {
                "m1": 0,
                "m2": 0,
                "chi2": 0.0,
                "p_value": 1.0,
                "taxon1": "A",
                "taxon2": "D",
            },
        ]
        assert rows[3]["m1"] == 2
        assert rows[3]["m2"] == 1
        assert rows[3]["chi2"] == pytest.approx(1 / 3)

    def test_ingroup_ascii_matrix_matches_vstack_reference(self):
        seq_dict = {
            "A": "acgt",
            "B": "t?ca",
            "C": "nn-x",
        }
        ingroup = ["A", "B", "C"]

        observed = RelativeRateTest._ingroup_ascii_matrix(seq_dict, ingroup, 4)
        expected = relative_rate_test_module.np.vstack(
            [
                relative_rate_test_module.np.frombuffer(
                    seq_dict[taxon].upper().encode("ascii"),
                    dtype=relative_rate_test_module.np.uint8,
                )
                for taxon in ingroup
            ]
        )

        relative_rate_test_module.np.testing.assert_array_equal(observed, expected)

    def test_vectorized_pairwise_builds_ingroup_matrix_without_vstack(
        self, monkeypatch
    ):
        svc = RelativeRateTest.__new__(RelativeRateTest)
        ingroup = ["A", "B", "C"]
        seq_dict = {
            "A": "ACGT",
            "B": "ACTT",
            "C": "TCGA",
        }
        seq_out = "ACGT"

        def fail_vstack(*_args, **_kwargs):
            raise AssertionError("vectorized setup should use joined byte buffer")

        monkeypatch.setattr(relative_rate_test_module.np, "vstack", fail_vstack)

        results = svc._run_pairwise_tests_vectorized(seq_dict, ingroup, seq_out)

        assert [(row["taxon1"], row["taxon2"]) for row in results] == [
            ("A", "B"),
            ("A", "C"),
            ("B", "C"),
        ]

    def test_multiple_testing_corrections_apply_to_pairwise_results(self):
        svc = RelativeRateTest.__new__(RelativeRateTest)
        ingroup = ["A", "B", "C"]
        seq_dict = {
            "A": "ACGTACGT",
            "B": "ACTTACGT",
            "C": "TCGTACGA",
        }
        seq_out = "ACGTACGT"

        results = svc._run_pairwise_tests(seq_dict, ingroup, seq_out)
        svc._add_multiple_testing_corrections(results)

        assert len(results) == 3
        assert all("p_bonf" in row for row in results)
        assert all("p_fdr" in row for row in results)

    def test_large_pairwise_builtin_rows_preserve_scalar_branch_results(
        self, monkeypatch
    ):
        ingroup = ["A", "B", "C", "D"]
        m1_matrix = np.array(
            [
                [0, 4, 2, 8],
                [3, 0, 5, 1],
                [7, 6, 0, 9],
                [2, 4, 3, 0],
            ],
            dtype=np.int32,
        )

        monkeypatch.setattr(
            relative_rate_test_module,
            "_PAIRWISE_BUILTIN_ROW_MIN_TAXA",
            999,
        )
        scalar_rows = RelativeRateTest._pairwise_results_large(m1_matrix, ingroup)

        monkeypatch.setattr(
            relative_rate_test_module,
            "_PAIRWISE_BUILTIN_ROW_MIN_TAXA",
            0,
        )
        builtin_rows = RelativeRateTest._pairwise_results_large(m1_matrix, ingroup)

        assert builtin_rows == scalar_rows
        assert [(row["taxon1"], row["taxon2"]) for row in builtin_rows] == [
            ("A", "B"),
            ("A", "C"),
            ("A", "D"),
            ("B", "C"),
            ("B", "D"),
            ("C", "D"),
        ]

    def test_multiple_testing_corrections_preserve_result_order(self, monkeypatch):
        svc = RelativeRateTest.__new__(RelativeRateTest)
        results = [
            {"taxon1": "A", "taxon2": "B", "p_value": 0.2},
            {"taxon1": "A", "taxon2": "C", "p_value": 0.01},
            {"taxon1": "B", "taxon2": "C", "p_value": 0.5},
        ]

        monkeypatch.setattr(
            svc,
            "_bonferroni",
            lambda p_values: np.array([0.6, 0.03, 1.0]),
        )
        monkeypatch.setattr(
            svc,
            "_fdr",
            lambda p_values: np.array([0.3, 0.03, 0.5]),
        )

        svc._add_multiple_testing_corrections(results)

        assert [
            (row["taxon1"], row["taxon2"], row["p_bonf"], row["p_fdr"])
            for row in results
        ] == [
            ("A", "B", 0.6, 0.3),
            ("A", "C", 0.03, 0.03),
            ("B", "C", 1.0, 0.5),
        ]


class TestTextOutput:
    def test_output_single_json_builds_rows_in_one_pass(self, mocker):
        class CountingResults(list):
            iter_calls = 0

            def __iter__(self):
                self.iter_calls += 1
                return super().__iter__()

        svc = RelativeRateTest(_make_args(json=True))
        print_json = mocker.patch.object(relative_rate_test_module, "print_json")
        results = CountingResults([
            {
                "taxon1": "A",
                "taxon2": "B",
                "m1": 1,
                "m2": 3,
                "chi2": 2.12345,
                "p_value": 0.0455,
                "p_bonf": 0.091,
                "p_fdr": 0.04,
            },
            {
                "taxon1": "A",
                "taxon2": "C",
                "m1": 2,
                "m2": 2,
                "chi2": 0.0,
                "p_value": 1.0,
                "p_bonf": 1.0,
                "p_fdr": 1.0,
            },
        ])

        svc._output_single("O", results)

        assert results.iter_calls == 1
        print_json.assert_called_once_with({
            "outgroup": "O",
            "n_ingroup_taxa": 3,
            "n_tests": 2,
            "results": [
                {
                    "taxon1": "A",
                    "taxon2": "B",
                    "m1": 1,
                    "m2": 3,
                    "chi2": 2.1235,
                    "p_value": 0.0455,
                    "p_bonferroni": 0.091,
                    "p_fdr": 0.04,
                },
                {
                    "taxon1": "A",
                    "taxon2": "C",
                    "m1": 2,
                    "m2": 2,
                    "chi2": 0.0,
                    "p_value": 1.0,
                    "p_bonferroni": 1.0,
                    "p_fdr": 1.0,
                },
            ],
        })

    def test_output_single_batches_text_rows(self, mocker):
        svc = RelativeRateTest(_make_args())
        printed = mocker.patch("builtins.print")
        results = [
            {
                "taxon1": "A",
                "taxon2": "B",
                "m1": 1,
                "m2": 3,
                "chi2": 2.0,
                "p_value": 0.0455,
                "p_bonf": 0.091,
                "p_fdr": 0.04,
            },
            {
                "taxon1": "A",
                "taxon2": "C",
                "m1": 2,
                "m2": 2,
                "chi2": 0.0,
                "p_value": 1.0,
                "p_bonf": 1.0,
                "p_fdr": 1.0,
            },
        ]

        svc._output_single("O", results)

        printed.assert_called_once_with(
            "Outgroup: O\n"
            "Number of ingroup taxa: 3\n"
            "Number of pairwise tests: 2\n"
            "---\n"
            "taxon1\ttaxon2\tm1\tm2\tchi2\tp_value\tp_bonf\tp_fdr\tsignificant\n"
            "A\tB\t1\t3\t2.0000\t4.5500e-02\t9.1000e-02\t4.0000e-02\t*\n"
            "A\tC\t2\t2\t0.0000\t1.0000e+00\t1.0000e+00\t1.0000e+00\t"
        )

    def test_output_batch_batches_text_rows(self, mocker):
        svc = RelativeRateTest(_make_args())
        printed = mocker.patch("builtins.print")
        all_results = {
            ("A", "B"): [
                {"p_value": 0.01, "chi2": 4.0},
                {"p_value": 0.50, "chi2": 2.0},
                {"p_value": 0.20, "chi2": 6.0},
            ],
            ("A", "C"): [
                {"p_value": 0.80, "chi2": 0.5},
            ],
        }

        svc._output_batch("O", 3, all_results)

        printed.assert_called_once_with(
            "Number of alignments: 3\n"
            "Outgroup: O\n"
            "---\n"
            "taxon1\ttaxon2\tn_reject\tn_total\tpct_reject\tmedian_chi2\n"
            "A\tB\t1\t3\t33.3\t4.0000\n"
            "A\tC\t0\t1\t0.0\t0.5000"
        )

    def test_output_batch_json_payload(self, mocker):
        svc = RelativeRateTest(_make_args(json=True))
        print_json = mocker.patch.object(relative_rate_test_module, "print_json")
        all_results = {
            ("A", "B"): [
                {"p_value": 0.01, "chi2": 4.0},
                {"p_value": 0.50, "chi2": 2.0},
                {"p_value": 0.20, "chi2": 6.0},
            ],
            ("A", "C"): [
                {"p_value": 0.80, "chi2": 0.5},
            ],
        }

        svc._output_batch("O", 3, all_results)

        print_json.assert_called_once_with({
            "outgroup": "O",
            "n_alignments": 3,
            "pairs": [
                {
                    "taxon1": "A",
                    "taxon2": "B",
                    "n_reject": 1,
                    "n_total": 3,
                    "pct_reject": 33.3,
                    "median_chi2": 4.0,
                },
                {
                    "taxon1": "A",
                    "taxon2": "C",
                    "n_reject": 0,
                    "n_total": 1,
                    "pct_reject": 0.0,
                    "median_chi2": 0.5,
                },
            ],
        })

    def test_summarize_batch_gene_results_matches_output_summary(self):
        gene_results = [
            {"p_value": 0.01, "chi2": 6.0},
            {"p_value": 0.20, "chi2": 2.0},
            {"p_value": 0.03, "chi2": 8.0},
            {"p_value": 0.50, "chi2": 4.0},
        ]

        assert RelativeRateTest._summarize_batch_gene_results(gene_results) == (
            2,
            4,
            50.0,
            6.0,
        )


class TestTajimaTest:
    def test_long_ascii_sequences_use_vectorized_path(self, monkeypatch):
        length = 3000
        seq_out = "A" * length
        seq_x = "C" * 100 + "A" * (length - 100)
        seq_y = "A" * length

        def fail_scalar(*_args, **_kwargs):
            raise AssertionError("long ASCII Tajima tests should use vector path")

        monkeypatch.setattr(
            RelativeRateTest,
            "_tajima_test_scalar",
            staticmethod(fail_scalar),
        )

        result = RelativeRateTest._tajima_test(seq_x, seq_y, seq_out)

        assert result["m1"] == 100
        assert result["m2"] == 0
        assert result["chi2"] == 100.0

    def test_mid_sized_ascii_sequences_use_vectorized_path(self, monkeypatch):
        length = 1500
        seq_out = "A" * length
        seq_x = "C" * 60 + "A" * (length - 60)
        seq_y = "A" * 90 + "G" * 30 + "A" * (length - 120)

        def fail_scalar(*_args, **_kwargs):
            raise AssertionError(
                "mid-sized ASCII Tajima tests should use vector path"
            )

        monkeypatch.setattr(
            RelativeRateTest,
            "_tajima_test_scalar",
            staticmethod(fail_scalar),
        )

        result = RelativeRateTest._tajima_test(seq_x, seq_y, seq_out)

        assert result["m1"] == 60
        assert result["m2"] == 30
        assert result["chi2"] == 10.0

    def test_vectorized_tajima_matches_scalar_with_lowercase_and_skips(self):
        seq_out = ("ACGT" * 750)
        seq_x = list(seq_out)
        seq_y = list(seq_out.lower())
        for idx in range(0, len(seq_out), 13):
            seq_x[idx] = "t" if seq_out[idx] != "T" else "a"
        for idx in range(5, len(seq_out), 17):
            seq_y[idx] = "c" if seq_out[idx] != "C" else "g"
        for idx in range(3, len(seq_out), 101):
            seq_x[idx] = "n"

        seq_x = "".join(seq_x)
        seq_y = "".join(seq_y)

        assert RelativeRateTest._tajima_test(
            seq_x,
            seq_y,
            seq_out,
        ) == RelativeRateTest._tajima_test_scalar(seq_x, seq_y, seq_out)

    def test_clean_ascii_tajima_skips_validity_mask(self, monkeypatch):
        length = 3000
        seq_out = "A" * length
        seq_x = "C" * 120 + "A" * (length - 120)
        seq_y = "A" * 200 + "G" * 80 + "A" * (length - 280)

        def fail_valid_mask(*_args, **_kwargs):
            raise AssertionError(
                "clean ASCII Tajima tests should not build a validity mask"
            )

        monkeypatch.setattr(
            RelativeRateTest,
            "_tajima_valid_mask",
            staticmethod(fail_valid_mask),
        )

        result = RelativeRateTest._tajima_test(seq_x, seq_y, seq_out)

        assert result["m1"] == 120
        assert result["m2"] == 80
        assert result["chi2"] == 8.0

    def test_equal_sequences(self):
        """Identical ingroup sequences: m1=m2=0, chi2=0, p=1.0."""
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGT",
            seq_y="ACGTACGT",
            seq_out="ACGTACGT",
        )
        assert result["m1"] == 0
        assert result["m2"] == 0
        assert result["chi2"] == 0.0
        assert result["p_value"] == 1.0

    def test_one_lineage_faster(self):
        """One lineage has unique changes, the other doesn't."""
        # x matches outgroup, y differs at 2 sites
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGT",
            seq_y="ACTTACTT",
            seq_out="ACGTACGT",
        )
        assert result["m1"] == 0
        assert result["m2"] == 2
        assert result["chi2"] == 2.0
        assert abs(result["p_value"] - 0.1573) < 0.001

    def test_symmetric_changes(self):
        """Equal unique changes on both lineages: chi2=0, p=1."""
        result = RelativeRateTest._tajima_test(
            seq_x="TCGTACGT",
            seq_y="ACGTACTT",
            seq_out="ACGTACGT",
        )
        assert result["m1"] == 1
        assert result["m2"] == 1
        assert result["chi2"] == 0.0
        assert result["p_value"] == 1.0

    def test_highly_unequal_rates(self):
        """Many changes on one lineage: should reject clock."""
        # seq_y differs from outgroup at 6 sites (positions 0,1,2,4,5,6)
        # positions 3,7 have T in all three sequences (no difference)
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGTACGT",
            seq_y="TTTTTTTTACGT",
            seq_out="ACGTACGTACGT",
        )
        assert result["m1"] == 0
        assert result["m2"] == 6
        assert result["chi2"] == 6.0
        assert result["p_value"] < 0.02

    def test_gaps_are_skipped(self):
        """Sites with gaps in any sequence should be excluded."""
        result = RelativeRateTest._tajima_test(
            seq_x="A-GTACGT",
            seq_y="ACTTACGT",
            seq_out="ACGTACGT",
        )
        # Position 1 has gap in x, skip it. Position 2: y=T, out=G, x=G -> m2.
        # Only non-gap informative site where y differs: pos 2 (y=T, out=G, x=G)
        assert result["m1"] == 0
        assert result["m2"] == 1

    def test_shared_changes_not_counted(self):
        """Sites where BOTH ingroup differ from outgroup are uninformative."""
        result = RelativeRateTest._tajima_test(
            seq_x="TCGTACGT",
            seq_y="TCGTACGT",
            seq_out="ACGTACGT",
        )
        # Position 0: both x and y differ from outgroup -> shared, skip
        assert result["m1"] == 0
        assert result["m2"] == 0
        assert result["chi2"] == 0.0


class TestIdentifyOutgroup:
    def test_simple_rooted_tree(self):
        """Outgroup should be the earliest-diverging taxon."""
        from io import StringIO
        from Bio import Phylo

        tree = Phylo.read(StringIO("(((A:0.1,B:0.1):0.05,C:0.15):0.1,O:0.2);"), "newick")
        outgroup = RelativeRateTest._identify_outgroup(tree)
        assert outgroup == "O"

    def test_two_taxa_at_root(self):
        """With two clades at root, outgroup is the smaller one."""
        from io import StringIO
        from Bio import Phylo

        tree = Phylo.read(StringIO("((A:0.1,B:0.1):0.1,O:0.2);"), "newick")
        outgroup = RelativeRateTest._identify_outgroup(tree)
        assert outgroup == "O"

    def test_identify_outgroup_uses_direct_singleton_scan(self, monkeypatch):
        from io import StringIO
        from Bio import Phylo
        from Bio.Phylo.BaseTree import TreeMixin

        tree = Phylo.read(
            StringIO("(((A:1,B:1):1,C:1):1,O:1);"),
            "newick",
        )

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("outgroup scan should not materialize terminals")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        outgroup = RelativeRateTest._identify_outgroup(tree)

        assert outgroup == "O"

    def test_identify_outgroup_generic_fallback_stops_after_second_tip(self):
        class Tip:
            def __init__(self, name):
                self.name = name

        class Child:
            def __init__(self, name, terminal_count):
                self.name = name
                self.terminal_count = terminal_count

            @property
            def clades(self):
                raise AttributeError

            def get_terminals(self):
                for idx in range(self.terminal_count):
                    if idx > 1:
                        raise AssertionError("fallback should stop after second tip")
                    yield Tip(f"{self.name}{idx}")

        class Tree:
            root = type(
                "Root",
                (),
                {"clades": [Child("large", 10), Child("out", 1)]},
            )()

        assert RelativeRateTest._identify_outgroup(Tree()) == "out0"

    def test_identify_outgroup_generic_fallback_raises_without_singleton(self):
        class Tip:
            name = "tip"

        class Child:
            @property
            def clades(self):
                raise AttributeError

            def get_terminals(self):
                yield Tip()
                yield Tip()
                raise AssertionError("fallback should stop after second tip")

        class Tree:
            root = type("Root", (), {"clades": [Child(), Child()]})()

        with pytest.raises(relative_rate_test_module.PhykitUserError):
            RelativeRateTest._identify_outgroup(Tree())


class TestMultipleTestingCorrection:
    def test_bonferroni(self):
        p_values = [0.01, 0.04, 0.06]
        corrected = RelativeRateTest._bonferroni(p_values)
        assert corrected[0] == 0.03  # 0.01 * 3
        assert corrected[1] == 0.12  # 0.04 * 3
        assert corrected[2] == 0.18  # 0.06 * 3

    def test_bonferroni_capped_at_one(self):
        p_values = [0.5, 0.8]
        corrected = RelativeRateTest._bonferroni(p_values)
        assert corrected[0] == 1.0
        assert corrected[1] == 1.0

    def test_bonferroni_empty_list(self):
        assert RelativeRateTest._bonferroni([]) == []

    def test_fdr(self):
        p_values = [0.01, 0.04, 0.06]
        corrected = RelativeRateTest._fdr(p_values)
        # BH-FDR: rank 1: 0.01*3/1=0.03, rank 2: 0.04*3/2=0.06, rank 3: 0.06*3/3=0.06
        assert abs(corrected[0] - 0.03) < 1e-10
        assert abs(corrected[1] - 0.06) < 1e-10
        assert abs(corrected[2] - 0.06) < 1e-10

    def test_fdr_matches_scalar_reference_with_ties(self):
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

        assert RelativeRateTest._fdr(p_values) == pytest.approx(expected)

    def test_medium_corrections_match_scalar_references(self):
        p_values = [((idx * 37) % 101) / 1000 for idx in range(32)]

        expected_bonf = [min(value * len(p_values), 1.0) for value in p_values]
        assert RelativeRateTest._bonferroni(p_values) == pytest.approx(expected_bonf)

        indexed = sorted(enumerate(p_values), key=lambda x: x[1])
        expected_fdr = [0.0] * len(p_values)
        previous = 1.0
        for rank_minus_1 in range(len(p_values) - 1, -1, -1):
            original_idx, p_value = indexed[rank_minus_1]
            rank = rank_minus_1 + 1
            adjusted = min(p_value * len(p_values) / rank, previous)
            adjusted = min(adjusted, 1.0)
            expected_fdr[original_idx] = adjusted
            previous = adjusted

        assert RelativeRateTest._fdr(p_values) == pytest.approx(expected_fdr)

    def test_vector_fdr_preserves_input_values(self):
        p_values = [((idx * 37) % 101) / 101 for idx in range(64)]
        original = list(p_values)

        RelativeRateTest._fdr(p_values)

        assert p_values == original

    def test_small_corrections_do_not_import_numpy(self):
        code = """
import sys
from phykit.services.tree.relative_rate_test import RelativeRateTest
results = [
    {"p_value": 0.20},
    {"p_value": 0.01},
    {"p_value": 0.03},
]
svc = RelativeRateTest.__new__(RelativeRateTest)
svc._add_multiple_testing_corrections(results)
assert [row["p_bonf"] for row in results] == [0.6000000000000001, 0.03, 0.09]
assert [round(row["p_fdr"], 10) for row in results] == [0.2, 0.03, 0.045]
assert "numpy" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)


class TestMatchesRReference:
    """Validate against R's pegas::rr.test() reference values.

    R code used to generate reference values:
        library(pegas); library(ape)
        aln <- read.FASTA("rrt_test.fa")
        rr.test(aln[["A"]], aln[["B"]], aln[["O"]])
    """

    def test_a_vs_b_outgroup_o(self):
        """R: Chi=4.0, Pval=0.0455002639."""
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGTACGTACGTACGT",
            seq_y="ACGTACGTACTAACTAACGT",
            seq_out="ACGTACGTACGTACGTACGT",
        )
        assert result["m1"] == 0
        assert result["m2"] == 4
        assert result["chi2"] == 4.0
        assert abs(result["p_value"] - 0.0455002639) < 1e-8

    def test_a_vs_c_outgroup_o(self):
        """R: Chi=2.0, Pval=0.1572992071."""
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGTACGTACGTACGT",
            seq_y="ACTAACGTACGTACGTACGT",
            seq_out="ACGTACGTACGTACGTACGT",
        )
        assert result["m1"] == 0
        assert result["m2"] == 2
        assert result["chi2"] == 2.0
        assert abs(result["p_value"] - 0.1572992071) < 1e-8

    def test_b_vs_c_outgroup_o(self):
        """R: Chi=0.6667, Pval=0.4142161782."""
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGTACTAACTAACGT",
            seq_y="ACTAACGTACGTACGTACGT",
            seq_out="ACGTACGTACGTACGTACGT",
        )
        assert result["m1"] == 4
        assert result["m2"] == 2
        assert abs(result["chi2"] - 0.6666666667) < 1e-8
        assert abs(result["p_value"] - 0.4142161782) < 1e-8

    def test_one_lineage_faster_matches_r(self):
        """R: Chi=2.0, Pval=0.1572992071."""
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGT",
            seq_y="ACTTACTT",
            seq_out="ACGTACGT",
        )
        assert result["m1"] == 0
        assert result["m2"] == 2
        assert result["chi2"] == 2.0
        assert abs(result["p_value"] - 0.1572992071) < 1e-8

    def test_tajima_scalar_p_value_does_not_import_scipy(self, monkeypatch):
        real_import = __import__

        def fail_scipy_import(name, *args, **kwargs):
            if name == "scipy" or name.startswith("scipy."):
                raise AssertionError("df=1 chi-square p-value should not import scipy")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr("builtins.__import__", fail_scipy_import)

        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGT",
            seq_y="ACTTACTT",
            seq_out="ACGTACGT",
        )

        assert abs(result["p_value"] - 0.1572992071) < 1e-8

    def test_equal_sequences_edge_case(self):
        """R returns NaN for m1+m2=0; we return chi2=0, p=1.0."""
        result = RelativeRateTest._tajima_test(
            seq_x="ACGTACGT",
            seq_y="ACGTACGT",
            seq_out="ACGTACGT",
        )
        assert result["chi2"] == 0.0
        assert result["p_value"] == 1.0
