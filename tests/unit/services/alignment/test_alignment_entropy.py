import pytest
import subprocess
import sys
import numpy as np
from argparse import Namespace
import builtins
from math import isclose
from types import SimpleNamespace

from phykit.services.alignment.alignment_entropy import AlignmentEntropy
import phykit.services.alignment.alignment_entropy as alignment_entropy_module


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa", verbose=False)
    return Namespace(**kwargs)


class TestAlignmentEntropy(object):
    def test_module_import_does_not_import_numpy_or_biopython_align(self):
        code = """
import sys
import phykit.services.alignment.alignment_entropy as module
assert hasattr(module.np, "__getattr__")
assert "json" not in sys.modules
assert "hashlib" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
assert "typing" not in sys.modules
assert "numpy" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_init_sets_alignment_file_path(self, args):
        entropy = AlignmentEntropy(args)
        assert entropy.alignment_file_path == args.alignment
        assert entropy.verbose is False

    def test_site_entropies(self, alignment_simple, args):
        entropy = AlignmentEntropy(args)
        entropies = entropy.calculate_site_entropies(alignment_simple, is_protein=False)

        assert len(entropies) == 6
        assert isclose(entropies[0], 0.0, rel_tol=0.001)
        assert isclose(entropies[1], 1.0, rel_tol=0.001)
        assert isclose(entropies[2], 0.97095, rel_tol=0.001)
        assert isclose(entropies[3], 0.0, rel_tol=0.001)
        assert isclose(entropies[4], 0.97095, rel_tol=0.001)
        assert isclose(entropies[5], 1.0, rel_tol=0.001)

    def test_process_args_defaults(self):
        parsed = AlignmentEntropy(Namespace(alignment="x.fa", verbose=True)).process_args(
            Namespace(alignment="x.fa", verbose=True)
        )
        assert parsed["json_output"] is False
        assert parsed["plot"] is False
        assert parsed["plot_output"] == "alignment_entropy_plot.png"

    def test_site_entropies_protein_invalid_characters(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        entropy = AlignmentEntropy(Namespace(alignment="x.fa", verbose=False))
        alignment = [
            SeqRecord(Seq("AX"), id="t1"),
            SeqRecord(Seq("AX"), id="t2"),
        ]
        entropies = entropy.calculate_site_entropies(alignment, is_protein=True)
        assert entropies[0] == 0.0
        assert entropies[1] == 0.0

    def test_site_entropies_ignores_ambiguous_sites_and_uppercases(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        entropy = AlignmentEntropy(Namespace(alignment="x.fa", verbose=False))
        alignment = [
            SeqRecord(Seq("aN-"), id="t1"),
            SeqRecord(Seq("AN?"), id="t2"),
            SeqRecord(Seq("tNX"), id="t3"),
        ]

        entropies = entropy.calculate_site_entropies(alignment, is_protein=False)

        assert len(entropies) == 3
        assert isclose(entropies[0], 0.918295, rel_tol=0.001)
        assert entropies[1] == 0.0
        assert entropies[2] == 0.0

    def test_site_entropies_ascii_matrix_path(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        entropy = AlignmentEntropy(Namespace(alignment="x.fa", verbose=False))
        alignment = [
            SeqRecord(Seq("aa"), id="t1"),
            SeqRecord(Seq("at"), id="t2"),
            SeqRecord(Seq("tt"), id="t3"),
        ]

        entropies = entropy.calculate_site_entropies(alignment, is_protein=False)

        assert len(entropies) == 2
        assert isclose(entropies[0], 0.918295, rel_tol=0.001)
        assert isclose(entropies[1], 0.918295, rel_tol=0.001)

    def test_site_entropies_ascii_matrix_path_avoids_isin(self, mocker):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        entropy = AlignmentEntropy(Namespace(alignment="x.fa", verbose=False))
        alignment = [
            SeqRecord(Seq("aa"), id="t1"),
            SeqRecord(Seq("a-"), id="t2"),
            SeqRecord(Seq("tt"), id="t3"),
        ]
        mocker.patch(
            "phykit.services.alignment.alignment_entropy.np.isin",
            side_effect=AssertionError("ASCII path should use the lookup table"),
        )

        entropies = entropy.calculate_site_entropies(alignment, is_protein=False)

        assert len(entropies) == 2
        assert isclose(entropies[0], 0.918295, rel_tol=0.001)
        assert isclose(entropies[1], 1.0, rel_tol=0.001)

    def test_site_entropies_ascii_path_uses_gap_code_reduction(self, mocker):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        entropy = AlignmentEntropy(Namespace(alignment="x.fa", verbose=False))
        alignment = [
            SeqRecord(Seq("aa"), id="t1"),
            SeqRecord(Seq("a-"), id="t2"),
            SeqRecord(Seq("tt"), id="t3"),
        ]
        gap_codes_spy = mocker.spy(alignment_entropy_module, "_get_gap_codes")

        entropies = entropy.calculate_site_entropies(alignment, is_protein=False)

        gap_codes_spy.assert_called_once_with(False)
        assert len(entropies) == 2
        assert isclose(entropies[0], 0.918295, rel_tol=0.001)
        assert isclose(entropies[1], 1.0, rel_tol=0.001)

    def test_site_entropies_ascii_path_skips_gap_char_set(self, mocker):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        entropy = AlignmentEntropy(Namespace(alignment="x.fa", verbose=False))
        alignment = [
            SeqRecord(Seq("aa"), id="t1"),
            SeqRecord(Seq("a-"), id="t2"),
            SeqRecord(Seq("tt"), id="t3"),
        ]
        mocker.patch.object(
            entropy,
            "get_gap_chars",
            side_effect=AssertionError("ASCII path should use byte gap codes"),
        )

        entropies = entropy.calculate_site_entropies(alignment, is_protein=False)

        assert len(entropies) == 2
        assert isclose(entropies[0], 0.918295, rel_tol=0.001)
        assert isclose(entropies[1], 1.0, rel_tol=0.001)

    def test_sparse_protein_entropy_avoids_where_log_terms(self, monkeypatch):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        entropy = AlignmentEntropy(Namespace(alignment="x.fa", verbose=False))
        alignment = [
            SeqRecord(Seq("ACDEFGHIK"), id="t1"),
            SeqRecord(Seq("AAAAAAAAA"), id="t2"),
            SeqRecord(Seq("AAAAAAAAA"), id="t3"),
        ]

        def fail_where(*args, **kwargs):
            raise AssertionError("entropy terms should use masked np.log2")

        monkeypatch.setattr(
            alignment_entropy_module.np, "where", fail_where, raising=False
        )

        entropies = entropy.calculate_site_entropies(alignment, is_protein=True)

        assert len(entropies) == 9
        assert entropies[0] == 0.0
        assert all(value >= 0.0 for value in entropies)

    def test_site_entropies_single_valid_symbol_returns_zeroes_directly(
        self, mocker
    ):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        entropy = AlignmentEntropy(Namespace(alignment="x.fa", verbose=False))
        alignment = [
            SeqRecord(Seq("AAAA"), id="t1"),
            SeqRecord(Seq("A-A?"), id="t2"),
            SeqRecord(Seq("AAAA"), id="t3"),
        ]
        mocker.patch(
            "phykit.services.alignment.alignment_entropy.np.vstack",
            side_effect=AssertionError("single-symbol entropy should skip counts"),
        )

        entropies = entropy.calculate_site_entropies(alignment, is_protein=False)

        assert entropies == [0.0, 0.0, 0.0, 0.0]

    def test_site_entropies_identical_sequences_return_zeroes_before_matrix(
        self, mocker
    ):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        entropy = AlignmentEntropy(Namespace(alignment="x.fa", verbose=False))
        alignment = [
            SeqRecord(Seq("AcGt"), id="t1"),
            SeqRecord(Seq("aCgT"), id="t2"),
            SeqRecord(Seq("ACGT"), id="t3"),
            SeqRecord(Seq("acgt"), id="t4"),
        ]
        mocker.patch(
            "phykit.services.alignment.alignment_entropy.np.frombuffer",
            side_effect=AssertionError("identical alignments should skip matrix setup"),
        )

        entropies = entropy.calculate_site_entropies(alignment, is_protein=False)

        assert entropies == [0.0, 0.0, 0.0, 0.0]

    def test_entropy_from_ascii_codes_matches_repeated_counts(self):
        alignment_array = np.array(
            [
                list(b"ACDEFX"),
                list(b"ACGHIX"),
                list(b"LMNPQ-"),
            ],
            dtype=np.uint8,
        )
        invalid_mask = np.zeros(alignment_array.shape, dtype=np.bool_)
        for gap_code in alignment_entropy_module._get_gap_codes(is_protein=True):
            invalid_mask |= alignment_array == gap_code
        valid_mask = ~invalid_mask
        valid_chars = np.unique(alignment_array[valid_mask])

        observed = alignment_entropy_module._entropy_from_ascii_codes(
            alignment_array,
            valid_mask,
            valid_chars,
            block_size=2,
        )

        counts = np.vstack(
            [
                np.sum(alignment_array == char, axis=0)
                for char in valid_chars
            ]
        ).astype(np.float64)
        totals = np.sum(counts, axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            probs = np.divide(
                counts,
                totals,
                out=np.zeros_like(counts, dtype=np.float64),
                where=totals > 0,
            )
            expected = -np.sum(
                np.where(probs > 0, probs * np.log2(probs), 0.0),
                axis=0,
            )

        np.testing.assert_allclose(observed, expected)

    def test_entropy_from_ascii_codes_all_valid_skips_mask_scan(self, mocker):
        alignment_array = np.array(
            [
                list(b"ACDEF"),
                list(b"ACGHI"),
                list(b"LMNPQ"),
            ],
            dtype=np.uint8,
        )
        valid_mask = np.ones_like(alignment_array, dtype=bool)
        valid_chars = np.unique(alignment_array)
        expected = alignment_entropy_module._entropy_from_ascii_codes(
            alignment_array,
            valid_mask,
            valid_chars,
            block_size=2,
        )
        mocker.patch(
            "phykit.services.alignment.alignment_entropy.np.any",
            side_effect=AssertionError("all-valid blocks should not scan masks"),
        )

        observed = alignment_entropy_module._entropy_from_ascii_codes(
            alignment_array,
            valid_mask,
            valid_chars,
            block_size=2,
            all_valid=True,
        )

        np.testing.assert_allclose(observed, expected)

    def test_entropy_values_to_list_uses_array_tolist(self):
        class NoIterArray:
            def __iter__(self):
                raise AssertionError("entropy conversion should use tolist")

            def tolist(self):
                return [0.0, 1.25, 2.5]

        assert alignment_entropy_module._entropy_values_to_list(NoIterArray()) == [
            0.0,
            1.25,
            2.5,
        ]

    def test_prepare_entropy_plot_series_matches_row_series(self):
        entropies = [0.123456, 0.5, 1.98765]
        rows = [
            {"site": idx, "entropy": round(entropy, 4)}
            for idx, entropy in enumerate(entropies, start=1)
        ]

        sites, values = alignment_entropy_module._prepare_entropy_plot_series(entropies)
        row_sites, row_values = alignment_entropy_module._prepare_entropy_row_plot_series(rows)

        np.testing.assert_array_equal(sites, row_sites)
        np.testing.assert_array_equal(values, row_values)

    def test_entropy_from_ascii_codes_avoids_tiled_site_offsets(self, mocker):
        alignment_array = np.array(
            [
                list(b"ACDEFX"),
                list(b"ACGHIX"),
                list(b"LMNPQ-"),
            ],
            dtype=np.uint8,
        )
        invalid_mask = np.zeros(alignment_array.shape, dtype=np.bool_)
        for gap_code in alignment_entropy_module._get_gap_codes(is_protein=True):
            invalid_mask |= alignment_array == gap_code
        valid_mask = ~invalid_mask
        valid_chars = np.unique(alignment_array[valid_mask])
        mocker.patch(
            "phykit.services.alignment.alignment_entropy.np.tile",
            side_effect=AssertionError("ASCII entropy path should avoid tiled offsets"),
        )

        observed = alignment_entropy_module._entropy_from_ascii_codes(
            alignment_array,
            valid_mask,
            valid_chars,
            block_size=2,
        )

        assert len(observed) == alignment_array.shape[1]
        assert observed[-1] == 0.0

    def test_site_entropies_protein_uses_chunked_ascii_counts(self, mocker):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        entropy = AlignmentEntropy(Namespace(alignment="x.fa", verbose=False))
        alignment = [
            SeqRecord(Seq("ACDEFGHIK"), id="t1"),
            SeqRecord(Seq("LMNPQRSTV"), id="t2"),
            SeqRecord(Seq("WYACDEFGH"), id="t3"),
        ]
        spy = mocker.spy(alignment_entropy_module, "_entropy_from_ascii_codes")

        entropies = entropy.calculate_site_entropies(alignment, is_protein=True)

        assert spy.call_count == 1
        assert spy.call_args.kwargs["all_valid"] is True
        assert len(entropies) == 9

    def test_site_entropies_unicode_fallback(self):
        entropy = AlignmentEntropy(Namespace(alignment="x.fa", verbose=False))
        alignment = [
            SimpleNamespace(seq="A\u00d1", id="t1"),
            SimpleNamespace(seq="A\u00d1", id="t2"),
            SimpleNamespace(seq="T\u00d1", id="t3"),
        ]

        entropies = entropy.calculate_site_entropies(alignment, is_protein=False)

        assert len(entropies) == 2
        assert isclose(entropies[0], 0.918295, rel_tol=0.001)
        assert entropies[1] == 0.0

    def test_plot_alignment_entropy_creates_file(self, tmp_path):
        pytest.importorskip("matplotlib")
        output = tmp_path / "entropy.png"
        svc = AlignmentEntropy(
            Namespace(alignment="x.fa", verbose=False, plot=True, plot_output=str(output))
        )
        svc._plot_alignment_entropy([{"site": 1, "entropy": 0.2}, {"site": 2, "entropy": 0.8}])
        assert output.exists()

    def test_plot_alignment_entropy_importerror_exits(self, monkeypatch, capsys):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name.startswith("matplotlib"):
                raise ImportError("no matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        svc = AlignmentEntropy(Namespace(alignment="x.fa", verbose=False))
        with pytest.raises(SystemExit) as exc:
            svc._plot_alignment_entropy([{"site": 1, "entropy": 0.2}])
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "matplotlib is required for --plot in alignment_entropy" in out

    def test_run_json_verbose_and_plot(self, tmp_path, monkeypatch):
        pytest.importorskip("matplotlib")
        output = tmp_path / "entropy.json.png"
        captured = {}
        svc = AlignmentEntropy(
            Namespace(
                alignment="x.fa",
                verbose=True,
                json=True,
                plot=True,
                plot_output=str(output),
            )
        )
        monkeypatch.setattr(svc, "get_alignment_and_format", lambda: ("aln", None, False))
        monkeypatch.setattr(svc, "calculate_site_entropies", lambda _a, _p: [0.25, 0.75])
        monkeypatch.setattr(svc, "_plot_alignment_entropy", lambda rows: captured.setdefault("rows", rows))
        monkeypatch.setattr(alignment_entropy_module, "print_json", lambda payload: captured.setdefault("payload", payload))

        svc.run()

        assert captured["rows"] == [{"site": 1, "entropy": 0.25}, {"site": 2, "entropy": 0.75}]
        assert captured["payload"]["verbose"] is True
        assert captured["payload"]["plot_output"] == str(output)
        assert captured["payload"]["rows"][0]["site"] == 1

    def test_run_nonverbose_terminal_output(self, monkeypatch, capsys):
        svc = AlignmentEntropy(Namespace(alignment="x.fa", verbose=False, json=False, plot=False))
        monkeypatch.setattr(svc, "get_alignment_and_format", lambda: ("aln", None, False))
        monkeypatch.setattr(svc, "calculate_site_entropies", lambda _a, _p: [0.5, 1.0])

        svc.run()

        out, _ = capsys.readouterr()
        assert out.strip() == "0.75"

    def test_run_nonverbose_plot_uses_value_series(self, monkeypatch, capsys):
        captured = {}
        svc = AlignmentEntropy(
            Namespace(alignment="x.fa", verbose=False, json=False, plot=True, plot_output="entropy.png")
        )
        monkeypatch.setattr(svc, "get_alignment_and_format", lambda: ("aln", None, False))
        monkeypatch.setattr(svc, "calculate_site_entropies", lambda _a, _p: [0.25, 0.75])
        monkeypatch.setattr(
            svc,
            "_plot_alignment_entropy_values",
            lambda values: captured.setdefault("values", values),
        )

        svc.run()

        out, _ = capsys.readouterr()
        assert captured["values"] == [0.25, 0.75]
        assert out == "0.5\nSaved alignment entropy plot: entropy.png\n"

    def test_run_nonverbose_json_plot_uses_value_series(self, monkeypatch):
        captured = {}
        svc = AlignmentEntropy(
            Namespace(alignment="x.fa", verbose=False, json=True, plot=True, plot_output="entropy.png")
        )
        monkeypatch.setattr(svc, "get_alignment_and_format", lambda: ("aln", None, False))
        monkeypatch.setattr(svc, "calculate_site_entropies", lambda _a, _p: [0.25, 0.75])
        monkeypatch.setattr(
            svc,
            "_plot_alignment_entropy_values",
            lambda values: captured.setdefault("values", values),
        )
        monkeypatch.setattr(
            alignment_entropy_module,
            "print_json",
            lambda payload: captured.setdefault("payload", payload),
        )

        svc.run()

        assert captured["values"] == [0.25, 0.75]
        assert captured["payload"] == {
            "verbose": False,
            "mean_entropy": 0.5,
            "plot_output": "entropy.png",
        }

    def test_run_nonverbose_summary_does_not_use_numpy_mean(
        self, monkeypatch, capsys
    ):
        class NoRowRound(float):
            def __round__(self, *_args, **_kwargs):
                raise AssertionError("summary-only output should not build rows")

        svc = AlignmentEntropy(Namespace(alignment="x.fa", verbose=False, json=False, plot=False))
        monkeypatch.setattr(svc, "get_alignment_and_format", lambda: ("aln", None, False))
        monkeypatch.setattr(
            svc,
            "calculate_site_entropies",
            lambda _a, _p: [NoRowRound(0.5), NoRowRound(1.0)],
        )
        monkeypatch.setattr(
            "phykit.services.alignment.alignment_entropy.np.mean",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("summary should use the existing Python list")
            ),
        )

        svc.run()

        out, _ = capsys.readouterr()
        assert out.strip() == "0.75"

    def test_run_nonverbose_json_summary_does_not_use_numpy_mean(
        self, monkeypatch, capsys
    ):
        svc = AlignmentEntropy(Namespace(alignment="x.fa", verbose=False, json=True, plot=False))
        monkeypatch.setattr(svc, "get_alignment_and_format", lambda: ("aln", None, False))
        monkeypatch.setattr(svc, "calculate_site_entropies", lambda _a, _p: [0.5, 1.0])
        monkeypatch.setattr(
            "phykit.services.alignment.alignment_entropy.np.mean",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("json summary should use the existing Python list")
            ),
        )

        svc.run()

        out, _ = capsys.readouterr()
        assert '"mean_entropy": 0.75' in out

    def test_run_verbose_terminal_rows(self, monkeypatch, capsys):
        svc = AlignmentEntropy(Namespace(alignment="x.fa", verbose=True, json=False, plot=False))
        monkeypatch.setattr(svc, "get_alignment_and_format", lambda: ("aln", None, True))
        monkeypatch.setattr(svc, "calculate_site_entropies", lambda _a, _p: [0.0, 1.23456])

        svc.run()

        out, _ = capsys.readouterr()
        assert out == "1\t0.0\n2\t1.2346\n"
