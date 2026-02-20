import pytest
from argparse import Namespace
import builtins
from math import isclose

from phykit.services.alignment.alignment_entropy import AlignmentEntropy
import phykit.services.alignment.alignment_entropy as alignment_entropy_module


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa", verbose=False)
    return Namespace(**kwargs)


class TestAlignmentEntropy(object):
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

    def test_run_verbose_terminal_rows(self, monkeypatch, capsys):
        svc = AlignmentEntropy(Namespace(alignment="x.fa", verbose=True, json=False, plot=False))
        monkeypatch.setattr(svc, "get_alignment_and_format", lambda: ("aln", None, True))
        monkeypatch.setattr(svc, "calculate_site_entropies", lambda _a, _p: [0.0, 1.23456])

        svc.run()

        out, _ = capsys.readouterr()
        assert "1\t0.0" in out
        assert "2\t1.2346" in out
