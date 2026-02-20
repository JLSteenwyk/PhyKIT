from argparse import Namespace

import pytest
import builtins
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.rcvt import RelativeCompositionVariabilityTaxon
import phykit.services.alignment.rcvt as rcvt_module


def _alignment():
    return MultipleSeqAlignment(
        [
            SeqRecord(Seq("ACGT"), id="t1"),
            SeqRecord(Seq("AC-T"), id="t2"),
            SeqRecord(Seq("TCGT"), id="t3"),
        ]
    )


@pytest.fixture
def args():
    return Namespace(alignment="/some/path/to/file.fa")


class TestRCVT:
    def test_init_sets_expected_attrs(self, args):
        service = RelativeCompositionVariabilityTaxon(args)
        assert service.alignment_file_path == args.alignment
        assert service.json_output is False
        assert service.plot is False
        assert service.plot_output == "rcvt_plot.png"

    def test_calculate_rows_nucleotide(self, args):
        service = RelativeCompositionVariabilityTaxon(args)
        rows = service.calculate_rows(_alignment(), is_protein=False)
        assert len(rows) == 3
        assert rows[0]["taxon"] == "t1"
        assert isinstance(rows[0]["rcvt"], float)

    def test_calculate_rows_empty_alignment(self, args):
        service = RelativeCompositionVariabilityTaxon(args)
        rows = service.calculate_rows(MultipleSeqAlignment([]), is_protein=False)
        assert rows == []

    def test_run_json_with_plot(self, mocker):
        args = Namespace(alignment="/some/path/to/file.fa", json=True, plot=True, plot_output="rcvt.png")
        service = RelativeCompositionVariabilityTaxon(args)
        alignment = _alignment()
        mocker.patch.object(
            RelativeCompositionVariabilityTaxon,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocked_plot = mocker.patch.object(RelativeCompositionVariabilityTaxon, "_plot_rcvt")
        mocked_json = mocker.patch("phykit.services.alignment.rcvt.print_json")

        service.run()

        mocked_plot.assert_called_once()
        payload = mocked_json.call_args.args[0]
        assert payload["rows"] == payload["taxa"]
        assert payload["plot_output"] == "rcvt.png"

    def test_process_args_defaults(self):
        parsed = RelativeCompositionVariabilityTaxon(Namespace(alignment="x.fa")).process_args(
            Namespace(alignment="x.fa")
        )
        assert parsed["json_output"] is False
        assert parsed["plot"] is False
        assert parsed["plot_output"] == "rcvt_plot.png"

    def test_calculate_rows_all_invalid_chars(self, args):
        service = RelativeCompositionVariabilityTaxon(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("----"), id="t1"),
                SeqRecord(Seq("XXXX"), id="t2"),
            ]
        )
        rows = service.calculate_rows(alignment, is_protein=True)
        assert rows == [{"taxon": "t1", "rcvt": 0.0}, {"taxon": "t2", "rcvt": 0.0}]

    def test_plot_rcvt_creates_file(self, tmp_path):
        pytest.importorskip("matplotlib")
        out = tmp_path / "rcvt.png"
        service = RelativeCompositionVariabilityTaxon(
            Namespace(alignment="x.fa", plot=True, plot_output=str(out))
        )
        service._plot_rcvt([{"taxon": "a", "rcvt": 0.1}, {"taxon": "b", "rcvt": 0.2}])
        assert out.exists()

    def test_plot_rcvt_importerror(self, monkeypatch, capsys):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name.startswith("matplotlib"):
                raise ImportError("no matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        service = RelativeCompositionVariabilityTaxon(Namespace(alignment="x.fa"))
        with pytest.raises(SystemExit) as exc:
            service._plot_rcvt([{"taxon": "a", "rcvt": 0.1}])
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "matplotlib is required for --plot in rcvt" in out

    def test_run_text_with_plot(self, mocker, capsys):
        service = RelativeCompositionVariabilityTaxon(
            Namespace(alignment="x.fa", json=False, plot=True, plot_output="rcvt_plot.png")
        )
        mocker.patch.object(
            RelativeCompositionVariabilityTaxon,
            "get_alignment_and_format",
            return_value=(_alignment(), "fasta", False),
        )
        mocker.patch.object(
            RelativeCompositionVariabilityTaxon,
            "calculate_rows",
            return_value=[{"taxon": "a", "rcvt": 0.1}, {"taxon": "b", "rcvt": 0.2}],
        )
        mocked_plot = mocker.patch.object(RelativeCompositionVariabilityTaxon, "_plot_rcvt")
        service.run()
        out, _ = capsys.readouterr()
        assert "a\t0.1" in out
        assert "b\t0.2" in out
        assert "Saved RCVT plot: rcvt_plot.png" in out
        mocked_plot.assert_called_once()
