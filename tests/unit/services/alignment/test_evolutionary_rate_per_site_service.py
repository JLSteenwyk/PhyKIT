from argparse import Namespace

import pytest
import builtins
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.evolutionary_rate_per_site import EvolutionaryRatePerSite
import phykit.services.alignment.evolutionary_rate_per_site as erps_module


def _alignment():
    return MultipleSeqAlignment(
        [
            SeqRecord(Seq("ACGT"), id="t1"),
            SeqRecord(Seq("A-GT"), id="t2"),
            SeqRecord(Seq("TCGT"), id="t3"),
        ]
    )


@pytest.fixture
def args():
    return Namespace(alignment="/some/path/to/file.fa")


class TestEvolutionaryRatePerSite:
    def test_init_sets_expected_attrs(self, args):
        service = EvolutionaryRatePerSite(args)
        assert service.alignment_file_path == args.alignment
        assert service.json_output is False
        assert service.plot is False
        assert service.plot_output == "evolutionary_rate_per_site_plot.png"

    def test_remove_gap_characters(self, args):
        service = EvolutionaryRatePerSite(args)
        assert service.remove_gap_characters("A-C?T", ["-", "?"]) == "ACT"

    def test_get_number_of_occurrences_per_character(self, args):
        service = EvolutionaryRatePerSite(args)
        alignment = _alignment()
        counts = service.get_number_of_occurrences_per_character(alignment, 0, ["-"])
        assert counts["A"] == 2
        assert counts["T"] == 1

    def test_calculate_pic(self, args):
        service = EvolutionaryRatePerSite(args)
        pic = service.calculate_pic({"A": 2, "T": 1})
        assert round(pic, 4) == 0.4444

    def test_calculate_evolutionary_rate_per_site(self, args):
        service = EvolutionaryRatePerSite(args)
        values = service.calculate_evolutionary_rate_per_site(_alignment(), is_protein=False)
        assert len(values) == 4
        assert all(0 <= v <= 1 for v in values)

    def test_run_json_with_plot(self, mocker):
        args = Namespace(alignment="/some/path/to/file.fa", json=True, plot=True, plot_output="erps.png")
        service = EvolutionaryRatePerSite(args)
        alignment = _alignment()
        mocker.patch.object(
            EvolutionaryRatePerSite,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocked_plot = mocker.patch.object(EvolutionaryRatePerSite, "_plot_evolutionary_rate_per_site")
        mocked_json = mocker.patch("phykit.services.alignment.evolutionary_rate_per_site.print_json")

        service.run()

        mocked_plot.assert_called_once()
        payload = mocked_json.call_args.args[0]
        assert payload["rows"] == payload["sites"]
        assert payload["plot_output"] == "erps.png"

    def test_process_args_defaults(self):
        parsed = EvolutionaryRatePerSite(Namespace(alignment="x.fa")).process_args(
            Namespace(alignment="x.fa")
        )
        assert parsed["json_output"] is False
        assert parsed["plot"] is False
        assert parsed["plot_output"] == "evolutionary_rate_per_site_plot.png"

    def test_calculate_evolutionary_rate_per_site_all_gaps(self, args):
        service = EvolutionaryRatePerSite(args)
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("----"), id="t1"),
                SeqRecord(Seq("????"), id="t2"),
            ]
        )
        values = service.calculate_evolutionary_rate_per_site(alignment, is_protein=False)
        assert values == [0, 0, 0, 0]

    def test_plot_evolutionary_rate_per_site_creates_file(self, tmp_path):
        pytest.importorskip("matplotlib")
        out = tmp_path / "erps.png"
        service = EvolutionaryRatePerSite(
            Namespace(alignment="x.fa", plot=True, plot_output=str(out))
        )
        service._plot_evolutionary_rate_per_site(
            [{"site": 1, "evolutionary_rate": 0.2}, {"site": 2, "evolutionary_rate": 0.8}]
        )
        assert out.exists()

    def test_plot_evolutionary_rate_per_site_empty_rows(self):
        pytest.importorskip("matplotlib")
        service = EvolutionaryRatePerSite(Namespace(alignment="x.fa", plot=True))
        service._plot_evolutionary_rate_per_site([])

    def test_plot_evolutionary_rate_per_site_importerror(self, monkeypatch, capsys):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name.startswith("matplotlib"):
                raise ImportError("no matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        service = EvolutionaryRatePerSite(Namespace(alignment="x.fa"))
        with pytest.raises(SystemExit) as exc:
            service._plot_evolutionary_rate_per_site([{"site": 1, "evolutionary_rate": 0.2}])
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "matplotlib is required for --plot in evolutionary_rate_per_site" in out

    def test_run_text_with_plot_message(self, mocker, capsys):
        service = EvolutionaryRatePerSite(
            Namespace(alignment="x.fa", json=False, plot=True, plot_output="erps.png")
        )
        mocker.patch.object(
            EvolutionaryRatePerSite,
            "get_alignment_and_format",
            return_value=(_alignment(), "fasta", False),
        )
        mocker.patch.object(
            EvolutionaryRatePerSite,
            "calculate_evolutionary_rate_per_site",
            return_value=[0.2, 0.4],
        )
        mocked_plot = mocker.patch.object(EvolutionaryRatePerSite, "_plot_evolutionary_rate_per_site")
        service.run()
        out, _ = capsys.readouterr()
        assert "1\t0.2" in out
        assert "2\t0.4" in out
        assert "Saved evolutionary-rate plot: erps.png" in out
        mocked_plot.assert_called_once()
