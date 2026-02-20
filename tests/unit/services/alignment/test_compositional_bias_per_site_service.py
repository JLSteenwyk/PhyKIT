from argparse import Namespace

import numpy as np
import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.compositional_bias_per_site import (
    CompositionalBiasPerSite,
)
import phykit.services.alignment.compositional_bias_per_site as cbps_module
import builtins


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


class TestCompositionalBiasPerSite:
    def test_init_sets_expected_attrs(self, args):
        service = CompositionalBiasPerSite(args)
        assert service.alignment_file_path == args.alignment
        assert service.json_output is False
        assert service.plot is False
        assert service.plot_output == "compositional_bias_per_site_plot.png"

    def test_get_number_of_occurrences_per_character(self, args):
        service = CompositionalBiasPerSite(args)
        counts = service.get_number_of_occurrences_per_character(_alignment(), 0, is_protein=False)
        assert sorted(counts) == [1, 2]

    def test_calculate_compositional_bias_per_site(self, args):
        service = CompositionalBiasPerSite(args)
        stat_res, corrected = service.calculate_compositional_bias_per_site(_alignment(), is_protein=False)
        assert len(stat_res) == 4
        assert len(corrected) == 4

    def test_build_rows_handles_nan(self, args):
        service = CompositionalBiasPerSite(args)
        stat_res = [
            type("R", (), {"statistic": 0.0, "pvalue": np.nan})(),
            type("R", (), {"statistic": 0.2, "pvalue": 0.6547})(),
        ]
        rows = service._build_rows(stat_res, ["nan", 1.0])
        assert rows[0] == {
            "site": 1,
            "chi_square": 0.0,
            "p_value_corrected": None,
            "p_value": None,
        }
        assert rows[1] == {
            "site": 2,
            "chi_square": 0.2,
            "p_value_corrected": 1.0,
            "p_value": 0.6547,
        }

    def test_run_json_with_plot(self, mocker):
        args = Namespace(alignment="/some/path/to/file.fa", json=True, plot=True, plot_output="cbps.png")
        service = CompositionalBiasPerSite(args)
        mocker.patch.object(
            CompositionalBiasPerSite,
            "get_alignment_and_format",
            return_value=(_alignment(), "fasta", False),
        )
        mocker.patch.object(
            CompositionalBiasPerSite,
            "calculate_compositional_bias_per_site",
            return_value=(
                [type("R", (), {"statistic": 0.2, "pvalue": 0.6547})()],
                [1.0],
            ),
        )
        mocked_plot = mocker.patch.object(CompositionalBiasPerSite, "_plot_compositional_bias_manhattan")
        mocked_json = mocker.patch("phykit.services.alignment.compositional_bias_per_site.print_json")

        service.run()

        mocked_plot.assert_called_once()
        payload = mocked_json.call_args.args[0]
        assert payload["rows"] == payload["sites"]
        assert payload["plot_output"] == "cbps.png"

    def test_process_args_defaults(self):
        parsed = CompositionalBiasPerSite(Namespace(alignment="x.fa")).process_args(
            Namespace(alignment="x.fa")
        )
        assert parsed["json_output"] is False
        assert parsed["plot"] is False
        assert parsed["plot_output"] == "compositional_bias_per_site_plot.png"

    def test_calculate_compositional_bias_per_site_all_gaps(self):
        service = CompositionalBiasPerSite(Namespace(alignment="x.fa"))
        alignment = MultipleSeqAlignment(
            [
                SeqRecord(Seq("----"), id="t1"),
                SeqRecord(Seq("????"), id="t2"),
            ]
        )
        stat_res, corrected = service.calculate_compositional_bias_per_site(alignment, is_protein=False)
        assert len(stat_res) == 4
        assert corrected == ["nan", "nan", "nan", "nan"]

    def test_plot_compositional_bias_manhattan_creates_file(self, tmp_path):
        pytest.importorskip("matplotlib")
        output = tmp_path / "cbps_plot.png"
        service = CompositionalBiasPerSite(
            Namespace(alignment="x.fa", plot=True, plot_output=str(output))
        )
        rows = [
            {"site": 1, "chi_square": 1.0, "p_value_corrected": 0.01, "p_value": 0.01},
            {"site": 2, "chi_square": 0.5, "p_value_corrected": 0.5, "p_value": 0.5},
            {"site": 3, "chi_square": 0.2, "p_value_corrected": None, "p_value": None},
        ]
        service._plot_compositional_bias_manhattan(rows)
        assert output.exists()

    def test_plot_compositional_bias_manhattan_importerror(self, monkeypatch, capsys):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name.startswith("matplotlib"):
                raise ImportError("no matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        service = CompositionalBiasPerSite(Namespace(alignment="x.fa"))
        with pytest.raises(SystemExit) as exc:
            service._plot_compositional_bias_manhattan([{"site": 1, "chi_square": 0.0, "p_value_corrected": 0.5, "p_value": 0.5}])
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "matplotlib is required for --plot in compositional_bias_per_site" in out

    def test_run_text_output_with_plot_message(self, mocker, capsys):
        service = CompositionalBiasPerSite(
            Namespace(alignment="x.fa", json=False, plot=True, plot_output="cbps.png")
        )
        mocker.patch.object(
            CompositionalBiasPerSite, "get_alignment_and_format", return_value=(_alignment(), "fasta", False)
        )
        mocker.patch.object(
            CompositionalBiasPerSite,
            "calculate_compositional_bias_per_site",
            return_value=(
                [type("R", (), {"statistic": 2.0, "pvalue": 0.0456})()],
                [0.1234],
            ),
        )
        mocker.patch.object(CompositionalBiasPerSite, "_plot_compositional_bias_manhattan")

        service.run()

        out, _ = capsys.readouterr()
        assert "1\t2.0\t0.1234\t0.0456" in out
        assert "Saved compositional bias plot: cbps.png" in out
