from argparse import Namespace

import numpy as np
import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.compositional_bias_per_site import (
    CompositionalBiasPerSite,
)


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
