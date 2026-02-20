from argparse import Namespace

import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.evolutionary_rate_per_site import EvolutionaryRatePerSite


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
