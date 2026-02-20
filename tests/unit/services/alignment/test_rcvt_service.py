from argparse import Namespace

import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.rcvt import RelativeCompositionVariabilityTaxon


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
