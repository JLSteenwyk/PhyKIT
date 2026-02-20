from argparse import Namespace

import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.gc_content import GCContent


def _alignment(records):
    return MultipleSeqAlignment(records)


@pytest.fixture
def args():
    return Namespace(fasta="/some/path/to/file.fa", verbose=False)


class TestGCContent:
    def test_init_sets_expected_attrs(self, args):
        service = GCContent(args)
        assert service.fasta == args.fasta
        assert service.verbose is False
        assert service.json_output is False

    def test_calculate_gc_total_value(self, args):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("ACGT"), id="a"),
                SeqRecord(Seq("GGTT"), id="b"),
            ]
        )
        assert service.calculate_gc_total_value(records, is_protein=False) == 0.5

    def test_calculate_gc_total_value_exits_for_empty_cleaned_seq(self, mocker, args):
        service = GCContent(args)
        records = _alignment([SeqRecord(Seq("----"), id="a")])
        mocked_print = mocker.patch("builtins.print")
        with pytest.raises(SystemExit) as excinfo:
            service.calculate_gc_total_value(records, is_protein=False)
        assert excinfo.value.code == 2
        mocked_print.assert_called_once()

    def test_calculate_gc_per_sequence_data(self, args):
        service = GCContent(args)
        records = _alignment(
            [
                SeqRecord(Seq("G-C"), id="a"),  # cleaned GC = 2/2
                SeqRecord(Seq("A-T"), id="b"),  # cleaned GC = 0/2
            ]
        )
        assert service.calculate_gc_per_sequence_data(records, is_protein=False) == [
            ("a", 1.0),
            ("b", 0.0),
        ]

    def test_run_json_summary(self, mocker):
        args = Namespace(fasta="/some/path/to/file.fa", verbose=False, json=True)
        records = _alignment([SeqRecord(Seq("ACGT"), id="a")])
        mocker.patch("phykit.services.alignment.gc_content.get_alignment_and_format", return_value=(records, "fasta", False))
        mocked_json = mocker.patch("phykit.services.alignment.gc_content.print_json")
        service = GCContent(args)
        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload == {"verbose": False, "gc_content": 0.5}

    def test_run_json_verbose(self, mocker):
        args = Namespace(fasta="/some/path/to/file.fa", verbose=True, json=True)
        records = _alignment([SeqRecord(Seq("G-C"), id="a")])
        mocker.patch("phykit.services.alignment.gc_content.get_alignment_and_format", return_value=(records, "fasta", False))
        mocked_json = mocker.patch("phykit.services.alignment.gc_content.print_json")
        service = GCContent(args)
        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload["verbose"] is True
        assert payload["rows"] == [{"taxon": "a", "gc_content": 1.0}]
        assert payload["sequences"] == payload["rows"]

    def test_run_exits_for_protein_alignment(self, mocker):
        args = Namespace(fasta="/some/path/to/file.fa", verbose=False, json=False)
        records = _alignment([SeqRecord(Seq("MSTV"), id="a")])
        mocker.patch("phykit.services.alignment.gc_content.get_alignment_and_format", return_value=(records, "fasta", True))
        mocked_print = mocker.patch("builtins.print")
        service = GCContent(args)
        with pytest.raises(SystemExit) as excinfo:
            service.run()
        assert excinfo.value.code == 2
        mocked_print.assert_called_with("GC content can't be calculated for protein sequences")
