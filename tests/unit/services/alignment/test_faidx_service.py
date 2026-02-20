from argparse import Namespace

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.faidx import Faidx


@pytest.fixture
def args():
    return Namespace(fasta="/some/path/to/file.fa", entry="1,2")


class TestFaidx:
    def test_init_sets_expected_attrs(self, args):
        service = Faidx(args)
        assert service.fasta == args.fasta
        assert service.entry == args.entry
        assert service.json_output is False

    def test_run_prints_requested_entries(self, mocker, args):
        records = {
            "1": SeqRecord(Seq("A-GTAT"), id="1", name="1"),
            "2": SeqRecord(Seq("A-G-AT"), id="2", name="2"),
        }
        mocker.patch("phykit.services.alignment.faidx.SeqIO.index", return_value=records)
        mocked_print = mocker.patch("builtins.print")

        service = Faidx(args)
        service.run()

        assert mocked_print.call_count == 2
        mocked_print.assert_any_call(">1\nA-GTAT")
        mocked_print.assert_any_call(">2\nA-G-AT")

    def test_run_json_prints_structured_payload(self, mocker):
        args = Namespace(fasta="/some/path/to/file.fa", entry="1, 2 ,", json=True)
        records = {
            "1": SeqRecord(Seq("A-GTAT"), id="1", name="1"),
            "2": SeqRecord(Seq("A-G-AT"), id="2", name="2"),
        }
        mocker.patch("phykit.services.alignment.faidx.SeqIO.index", return_value=records)
        mocked_json = mocker.patch("phykit.services.alignment.faidx.print_json")

        service = Faidx(args)
        service.run()

        payload = mocked_json.call_args.args[0]
        assert payload["rows"][0] == {"entry": "1", "name": "1", "sequence": "A-GTAT"}
        assert payload["rows"][1] == {"entry": "2", "name": "2", "sequence": "A-G-AT"}
        assert payload["entries"] == payload["rows"]
