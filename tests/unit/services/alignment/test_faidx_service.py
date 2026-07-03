from argparse import Namespace
import subprocess
import sys

import pytest
from phykit.services.alignment.faidx import Faidx


@pytest.fixture
def args():
    return Namespace(fasta="/some/path/to/file.fa", entry="1,2")


class TestFaidx:
    def test_module_import_does_not_import_biopython_fasta_parser(self):
        code = """
import sys
import phykit.services.alignment.faidx
assert "typing" not in sys.modules
assert "Bio.SeqIO" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_init_sets_expected_attrs(self, args):
        service = Faidx(args)
        assert service.fasta == args.fasta
        assert service.entry == args.entry
        assert service.json_output is False

    def test_run_prints_requested_entries(self, mocker, tmp_path):
        path = tmp_path / "alignment.fa"
        path.write_text(">1 description\nA- GT\nAT\n>2 second\nA-G-AT\n")
        args = Namespace(fasta=str(path), entry="1,2")
        mocked_print = mocker.patch("builtins.print")

        service = Faidx(args)
        service.run()

        mocked_print.assert_called_once_with(">1\nA-GTAT\n>2\nA-G-AT")

    def test_run_text_output_preserves_requested_entry_order(self, mocker):
        args = Namespace(fasta="/some/path/to/file.fa", entry="2,1")
        mocker.patch.object(
            Faidx,
            "_fetch_entries",
            return_value={"1": "AAAA", "2": "CCCC"},
        )
        mocked_print = mocker.patch("builtins.print")

        service = Faidx(args)
        service.run()

        mocked_print.assert_called_once_with(">2\nCCCC\n>1\nAAAA")

    def test_run_json_prints_structured_payload(self, mocker, tmp_path):
        path = tmp_path / "alignment.fa"
        path.write_text(">1 description\nA-GTAT\n>2 second\nA-G-AT\n")
        args = Namespace(fasta=str(path), entry="1, 2 ,", json=True)
        mocked_json = mocker.patch("phykit.services.alignment.faidx.print_json")

        service = Faidx(args)
        service.run()

        payload = mocked_json.call_args.args[0]
        assert payload["rows"][0] == {"entry": "1", "name": "1", "sequence": "A-GTAT"}
        assert payload["rows"][1] == {"entry": "2", "name": "2", "sequence": "A-G-AT"}
        assert payload["entries"] == payload["rows"]

    def test_fetch_entries_returns_plain_sequences(self, tmp_path):
        path = tmp_path / "alignment.fa"
        path.write_text(">1 description\nA-GTAT\n>2 second\nA-G-AT\n")

        records = Faidx._fetch_entries(str(path), ["1", "2"])

        assert records == {"1": "A-GTAT", "2": "A-G-AT"}

    def test_fetch_entries_raises_for_missing_entry(self, tmp_path):
        path = tmp_path / "alignment.fa"
        path.write_text(">1\nA-GTAT\n")

        with pytest.raises(KeyError):
            Faidx._fetch_entries(str(path), ["2"])

    def test_fetch_entries_raises_for_duplicate_ids(self, tmp_path):
        path = tmp_path / "alignment.fa"
        path.write_text(">1\nA-GTAT\n>1 duplicate\nA-G-AT\n")

        with pytest.raises(ValueError, match="Duplicate key"):
            Faidx._fetch_entries(str(path), ["1"])

    def test_fetch_entries_scans_later_duplicate_after_requested_entry(self, tmp_path):
        path = tmp_path / "alignment.fa"
        path.write_text(">1\nA-GTAT\n>2\nA-G-AT\n>1 duplicate\nTTTT\n")

        with pytest.raises(ValueError, match="Duplicate key"):
            Faidx._fetch_entries(str(path), ["1"])
