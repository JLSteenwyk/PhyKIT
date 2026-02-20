from argparse import Namespace

import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.alignment_recoding import AlignmentRecoding


@pytest.fixture
def args():
    return Namespace(alignment="/some/path/to/file.fa", code="RY-nucleotide")


class TestAlignmentRecoding:
    def test_init_sets_expected_attrs(self, args):
        service = AlignmentRecoding(args)
        assert service.alignment_file_path == args.alignment
        assert service.code == "RY-nucleotide"
        assert service.json_output is False

    def test_read_recoding_table_requires_code(self, args):
        service = AlignmentRecoding(args)
        with pytest.raises(SystemExit) as excinfo:
            service.read_recoding_table(None)
        assert excinfo.value.code == 2

    def test_read_recoding_table_from_custom_file(self, tmp_path, args):
        table = tmp_path / "custom.txt"
        table.write_text("R A\nY C\n")
        service = AlignmentRecoding(args)
        recoding = service.read_recoding_table(str(table))
        assert recoding == {"A": "R", "C": "Y"}

    def test_recode_alignment_preserves_gap_chars(self, args):
        service = AlignmentRecoding(args)
        alignment = MultipleSeqAlignment(
            [SeqRecord(Seq("A-C"), id="t1"), SeqRecord(Seq("C-A"), id="t2")]
        )
        recoded = service.recode_alignment(
            alignment,
            recoding_table={"A": "R", "C": "Y"},
            is_protein=False,
        )
        assert recoded["t1"] == ["R", "-", "Y"]
        assert recoded["t2"] == ["Y", "-", "R"]

    def test_run_json_output(self, mocker):
        args = Namespace(alignment="/some/path/to/file.fa", code="RY-nucleotide", json=True)
        service = AlignmentRecoding(args)
        alignment = MultipleSeqAlignment([SeqRecord(Seq("AC"), id="t1")])
        mocker.patch.object(
            AlignmentRecoding,
            "get_alignment_and_format",
            return_value=(alignment, "fasta", False),
        )
        mocker.patch.object(
            AlignmentRecoding,
            "read_recoding_table",
            return_value={"A": "R", "C": "Y"},
        )
        mocked_json = mocker.patch(
            "phykit.services.alignment.alignment_recoding.print_json"
        )

        service.run()

        payload = mocked_json.call_args.args[0]
        assert payload["code"] == "RY-nucleotide"
        assert payload["taxa"] == [{"taxon": "t1", "sequence": "RY"}]
