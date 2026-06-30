from argparse import Namespace
import subprocess
import sys

import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.alignment_recoding import AlignmentRecoding


def test_module_import_does_not_import_bio_align_or_json_helpers():
    code = """
import sys
import phykit.services.alignment.alignment_recoding as module
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


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

    def test_read_recoding_table_ignores_extra_columns_and_whitespace(self, tmp_path, args):
        table = tmp_path / "custom.txt"
        table.write_text(" r   a ignored\nY\tc\ttrailing\n")
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

    def test_recode_alignment_recodes_lowercase_and_preserves_ambiguous_gaps(self, args):
        service = AlignmentRecoding(args)
        alignment = MultipleSeqAlignment(
            [SeqRecord(Seq("acgtnNxX?-"), id="t1")]
        )
        recoded = service.recode_alignment(
            alignment,
            recoding_table={"A": "R", "C": "Y", "G": "R", "T": "Y"},
            is_protein=False,
        )
        assert recoded["t1"] == ["R", "Y", "R", "Y", "n", "N", "x", "X", "?", "-"]

    def test_build_translation_table_skips_uppercase_and_lowercase_gaps(self):
        table = AlignmentRecoding._build_translation_table(
            {"A": "R", "N": "Y", "X": "Z"},
            ["N", "n", "X", "x"],
        )

        assert table == {ord("A"): "R", ord("a"): "R"}

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

    def test_run_text_output(self, mocker):
        args = Namespace(alignment="/some/path/to/file.fa", code="RY-nucleotide")
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
        mocked_print = mocker.patch("builtins.print")

        service.run()

        mocked_print.assert_called_once_with(">t1\nRY")

    def test_run_uses_string_recoding_path(self, mocker):
        args = Namespace(alignment="/some/path/to/file.fa", code="RY-nucleotide")
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
        mocker.patch.object(
            service,
            "recode_alignment",
            side_effect=AssertionError("run should avoid list materialization"),
        )
        mocked_print = mocker.patch("builtins.print")

        service.run()

        mocked_print.assert_called_once_with(">t1\nRY")
