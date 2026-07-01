import pytest
import subprocess
import sys
from argparse import Namespace
from mock import call

from phykit.services.alignment.alignment_length import AlignmentLength


def test_fasta_lightweight_modules_defer_heavy_imports():
    modules = [
        "phykit.services.alignment.alignment_length",
        "phykit.services.alignment.faidx",
        "phykit.services.alignment.rename_fasta_entries",
        "phykit.services.alignment.alignment_subsample",
    ]
    code = f"""
import sys
modules = {modules!r}
for module_name in modules:
    module = __import__(module_name, fromlist=["*"])
    assert callable(module.print_json)

assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "Bio.SeqIO" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_alignment_length_import_does_not_import_typing():
    code = """
import sys
import phykit.services.alignment.alignment_length

assert "typing" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def args():
    kwargs = dict(alignment="/some/path/to/file.fa")
    return Namespace(**kwargs)


class TestAlignmentLength(object):
    def test_init_sets_alignment_file_path(self, args):
        aln = AlignmentLength(args)
        assert aln.alignment_file_path == args.alignment
        assert aln.output_file_path is None

    def test_alignment_length_is_printed(self, mocker, args):
        expected_length = "6"
        aln = mocker.MagicMock(
            get_alignment_length=mocker.MagicMock(return_value=expected_length)
        )
        mocked_print = mocker.patch("builtins.print")
        mocked_get_alignment_and_format = mocker.patch(
            "phykit.services.alignment.alignment_length.AlignmentLength.get_alignment_and_format",
            return_value=(aln, '', '')
        )
        aln_len = AlignmentLength(args)
        _ = aln_len.run()

        assert mocked_get_alignment_and_format.called
        assert mocked_print.mock_calls == [
            call(expected_length)
        ]

    def test_alignment_length_fast_fasta_path_skips_full_parser(self, mocker, tmp_path):
        path = tmp_path / "alignment.fa"
        path.write_text(">a\nACGT\n>b\nTGCA\n")
        args = Namespace(alignment=str(path))
        mocked_get_alignment_and_format = mocker.patch(
            "phykit.services.alignment.alignment_length.AlignmentLength.get_alignment_and_format",
            side_effect=AssertionError("valid FASTA should use the length fast path"),
        )
        mocked_print = mocker.patch("builtins.print")

        AlignmentLength(args).run()

        mocked_get_alignment_and_format.assert_not_called()
        mocked_print.assert_called_once_with(4)

    def test_alignment_length_fast_fasta_path_counts_wrapped_spaced_sequences(
        self, mocker, tmp_path
    ):
        path = tmp_path / "alignment.fa"
        path.write_text(">a description\nAC GT\nA\n>b\nTG CA\nT\n")
        args = Namespace(alignment=str(path))
        mocked_get_alignment_and_format = mocker.patch(
            "phykit.services.alignment.alignment_length.AlignmentLength.get_alignment_and_format",
            side_effect=AssertionError("valid FASTA should use the length fast path"),
        )
        mocked_print = mocker.patch("builtins.print")

        AlignmentLength(args).run()

        mocked_get_alignment_and_format.assert_not_called()
        mocked_print.assert_called_once_with(5)

    def test_alignment_length_fast_fasta_path_strips_trailing_whitespace(
        self, mocker, tmp_path
    ):
        path = tmp_path / "alignment.fa"
        path.write_text(">a\nACGT \t\n>b\nTGCA\t\n")
        args = Namespace(alignment=str(path))
        mocked_get_alignment_and_format = mocker.patch(
            "phykit.services.alignment.alignment_length.AlignmentLength.get_alignment_and_format",
            side_effect=AssertionError("valid FASTA should use the length fast path"),
        )
        mocked_print = mocker.patch("builtins.print")

        AlignmentLength(args).run()

        mocked_get_alignment_and_format.assert_not_called()
        mocked_print.assert_called_once_with(4)

    def test_alignment_length_inconsistent_fasta_falls_back(self, mocker, tmp_path):
        path = tmp_path / "alignment.fa"
        path.write_text(">a\nACGT\n>b\nTGCAA\n")
        args = Namespace(alignment=str(path))
        aln = mocker.MagicMock(
            get_alignment_length=mocker.MagicMock(return_value=5)
        )
        mocked_get_alignment_and_format = mocker.patch(
            "phykit.services.alignment.alignment_length.AlignmentLength.get_alignment_and_format",
            return_value=(aln, "fasta", False),
        )
        mocked_print = mocker.patch("builtins.print")

        AlignmentLength(args).run()

        mocked_get_alignment_and_format.assert_called_once()
        mocked_print.assert_called_once_with(5)

    def test_alignment_length_non_ascii_fasta_falls_back(self, mocker, tmp_path):
        path = tmp_path / "alignment.fa"
        path.write_text(">a\nAα\n>b\nTβ\n")
        args = Namespace(alignment=str(path))
        aln = mocker.MagicMock(
            get_alignment_length=mocker.MagicMock(return_value=2)
        )
        mocked_get_alignment_and_format = mocker.patch(
            "phykit.services.alignment.alignment_length.AlignmentLength.get_alignment_and_format",
            return_value=(aln, "fasta", False),
        )
        mocked_print = mocker.patch("builtins.print")

        AlignmentLength(args).run()

        mocked_get_alignment_and_format.assert_called_once()
        mocked_print.assert_called_once_with(2)
