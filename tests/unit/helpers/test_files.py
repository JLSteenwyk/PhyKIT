import time
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from phykit.helpers.files import (
    _FILE_FORMAT_VALUES,
    _cached_alignment_read,
    _cached_detect_format_by_content,
    _detect_format_by_content,
    FileFormat,
    _get_file_hash,
    get_alignment_and_format,
    is_protein_alignment,
    read_single_column_file_to_list,
)
from phykit.errors import PhykitUserError


def test_module_import_does_not_import_biopython_alignment_reader():
    code = """
import sys
import phykit.helpers.files
assert "typing" not in sys.modules
assert "hashlib" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
assert "Bio.Align" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


class TestFileErrorHandling:
    def test_get_alignment_and_format_error_handling(self):
        file_path = "not_real_file"
        with pytest.raises(PhykitUserError) as excinfo:
            get_alignment_and_format(file_path)
        assert excinfo.value.code == 2
        assert "corresponds to no such file" in excinfo.value.messages[0]

    def test_get_read_single_column_file_to_list_error_handling(self):
        file_path = "not_real_file"
        with pytest.raises(PhykitUserError) as excinfo:
            read_single_column_file_to_list(file_path)
        assert excinfo.value.code == 2
        assert "corresponds to no such file or directory" in excinfo.value.messages[0]


class TestFormatDetection:
    def test_cached_file_format_values_match_file_format_enum(self):
        assert _FILE_FORMAT_VALUES == tuple(
            file_format.value for file_format in FileFormat
        )

    def test_detect_format_fasta(self, tmp_path: Path):
        aln = tmp_path / "test.fa"
        aln.write_text(">a\nACGT\n>b\nACGT\n")
        assert _detect_format_by_content(str(aln)) == "fasta"

    def test_detect_format_fasta_with_leading_whitespace(self, tmp_path: Path):
        aln = tmp_path / "test.fa"
        aln.write_text("  >a\nACGT\n")
        assert _detect_format_by_content(str(aln)) == "fasta"

    def test_detect_format_clustal(self, tmp_path: Path):
        aln = tmp_path / "test.aln"
        aln.write_text("CLUSTAL W(1.82)\n\n")
        assert _detect_format_by_content(str(aln)) == "clustal"

    def test_detect_format_stockholm(self, tmp_path: Path):
        aln = tmp_path / "test.sto"
        aln.write_text("# STOCKHOLM 1.0\n")
        assert _detect_format_by_content(str(aln)) == "stockholm"

    def test_detect_format_phylip_like_header(self, tmp_path: Path):
        aln = tmp_path / "test.phy"
        aln.write_text("2 4\n")
        assert _detect_format_by_content(str(aln)) == "phylip"

    def test_detect_format_rejects_phylip_header_with_extra_columns(
        self,
        tmp_path: Path,
    ):
        aln = tmp_path / "test.phy"
        aln.write_text("2 4 extra\n")
        assert _detect_format_by_content(str(aln)) is None

    def test_detect_format_digit_start_non_phylip_header(self, tmp_path: Path):
        aln = tmp_path / "test.phy"
        aln.write_text("2 taxa 4 sites\n")
        assert _detect_format_by_content(str(aln)) is None

    def test_detect_format_unknown(self, tmp_path: Path):
        aln = tmp_path / "test.txt"
        aln.write_text("hello world\n")
        assert _detect_format_by_content(str(aln)) is None

    def test_detect_format_unknown_non_digit_header_avoids_split(self, mocker):
        class NoSplitLine(str):
            def strip(self):
                return self

            def split(self, *_args, **_kwargs):
                raise AssertionError("non-digit unknown headers should not split")

        class FakeFile:
            def __enter__(self):
                return self

            def __exit__(self, *_args):
                return None

            def readline(self):
                return NoSplitLine("hello world\n")

        mocker.patch("builtins.open", return_value=FakeFile())

        assert _detect_format_by_content("unknown.txt") is None


class TestAlignmentReadAndType:
    def test_get_alignment_and_format_reads_fasta_nucleotide(self, tmp_path: Path):
        aln = tmp_path / "nucl.fa"
        aln.write_text(">a\nACGTN-\n>b\nACGTN-\n")
        alignment, fmt, is_protein = get_alignment_and_format(str(aln))
        assert fmt == "fasta"
        assert len(alignment) == 2
        assert alignment.get_alignment_length() == 6
        assert is_protein is False

    def test_get_alignment_and_format_reads_fasta_protein(self, tmp_path: Path):
        aln = tmp_path / "prot.fa"
        aln.write_text(">a\nMSTV\n>b\nMSTL\n")
        alignment, fmt, is_protein = get_alignment_and_format(str(aln))
        assert fmt == "fasta"
        assert len(alignment) == 2
        assert alignment.get_alignment_length() == 4
        assert is_protein is True

    def test_get_alignment_and_format_reuses_cached_format_detection(
        self,
        tmp_path: Path,
        mocker,
    ):
        aln = tmp_path / "test.fa"
        aln.write_text(">a\nACGT\n>b\nACGT\n")
        _cached_alignment_read.cache_clear()
        _cached_detect_format_by_content.cache_clear()
        detect = mocker.spy(
            sys.modules["phykit.helpers.files"],
            "_detect_format_by_content",
        )

        first_alignment, first_format, first_is_protein = get_alignment_and_format(
            str(aln)
        )
        second_alignment, second_format, second_is_protein = get_alignment_and_format(
            str(aln)
        )

        assert first_alignment is second_alignment
        assert first_format == second_format == "fasta"
        assert first_is_protein is second_is_protein is False
        detect.assert_called_once_with(str(aln))

    def test_get_alignment_and_format_unknown_format_exits(self, tmp_path: Path):
        bad = tmp_path / "bad.aln"
        bad.write_text("not-an-alignment\nstill-not-an-alignment\n")
        with pytest.raises(PhykitUserError) as excinfo:
            get_alignment_and_format(str(bad))
        assert excinfo.value.code == 2

    def test_is_protein_alignment_handles_nucl_and_protein(self, tmp_path: Path):
        nucl = tmp_path / "n.fa"
        nucl.write_text(">a\nACGTN\n>b\nACGT-\n")
        prot = tmp_path / "p.fa"
        prot.write_text(">a\nMSTV\n>b\nMSTL\n")
        nucl_alignment, _, _ = get_alignment_and_format(str(nucl))
        prot_alignment, _, _ = get_alignment_and_format(str(prot))
        assert is_protein_alignment(nucl_alignment) is False
        assert is_protein_alignment(prot_alignment) is True

    def test_is_protein_alignment_ascii_fast_path_handles_lowercase(self):
        alignment = [
            SimpleNamespace(seq="acgtun?-*"),
            SimpleNamespace(seq="ACGTUN?-*"),
        ]

        assert is_protein_alignment(alignment) is False

        protein_alignment = alignment + [SimpleNamespace(seq="acgtM")]

        assert is_protein_alignment(protein_alignment) is True

    def test_is_protein_alignment_unicode_fallback(self):
        alignment = [SimpleNamespace(seq="ACGT\u03a9")]

        assert is_protein_alignment(alignment) is True


class TestHelpers:
    def test_get_file_hash_changes_when_file_changes(self, tmp_path: Path):
        fp = tmp_path / "x.fa"
        fp.write_text(">a\nACGT\n")
        first_hash = _get_file_hash(str(fp))
        first_stat = fp.stat()
        assert first_hash == f"{fp}_{first_stat.st_size}_{first_stat.st_mtime_ns}"

        # Ensure mtime changes on fast filesystems.
        time.sleep(1.1)
        fp.write_text(">a\nACGTA\n")
        second_hash = _get_file_hash(str(fp))
        assert first_hash != second_hash

    def test_read_single_column_file_to_list_trims_whitespace(self, tmp_path: Path):
        fp = tmp_path / "items.txt"
        fp.write_text("  a  \n b\nc \n")
        assert read_single_column_file_to_list(str(fp)) == ["a", "b", "c"]

    def test_read_single_column_file_to_list_strips_tabs_and_missing_final_newline(
        self,
        tmp_path: Path,
    ):
        fp = tmp_path / "items.txt"
        fp.write_text("\ta\t\n b \n\tc")
        assert read_single_column_file_to_list(str(fp)) == ["a", "b", "c"]

    def test_read_single_column_file_to_list_streams_and_preserves_blank_lines(
        self,
        mocker,
    ):
        class StreamingOnlyFile:
            def __enter__(self):
                return self

            def __exit__(self, *_args):
                return None

            def __iter__(self):
                return iter(["  a  \n", "\n", "\tb\t"])

        mocker.patch("builtins.open", return_value=StreamingOnlyFile())

        assert read_single_column_file_to_list("items.txt") == ["a", "", "b"]
