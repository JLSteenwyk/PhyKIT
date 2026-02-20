import time
from pathlib import Path

import pytest

from phykit.helpers.files import (
    _detect_format_by_content,
    _get_file_hash,
    get_alignment_and_format,
    is_protein_alignment,
    read_single_column_file_to_list,
)


class TestFileErrorHandling:
    def test_get_alignment_and_format_error_handling(self):
        file_path = "not_real_file"
        with pytest.raises(SystemExit) as excinfo:
            get_alignment_and_format(file_path)
        assert excinfo.type is SystemExit

    def test_get_read_single_column_file_to_list_error_handling(self):
        file_path = "not_real_file"
        with pytest.raises(SystemExit) as excinfo:
            read_single_column_file_to_list(file_path)
        assert excinfo.type is SystemExit


class TestFormatDetection:
    def test_detect_format_fasta(self, tmp_path: Path):
        aln = tmp_path / "test.fa"
        aln.write_text(">a\nACGT\n>b\nACGT\n")
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

    def test_detect_format_unknown(self, tmp_path: Path):
        aln = tmp_path / "test.txt"
        aln.write_text("hello world\n")
        assert _detect_format_by_content(str(aln)) is None


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

    def test_get_alignment_and_format_unknown_format_exits(self, tmp_path: Path):
        bad = tmp_path / "bad.aln"
        bad.write_text("not-an-alignment\nstill-not-an-alignment\n")
        with pytest.raises(SystemExit) as excinfo:
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


class TestHelpers:
    def test_get_file_hash_changes_when_file_changes(self, tmp_path: Path):
        fp = tmp_path / "x.fa"
        fp.write_text(">a\nACGT\n")
        first_hash = _get_file_hash(str(fp))

        # Ensure mtime changes on fast filesystems.
        time.sleep(1.1)
        fp.write_text(">a\nACGTA\n")
        second_hash = _get_file_hash(str(fp))
        assert first_hash != second_hash

    def test_read_single_column_file_to_list_trims_whitespace(self, tmp_path: Path):
        fp = tmp_path / "items.txt"
        fp.write_text("  a  \n b\nc \n")
        assert read_single_column_file_to_list(str(fp)) == ["a", "b", "c"]
