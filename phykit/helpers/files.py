from __future__ import annotations

from enum import Enum
from functools import lru_cache
import os

from ..errors import PhykitUserError


_NUCLEOTIDE_CHARS = {
    "A", "C", "G", "T", "U", "-", "N", "?", "*"
}
_NUCLEOTIDE_BYTES = b"ACGTUacgtu-Nn?*"


class FileFormat(Enum):
    fasta = "fasta"
    clustal = "clustal"
    maf = "maf"
    mauve = "mauve"
    phylip = "phylip"
    phylip_seq = "phylip-sequential"
    phylip_rel = "phylip-relaxed"
    stockholm = "stockholm"


_FILE_FORMAT_VALUES = tuple(file_format.value for file_format in FileFormat)


def _get_file_hash(file_path: str) -> str:
    """Calculate a hash for file content to use as cache key."""
    # Use file path, size, and modification time for cache key
    # This is faster than hashing file contents
    stat = os.stat(file_path)
    return f"{file_path}_{stat.st_size}_{stat.st_mtime_ns}"

def _detect_format_by_content(file_path: str) -> str | None:
    """Attempt to detect file format by examining file content."""
    with open(file_path) as f:
        first_line = f.readline().strip()
        first_char = first_line[:1]

        # Quick format detection based on first line
        if first_char == '>':
            return 'fasta'
        elif first_line.startswith('CLUSTAL'):
            return 'clustal'
        elif first_char == '#':
            # Could be Stockholm
            if 'STOCKHOLM' in first_line:
                return 'stockholm'
        elif first_line.isdigit():
            return 'phylip'
        elif first_char and "0" <= first_char <= "9":
            parts = first_line.split(None, 2)
            if len(parts) == 2 and parts[0].isdigit():
                return 'phylip'

    return None

@lru_cache(maxsize=32)
def _cached_detect_format_by_content(
    file_hash: str,
    file_path: str,
) -> str | None:
    return _detect_format_by_content(file_path)

@lru_cache(maxsize=32)
def _cached_alignment_read(
    file_hash: str,
    file_path: str,
    file_format: str,
) -> tuple["MultipleSeqAlignment", bool]:
    """Cached reading of alignment files."""
    from Bio import AlignIO

    with open(file_path) as f:
        alignment = AlignIO.read(f, file_format)
    return alignment, is_protein_alignment(alignment)

def get_alignment_and_format(
    alignment_file_path: str
) -> tuple["MultipleSeqAlignment", str, bool]:
    try:
        file_hash = _get_file_hash(alignment_file_path)
    except FileNotFoundError:
        raise PhykitUserError(
            [
                f"{alignment_file_path} corresponds to no such file.",
                "Please check file name and pathing",
            ],
            code=2,
        )

    # Try to detect format by content first
    detected_format = _cached_detect_format_by_content(
        file_hash,
        alignment_file_path,
    )

    # If format was detected, try it first
    if detected_format:
        try:
            alignment, is_protein = _cached_alignment_read(
                file_hash, alignment_file_path, detected_format
            )
            return alignment, detected_format, is_protein
        except (ValueError, AssertionError):
            pass

    # Fall back to trying all formats
    for file_format in _FILE_FORMAT_VALUES:
        # Skip the already tried format
        if detected_format and file_format == detected_format:
            continue

        try:
            alignment, is_protein = _cached_alignment_read(
                file_hash, alignment_file_path, file_format
            )
            return alignment, file_format, is_protein
        except (ValueError, AssertionError):
            continue

    # If we get here, no format worked
    raise PhykitUserError(
        [
            f"Could not determine format for {alignment_file_path}",
            "Please ensure the file is in a supported format",
        ],
        code=2,
    )


def is_protein_alignment(alignment: "MultipleSeqAlignment") -> bool:
    for record in alignment:
        sequence = str(record.seq)
        try:
            has_non_nucleotide = bool(
                sequence.encode("ascii").translate(None, _NUCLEOTIDE_BYTES)
            )
        except UnicodeEncodeError:
            seq_set = set(sequence.upper())
            has_non_nucleotide = bool(seq_set - _NUCLEOTIDE_CHARS)
        if has_non_nucleotide:
            return True

    return False


def read_single_column_file_to_list(single_col_file_path: str) -> list:
    try:
        with open(single_col_file_path) as f:
            return [line.strip() for line in f.read().splitlines()]
    except FileNotFoundError:
        raise PhykitUserError(
            [
                f"{single_col_file_path} corresponds to no such file or directory.",
                "Please check file name and pathing",
            ],
            code=2,
        )
