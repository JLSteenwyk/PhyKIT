from enum import Enum
import sys
from typing import Tuple, Optional
from functools import lru_cache
import hashlib
import os

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


class FileFormat(Enum):
    fasta = "fasta"
    clustal = "clustal"
    maf = "maf"
    mauve = "mauve"
    phylip = "phylip"
    phylip_seq = "phylip-sequential"
    phylip_rel = "phylip-relaxed"
    stockholm = "stockholm"


def _get_file_hash(file_path: str) -> str:
    """Calculate a hash for file content to use as cache key."""
    # Use file path, size, and modification time for cache key
    # This is faster than hashing file contents
    stat = os.stat(file_path)
    cache_key = f"{file_path}_{stat.st_size}_{stat.st_mtime}"
    return hashlib.md5(cache_key.encode()).hexdigest()

def _detect_format_by_content(file_path: str) -> Optional[str]:
    """Attempt to detect file format by examining file content."""
    with open(file_path, 'r') as f:
        first_line = f.readline().strip()

        # Quick format detection based on first line
        if first_line.startswith('>'):
            return 'fasta'
        elif first_line.startswith('CLUSTAL'):
            return 'clustal'
        elif first_line.startswith('#'):
            # Could be Stockholm
            if 'STOCKHOLM' in first_line:
                return 'stockholm'
        elif first_line.isdigit() or (len(first_line.split()) == 2 and
                                      first_line.split()[0].isdigit()):
            return 'phylip'

    return None

@lru_cache(maxsize=32)
def _cached_alignment_read(file_hash: str, file_path: str, file_format: str) -> Tuple[MultipleSeqAlignment, bool]:
    """Cached reading of alignment files."""
    with open(file_path) as f:
        alignment = AlignIO.read(f, file_format)
    return alignment, is_protein_alignment(alignment)

def get_alignment_and_format(
    alignment_file_path: str
) -> Tuple[MultipleSeqAlignment, str, bool]:
    # Check if file exists first
    if not os.path.exists(alignment_file_path):
        print(f"{alignment_file_path} corresponds to no such file.")
        print("Please check file name and pathing")
        sys.exit(2)

    # Try to detect format by content first
    detected_format = _detect_format_by_content(alignment_file_path)

    # Get file hash for caching
    file_hash = _get_file_hash(alignment_file_path)

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
    for fileFormat in FileFormat:
        # Skip the already tried format
        if detected_format and fileFormat.value == detected_format:
            continue

        try:
            alignment, is_protein = _cached_alignment_read(
                file_hash, alignment_file_path, fileFormat.value
            )
            return alignment, fileFormat.value, is_protein
        except (ValueError, AssertionError):
            continue

    # If we get here, no format worked
    print(f"Could not determine format for {alignment_file_path}")
    print("Please ensure the file is in a supported format")
    sys.exit(2)


def is_protein_alignment(alignment: MultipleSeqAlignment) -> bool:
    nucleotide_set = {
        "A", "C", "G", "T", "U", "-", "N", "?", "*"
    }

    for record in alignment:
        seq_set = set(record.seq.upper())
        if seq_set - nucleotide_set:
            # if there are chars that are not in the nucl set,
            # it's likely a protein sequence
            return True

    return False


def read_single_column_file_to_list(single_col_file_path: str) -> list:
    try:
        with open(single_col_file_path) as f:
            return [line.rstrip("\n").strip() for line in f]
    except FileNotFoundError:
        print(f"{single_col_file_path} corresponds to no such file or directory.")
        print("Please check file name and pathing")
        sys.exit(2)
