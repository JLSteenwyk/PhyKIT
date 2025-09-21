"""
Streaming utilities for memory-efficient processing of large files
"""

from typing import Iterator, Tuple, Optional
import mmap
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class StreamingFastaReader:
    """
    Memory-efficient streaming reader for large FASTA files.
    Uses memory mapping to avoid loading entire file into memory.
    """

    def __init__(self, file_path: str, chunk_size: int = 1000):
        """
        Initialize streaming reader.

        Args:
            file_path: Path to FASTA file
            chunk_size: Number of sequences to yield at once
        """
        self.file_path = file_path
        self.chunk_size = chunk_size
        self.file_size = os.path.getsize(file_path)

    def stream_sequences(self) -> Iterator[SeqRecord]:
        """
        Stream sequences one at a time.
        """
        with open(self.file_path, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                yield record

    def stream_chunks(self) -> Iterator[list]:
        """
        Stream sequences in chunks for batch processing.
        """
        chunk = []
        for record in self.stream_sequences():
            chunk.append(record)
            if len(chunk) >= self.chunk_size:
                yield chunk
                chunk = []

        # Yield remaining sequences
        if chunk:
            yield chunk

    def get_sequence_count(self) -> int:
        """
        Count sequences without loading entire file.
        """
        count = 0
        with open(self.file_path, 'rb') as f:
            with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mmapped_file:
                for line in iter(mmapped_file.readline, b""):
                    if line.startswith(b'>'):
                        count += 1
        return count

    def get_sequence_at_position(self, position: int) -> Optional[SeqRecord]:
        """
        Get a specific sequence by position without loading entire file.
        """
        for i, record in enumerate(self.stream_sequences()):
            if i == position:
                return record
        return None


class MemoryEfficientAlignmentProcessor:
    """
    Process large alignments with minimal memory footprint.
    """

    @staticmethod
    def calculate_column_stats_streaming(file_path: str) -> dict:
        """
        Calculate column statistics without loading entire alignment.

        Returns:
            Dictionary with column statistics
        """
        reader = StreamingFastaReader(file_path)

        # First pass: get dimensions
        num_seqs = 0
        seq_length = None

        for record in reader.stream_sequences():
            if seq_length is None:
                seq_length = len(record.seq)
            num_seqs += 1

        # Initialize column stats
        column_stats = {
            'variable_sites': [],
            'gap_counts': [],
            'conservation': []
        }

        # Process in chunks to maintain memory efficiency
        for col_idx in range(seq_length):
            column_chars = []
            gap_count = 0

            for record in reader.stream_sequences():
                char = str(record.seq[col_idx]).upper()
                column_chars.append(char)
                if char in ['-', '?', 'X', 'N']:
                    gap_count += 1

            # Calculate statistics
            unique_chars = set(column_chars)
            is_variable = len(unique_chars) > 1
            conservation_score = 1 - (len(unique_chars) / len(column_chars))

            column_stats['variable_sites'].append(is_variable)
            column_stats['gap_counts'].append(gap_count)
            column_stats['conservation'].append(conservation_score)

        return column_stats

    @staticmethod
    def process_large_alignment_in_batches(
        file_path: str,
        processing_func,
        batch_size: int = 100
    ):
        """
        Process large alignment in batches to avoid memory issues.

        Args:
            file_path: Path to alignment file
            processing_func: Function to apply to each batch
            batch_size: Number of sequences per batch

        Returns:
            Aggregated results from processing function
        """
        reader = StreamingFastaReader(file_path, chunk_size=batch_size)
        results = []

        for batch_num, batch in enumerate(reader.stream_chunks()):
            batch_result = processing_func(batch)
            results.append(batch_result)

        return results