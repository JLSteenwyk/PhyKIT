"""
Streaming utilities for memory-efficient processing of large files
"""

from __future__ import annotations

import os


_COUNT_CHUNK_SIZE = 1024 * 1024


class StreamingFastaReader:
    """
    Memory-efficient streaming reader for large FASTA files.
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

    def stream_sequences(self) -> object:
        """
        Stream sequences one at a time.
        """
        from Bio import SeqIO

        with open(self.file_path) as handle:
            yield from SeqIO.parse(handle, "fasta")

    def stream_chunks(self) -> object:
        """
        Stream sequences in chunks for batch processing.
        """
        from itertools import islice

        sequence_iter = iter(self.stream_sequences())
        chunk_size = self.chunk_size
        while True:
            chunk = list(islice(sequence_iter, chunk_size))
            if not chunk:
                break
            yield chunk

    def get_sequence_count(self) -> int:
        """
        Count sequences without loading entire file.
        """
        count = 0
        saw_data = False
        previous = b""
        with open(self.file_path, 'rb') as f:
            while True:
                chunk = f.read(_COUNT_CHUNK_SIZE)
                if not chunk:
                    break
                if not saw_data:
                    saw_data = True
                    if chunk[:1] == b'>':
                        count = 1
                elif previous == b"\n" and chunk[:1] == b">":
                    count += 1
                count += chunk.count(b"\n>")
                previous = chunk[-1:]

        if not saw_data:
            raise ValueError("cannot mmap an empty file")

        return count

    def get_sequence_at_position(self, position: int) -> object | None:
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
        seq_length = None
        num_seqs = 0
        unique_chars_by_col = None
        gap_counts = None
        gap_chars = {"-", "?", "X", "N"}

        for record in reader.stream_sequences():
            seq = str(record.seq).upper()
            if seq_length is None:
                seq_length = len(seq)
                unique_chars_by_col = [set() for _ in range(seq_length)]
                gap_counts = [0] * seq_length
            num_seqs += 1

            for col_idx, char in enumerate(seq):
                unique_chars_by_col[col_idx].add(char)
                if char in gap_chars:
                    gap_counts[col_idx] += 1

        if seq_length is None:
            return {
                'variable_sites': [],
                'gap_counts': [],
                'conservation': []
            }

        return {
            'variable_sites': [
                len(unique_chars) > 1 for unique_chars in unique_chars_by_col
            ],
            'gap_counts': gap_counts,
            'conservation': [
                1 - (len(unique_chars) / num_seqs)
                for unique_chars in unique_chars_by_col
            ],
        }

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

        for batch in reader.stream_chunks():
            batch_result = processing_func(batch)
            results.append(batch_result)

        return results
