"""
Unit tests for streaming utilities
"""

import unittest
from unittest.mock import patch, MagicMock
import tempfile
import os
import subprocess
import sys
import itertools
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.helpers.streaming import (
    StreamingFastaReader,
    MemoryEfficientAlignmentProcessor
)


def test_module_import_does_not_import_biopython():
    code = """
import sys
import phykit.helpers.streaming

assert "typing" not in sys.modules
assert "Bio" not in sys.modules
assert "Bio.SeqIO" not in sys.modules
assert "Bio.SeqRecord" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


class TestStreamingFastaReader(unittest.TestCase):
    """Test StreamingFastaReader class"""

    def setUp(self):
        """Set up test fixtures"""
        # Create temporary FASTA file
        self.temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta')
        self.temp_file.write(">seq1\n")
        self.temp_file.write("ATCGATCG\n")
        self.temp_file.write(">seq2\n")
        self.temp_file.write("GCTAGCTA\n")
        self.temp_file.write(">seq3\n")
        self.temp_file.write("TTAATTAA\n")
        self.temp_file.close()

        self.reader = StreamingFastaReader(self.temp_file.name, chunk_size=2)

    def tearDown(self):
        """Clean up test fixtures"""
        if os.path.exists(self.temp_file.name):
            os.unlink(self.temp_file.name)

    def test_init(self):
        """Test reader initialization"""
        self.assertEqual(self.reader.file_path, self.temp_file.name)
        self.assertEqual(self.reader.chunk_size, 2)
        self.assertGreater(self.reader.file_size, 0)

    def test_stream_sequences(self):
        """Test streaming sequences one at a time"""
        sequences = list(self.reader.stream_sequences())

        self.assertEqual(len(sequences), 3)
        self.assertEqual(sequences[0].id, "seq1")
        self.assertEqual(str(sequences[0].seq), "ATCGATCG")
        self.assertEqual(sequences[1].id, "seq2")
        self.assertEqual(str(sequences[1].seq), "GCTAGCTA")
        self.assertEqual(sequences[2].id, "seq3")
        self.assertEqual(str(sequences[2].seq), "TTAATTAA")

    def test_stream_chunks(self):
        """Test streaming sequences in chunks"""
        chunks = list(self.reader.stream_chunks())

        # With chunk_size=2, we should get 2 chunks
        self.assertEqual(len(chunks), 2)
        self.assertEqual(len(chunks[0]), 2)  # First chunk has 2 sequences
        self.assertEqual(len(chunks[1]), 1)  # Second chunk has remaining 1 sequence

        # Check first chunk
        self.assertEqual(chunks[0][0].id, "seq1")
        self.assertEqual(chunks[0][1].id, "seq2")

        # Check second chunk
        self.assertEqual(chunks[1][0].id, "seq3")

    def test_stream_chunks_exact_division(self):
        """Test streaming when sequences divide evenly into chunks"""
        # Create reader with chunk_size=3 (exact division for 3 sequences)
        reader = StreamingFastaReader(self.temp_file.name, chunk_size=3)
        chunks = list(reader.stream_chunks())

        self.assertEqual(len(chunks), 1)
        self.assertEqual(len(chunks[0]), 3)

    def test_stream_chunks_batches_with_islice(self):
        """Test chunking uses iterator batches while preserving records."""
        original_islice = itertools.islice
        calls = []

        def tracking_islice(iterator, chunk_size):
            calls.append(chunk_size)
            return original_islice(iterator, chunk_size)

        with patch("itertools.islice", tracking_islice):
            chunks = list(self.reader.stream_chunks())

        self.assertEqual([[record.id for record in chunk] for chunk in chunks], [
            ["seq1", "seq2"],
            ["seq3"],
        ])
        self.assertEqual(calls, [2, 2, 2])

    def test_get_sequence_count(self):
        """Test counting sequences without loading entire file"""
        count = self.reader.get_sequence_count()
        self.assertEqual(count, 3)

    def test_get_sequence_count_across_chunk_boundary(self):
        """Test sequence counting when a header marker spans chunks"""
        with patch("phykit.helpers.streaming._COUNT_CHUNK_SIZE", 9):
            reader = StreamingFastaReader(self.temp_file.name)
            count = reader.get_sequence_count()

        self.assertEqual(count, 3)

    def test_get_sequence_count_when_chunk_starts_with_header(self):
        """Test sequence counting when a chunk begins with a header."""
        with patch("phykit.helpers.streaming._COUNT_CHUNK_SIZE", 15):
            reader = StreamingFastaReader(self.temp_file.name)
            count = reader.get_sequence_count()

        self.assertEqual(count, 3)

    def test_get_sequence_at_position(self):
        """Test getting specific sequence by position"""
        # Get first sequence
        seq = self.reader.get_sequence_at_position(0)
        self.assertIsNotNone(seq)
        self.assertEqual(seq.id, "seq1")

        # Get middle sequence
        seq = self.reader.get_sequence_at_position(1)
        self.assertIsNotNone(seq)
        self.assertEqual(seq.id, "seq2")

        # Get last sequence
        seq = self.reader.get_sequence_at_position(2)
        self.assertIsNotNone(seq)
        self.assertEqual(seq.id, "seq3")

        # Get non-existent sequence
        seq = self.reader.get_sequence_at_position(10)
        self.assertIsNone(seq)

    def test_empty_file(self):
        """Test handling of empty FASTA file"""
        empty_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta')
        empty_file.close()

        try:
            reader = StreamingFastaReader(empty_file.name)
            sequences = list(reader.stream_sequences())
            self.assertEqual(len(sequences), 0)

            # Empty file will raise ValueError with mmap
            with self.assertRaises(ValueError):
                reader.get_sequence_count()
        finally:
            os.unlink(empty_file.name)


class TestMemoryEfficientAlignmentProcessor(unittest.TestCase):
    """Test MemoryEfficientAlignmentProcessor class"""

    def setUp(self):
        """Set up test fixtures"""
        # Create temporary alignment file
        self.temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta')
        self.temp_file.write(">seq1\n")
        self.temp_file.write("ATCG-TCG\n")
        self.temp_file.write(">seq2\n")
        self.temp_file.write("ATCGATCG\n")
        self.temp_file.write(">seq3\n")
        self.temp_file.write("ANCG-TCG\n")
        self.temp_file.close()

    def tearDown(self):
        """Clean up test fixtures"""
        if os.path.exists(self.temp_file.name):
            os.unlink(self.temp_file.name)

    @patch('phykit.helpers.streaming.StreamingFastaReader')
    def test_calculate_column_stats_streaming(self, mock_reader_class):
        """Test calculating column statistics without loading entire alignment"""
        # Setup mock reader
        mock_reader = MagicMock()
        mock_reader_class.return_value = mock_reader

        # Create mock sequences
        seq1 = SeqRecord(Seq("ATCG-TCG"), id="seq1")
        seq2 = SeqRecord(Seq("ATCGATCG"), id="seq2")
        seq3 = SeqRecord(Seq("ANCG-TCG"), id="seq3")

        mock_reader.stream_sequences.return_value = iter([seq1, seq2, seq3])

        stats = MemoryEfficientAlignmentProcessor.calculate_column_stats_streaming(
            "dummy.fasta"
        )

        # Check dimensions
        self.assertIn('variable_sites', stats)
        self.assertIn('gap_counts', stats)
        self.assertIn('conservation', stats)

        # Check column stats length
        self.assertEqual(len(stats['variable_sites']), 8)
        self.assertEqual(len(stats['gap_counts']), 8)
        self.assertEqual(len(stats['conservation']), 8)

        # Check specific columns
        # Column 0: A, A, A - not variable
        self.assertFalse(stats['variable_sites'][0])
        # Column 1: T, T, N - variable (N counts as different)
        self.assertTrue(stats['variable_sites'][1])
        # Column 4: -, A, - - has gaps
        self.assertGreater(stats['gap_counts'][4], 0)
        mock_reader.stream_sequences.assert_called_once()

    def test_calculate_column_stats_streaming_real_file(self):
        """Test column stats calculation with real file"""
        stats = MemoryEfficientAlignmentProcessor.calculate_column_stats_streaming(
            self.temp_file.name
        )

        # Check that we got results
        self.assertEqual(len(stats['variable_sites']), 8)
        self.assertEqual(len(stats['gap_counts']), 8)
        self.assertEqual(len(stats['conservation']), 8)

        # Check specific values
        # Column 4 has gaps in seq1 and seq3
        self.assertEqual(stats['gap_counts'][4], 2)

    @patch('phykit.helpers.streaming.StreamingFastaReader')
    def test_process_large_alignment_in_batches(self, mock_reader_class):
        """Test processing large alignment in batches"""
        # Setup mock reader
        mock_reader = MagicMock()
        mock_reader_class.return_value = mock_reader

        # Create mock sequences
        batch1 = [
            SeqRecord(Seq("ATCG"), id="seq1"),
            SeqRecord(Seq("ATCG"), id="seq2")
        ]
        batch2 = [
            SeqRecord(Seq("GCTA"), id="seq3"),
            SeqRecord(Seq("GCTA"), id="seq4")
        ]

        mock_reader.stream_chunks.return_value = iter([batch1, batch2])

        # Define processing function
        def process_batch(batch):
            return len(batch)

        results = MemoryEfficientAlignmentProcessor.process_large_alignment_in_batches(
            "dummy.fasta", process_batch, batch_size=2
        )

        self.assertEqual(len(results), 2)
        self.assertEqual(results[0], 2)  # First batch has 2 sequences
        self.assertEqual(results[1], 2)  # Second batch has 2 sequences

    @patch('phykit.helpers.streaming.StreamingFastaReader')
    def test_process_large_alignment_in_batches_does_not_enumerate_chunks(
        self, mock_reader_class
    ):
        """Test batch processing does not build unused chunk indices."""
        mock_reader = MagicMock()
        mock_reader_class.return_value = mock_reader
        mock_reader.stream_chunks.return_value = iter([[1], [2]])

        with patch(
            "builtins.enumerate",
            side_effect=AssertionError("batch processor should not enumerate"),
        ):
            results = MemoryEfficientAlignmentProcessor.process_large_alignment_in_batches(
                "dummy.fasta",
                len,
                batch_size=1,
            )

        self.assertEqual(results, [1, 1])

    def test_process_large_alignment_with_aggregation(self):
        """Test processing with custom aggregation function"""
        # Define processing function that returns sequence lengths
        def get_sequence_lengths(batch):
            return [len(record.seq) for record in batch]

        results = MemoryEfficientAlignmentProcessor.process_large_alignment_in_batches(
            self.temp_file.name, get_sequence_lengths, batch_size=2
        )

        # Flatten results
        all_lengths = []
        for batch_result in results:
            all_lengths.extend(batch_result)

        # All sequences should be 8 characters long
        self.assertEqual(all_lengths, [8, 8, 8])

    def test_empty_alignment(self):
        """Test processing empty alignment"""
        empty_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta')
        empty_file.close()

        try:
            def process_batch(batch):
                return len(batch)

            results = MemoryEfficientAlignmentProcessor.process_large_alignment_in_batches(
                empty_file.name, process_batch, batch_size=10
            )

            self.assertEqual(len(results), 0)
        finally:
            os.unlink(empty_file.name)
