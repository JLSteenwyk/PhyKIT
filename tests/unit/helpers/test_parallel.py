"""
Unit tests for parallel processing utilities
"""

import unittest
from unittest.mock import Mock, patch, MagicMock
from typing import Tuple
import multiprocessing as mp
import numpy as np
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

from phykit.helpers.parallel import (
    ParallelProcessor,
    BatchFileProcessor,
    NumpyParallel
)


class TestParallelProcessor(unittest.TestCase):
    """Test ParallelProcessor class"""

    def test_get_optimal_workers(self):
        """Test optimal worker calculation"""
        # Small dataset
        workers = ParallelProcessor.get_optimal_workers(5, min_chunk_size=10)
        self.assertEqual(workers, 1)  # Too small for parallelization

        # Medium dataset
        workers = ParallelProcessor.get_optimal_workers(100, min_chunk_size=10)
        self.assertGreaterEqual(workers, 1)
        self.assertLessEqual(workers, 8)  # Capped at 8

        # Large dataset
        workers = ParallelProcessor.get_optimal_workers(1000, min_chunk_size=10)
        self.assertGreaterEqual(workers, 1)
        self.assertLessEqual(workers, 8)

        # Very large dataset (should still be capped at 8)
        workers = ParallelProcessor.get_optimal_workers(10000, min_chunk_size=10)
        self.assertEqual(workers, 8)

    def test_chunk_data(self):
        """Test data chunking"""
        # Even division
        data = list(range(10))
        chunks = ParallelProcessor.chunk_data(data, 2)
        self.assertEqual(len(chunks), 2)
        self.assertEqual(chunks[0], [0, 1, 2, 3, 4])
        self.assertEqual(chunks[1], [5, 6, 7, 8, 9])

        # Uneven division
        data = list(range(11))
        chunks = ParallelProcessor.chunk_data(data, 3)
        self.assertGreaterEqual(len(chunks), 3)
        self.assertEqual(sum(len(chunk) for chunk in chunks), 11)

        # Single chunk
        data = list(range(5))
        chunks = ParallelProcessor.chunk_data(data, 1)
        self.assertEqual(len(chunks), 1)
        self.assertEqual(chunks[0], data)

        # Empty data
        chunks = ParallelProcessor.chunk_data([], 5)
        self.assertEqual(len(chunks), 0)

    def test_parallel_map_empty_data(self):
        """Test parallel_map with empty data"""
        result = ParallelProcessor.parallel_map(lambda x: x * 2, [])
        self.assertEqual(result, [])

    def test_parallel_map_small_dataset(self):
        """Test parallel_map with small dataset (sequential processing)"""
        def square(x):
            return x ** 2

        data = [1, 2, 3, 4]
        result = ParallelProcessor.parallel_map(square, data)
        self.assertEqual(result, [1, 4, 9, 16])

    @patch('phykit.helpers.parallel.ProcessPoolExecutor')
    def test_parallel_map_large_dataset_processes(self, mock_executor_class):
        """Test parallel_map with large dataset using processes"""
        # Setup mock
        mock_executor = MagicMock()
        mock_executor_class.return_value.__enter__.return_value = mock_executor
        mock_executor.map.return_value = [2, 4, 6, 8, 10]

        def double(x):
            return x * 2

        data = list(range(1, 6))
        # Force parallel processing by using large dataset
        data = data * 10  # 50 elements

        result = ParallelProcessor.parallel_map(double, data, num_workers=2)

        # Should use ProcessPoolExecutor
        mock_executor_class.assert_called_once_with(max_workers=2)
        mock_executor.map.assert_called_once()

    @patch('phykit.helpers.parallel.ThreadPoolExecutor')
    def test_parallel_map_with_threads(self, mock_executor_class):
        """Test parallel_map with threads"""
        # Setup mock
        mock_executor = MagicMock()
        mock_executor_class.return_value.__enter__.return_value = mock_executor
        mock_executor.map.return_value = [2, 4, 6]

        def double(x):
            return x * 2

        data = list(range(1, 4)) * 10  # Make it large enough for parallel

        result = ParallelProcessor.parallel_map(
            double, data, num_workers=2, use_threads=True
        )

        # Should use ThreadPoolExecutor
        mock_executor_class.assert_called_once_with(max_workers=2)

    @patch('phykit.helpers.parallel.sys.stderr.isatty')
    @patch('phykit.helpers.parallel.ProcessPoolExecutor')
    def test_parallel_map_with_progress(self, mock_executor_class, mock_isatty):
        """Test parallel_map with progress bar disabled (non-TTY)"""
        # Disable progress bar by making stderr non-TTY
        # This avoids issues with tqdm and as_completed in tests
        mock_isatty.return_value = False

        # Setup mock executor
        mock_executor = MagicMock()
        mock_executor_class.return_value.__enter__.return_value = mock_executor
        mock_executor.map.return_value = [2, 4, 6, 8, 10, 12] * 5  # Results for 30 items

        def double(x):
            return x * 2

        data = list(range(1, 4)) * 10  # 30 items

        # Should work without progress bar
        result = ParallelProcessor.parallel_map(
            double, data, show_progress=True, num_workers=2
        )

        # Map should be called since we disabled TTY
        mock_executor.map.assert_called_once()

    @patch('phykit.helpers.parallel.sys')
    @patch('phykit.helpers.parallel.ProcessPoolExecutor')
    def test_parallel_map_with_tqdm_import_error(self, mock_executor_class, mock_sys):
        """Test parallel_map when tqdm is not available"""
        # Make stderr appear as TTY
        mock_sys.stderr.isatty.return_value = True

        # Setup mock executor
        mock_executor = MagicMock()
        mock_executor_class.return_value.__enter__.return_value = mock_executor
        mock_executor.map.return_value = [2, 4, 6] * 10

        def double(x):
            return x * 2

        data = list(range(1, 4)) * 10  # Large dataset

        # Mock tqdm import to fail
        with patch.dict('sys.modules', {'tqdm': None}):
            result = ParallelProcessor.parallel_map(
                double, data, show_progress=True, num_workers=2
            )

        # Should fall back to regular map
        mock_executor.map.assert_called_once()

    def test_parallel_reduce(self):
        """Test parallel_reduce function"""
        def square(x):
            return x ** 2

        def add(x, y):
            return x + y

        data = [1, 2, 3, 4]
        result = ParallelProcessor.parallel_reduce(
            square, data, add, initial_value=0
        )
        # 1^2 + 2^2 + 3^2 + 4^2 = 1 + 4 + 9 + 16 = 30
        self.assertEqual(result, 30)

    def test_parallel_reduce_without_initial(self):
        """Test parallel_reduce without initial value"""
        def identity(x):
            return x

        def multiply(x, y):
            return x * y

        data = [2, 3, 4]
        result = ParallelProcessor.parallel_reduce(
            identity, data, multiply
        )
        # 2 * 3 * 4 = 24
        self.assertEqual(result, 24)

    def test_parallel_reduce_empty_data(self):
        """Test parallel_reduce with empty data"""
        result = ParallelProcessor.parallel_reduce(
            lambda x: x, [], lambda x, y: x + y
        )
        self.assertIsNone(result)


class TestBatchFileProcessor(unittest.TestCase):
    """Test BatchFileProcessor class"""

    @patch('phykit.helpers.parallel.ParallelProcessor.parallel_map')
    def test_process_files(self, mock_parallel_map):
        """Test processing multiple files"""
        mock_parallel_map.return_value = ["result1", "result2", "result3"]

        def process_func(file_path):
            return f"processed_{file_path}"

        file_paths = ["file1.txt", "file2.txt", "file3.txt"]
        result = BatchFileProcessor.process_files(file_paths, process_func)

        self.assertEqual(result, ["result1", "result2", "result3"])
        mock_parallel_map.assert_called_once()

    def test_process_files_empty_list(self):
        """Test processing empty file list"""
        result = BatchFileProcessor.process_files([], lambda x: x)
        self.assertEqual(result, [])

    @patch('phykit.helpers.parallel.ParallelProcessor.parallel_map')
    def test_process_files_with_aggregation(self, mock_parallel_map):
        """Test processing files with aggregation function"""
        mock_parallel_map.return_value = [1, 2, 3, 4]

        def process_func(file_path):
            return len(file_path)

        def aggregate_func(results):
            return sum(results)

        file_paths = ["a.txt", "bb.txt", "ccc.txt", "dddd.txt"]
        result = BatchFileProcessor.process_files(
            file_paths, process_func, aggregate_func=aggregate_func
        )

        self.assertEqual(result, 10)  # Sum of [1, 2, 3, 4]

    @patch('phykit.helpers.parallel.ParallelProcessor.parallel_map')
    def test_process_file_pairs(self, mock_parallel_map):
        """Test processing file pairs"""
        mock_parallel_map.return_value = ["pair1_result", "pair2_result"]

        def process_func(file1, file2):
            return f"{file1}_{file2}"

        file_pairs = [("a.txt", "b.txt"), ("c.txt", "d.txt")]
        result = BatchFileProcessor.process_file_pairs(file_pairs, process_func)

        self.assertEqual(result, ["pair1_result", "pair2_result"])
        mock_parallel_map.assert_called_once()

        # Check that the wrapper function was created correctly
        wrapper_func = mock_parallel_map.call_args[0][0]
        test_result = wrapper_func(("test1", "test2"))
        # The wrapper should call process_func with unpacked tuple
        self.assertEqual(test_result, "test1_test2")


class TestNumpyParallel(unittest.TestCase):
    """Test NumpyParallel class"""

    def test_parallel_apply_along_axis_columns(self):
        """Test applying function along columns (axis=0)"""
        array = np.array([[1, 2, 3],
                          [4, 5, 6],
                          [7, 8, 9]])

        def col_sum(col):
            return np.sum(col)

        with patch('phykit.helpers.parallel.ParallelProcessor.parallel_map') as mock_map:
            mock_map.return_value = [12, 15, 18]  # Column sums

            result = NumpyParallel.parallel_apply_along_axis(
                col_sum, 0, array
            )

            # Result should be transposed back
            expected = np.array([12, 15, 18])
            np.testing.assert_array_equal(result, expected)

    def test_parallel_apply_along_axis_rows(self):
        """Test applying function along rows (axis=1)"""
        array = np.array([[1, 2, 3],
                          [4, 5, 6],
                          [7, 8, 9]])

        def row_sum(row):
            return np.sum(row)

        with patch('phykit.helpers.parallel.ParallelProcessor.parallel_map') as mock_map:
            mock_map.return_value = [6, 15, 24]  # Row sums

            result = NumpyParallel.parallel_apply_along_axis(
                row_sum, 1, array
            )

            expected = np.array([6, 15, 24])
            np.testing.assert_array_equal(result, expected)

    def test_parallel_apply_along_axis_invalid_axis(self):
        """Test applying function with invalid axis"""
        array = np.array([[1, 2], [3, 4]])

        with self.assertRaises(ValueError) as context:
            NumpyParallel.parallel_apply_along_axis(
                lambda x: x, 2, array
            )

        self.assertIn("Unsupported axis: 2", str(context.exception))

    def test_parallel_pairwise_operation_symmetric(self):
        """Test pairwise operations with symmetric matrix"""
        items = ["A", "B", "C"]

        def distance_func(item1, item2):
            if item1 == item2:
                return 0
            return abs(ord(item1) - ord(item2))

        with patch('phykit.helpers.parallel.ParallelProcessor.parallel_map') as mock_map:
            # Mock results for upper triangle pairs
            mock_map.return_value = [
                (0, 1, 1),  # A-B distance
                (0, 2, 2),  # A-C distance
                (1, 2, 1),  # B-C distance
            ]

            result = NumpyParallel.parallel_pairwise_operation(
                items, distance_func, symmetric=True
            )

            # Check that result is symmetric
            expected = np.array([[0, 1, 2],
                                [1, 0, 1],
                                [2, 1, 0]])
            np.testing.assert_array_equal(result, expected)

    def test_parallel_pairwise_operation_non_symmetric(self):
        """Test pairwise operations with non-symmetric matrix"""
        items = [1, 2, 3]

        def subtract_func(item1, item2):
            return item1 - item2

        with patch('phykit.helpers.parallel.ParallelProcessor.parallel_map') as mock_map:
            # Mock results for all pairs
            mock_map.return_value = [
                (0, 0, 0),   # 1-1
                (0, 1, -1),  # 1-2
                (0, 2, -2),  # 1-3
                (1, 0, 1),   # 2-1
                (1, 1, 0),   # 2-2
                (1, 2, -1),  # 2-3
                (2, 0, 2),   # 3-1
                (2, 1, 1),   # 3-2
                (2, 2, 0),   # 3-3
            ]

            result = NumpyParallel.parallel_pairwise_operation(
                items, subtract_func, symmetric=False
            )

            expected = np.array([[0, -1, -2],
                                [1, 0, -1],
                                [2, 1, 0]])
            np.testing.assert_array_equal(result, expected)