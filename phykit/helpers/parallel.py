"""
Parallel processing utilities for batch operations
"""

import multiprocessing as mp
from functools import partial
from typing import List, Any, Callable, Optional, Tuple
import numpy as np
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import sys


class ParallelProcessor:
    """
    Utility class for parallel processing of batch operations.
    """

    @staticmethod
    def get_optimal_workers(data_size: int, min_chunk_size: int = 10) -> int:
        """
        Determine optimal number of workers based on data size.

        Args:
            data_size: Size of data to process
            min_chunk_size: Minimum size per chunk

        Returns:
            Optimal number of workers
        """
        max_workers = mp.cpu_count()
        optimal_workers = min(max_workers, max(1, data_size // min_chunk_size))
        return min(optimal_workers, 8)  # Cap at 8 to avoid overhead

    @staticmethod
    def chunk_data(data: List[Any], num_chunks: int) -> List[List[Any]]:
        """
        Split data into chunks for parallel processing.

        Args:
            data: Data to split
            num_chunks: Number of chunks

        Returns:
            List of data chunks
        """
        chunk_size = max(1, len(data) // num_chunks)
        chunks = []

        for i in range(0, len(data), chunk_size):
            chunk = data[i:i + chunk_size]
            if chunk:
                chunks.append(chunk)

        return chunks

    @staticmethod
    def parallel_map(
        func: Callable,
        data: List[Any],
        num_workers: Optional[int] = None,
        use_threads: bool = False,
        show_progress: bool = False
    ) -> List[Any]:
        """
        Apply function to data in parallel.

        Args:
            func: Function to apply
            data: Data to process
            num_workers: Number of workers (auto-determined if None)
            use_threads: Use threads instead of processes
            show_progress: Show progress bar

        Returns:
            List of results
        """
        if not data:
            return []

        # Determine number of workers
        if num_workers is None:
            num_workers = ParallelProcessor.get_optimal_workers(len(data))

        # For small datasets, use sequential processing
        if len(data) < 20 or num_workers == 1:
            return [func(item) for item in data]

        # Choose executor type
        executor_class = ThreadPoolExecutor if use_threads else ProcessPoolExecutor

        results = []
        with executor_class(max_workers=num_workers) as executor:
            if show_progress and sys.stderr.isatty():
                try:
                    from tqdm import tqdm
                    futures = {executor.submit(func, item): i for i, item in enumerate(data)}
                    results = [None] * len(data)

                    for future in tqdm(as_completed(futures), total=len(data), desc="Processing"):
                        idx = futures[future]
                        results[idx] = future.result()
                except ImportError:
                    # Fallback without progress bar
                    results = list(executor.map(func, data))
            else:
                results = list(executor.map(func, data))

        return results

    @staticmethod
    def parallel_reduce(
        func: Callable,
        data: List[Any],
        reduce_func: Callable,
        initial_value: Any = None,
        num_workers: Optional[int] = None
    ) -> Any:
        """
        Apply function to data in parallel and reduce results.

        Args:
            func: Function to apply to each item
            data: Data to process
            reduce_func: Function to reduce results
            initial_value: Initial value for reduction
            num_workers: Number of workers

        Returns:
            Reduced result
        """
        # Apply function in parallel
        results = ParallelProcessor.parallel_map(func, data, num_workers)

        # Reduce results
        if initial_value is not None:
            result = initial_value
            for item in results:
                result = reduce_func(result, item)
        else:
            if not results:
                return None
            result = results[0]
            for item in results[1:]:
                result = reduce_func(result, item)

        return result


class BatchFileProcessor:
    """
    Process multiple files in parallel.
    """

    @staticmethod
    def process_files(
        file_paths: List[str],
        processing_func: Callable,
        num_workers: Optional[int] = None,
        aggregate_func: Optional[Callable] = None
    ) -> Any:
        """
        Process multiple files in parallel.

        Args:
            file_paths: List of file paths
            processing_func: Function to process each file
            num_workers: Number of workers
            aggregate_func: Function to aggregate results

        Returns:
            Processed results or aggregated result
        """
        if not file_paths:
            return []

        # Process files in parallel
        results = ParallelProcessor.parallel_map(
            processing_func,
            file_paths,
            num_workers,
            show_progress=True
        )

        # Aggregate results if function provided
        if aggregate_func:
            return aggregate_func(results)

        return results

    @staticmethod
    def process_file_pairs(
        file_pairs: List[Tuple[str, str]],
        processing_func: Callable,
        num_workers: Optional[int] = None
    ) -> List[Any]:
        """
        Process pairs of files in parallel.

        Args:
            file_pairs: List of file path pairs
            processing_func: Function to process each pair
            num_workers: Number of workers

        Returns:
            List of results
        """
        def process_pair(pair):
            return processing_func(pair[0], pair[1])

        return ParallelProcessor.parallel_map(
            process_pair,
            file_pairs,
            num_workers
        )


class NumpyParallel:
    """
    Utilities for parallel NumPy operations.
    """

    @staticmethod
    def parallel_apply_along_axis(
        func: Callable,
        axis: int,
        array: np.ndarray,
        num_workers: Optional[int] = None
    ) -> np.ndarray:
        """
        Apply function along axis in parallel.

        Args:
            func: Function to apply
            axis: Axis along which to apply function
            array: NumPy array
            num_workers: Number of workers

        Returns:
            Result array
        """
        if axis == 0:
            # Process columns
            results = ParallelProcessor.parallel_map(
                lambda col: func(array[:, col]),
                list(range(array.shape[1])),
                num_workers
            )
            return np.array(results).T
        elif axis == 1:
            # Process rows
            results = ParallelProcessor.parallel_map(
                lambda row: func(array[row, :]),
                list(range(array.shape[0])),
                num_workers
            )
            return np.array(results)
        else:
            raise ValueError(f"Unsupported axis: {axis}")

    @staticmethod
    def parallel_pairwise_operation(
        items: List[Any],
        operation_func: Callable,
        num_workers: Optional[int] = None,
        symmetric: bool = True
    ) -> np.ndarray:
        """
        Perform pairwise operations in parallel.

        Args:
            items: List of items
            operation_func: Function to apply to pairs
            num_workers: Number of workers
            symmetric: Whether operation is symmetric

        Returns:
            Matrix of pairwise results
        """
        n = len(items)
        result_matrix = np.zeros((n, n))

        # Generate pairs
        pairs = []
        for i in range(n):
            for j in range(i + 1 if symmetric else 0, n):
                pairs.append((i, j, items[i], items[j]))

        # Process pairs in parallel
        def process_pair(pair_data):
            i, j, item1, item2 = pair_data
            return i, j, operation_func(item1, item2)

        results = ParallelProcessor.parallel_map(
            process_pair,
            pairs,
            num_workers
        )

        # Fill result matrix
        for i, j, value in results:
            result_matrix[i, j] = value
            if symmetric:
                result_matrix[j, i] = value

        return result_matrix