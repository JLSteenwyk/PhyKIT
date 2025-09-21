"""
Caching utilities for expensive computations
"""

import hashlib
import pickle
import os
import tempfile
from functools import wraps, lru_cache
from typing import Any, Callable
import json


class ResultCache:
    """
    File-based cache for expensive computation results.
    """

    def __init__(self, cache_dir: str = None):
        """
        Initialize cache.

        Args:
            cache_dir: Directory for cache files (uses temp dir if None)
        """
        if cache_dir is None:
            cache_dir = os.path.join(tempfile.gettempdir(), 'phykit_cache')

        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)

    def _get_cache_key(self, *args, **kwargs) -> str:
        """Generate a unique cache key from function arguments."""
        # Create a string representation of arguments
        key_parts = []

        for arg in args:
            if isinstance(arg, (str, int, float, bool)):
                key_parts.append(str(arg))
            elif hasattr(arg, '__dict__'):
                # For objects, use their attributes
                key_parts.append(json.dumps(vars(arg), sort_keys=True, default=str))
            else:
                key_parts.append(str(arg))

        for k, v in sorted(kwargs.items()):
            key_parts.append(f"{k}={v}")

        key_string = "_".join(key_parts)
        return hashlib.md5(key_string.encode()).hexdigest()

    def get(self, cache_key: str) -> Any:
        """Retrieve cached result."""
        cache_file = os.path.join(self.cache_dir, f"{cache_key}.pkl")

        if os.path.exists(cache_file):
            try:
                with open(cache_file, 'rb') as f:
                    return pickle.load(f)
            except:
                # Cache corrupted, remove it
                os.remove(cache_file)

        return None

    def set(self, cache_key: str, value: Any) -> None:
        """Store result in cache."""
        cache_file = os.path.join(self.cache_dir, f"{cache_key}.pkl")

        try:
            with open(cache_file, 'wb') as f:
                pickle.dump(value, f)
        except:
            # Caching failed, continue without caching
            pass

    def clear(self) -> None:
        """Clear all cached results."""
        for file in os.listdir(self.cache_dir):
            if file.endswith('.pkl'):
                os.remove(os.path.join(self.cache_dir, file))


def cached_computation(cache_instance: ResultCache = None):
    """
    Decorator for caching expensive computation results.

    Usage:
        @cached_computation()
        def expensive_function(param1, param2):
            # Expensive computation
            return result
    """
    if cache_instance is None:
        cache_instance = ResultCache()

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Generate cache key
            cache_key = cache_instance._get_cache_key(func.__name__, *args, **kwargs)

            # Try to get cached result
            cached_result = cache_instance.get(cache_key)
            if cached_result is not None:
                return cached_result

            # Compute result
            result = func(*args, **kwargs)

            # Cache result
            cache_instance.set(cache_key, result)

            return result

        # Add method to clear cache for this function
        wrapper.clear_cache = cache_instance.clear

        return wrapper

    return decorator


# Specialized caching for tree operations
@lru_cache(maxsize=128)
def cached_tree_distance(tree_pickle: bytes, tip1: str, tip2: str) -> float:
    """
    Cache tree distance calculations.

    Args:
        tree_pickle: Pickled tree object
        tip1: First tip name
        tip2: Second tip name

    Returns:
        Distance between tips
    """
    import pickle
    tree = pickle.loads(tree_pickle)
    return tree.distance(tip1, tip2)


# Specialized caching for alignment operations
class AlignmentCache:
    """
    Specialized cache for alignment operations.
    """

    def __init__(self):
        self._column_cache = {}
        self._stats_cache = {}

    @lru_cache(maxsize=1024)
    def get_column(self, alignment_hash: str, column_idx: int) -> str:
        """
        Get cached alignment column.
        """
        return self._column_cache.get(f"{alignment_hash}_{column_idx}")

    def set_column(self, alignment_hash: str, column_idx: int, column: str) -> None:
        """
        Cache alignment column.
        """
        self._column_cache[f"{alignment_hash}_{column_idx}"] = column

    @lru_cache(maxsize=128)
    def get_stats(self, alignment_hash: str, stat_type: str) -> Any:
        """
        Get cached alignment statistics.
        """
        return self._stats_cache.get(f"{alignment_hash}_{stat_type}")

    def set_stats(self, alignment_hash: str, stat_type: str, stats: Any) -> None:
        """
        Cache alignment statistics.
        """
        self._stats_cache[f"{alignment_hash}_{stat_type}"] = stats

    def clear(self) -> None:
        """
        Clear all caches.
        """
        self._column_cache.clear()
        self._stats_cache.clear()
        self.get_column.cache_clear()
        self.get_stats.cache_clear()


# Global cache instances
_result_cache = ResultCache()
_alignment_cache = AlignmentCache()


def get_result_cache() -> ResultCache:
    """Get global result cache instance."""
    return _result_cache


def get_alignment_cache() -> AlignmentCache:
    """Get global alignment cache instance."""
    return _alignment_cache