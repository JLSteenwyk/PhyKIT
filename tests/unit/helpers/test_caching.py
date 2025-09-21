"""
Unit tests for caching utilities
"""

import os
import tempfile
import shutil
import pickle
from unittest import TestCase
from unittest.mock import Mock, patch, MagicMock
import json

from phykit.helpers.caching import (
    ResultCache,
    cached_computation,
    cached_tree_distance,
    AlignmentCache,
    get_result_cache,
    get_alignment_cache
)


class TestResultCache(TestCase):
    """Test ResultCache class"""

    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = tempfile.mkdtemp()
        self.cache = ResultCache(cache_dir=self.temp_dir)

    def tearDown(self):
        """Clean up test fixtures"""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    def test_init_with_default_dir(self):
        """Test initialization with default directory"""
        cache = ResultCache()
        self.assertTrue(os.path.exists(cache.cache_dir))
        self.assertIn('phykit_cache', cache.cache_dir)
        # Clean up
        if os.path.exists(cache.cache_dir):
            shutil.rmtree(cache.cache_dir)

    def test_init_with_custom_dir(self):
        """Test initialization with custom directory"""
        custom_dir = os.path.join(self.temp_dir, 'custom_cache')
        cache = ResultCache(cache_dir=custom_dir)
        self.assertEqual(cache.cache_dir, custom_dir)
        self.assertTrue(os.path.exists(custom_dir))

    def test_get_cache_key_basic_types(self):
        """Test cache key generation with basic types"""
        key1 = self.cache._get_cache_key("test", 123, 45.6, True)
        key2 = self.cache._get_cache_key("test", 123, 45.6, True)
        key3 = self.cache._get_cache_key("test", 123, 45.6, False)

        # Same arguments should produce same key
        self.assertEqual(key1, key2)
        # Different arguments should produce different key
        self.assertNotEqual(key1, key3)

    def test_get_cache_key_with_objects(self):
        """Test cache key generation with objects"""
        class TestObj:
            def __init__(self, value):
                self.value = value

        obj1 = TestObj(10)
        obj2 = TestObj(10)
        obj3 = TestObj(20)

        key1 = self.cache._get_cache_key(obj1)
        key2 = self.cache._get_cache_key(obj2)
        key3 = self.cache._get_cache_key(obj3)

        # Objects with same attributes should produce same key
        self.assertEqual(key1, key2)
        # Objects with different attributes should produce different key
        self.assertNotEqual(key1, key3)

    def test_get_cache_key_with_kwargs(self):
        """Test cache key generation with keyword arguments"""
        key1 = self.cache._get_cache_key("test", foo="bar", baz=123)
        key2 = self.cache._get_cache_key("test", baz=123, foo="bar")  # Different order
        key3 = self.cache._get_cache_key("test", foo="bar", baz=456)

        # Order shouldn't matter for kwargs
        self.assertEqual(key1, key2)
        # Different values should produce different key
        self.assertNotEqual(key1, key3)

    def test_set_and_get(self):
        """Test setting and getting cached values"""
        test_key = "test_key_123"
        test_value = {"data": [1, 2, 3], "name": "test"}

        # Initially should return None
        self.assertIsNone(self.cache.get(test_key))

        # Set value
        self.cache.set(test_key, test_value)

        # Should retrieve same value
        retrieved = self.cache.get(test_key)
        self.assertEqual(retrieved, test_value)

    def test_get_corrupted_cache(self):
        """Test handling of corrupted cache files"""
        test_key = "corrupted_key"
        cache_file = os.path.join(self.cache.cache_dir, f"{test_key}.pkl")

        # Create corrupted cache file
        with open(cache_file, 'wb') as f:
            f.write(b"corrupted data")

        # Should return None and remove corrupted file
        result = self.cache.get(test_key)
        self.assertIsNone(result)
        self.assertFalse(os.path.exists(cache_file))

    def test_set_with_exception(self):
        """Test set method handles exceptions gracefully"""
        test_key = "test_key"

        # Make cache_dir read-only to cause exception
        os.chmod(self.cache.cache_dir, 0o444)

        try:
            # Should not raise exception
            self.cache.set(test_key, "test_value")
        except:
            self.fail("set() raised exception unexpectedly")
        finally:
            # Restore permissions
            os.chmod(self.cache.cache_dir, 0o755)

    def test_clear(self):
        """Test clearing cache"""
        # Add some cache files
        for i in range(3):
            self.cache.set(f"key_{i}", f"value_{i}")

        # Verify files exist
        cache_files = [f for f in os.listdir(self.cache.cache_dir) if f.endswith('.pkl')]
        self.assertEqual(len(cache_files), 3)

        # Clear cache
        self.cache.clear()

        # Verify files are removed
        cache_files = [f for f in os.listdir(self.cache.cache_dir) if f.endswith('.pkl')]
        self.assertEqual(len(cache_files), 0)


class TestCachedComputation(TestCase):
    """Test cached_computation decorator"""

    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = tempfile.mkdtemp()
        self.cache = ResultCache(cache_dir=self.temp_dir)

    def tearDown(self):
        """Clean up test fixtures"""
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    def test_decorator_caches_results(self):
        """Test that decorator caches function results"""
        call_count = 0

        @cached_computation(self.cache)
        def expensive_function(x, y):
            nonlocal call_count
            call_count += 1
            return x + y

        # First call
        result1 = expensive_function(2, 3)
        self.assertEqual(result1, 5)
        self.assertEqual(call_count, 1)

        # Second call with same arguments (should use cache)
        result2 = expensive_function(2, 3)
        self.assertEqual(result2, 5)
        self.assertEqual(call_count, 1)  # Function not called again

        # Call with different arguments
        result3 = expensive_function(3, 4)
        self.assertEqual(result3, 7)
        self.assertEqual(call_count, 2)

    def test_decorator_with_default_cache(self):
        """Test decorator with default cache instance"""
        @cached_computation()
        def test_function(x):
            return x * 2

        result1 = test_function(5)
        self.assertEqual(result1, 10)

        # Should have clear_cache method
        self.assertTrue(hasattr(test_function, 'clear_cache'))
        test_function.clear_cache()

    def test_clear_cache_method(self):
        """Test clear_cache method added by decorator"""
        @cached_computation(self.cache)
        def test_function(x):
            return x ** 2

        # Add some cached results
        result1 = test_function(2)
        result2 = test_function(3)

        # Clear cache
        test_function.clear_cache()

        # Verify cache is cleared
        cache_files = [f for f in os.listdir(self.cache.cache_dir) if f.endswith('.pkl')]
        self.assertEqual(len(cache_files), 0)


class TestCachedTreeDistance(TestCase):
    """Test cached_tree_distance function"""

    @patch('phykit.helpers.caching.pickle.loads')
    def test_cached_tree_distance(self, mock_loads):
        """Test tree distance caching"""
        # Create mock tree
        mock_tree = Mock()
        mock_tree.distance.return_value = 1.5
        mock_loads.return_value = mock_tree

        # Create fake pickle data
        tree_pickle = b"fake_pickle_data"

        # Call function
        distance = cached_tree_distance(tree_pickle, "tip1", "tip2")
        self.assertEqual(distance, 1.5)
        mock_tree.distance.assert_called_once_with("tip1", "tip2")

        # Call again (should use cache due to LRU)
        distance2 = cached_tree_distance(tree_pickle, "tip1", "tip2")
        self.assertEqual(distance2, 1.5)
        # Still only called once due to caching
        mock_tree.distance.assert_called_once()

        # Clear cache
        cached_tree_distance.cache_clear()


class TestAlignmentCache(TestCase):
    """Test AlignmentCache class"""

    def setUp(self):
        """Set up test fixtures"""
        self.cache = AlignmentCache()

    def test_column_cache(self):
        """Test column caching"""
        # Set column
        self.cache.set_column("hash1", 0, "ATCG")

        # Clear LRU cache to ensure fresh lookup
        self.cache.get_column.cache_clear()

        # Should retrieve column
        result = self.cache.get_column("hash1", 0)
        self.assertEqual(result, "ATCG")

        # Test with different key
        result2 = self.cache.get_column("hash2", 0)
        self.assertIsNone(result2)

    def test_stats_cache(self):
        """Test statistics caching"""
        # Set stats
        stats = {"mean": 0.5, "std": 0.1}
        self.cache.set_stats("hash1", "mean", stats)

        # Clear LRU cache to ensure fresh lookup
        self.cache.get_stats.cache_clear()

        # Should retrieve stats
        result = self.cache.get_stats("hash1", "mean")
        self.assertEqual(result, stats)

        # Test with different key
        result2 = self.cache.get_stats("hash2", "mean")
        self.assertIsNone(result2)

    def test_clear(self):
        """Test clearing all caches"""
        # Add some cached data
        self.cache.set_column("hash1", 0, "ATCG")
        self.cache.set_stats("hash1", "mean", {"mean": 0.5})

        # Clear caches
        self.cache.clear()

        # Verify caches are cleared
        self.assertEqual(len(self.cache._column_cache), 0)
        self.assertEqual(len(self.cache._stats_cache), 0)


class TestGlobalCaches(TestCase):
    """Test global cache instances"""

    def test_get_result_cache(self):
        """Test getting global result cache"""
        cache1 = get_result_cache()
        cache2 = get_result_cache()

        # Should return same instance
        self.assertIs(cache1, cache2)
        self.assertIsInstance(cache1, ResultCache)

    def test_get_alignment_cache(self):
        """Test getting global alignment cache"""
        cache1 = get_alignment_cache()
        cache2 = get_alignment_cache()

        # Should return same instance
        self.assertIs(cache1, cache2)
        self.assertIsInstance(cache1, AlignmentCache)