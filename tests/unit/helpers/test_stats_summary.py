"""
Unit tests for stats summary utilities
"""

import unittest
from unittest.mock import patch, Mock
import statistics as stat
import numpy as np
from io import StringIO

from phykit.helpers.stats_summary import (
    calculate_summary_statistics_from_arr,
    calculate_summary_statistics_from_dict,
    print_summary_statistics
)


class TestCalculateSummaryStatisticsFromArr(unittest.TestCase):
    """Test calculate_summary_statistics_from_arr function"""

    def test_basic_statistics(self):
        """Test calculation of basic statistics from array"""
        data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        stats = calculate_summary_statistics_from_arr(data)

        self.assertAlmostEqual(stats['mean'], 5.5)
        self.assertEqual(stats['median'], 5.5)
        self.assertAlmostEqual(stats['twenty_fifth'], 3.25)
        self.assertAlmostEqual(stats['seventy_fifth'], 7.75)
        self.assertEqual(stats['minimum'], 1)
        self.assertEqual(stats['maximum'], 10)
        self.assertAlmostEqual(stats['standard_deviation'], stat.stdev(data))
        self.assertAlmostEqual(stats['variance'], stat.variance(data))

    def test_single_value(self):
        """Test statistics with single value array"""
        data = [5]

        # Should raise StatisticsError for stdev and variance with single value
        with patch('sys.stdout', new_callable=StringIO) as mock_stdout:
            stats = calculate_summary_statistics_from_arr(data)
            output = mock_stdout.getvalue()

        # Check that error message was printed
        self.assertIn("There are no values to calculate summary statistics", output)
        # stats should be None since it wasn't set
        self.assertIsNone(stats)

    def test_identical_values(self):
        """Test statistics with identical values"""
        data = [3, 3, 3, 3, 3]
        stats = calculate_summary_statistics_from_arr(data)

        self.assertEqual(stats['mean'], 3)
        self.assertEqual(stats['median'], 3)
        self.assertEqual(stats['twenty_fifth'], 3)
        self.assertEqual(stats['seventy_fifth'], 3)
        self.assertEqual(stats['minimum'], 3)
        self.assertEqual(stats['maximum'], 3)
        self.assertEqual(stats['standard_deviation'], 0)
        self.assertEqual(stats['variance'], 0)

    def test_negative_values(self):
        """Test statistics with negative values"""
        data = [-5, -3, -1, 0, 1, 3, 5]
        stats = calculate_summary_statistics_from_arr(data)

        self.assertAlmostEqual(stats['mean'], 0)
        self.assertEqual(stats['median'], 0)
        self.assertEqual(stats['minimum'], -5)
        self.assertEqual(stats['maximum'], 5)

    def test_float_values(self):
        """Test statistics with float values"""
        data = [1.5, 2.5, 3.5, 4.5, 5.5]
        stats = calculate_summary_statistics_from_arr(data)

        self.assertAlmostEqual(stats['mean'], 3.5)
        self.assertEqual(stats['median'], 3.5)
        self.assertAlmostEqual(stats['twenty_fifth'], 2.5)
        self.assertAlmostEqual(stats['seventy_fifth'], 4.5)
        self.assertEqual(stats['minimum'], 1.5)
        self.assertEqual(stats['maximum'], 5.5)

    def test_percentiles(self):
        """Test percentile calculations"""
        data = list(range(1, 101))  # 1 to 100
        stats = calculate_summary_statistics_from_arr(data)

        self.assertAlmostEqual(stats['twenty_fifth'], 25.75)
        self.assertAlmostEqual(stats['seventy_fifth'], 75.25)
        self.assertEqual(stats['median'], 50.5)


class TestCalculateSummaryStatisticsFromDict(unittest.TestCase):
    """Test calculate_summary_statistics_from_dict function"""

    def test_basic_statistics(self):
        """Test calculation of basic statistics from dictionary"""
        data = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5}
        stats = calculate_summary_statistics_from_dict(data)

        self.assertEqual(stats['mean'], 3)
        self.assertEqual(stats['median'], 3)
        self.assertAlmostEqual(stats['twenty_fifth'], 2)
        self.assertAlmostEqual(stats['seventy_fifth'], 4)
        self.assertEqual(stats['minimum'], 1)
        self.assertEqual(stats['maximum'], 5)
        self.assertAlmostEqual(stats['standard_deviation'], stat.stdev([1, 2, 3, 4, 5]))
        self.assertAlmostEqual(stats['variance'], stat.variance([1, 2, 3, 4, 5]))

    def test_single_item_dict(self):
        """Test statistics with single item dictionary"""
        data = {'only': 42}

        with patch('sys.stdout', new_callable=StringIO) as mock_stdout:
            stats = calculate_summary_statistics_from_dict(data)
            output = mock_stdout.getvalue()

        # Check that error message was printed
        self.assertIn("There are no values to calculate summary statistics", output)
        # stats should be None since it wasn't set
        self.assertIsNone(stats)

    def test_empty_dict(self):
        """Test statistics with empty dictionary"""
        data = {}

        with patch('sys.stdout', new_callable=StringIO) as mock_stdout:
            stats = calculate_summary_statistics_from_dict(data)
            output = mock_stdout.getvalue()

        # Should print error message and return None
        self.assertIn("There are no values to calculate summary statistics", output)
        self.assertIsNone(stats)

    def test_mixed_values(self):
        """Test dictionary with mixed positive and negative values"""
        data = {'pos1': 10, 'neg1': -5, 'zero': 0, 'pos2': 15, 'neg2': -10}
        stats = calculate_summary_statistics_from_dict(data)

        self.assertEqual(stats['mean'], 2)
        self.assertEqual(stats['median'], 0)
        self.assertEqual(stats['minimum'], -10)
        self.assertEqual(stats['maximum'], 15)

    def test_float_values_in_dict(self):
        """Test dictionary with float values"""
        data = {'a': 1.1, 'b': 2.2, 'c': 3.3, 'd': 4.4, 'e': 5.5}
        stats = calculate_summary_statistics_from_dict(data)

        self.assertAlmostEqual(stats['mean'], 3.3)
        self.assertEqual(stats['median'], 3.3)
        self.assertAlmostEqual(stats['minimum'], 1.1)
        self.assertAlmostEqual(stats['maximum'], 5.5)


class TestPrintSummaryStatistics(unittest.TestCase):
    """Test print_summary_statistics function"""

    def test_print_normal_stats(self):
        """Test printing normal statistics"""
        stats = {
            'mean': 5.5555,
            'median': 5.4444,
            'twenty_fifth': 3.3333,
            'seventy_fifth': 7.7777,
            'minimum': 1.1111,
            'maximum': 10.9999,
            'standard_deviation': 2.2222,
            'variance': 4.9382
        }

        with patch('sys.stdout', new_callable=StringIO) as mock_stdout:
            print_summary_statistics(stats)
            output = mock_stdout.getvalue()

        # Check that all statistics are printed with proper rounding
        self.assertIn("mean: 5.5555", output)
        self.assertIn("median: 5.4444", output)
        self.assertIn("25th percentile: 3.3333", output)
        self.assertIn("75th percentile: 7.7777", output)
        self.assertIn("minimum: 1.1111", output)
        self.assertIn("maximum: 10.9999", output)
        self.assertIn("standard deviation: 2.2222", output)
        self.assertIn("variance: 4.9382", output)

    def test_print_integer_stats(self):
        """Test printing integer statistics"""
        stats = {
            'mean': 5,
            'median': 5,
            'twenty_fifth': 3,
            'seventy_fifth': 7,
            'minimum': 1,
            'maximum': 10,
            'standard_deviation': 2,
            'variance': 4
        }

        with patch('sys.stdout', new_callable=StringIO) as mock_stdout:
            print_summary_statistics(stats)
            output = mock_stdout.getvalue()

        # Check that integers are displayed with decimal places
        self.assertIn("mean: 5", output)
        self.assertIn("median: 5", output)
        self.assertIn("minimum: 1", output)
        self.assertIn("maximum: 10", output)

    def test_print_with_broken_pipe(self):
        """Test handling of BrokenPipeError"""
        stats = {
            'mean': 5.5,
            'median': 5.5,
            'twenty_fifth': 3.25,
            'seventy_fifth': 7.75,
            'minimum': 1,
            'maximum': 10,
            'standard_deviation': 2.87,
            'variance': 8.25
        }

        with patch('builtins.print', side_effect=BrokenPipeError):
            # Should not raise exception
            try:
                print_summary_statistics(stats)
            except BrokenPipeError:
                self.fail("BrokenPipeError was not caught")

    def test_print_extreme_values(self):
        """Test printing with extreme values"""
        stats = {
            'mean': 1e10,
            'median': 1e-10,
            'twenty_fifth': -1e10,
            'seventy_fifth': 0,
            'minimum': -1e20,
            'maximum': 1e20,
            'standard_deviation': 1e15,
            'variance': 1e30
        }

        with patch('sys.stdout', new_callable=StringIO) as mock_stdout:
            print_summary_statistics(stats)
            output = mock_stdout.getvalue()

        # Should handle scientific notation properly
        self.assertIn("mean:", output)
        self.assertIn("median:", output)
        self.assertIn("minimum:", output)
        self.assertIn("maximum:", output)

    def test_print_missing_keys(self):
        """Test printing with missing keys (should raise KeyError)"""
        stats = {
            'mean': 5.5,
            'median': 5.5
            # Missing other keys
        }

        with self.assertRaises(KeyError):
            print_summary_statistics(stats)