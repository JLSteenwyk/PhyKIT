"""
Unit tests for stats summary utilities
"""

import unittest
from unittest.mock import patch
import statistics as stat
from io import StringIO
import subprocess
import sys

import numpy as np

from phykit.helpers.stats_summary import (
    calculate_summary_statistics_from_arr,
    calculate_summary_statistics_from_dict,
    print_summary_statistics
)


class TestCalculateSummaryStatisticsFromArr(unittest.TestCase):
    """Test calculate_summary_statistics_from_arr function"""

    def test_module_import_does_not_import_numpy(self):
        code = (
            "import sys; "
            "import phykit.helpers.stats_summary; "
            "assert 'numpy' not in sys.modules"
        )

        subprocess.run([sys.executable, "-c", code], check=True)

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

    def test_no_values_message_is_batched_and_preserves_exact_text(self):
        """Test no-values diagnostics keep the legacy text in one print call."""
        data = [5]

        with patch('builtins.print') as mocked_print:
            stats = calculate_summary_statistics_from_arr(data)

        self.assertIsNone(stats)
        mocked_print.assert_called_once_with(
            "There are no values to calculate summary statistics for.\n\n"
            "Double check that the input alignment/phylogeny contains\n"
            "the properties you want to calculate summary statistics for."
        )

    def test_sized_no_value_array_input_does_not_import_numpy(self):
        code = (
            "import sys; "
            "import phykit.helpers.stats_summary as stats_summary; "
            "stats_summary.calculate_summary_statistics_from_arr([5]); "
            "assert 'numpy' not in sys.modules"
        )

        subprocess.run([sys.executable, "-c", code], check=True)

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

    def test_identical_values_skip_percentile_and_variance_reductions(self):
        """Test constant arrays avoid unnecessary summary reductions."""
        data = np.array([3.5, 3.5, 3.5, 3.5])
        with patch(
            'phykit.helpers.stats_summary.np.percentile',
            side_effect=AssertionError("constant values should skip percentile"),
        ), patch(
            'phykit.helpers.stats_summary.np.var',
            side_effect=AssertionError("constant values should skip variance"),
        ):
            stats = calculate_summary_statistics_from_arr(data)

        self.assertEqual(stats['mean'], 3.5)
        self.assertEqual(stats['median'], 3.5)
        self.assertEqual(stats['twenty_fifth'], 3.5)
        self.assertEqual(stats['seventy_fifth'], 3.5)
        self.assertEqual(stats['minimum'], 3.5)
        self.assertEqual(stats['maximum'], 3.5)
        self.assertEqual(stats['standard_deviation'], 0.0)
        self.assertEqual(stats['variance'], 0.0)

    def test_nonconstant_mean_uses_array_reduction(self):
        data = np.array([1.0, 2.0, 4.0, 8.0])
        with patch(
            'phykit.helpers.stats_summary.np.mean',
            side_effect=AssertionError("mean should use ndarray.mean"),
        ):
            stats = calculate_summary_statistics_from_arr(data)

        self.assertEqual(stats['mean'], 3.75)

    def test_exact_integer_mean_and_median_preserve_integer_type(self):
        """Test exact integer mean/median keep legacy CLI formatting."""
        data = [85, 85, 100, 100, 100, 100, 100]
        stats = calculate_summary_statistics_from_arr(data)

        self.assertEqual(stats['mean'], 95.71428571428571)
        self.assertEqual(stats['median'], 100)
        self.assertIs(type(stats['median']), int)
        self.assertEqual(stats['seventy_fifth'], 100.0)

        exact_stats = calculate_summary_statistics_from_arr([100, 100, 100])
        self.assertEqual(exact_stats['mean'], 100)
        self.assertIs(type(exact_stats['mean']), int)
        self.assertEqual(exact_stats['median'], 100)
        self.assertIs(type(exact_stats['median']), int)

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

    def test_numpy_array_input(self):
        """Test summary statistics from an existing NumPy array."""
        data = np.array([1.0, 2.0, 3.0, 4.0])
        stats = calculate_summary_statistics_from_arr(data)

        self.assertAlmostEqual(stats['mean'], 2.5)
        self.assertAlmostEqual(stats['median'], 2.5)
        self.assertAlmostEqual(stats['standard_deviation'], stat.stdev(data))
        self.assertAlmostEqual(stats['variance'], stat.variance(data))

    def test_percentiles_and_median_calculated_together(self):
        """Test percentile setup asks NumPy for quartiles and median in one pass."""
        data = np.array([1.0, 2.0, 3.0, 4.0])
        with patch(
            'phykit.helpers.stats_summary.np.percentile',
            return_value=np.array([1.75, 2.5, 3.25]),
        ) as mock_percentile:
            stats = calculate_summary_statistics_from_arr(data)

        mock_percentile.assert_called_once()
        np.testing.assert_array_equal(mock_percentile.call_args.args[1], [25, 50, 75])
        self.assertEqual(stats['twenty_fifth'], 1.75)
        self.assertEqual(stats['median'], 2.5)
        self.assertEqual(stats['seventy_fifth'], 3.25)

    def test_standard_deviation_reuses_variance_reduction(self):
        """Test standard deviation is derived without a second NumPy reduction."""
        data = np.array([1.0, 2.0, 3.0, 4.0])
        with patch(
            'phykit.helpers.stats_summary.np.std',
            side_effect=AssertionError("standard deviation should reuse variance"),
        ):
            stats = calculate_summary_statistics_from_arr(data)

        self.assertAlmostEqual(stats['variance'], stat.variance(data))
        self.assertAlmostEqual(stats['standard_deviation'], stat.stdev(data))


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

    def test_sized_no_value_dict_input_does_not_import_numpy(self):
        code = (
            "import sys; "
            "import phykit.helpers.stats_summary as stats_summary; "
            "stats_summary.calculate_summary_statistics_from_dict({}); "
            "assert 'numpy' not in sys.modules"
        )

        subprocess.run([sys.executable, "-c", code], check=True)

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

    def test_numeric_dict_uses_fromiter(self):
        """Test numeric dictionary summaries use the one-pass NumPy conversion."""
        data = {'a': 1.0, 'b': 2.0, 'c': 3.0}
        with patch('phykit.helpers.stats_summary.np.fromiter') as mock_fromiter:
            mock_fromiter.return_value = np.array([1.0, 2.0, 3.0])
            stats = calculate_summary_statistics_from_dict(data)

        mock_fromiter.assert_called_once()
        self.assertEqual(stats['mean'], 2.0)
        self.assertEqual(stats['minimum'], 1.0)
        self.assertEqual(stats['maximum'], 3.0)

    def test_dict_fromiter_fallback_preserves_generic_values(self):
        """Test fallback path when one-pass numeric conversion is unavailable."""
        data = {'a': 1.0, 'b': 2.0, 'c': 3.0}
        with patch(
            'phykit.helpers.stats_summary.np.fromiter',
            side_effect=TypeError,
        ) as mock_fromiter:
            stats = calculate_summary_statistics_from_dict(data)

        mock_fromiter.assert_called_once()
        self.assertEqual(stats['mean'], 2.0)
        self.assertEqual(stats['minimum'], 1.0)
        self.assertEqual(stats['maximum'], 3.0)


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

    def test_print_batches_summary_output(self):
        """Test summary statistics are emitted with one print call."""
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

        with patch('builtins.print') as mocked_print:
            print_summary_statistics(stats)

        mocked_print.assert_called_once_with(
            "mean: 5.5555\n"
            "median: 5.4444\n"
            "25th percentile: 3.3333\n"
            "75th percentile: 7.7777\n"
            "minimum: 1.1111\n"
            "maximum: 10.9999\n"
            "standard deviation: 2.2222\n"
            "variance: 4.9382"
        )

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
