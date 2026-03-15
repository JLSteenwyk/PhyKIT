"""
Unit tests for CovaryingEvolutionaryRates class
"""

import unittest
import builtins
from unittest.mock import Mock, MagicMock, patch
from concurrent.futures import Future
from argparse import Namespace
import pytest
import types
import sys

from phykit.services.tree.covarying_evolutionary_rates import CovaryingEvolutionaryRates


class TestCovaryingEvolutionaryRates(unittest.TestCase):
    """Test CovaryingEvolutionaryRates class"""

    def setUp(self):
        """Set up test fixtures"""
        self.args = Namespace(
            tree_zero="tree0.tre",
            tree_one="tree1.tre",
            reference="ref.tre",
            verbose=False
        )
        self.cov_rates = CovaryingEvolutionaryRates(self.args)

    def test_init(self):
        """Test initialization"""
        self.assertEqual(self.cov_rates.tree_file_path, "tree0.tre")
        self.assertEqual(self.cov_rates.tree1_file_path, "tree1.tre")
        self.assertEqual(self.cov_rates.reference, "ref.tre")
        self.assertFalse(self.cov_rates.verbose)

    def test_process_args(self):
        """Test argument processing"""
        args = Namespace(
            tree_zero="test0.tre",
            tree_one="test1.tre",
            reference="test_ref.tre",
            verbose=True
        )
        processed = self.cov_rates.process_args(args)

        self.assertEqual(processed["tree_file_path"], "test0.tre")
        self.assertEqual(processed["tree1_file_path"], "test1.tre")
        self.assertEqual(processed["reference"], "test_ref.tre")
        self.assertTrue(processed["verbose"])

    def test_get_indices_of_outlier_branch_lengths(self):
        """Test outlier detection in branch lengths"""
        # Test with normal values
        corr_branch_lengths = [1.0, 2.0, 3.0, 4.0]
        outlier_indices = []
        result = self.cov_rates.get_indices_of_outlier_branch_lengths(
            corr_branch_lengths, outlier_indices
        )
        self.assertEqual(result, [])

        # Test with outliers (>5)
        corr_branch_lengths = [1.0, 6.0, 3.0, -8.0, 4.0]
        outlier_indices = []
        result = self.cov_rates.get_indices_of_outlier_branch_lengths(
            corr_branch_lengths, outlier_indices
        )
        self.assertIn(1, result)  # 6.0 > 5
        self.assertIn(3, result)  # abs(-8.0) > 5

        # Test with NaN values
        corr_branch_lengths = [1.0, float('nan'), 3.0, 4.0]
        outlier_indices = []
        result = self.cov_rates.get_indices_of_outlier_branch_lengths(
            corr_branch_lengths, outlier_indices
        )
        self.assertIn(1, result)  # NaN value

        # Test with existing outlier indices
        corr_branch_lengths = [1.0, 2.0, 7.0, 4.0]
        outlier_indices = [0, 4]
        result = self.cov_rates.get_indices_of_outlier_branch_lengths(
            corr_branch_lengths, outlier_indices
        )
        self.assertIn(0, result)  # Existing
        self.assertIn(2, result)  # New (7.0 > 5)
        self.assertIn(4, result)  # Existing

    def test_remove_outliers_based_on_indices(self):
        """Test outlier removal"""
        # Test with numeric list
        corr_branch_lengths = [1.0, 2.0, 3.0, 4.0, 5.0]
        outlier_indices = [1, 3]
        result = self.cov_rates.remove_outliers_based_on_indices(
            corr_branch_lengths, outlier_indices
        )
        self.assertEqual(result, [1.0, 3.0, 5.0])

        # Test with empty outlier indices
        corr_branch_lengths = [1.0, 2.0, 3.0]
        outlier_indices = []
        result = self.cov_rates.remove_outliers_based_on_indices(
            corr_branch_lengths, outlier_indices
        )
        self.assertEqual(result, [1.0, 2.0, 3.0])

        # Test with list of lists (tip names)
        corr_branch_lengths = [['a'], ['b', 'c'], ['d'], ['e', 'f']]
        outlier_indices = [1, 2]
        result = self.cov_rates.remove_outliers_based_on_indices(
            corr_branch_lengths, outlier_indices
        )
        self.assertEqual(result, [['a'], ['e', 'f']])

    def test_prune_tips(self):
        """Test tree pruning"""
        # Create mock tree
        mock_tree = Mock()
        tips_to_prune = ['tip1', 'tip2', 'tip3']

        result = self.cov_rates.prune_tips(mock_tree, tips_to_prune)

        # Check that prune was called for each tip
        self.assertEqual(mock_tree.prune.call_count, 3)
        mock_tree.prune.assert_any_call('tip1')
        mock_tree.prune.assert_any_call('tip2')
        mock_tree.prune.assert_any_call('tip3')
        self.assertEqual(result, mock_tree)

    @patch('pickle.loads')
    def test_process_terminal_batch(self, mock_loads):
        """Test static method for processing terminal batch"""
        # Create mock trees
        mock_tree0 = Mock()
        mock_tree1 = Mock()

        # Setup common_ancestor returns
        mock_newtree0 = Mock()
        mock_newtree0.branch_length = 2.0
        mock_newtree1 = Mock()
        mock_newtree1.branch_length = 4.0

        mock_tree0.common_ancestor.return_value = mock_newtree0
        mock_tree1.common_ancestor.return_value = mock_newtree1

        # Setup mock pickle loads
        mock_loads.side_effect = [mock_tree0, mock_tree1]

        # Create fake pickle data
        tree0_pickle = b'fake_tree0_pickle'
        tree1_pickle = b'fake_tree1_pickle'

        # Create terminal data
        terminals_data = [
            ('terminal1', 1.0, ['tip1', 'tip2']),
            ('terminal2', 2.0, ['tip3'])
        ]

        # Call the method
        results = CovaryingEvolutionaryRates._process_terminal_batch(
            tree0_pickle, tree1_pickle, terminals_data
        )

        # Check results
        self.assertEqual(len(results), 2)
        self.assertEqual(results[0][0], 2.0)  # bl0 = 2.0 / 1.0
        self.assertEqual(results[0][1], 4.0)  # bl1 = 4.0 / 1.0
        self.assertEqual(results[0][2], ['tip1', 'tip2'])

        self.assertEqual(results[1][0], 1.0)  # bl0 = 2.0 / 2.0
        self.assertEqual(results[1][1], 2.0)  # bl1 = 4.0 / 2.0
        self.assertEqual(results[1][2], ['tip3'])

    @patch('pickle.loads')
    def test_process_terminal_batch_with_none_branch_length(self, mock_loads):
        """Test terminal batch processing with None branch lengths"""
        mock_tree0 = Mock()
        mock_tree1 = Mock()

        mock_newtree0 = Mock()
        mock_newtree0.branch_length = None  # None branch length
        mock_newtree1 = Mock()
        mock_newtree1.branch_length = 4.0

        mock_tree0.common_ancestor.return_value = mock_newtree0
        mock_tree1.common_ancestor.return_value = mock_newtree1

        mock_loads.side_effect = [mock_tree0, mock_tree1]

        tree0_pickle = b'fake_tree0_pickle'
        tree1_pickle = b'fake_tree1_pickle'

        terminals_data = [('terminal1', 1.0, ['tip1'])]

        results = CovaryingEvolutionaryRates._process_terminal_batch(
            tree0_pickle, tree1_pickle, terminals_data
        )

        # Should return empty results since one branch length is None
        self.assertEqual(len(results), 0)

    @patch('pickle.loads')
    def test_process_terminal_batch_with_exception(self, mock_loads):
        """Test terminal batch processing with exception"""
        mock_tree0 = Mock()
        mock_tree1 = Mock()

        # Make common_ancestor raise exception
        mock_tree0.common_ancestor.side_effect = Exception("Test error")

        mock_loads.side_effect = [mock_tree0, mock_tree1]

        tree0_pickle = b'fake_tree0_pickle'
        tree1_pickle = b'fake_tree1_pickle'

        terminals_data = [('terminal1', 1.0, ['tip1'])]

        results = CovaryingEvolutionaryRates._process_terminal_batch(
            tree0_pickle, tree1_pickle, terminals_data
        )

        # Should return empty results due to exception
        self.assertEqual(len(results), 0)

    @patch('pickle.loads')
    def test_process_nonterminal_batch(self, mock_loads):
        """Test static method for processing nonterminal batch"""
        mock_tree0 = Mock()
        mock_tree1 = Mock()

        mock_newtree0 = Mock()
        mock_newtree0.branch_length = 3.0
        mock_newtree1 = Mock()
        mock_newtree1.branch_length = 6.0

        mock_tree0.common_ancestor.return_value = mock_newtree0
        mock_tree1.common_ancestor.return_value = mock_newtree1

        mock_loads.side_effect = [mock_tree0, mock_tree1]

        tree0_pickle = b'fake_tree0_pickle'
        tree1_pickle = b'fake_tree1_pickle'

        nonterminals_data = [
            (['tip1', 'tip2'], 1.5),
            (['tip3', 'tip4'], 2.0)
        ]

        results = CovaryingEvolutionaryRates._process_nonterminal_batch(
            tree0_pickle, tree1_pickle, nonterminals_data
        )

        # Check results
        self.assertEqual(len(results), 2)
        self.assertEqual(results[0][0], 2.0)  # bl0 = 3.0 / 1.5
        self.assertEqual(results[0][1], 4.0)  # bl1 = 6.0 / 1.5
        self.assertEqual(results[0][2], ['tip1', 'tip2'])

    @patch('pickle.loads')
    def test_process_nonterminal_batch_with_none_branch_length(self, mock_loads):
        """Test nonterminal batch processing with None branch lengths"""
        mock_tree0 = Mock()
        mock_tree1 = Mock()

        mock_newtree0 = Mock()
        mock_newtree0.branch_length = 3.0
        mock_newtree1 = Mock()
        mock_newtree1.branch_length = None  # None branch length

        mock_tree0.common_ancestor.return_value = mock_newtree0
        mock_tree1.common_ancestor.return_value = mock_newtree1

        mock_loads.side_effect = [mock_tree0, mock_tree1]

        tree0_pickle = b'fake_tree0_pickle'
        tree1_pickle = b'fake_tree1_pickle'

        nonterminals_data = [(['tip1', 'tip2'], 1.5)]

        results = CovaryingEvolutionaryRates._process_nonterminal_batch(
            tree0_pickle, tree1_pickle, nonterminals_data
        )

        # Should return empty results since one branch length is None
        self.assertEqual(len(results), 0)

    @patch('pickle.dumps')
    @patch('phykit.services.tree.covarying_evolutionary_rates.ProcessPoolExecutor')
    def test_correct_branch_lengths_parallel(self, mock_executor_class, mock_pickle_dumps):
        """Test correct_branch_lengths with parallel processing (large dataset)"""
        # Create mock trees
        mock_t0 = Mock()
        mock_t1 = Mock()
        mock_sp = Mock()

        # Create many terminals and nonterminals to trigger parallel processing
        terminals = []
        for i in range(30):
            term = Mock()
            term.name = f'terminal{i}'
            term.branch_length = 1.0 + i * 0.1
            terminals.append(term)

        nonterminals = []
        for i in range(25):
            nonterm = Mock()
            nonterm.branch_length = 2.0 + i * 0.1
            nonterminals.append(nonterm)

        mock_sp.get_terminals.return_value = terminals
        mock_sp.get_nonterminals.return_value = nonterminals

        # Mock get_tip_names_from_tree
        self.cov_rates.get_tip_names_from_tree = Mock(side_effect=lambda x: [f'tip_{x.name}'] if hasattr(x, 'name') else ['tip_x'])

        # Mock pickle dumps
        mock_pickle_dumps.return_value = b'pickled_tree'

        # Mock executor
        mock_executor = MagicMock()
        mock_executor_class.return_value.__enter__.return_value = mock_executor

        # Create mock futures
        future1 = Mock(spec=Future)
        future1.result.return_value = [(1.0, 2.0, ['tip1']), (1.5, 2.5, ['tip2'])]
        future2 = Mock(spec=Future)
        future2.result.return_value = [(2.0, 3.0, ['tip3'])]

        # Create enough futures for all submit calls
        mock_futures = [future1, future2] + [Mock(spec=Future) for _ in range(10)]
        for f in mock_futures[2:]:
            f.result.return_value = []
        mock_executor.submit.side_effect = mock_futures

        # Mock as_completed
        with patch('phykit.services.tree.covarying_evolutionary_rates.as_completed') as mock_as_completed:
            mock_as_completed.return_value = [future1, future2]

            l0, l1, tip_names = self.cov_rates.correct_branch_lengths(mock_t0, mock_t1, mock_sp)

        # Check results
        self.assertEqual(len(l0), 3)
        self.assertEqual(len(l1), 3)
        self.assertEqual(len(tip_names), 3)

        self.assertEqual(l0, [1.0, 1.5, 2.0])
        self.assertEqual(l1, [2.0, 2.5, 3.0])
        self.assertEqual(tip_names, [['tip1'], ['tip2'], ['tip3']])

    def test_correct_branch_lengths_sequential(self):
        """Test correct_branch_lengths with sequential processing (small dataset)"""
        # Create mock trees
        mock_t0 = Mock()
        mock_t1 = Mock()
        mock_sp = Mock()

        # Create few terminals (< 50 total) to trigger sequential processing
        terminal1 = Mock()
        terminal1.name = 'terminal1'
        terminal1.branch_length = 2.0

        terminal2 = Mock()
        terminal2.name = 'terminal2'
        terminal2.branch_length = 1.0

        nonterminal1 = Mock()
        nonterminal1.branch_length = 3.0

        mock_sp.get_terminals.return_value = [terminal1, terminal2]
        mock_sp.get_nonterminals.return_value = [nonterminal1]

        # Mock get_tip_names_from_tree
        def mock_get_tips(tree):
            if hasattr(tree, 'name'):
                return [tree.name]
            else:
                return ['nonterminal_tips']

        self.cov_rates.get_tip_names_from_tree = Mock(side_effect=mock_get_tips)

        # Mock common_ancestor returns
        mock_newtree0 = Mock()
        mock_newtree0.branch_length = 4.0
        mock_newtree1 = Mock()
        mock_newtree1.branch_length = 6.0

        mock_t0.common_ancestor.return_value = mock_newtree0
        mock_t1.common_ancestor.return_value = mock_newtree1

        l0, l1, tip_names = self.cov_rates.correct_branch_lengths(mock_t0, mock_t1, mock_sp)

        # Check results - terminals
        self.assertIn(2.0, l0)  # 4.0 / 2.0
        self.assertIn(3.0, l1)  # 6.0 / 2.0
        self.assertIn(4.0, l0)  # 4.0 / 1.0
        self.assertIn(6.0, l1)  # 6.0 / 1.0

        # Check that nonterminal was processed
        self.assertEqual(len(l0), 3)
        self.assertEqual(len(l1), 3)
        self.assertEqual(len(tip_names), 3)

    def test_correct_branch_lengths_sequential_with_exception(self):
        """Test correct_branch_lengths sequential with exception handling"""
        mock_t0 = Mock()
        mock_t1 = Mock()
        mock_sp = Mock()

        terminal1 = Mock()
        terminal1.name = 'terminal1'
        terminal1.branch_length = 2.0

        mock_sp.get_terminals.return_value = [terminal1]
        mock_sp.get_nonterminals.return_value = []

        self.cov_rates.get_tip_names_from_tree = Mock(return_value=['tip1'])

        # Make common_ancestor raise exception
        mock_t0.common_ancestor.side_effect = Exception("Test error")

        l0, l1, tip_names = self.cov_rates.correct_branch_lengths(mock_t0, mock_t1, mock_sp)

        # The function appends tip_names before processing, so it will have the tip names but no branch lengths
        self.assertEqual(len(l0), 0)  # No branch lengths due to exception
        self.assertEqual(len(l1), 0)
        self.assertEqual(len(tip_names), 1)  # tip_names was appended before exception

    @patch('builtins.print')
    def test_run_non_verbose(self, mock_print):
        """Test run method in non-verbose mode"""
        # Create mock trees
        mock_tree0 = Mock()
        mock_tree1 = Mock()
        mock_ref_tree = Mock()

        # Mock tree reading methods
        self.cov_rates.read_tree_file = Mock(return_value=mock_tree0)
        self.cov_rates.read_tree1_file = Mock(return_value=mock_tree1)
        self.cov_rates.read_reference_tree_file = Mock(return_value=mock_ref_tree)

        # Mock tip names
        self.cov_rates.get_tip_names_from_tree = Mock(side_effect=[
            ['tip1', 'tip2', 'tip3'],  # tree_zero
            ['tip1', 'tip2', 'tip4'],  # tree_one
            ['tip1', 'tip2', 'tip5'],  # tree_ref
        ])

        # Mock shared_tips
        self.cov_rates.shared_tips = Mock(return_value=['tip1', 'tip2'])

        # Mock prune_tips
        self.cov_rates.prune_tips = Mock(side_effect=lambda tree, tips: tree)

        # Mock correct_branch_lengths
        self.cov_rates.correct_branch_lengths = Mock(return_value=(
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 3.0, 4.0, 5.0],
            [['tip1'], ['tip2'], ['tip1', 'tip2'], ['tip3']]
        ))

        # Mock outlier methods
        self.cov_rates.get_indices_of_outlier_branch_lengths = Mock(return_value=[])
        self.cov_rates.remove_outliers_based_on_indices = Mock(side_effect=lambda x, y: x)

        # Run the method
        self.cov_rates.run()

        # Check that the correlation was printed
        mock_print.assert_called_once()
        printed = mock_print.call_args[0][0]
        # Should print correlation and p-value
        self.assertIn('\t', printed)
        parts = printed.split('\t')
        self.assertEqual(len(parts), 2)

    @patch('builtins.print')
    def test_run_verbose(self, mock_print):
        """Test run method in verbose mode"""
        self.cov_rates.verbose = True

        # Create mock trees
        mock_tree0 = Mock()
        mock_tree1 = Mock()
        mock_ref_tree = Mock()

        # Mock tree reading methods
        self.cov_rates.read_tree_file = Mock(return_value=mock_tree0)
        self.cov_rates.read_tree1_file = Mock(return_value=mock_tree1)
        self.cov_rates.read_reference_tree_file = Mock(return_value=mock_ref_tree)

        # Mock tip names
        self.cov_rates.get_tip_names_from_tree = Mock(side_effect=[
            ['tip1', 'tip2', 'tip3'],  # tree_zero
            ['tip1', 'tip2', 'tip4'],  # tree_one
            ['tip1', 'tip2', 'tip5'],  # tree_ref
        ])

        # Mock shared_tips
        self.cov_rates.shared_tips = Mock(return_value=['tip1', 'tip2'])

        # Mock prune_tips
        self.cov_rates.prune_tips = Mock(side_effect=lambda tree, tips: tree)

        # Mock correct_branch_lengths with standardized values
        self.cov_rates.correct_branch_lengths = Mock(return_value=(
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [['tip1'], ['tip2'], ['tip1', 'tip2']]
        ))

        # Mock outlier methods
        self.cov_rates.get_indices_of_outlier_branch_lengths = Mock(return_value=[])
        self.cov_rates.remove_outliers_based_on_indices = Mock(side_effect=lambda x, y: x)

        # Run the method
        self.cov_rates.run()

        # In verbose mode, should print each value pair
        self.assertEqual(mock_print.call_count, 3)  # One for each pair

        # Check that semicolon-separated tips are printed
        for call_obj in mock_print.call_args_list:
            printed = call_obj[0][0]
            self.assertIn('\t', printed)  # Should have tab separators

    @patch('builtins.print')
    def test_run_with_outliers(self, mock_print):
        """Test run method with outlier removal"""
        # Create mock trees
        mock_tree0 = Mock()
        mock_tree1 = Mock()
        mock_ref_tree = Mock()

        # Mock tree reading methods
        self.cov_rates.read_tree_file = Mock(return_value=mock_tree0)
        self.cov_rates.read_tree1_file = Mock(return_value=mock_tree1)
        self.cov_rates.read_reference_tree_file = Mock(return_value=mock_ref_tree)

        # Mock tip names
        self.cov_rates.get_tip_names_from_tree = Mock(side_effect=[
            ['tip1', 'tip2'],
            ['tip1', 'tip2'],
            ['tip1', 'tip2'],
        ])

        # Mock shared_tips
        self.cov_rates.shared_tips = Mock(return_value=['tip1', 'tip2'])

        # Mock prune_tips
        self.cov_rates.prune_tips = Mock(side_effect=lambda tree, tips: tree)

        # Mock correct_branch_lengths with values that will trigger outlier removal
        self.cov_rates.correct_branch_lengths = Mock(return_value=(
            [1.0, 10.0, 3.0, 4.0],  # 10.0 will be an outlier
            [2.0, 15.0, 4.0, 5.0],  # 15.0 will be an outlier
            [['tip1'], ['outlier'], ['tip2'], ['tip3']]
        ))

        # Mock outlier detection to return index 1 (the outlier values)
        def mock_get_outliers(values, existing):
            if any(abs(v) > 5 for v in values):
                return [1]  # Index of outlier
            return existing

        self.cov_rates.get_indices_of_outlier_branch_lengths = Mock(side_effect=mock_get_outliers)

        # Mock outlier removal
        def mock_remove_outliers(values, indices):
            if not indices:
                return values
            return [v for i, v in enumerate(values) if i not in indices]

        self.cov_rates.remove_outliers_based_on_indices = Mock(side_effect=mock_remove_outliers)

        # Run the method
        self.cov_rates.run()

        # Check that outlier methods were called
        self.assertEqual(self.cov_rates.get_indices_of_outlier_branch_lengths.call_count, 2)
        self.assertEqual(self.cov_rates.remove_outliers_based_on_indices.call_count, 3)

        # Check that correlation was calculated and printed
        mock_print.assert_called_once()

    @patch('builtins.print')
    def test_run_with_broken_pipe_error(self, mock_print):
        """Test run method handling BrokenPipeError"""
        # Setup mocks similar to test_run_non_verbose
        mock_tree0 = Mock()
        mock_tree1 = Mock()
        mock_ref_tree = Mock()

        self.cov_rates.read_tree_file = Mock(return_value=mock_tree0)
        self.cov_rates.read_tree1_file = Mock(return_value=mock_tree1)
        self.cov_rates.read_reference_tree_file = Mock(return_value=mock_ref_tree)

        self.cov_rates.get_tip_names_from_tree = Mock(side_effect=[
            ['tip1'], ['tip1'], ['tip1']
        ])

        self.cov_rates.shared_tips = Mock(return_value=['tip1'])
        self.cov_rates.prune_tips = Mock(side_effect=lambda tree, tips: tree)

        self.cov_rates.correct_branch_lengths = Mock(return_value=(
            [1.0, 2.0], [2.0, 3.0], [['tip1'], ['tip2']]
        ))

        self.cov_rates.get_indices_of_outlier_branch_lengths = Mock(return_value=[])
        self.cov_rates.remove_outliers_based_on_indices = Mock(side_effect=lambda x, y: x)

        # Make print raise BrokenPipeError
        mock_print.side_effect = BrokenPipeError()

        # Should not raise exception
        try:
            self.cov_rates.run()
        except BrokenPipeError:
            self.fail("BrokenPipeError was not caught")

    def test_plot_covarying_rates_scatter_importerror(self):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name.startswith("matplotlib"):
                raise ImportError("no matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        with patch("builtins.__import__", side_effect=fake_import):
            with self.assertRaises(SystemExit) as exc:
                self.cov_rates._plot_covarying_rates_scatter([0.1, 0.2], [0.3, 0.4], [0.5, 0.01])
        self.assertEqual(exc.exception.code, 2)

    @patch("phykit.services.tree.covarying_evolutionary_rates.print_json")
    def test_run_json_non_verbose(self, mock_json):
        self.cov_rates.json_output = True
        self.cov_rates.verbose = False
        self.cov_rates.plot = False

        self.cov_rates.read_tree_file = Mock(return_value=Mock())
        self.cov_rates.read_tree1_file = Mock(return_value=Mock())
        self.cov_rates.read_reference_tree_file = Mock(return_value=Mock())
        self.cov_rates.get_tip_names_from_tree = Mock(side_effect=[["a", "b"], ["a", "b"], ["a", "b"]])
        self.cov_rates.shared_tips = Mock(return_value=["a", "b"])
        self.cov_rates.prune_tips = Mock(side_effect=lambda tree, tips: tree)
        self.cov_rates.correct_branch_lengths = Mock(return_value=([1.0, 2.0], [2.0, 3.0], [["a"], ["b"]]))
        self.cov_rates.get_indices_of_outlier_branch_lengths = Mock(return_value=[])
        self.cov_rates.remove_outliers_based_on_indices = Mock(side_effect=lambda x, y: x)

        self.cov_rates.run()
        payload = mock_json.call_args.args[0]
        self.assertFalse(payload["verbose"])
        self.assertIn("correlation", payload)
        self.assertIn("p_value", payload)

    @patch("phykit.services.tree.covarying_evolutionary_rates.print_json")
    def test_run_json_verbose_with_plot(self, mock_json):
        self.cov_rates.json_output = True
        self.cov_rates.verbose = True
        self.cov_rates.plot = True
        self.cov_rates.plot_output = "cover.png"
        self.cov_rates._plot_covarying_rates_scatter = Mock()

        self.cov_rates.read_tree_file = Mock(return_value=Mock())
        self.cov_rates.read_tree1_file = Mock(return_value=Mock())
        self.cov_rates.read_reference_tree_file = Mock(return_value=Mock())
        self.cov_rates.get_tip_names_from_tree = Mock(side_effect=[["a", "b"], ["a", "b"], ["a", "b"]])
        self.cov_rates.shared_tips = Mock(return_value=["a", "b"])
        self.cov_rates.prune_tips = Mock(side_effect=lambda tree, tips: tree)
        self.cov_rates.correct_branch_lengths = Mock(return_value=([1.0, 2.0], [2.0, 3.0], [["a"], ["b"]]))
        self.cov_rates.get_indices_of_outlier_branch_lengths = Mock(return_value=[])
        self.cov_rates.remove_outliers_based_on_indices = Mock(side_effect=lambda x, y: x)

        self.cov_rates.run()
        self.cov_rates._plot_covarying_rates_scatter.assert_called_once()
        payload = mock_json.call_args.args[0]
        self.assertTrue(payload["verbose"])
        self.assertEqual(payload["rows"], payload["branches"])
        self.assertEqual(payload["plot_output"], "cover.png")

    @patch("phykit.services.tree.covarying_evolutionary_rates.print_json")
    def test_run_json_non_verbose_with_plot_output(self, mock_json):
        self.cov_rates.json_output = True
        self.cov_rates.verbose = False
        self.cov_rates.plot = True
        self.cov_rates.plot_output = "plot-file.png"
        self.cov_rates._plot_covarying_rates_scatter = Mock()

        self.cov_rates.read_tree_file = Mock(return_value=Mock())
        self.cov_rates.read_tree1_file = Mock(return_value=Mock())
        self.cov_rates.read_reference_tree_file = Mock(return_value=Mock())
        self.cov_rates.get_tip_names_from_tree = Mock(side_effect=[["a", "b"], ["a", "b"], ["a", "b"]])
        self.cov_rates.shared_tips = Mock(return_value=["a", "b"])
        self.cov_rates.prune_tips = Mock(side_effect=lambda tree, tips: tree)
        self.cov_rates.correct_branch_lengths = Mock(return_value=([1.0, 2.0], [2.0, 3.0], [["a"], ["b"]]))
        self.cov_rates.get_indices_of_outlier_branch_lengths = Mock(return_value=[])
        self.cov_rates.remove_outliers_based_on_indices = Mock(side_effect=lambda x, y: x)

        self.cov_rates.run()
        payload = mock_json.call_args.args[0]
        self.assertEqual(payload["plot_output"], "plot-file.png")

    @patch("builtins.print")
    def test_run_non_verbose_with_plot_prints_saved_path(self, mock_print):
        self.cov_rates.plot = True
        self.cov_rates.plot_output = "saved.png"
        self.cov_rates._plot_covarying_rates_scatter = Mock()

        self.cov_rates.read_tree_file = Mock(return_value=Mock())
        self.cov_rates.read_tree1_file = Mock(return_value=Mock())
        self.cov_rates.read_reference_tree_file = Mock(return_value=Mock())
        self.cov_rates.get_tip_names_from_tree = Mock(side_effect=[["a", "b"], ["a", "b"], ["a", "b"]])
        self.cov_rates.shared_tips = Mock(return_value=["a", "b"])
        self.cov_rates.prune_tips = Mock(side_effect=lambda tree, tips: tree)
        self.cov_rates.correct_branch_lengths = Mock(return_value=([1.0, 2.0], [2.0, 3.0], [["a"], ["b"]]))
        self.cov_rates.get_indices_of_outlier_branch_lengths = Mock(return_value=[])
        self.cov_rates.remove_outliers_based_on_indices = Mock(side_effect=lambda x, y: x)

        self.cov_rates.run()
        printed_lines = [c.args[0] for c in mock_print.call_args_list]
        self.assertEqual(len(printed_lines), 2)
        self.assertTrue(any("Saved covarying rates plot: saved.png" == line for line in printed_lines))

    @patch("pickle.loads")
    def test_process_terminal_batch_reraises_system_exit(self, mock_loads):
        tree0 = Mock()
        tree1 = Mock()
        tree0.common_ancestor.side_effect = SystemExit(2)
        mock_loads.side_effect = [tree0, tree1]

        with self.assertRaises(SystemExit):
            CovaryingEvolutionaryRates._process_terminal_batch(
                b"t0", b"t1", [("terminal", 1.0, ["a"])]
            )

    @patch("pickle.loads")
    def test_process_nonterminal_batch_reraises_keyboard_interrupt(self, mock_loads):
        tree0 = Mock()
        tree1 = Mock()
        tree0.common_ancestor.side_effect = KeyboardInterrupt()
        mock_loads.side_effect = [tree0, tree1]

        with self.assertRaises(KeyboardInterrupt):
            CovaryingEvolutionaryRates._process_nonterminal_batch(
                b"t0", b"t1", [(["a", "b"], 1.0)]
            )

    def test_correct_branch_lengths_parallel_fallback(self):
        mock_t0 = Mock()
        mock_t1 = Mock()
        mock_sp = Mock()

        terminals = []
        for i in range(30):
            term = Mock()
            term.name = f"t{i}"
            term.branch_length = 1.0
            terminals.append(term)
        nonterminals = []
        for _ in range(25):
            nonterm = Mock()
            nonterm.branch_length = 1.0
            nonterminals.append(nonterm)

        mock_sp.get_terminals.return_value = terminals
        mock_sp.get_nonterminals.return_value = nonterminals
        self.cov_rates.get_tip_names_from_tree = Mock(side_effect=lambda node: [getattr(node, "name", "n")])

        node0 = Mock()
        node0.branch_length = 2.0
        node1 = Mock()
        node1.branch_length = 4.0
        mock_t0.common_ancestor.return_value = node0
        mock_t1.common_ancestor.return_value = node1

        with patch("phykit.services.tree.covarying_evolutionary_rates.pickle.dumps", return_value=b"pickled"), \
             patch("phykit.services.tree.covarying_evolutionary_rates.ProcessPoolExecutor", side_effect=OSError()):
            l0, l1, tip_names = self.cov_rates.correct_branch_lengths(mock_t0, mock_t1, mock_sp)

        self.assertTrue(len(l0) > 0)
        self.assertEqual(len(l0), len(l1))
        self.assertTrue(len(tip_names) > 0)

    def test_plot_covarying_rates_scatter_success(self):
        fake_pyplot = types.ModuleType("matplotlib.pyplot")

        class DummyLabel:
            def set_fontsize(self, *args, **kwargs):
                pass

        class DummyAxis:
            label = DummyLabel()

        class DummyAx:
            def __init__(self):
                self.scatter_called = False
                self.plot_called = False
                self.xaxis = DummyAxis()
                self.yaxis = DummyAxis()

            def scatter(self, *args, **kwargs):
                self.scatter_called = True

            def plot(self, *args, **kwargs):
                self.plot_called = True

            def legend(self, *args, **kwargs):
                pass

            def set_title(self, *args, **kwargs):
                pass

            def set_xlabel(self, *args, **kwargs):
                pass

            def set_ylabel(self, *args, **kwargs):
                pass

            def text(self, *args, **kwargs):
                pass

            @property
            def transAxes(self):
                return object()

        class DummyFig:
            def __init__(self):
                self.saved = False

            def tight_layout(self):
                pass

            def savefig(self, *args, **kwargs):
                self.saved = True

        dummy_ax = DummyAx()
        dummy_fig = DummyFig()
        fake_pyplot.subplots = lambda figsize: (dummy_fig, dummy_ax)
        fake_pyplot.close = lambda fig: None

        fake_matplotlib = types.ModuleType("matplotlib")
        fake_matplotlib.use = lambda backend: None
        fake_matplotlib.pyplot = fake_pyplot
        old_matplotlib = sys.modules.get("matplotlib")
        old_pyplot = sys.modules.get("matplotlib.pyplot")
        sys.modules["matplotlib"] = fake_matplotlib
        sys.modules["matplotlib.pyplot"] = fake_pyplot
        try:
            self.cov_rates._plot_covarying_rates_scatter([0.1, 0.2], [0.3, 0.4], [0.5, 0.01])
        finally:
            if old_matplotlib is not None:
                sys.modules["matplotlib"] = old_matplotlib
            else:
                sys.modules.pop("matplotlib", None)
            if old_pyplot is not None:
                sys.modules["matplotlib.pyplot"] = old_pyplot
            else:
                sys.modules.pop("matplotlib.pyplot", None)

        self.assertTrue(dummy_ax.scatter_called)
        self.assertTrue(dummy_ax.plot_called)
        self.assertTrue(dummy_fig.saved)


if __name__ == '__main__':
    unittest.main()
