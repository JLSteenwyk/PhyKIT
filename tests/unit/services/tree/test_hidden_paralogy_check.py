"""
Unit tests for HiddenParalogyCheck class
"""

import unittest
from unittest.mock import Mock, MagicMock, patch, call, mock_open
import sys
import tempfile
import os
from argparse import Namespace
from Bio import Phylo
from io import StringIO
import multiprocessing as mp

from phykit.services.tree.hidden_paralogy_check import HiddenParalogyCheck


class TestHiddenParalogyCheck(unittest.TestCase):
    """Test HiddenParalogyCheck class"""

    def setUp(self):
        """Set up test fixtures"""
        self.args = Namespace(
            tree="test_tree.tre",
            clade="test_clades.txt"
        )
        self.checker = HiddenParalogyCheck(self.args)

    def test_init(self):
        """Test initialization"""
        self.assertEqual(self.checker.tree_file_path, "test_tree.tre")
        self.assertEqual(self.checker.clade, "test_clades.txt")

    def test_process_args(self):
        """Test argument processing"""
        args = Namespace(
            tree="my_tree.tre",
            clade="my_clades.txt"
        )
        processed = self.checker.process_args(args)

        self.assertEqual(processed["tree_file_path"], "my_tree.tre")
        self.assertEqual(processed["clade"], "my_clades.txt")

    def test_read_clades_file_success(self):
        """Test reading clades file successfully"""
        test_content = "taxa1 taxa2 taxa3\ntaxa4 taxa5\ntaxa6\n"

        with patch("builtins.open", mock_open(read_data=test_content)):
            clades = self.checker.read_clades_file("test_clades.txt")

        self.assertEqual(len(clades), 3)
        self.assertEqual(clades[0], ["taxa1", "taxa2", "taxa3"])
        self.assertEqual(clades[1], ["taxa4", "taxa5"])
        self.assertEqual(clades[2], ["taxa6"])

    @patch('sys.exit')
    @patch('builtins.print')
    def test_read_clades_file_not_found(self, mock_print, mock_exit):
        """Test handling of missing clades file"""
        with patch("builtins.open", side_effect=FileNotFoundError):
            self.checker.read_clades_file("nonexistent.txt")

        mock_print.assert_called_once_with("Clade file not found. Please check the path.")
        mock_exit.assert_called_once_with(2)

    @patch('builtins.print')
    def test_print_results(self, mock_print):
        """Test printing results"""
        res_arr = [
            ["monophyletic", []],
            ["not_monophyletic", ["taxa1", "taxa2"]],
            ["insufficient_taxon_representation"]
        ]

        self.checker.print_results(res_arr)

        mock_print.assert_has_calls([
            call("monophyletic"),
            call("not_monophyletic"),
            call("insufficient_taxon_representation")
        ])

    @patch('phykit.services.tree.hidden_paralogy_check.Phylo.read')
    def test_process_clade_batch_monophyletic(self, mock_phylo_read):
        """Test processing a batch of clades - monophyletic case"""
        # Create mock tree
        mock_tree = Mock()
        mock_subtree = Mock()

        # Setup mock terminals for subtree
        mock_terminal1 = Mock()
        mock_terminal1.name = "taxa1"
        mock_terminal2 = Mock()
        mock_terminal2.name = "taxa2"

        mock_subtree.get_terminals.return_value = [mock_terminal1, mock_terminal2]
        mock_tree.common_ancestor.return_value = mock_subtree
        mock_tree.root_with_outgroup.return_value = None

        mock_phylo_read.return_value = mock_tree

        # Test data
        clade_batch = [["taxa1", "taxa2"]]
        master_tree_tips = frozenset(["taxa1", "taxa2", "taxa3", "taxa4"])

        results = HiddenParalogyCheck._process_clade_batch(
            clade_batch, "test.tre", master_tree_tips
        )

        self.assertEqual(len(results), 1)
        self.assertEqual(results[0][0], "monophyletic")
        self.assertEqual(results[0][1], [])

    @patch('phykit.services.tree.hidden_paralogy_check.Phylo.read')
    def test_process_clade_batch_not_monophyletic(self, mock_phylo_read):
        """Test processing a batch of clades - not monophyletic case"""
        # Create mock tree
        mock_tree = Mock()
        mock_subtree = Mock()

        # Setup mock terminals - subtree has extra taxa making it not monophyletic
        mock_terminal1 = Mock()
        mock_terminal1.name = "taxa1"
        mock_terminal2 = Mock()
        mock_terminal2.name = "taxa2"
        mock_terminal3 = Mock()
        mock_terminal3.name = "taxa3"  # Extra taxon

        mock_subtree.get_terminals.return_value = [mock_terminal1, mock_terminal2, mock_terminal3]
        mock_tree.common_ancestor.return_value = mock_subtree
        mock_tree.root_with_outgroup.return_value = None

        mock_phylo_read.return_value = mock_tree

        # Test data - clade only asks for taxa1 and taxa2
        clade_batch = [["taxa1", "taxa2"]]
        master_tree_tips = frozenset(["taxa1", "taxa2", "taxa3", "taxa4"])

        results = HiddenParalogyCheck._process_clade_batch(
            clade_batch, "test.tre", master_tree_tips
        )

        self.assertEqual(len(results), 1)
        self.assertEqual(results[0][0], "not_monophyletic")
        # taxa3 is the difference
        self.assertIn("taxa3", results[0][1])

    @patch('phykit.services.tree.hidden_paralogy_check.Phylo.read')
    def test_process_clade_batch_insufficient_taxa(self, mock_phylo_read):
        """Test processing a batch with insufficient taxa"""
        # Mock tree (will be read but not used)
        mock_tree = Mock()
        mock_phylo_read.return_value = mock_tree

        # Test data - only one taxon in common with master tree
        clade_batch = [["taxa1", "taxa_not_in_tree"]]
        master_tree_tips = frozenset(["taxa1", "taxa2", "taxa3"])

        results = HiddenParalogyCheck._process_clade_batch(
            clade_batch, "test.tre", master_tree_tips
        )

        self.assertEqual(len(results), 1)
        self.assertEqual(results[0][0], "insufficient_taxon_representation")

    @patch('phykit.services.tree.hidden_paralogy_check.Phylo.read')
    def test_process_clade_batch_processing_error(self, mock_phylo_read):
        """Test handling of processing errors"""
        # Create mock tree that raises exception
        mock_tree = Mock()
        mock_tree.root_with_outgroup.side_effect = ValueError("Cannot root")

        mock_phylo_read.return_value = mock_tree

        clade_batch = [["taxa1", "taxa2"]]
        master_tree_tips = frozenset(["taxa1", "taxa2", "taxa3"])

        results = HiddenParalogyCheck._process_clade_batch(
            clade_batch, "test.tre", master_tree_tips
        )

        self.assertEqual(len(results), 1)
        self.assertEqual(results[0][0], "processing_error")

    @patch('phykit.services.tree.hidden_paralogy_check.Phylo.read')
    def test_process_clade_batch_multiple_clades(self, mock_phylo_read):
        """Test processing multiple clades in a batch"""
        # First tree - monophyletic
        mock_tree1 = Mock()
        mock_subtree1 = Mock()
        mock_terminal1_1 = Mock()
        mock_terminal1_1.name = "taxa1"
        mock_terminal1_2 = Mock()
        mock_terminal1_2.name = "taxa2"
        mock_subtree1.get_terminals.return_value = [mock_terminal1_1, mock_terminal1_2]
        mock_tree1.common_ancestor.return_value = mock_subtree1
        mock_tree1.root_with_outgroup.return_value = None

        # Second tree for third clade - not monophyletic
        mock_tree2 = Mock()
        mock_subtree2 = Mock()
        mock_terminal2_1 = Mock()
        mock_terminal2_1.name = "taxa3"
        mock_terminal2_2 = Mock()
        mock_terminal2_2.name = "taxa4"
        mock_terminal2_3 = Mock()
        mock_terminal2_3.name = "taxa5"
        mock_subtree2.get_terminals.return_value = [mock_terminal2_1, mock_terminal2_2, mock_terminal2_3]
        mock_tree2.common_ancestor.return_value = mock_subtree2
        mock_tree2.root_with_outgroup.return_value = None

        # Set up side_effect to return trees for each Phylo.read call
        # Actually, the code reads tree for each clade before checking
        mock_tree_insufficient = Mock()  # Tree for insufficient taxa clade
        mock_phylo_read.side_effect = [mock_tree1, mock_tree_insufficient, mock_tree2]

        clade_batch = [
            ["taxa1", "taxa2"],  # Will be monophyletic
            ["taxa_not_exist"],  # Insufficient taxa (won't read tree)
            ["taxa3", "taxa4"]   # Will be not monophyletic (has taxa5)
        ]
        master_tree_tips = frozenset(["taxa1", "taxa2", "taxa3", "taxa4", "taxa5"])

        results = HiddenParalogyCheck._process_clade_batch(
            clade_batch, "test.tre", master_tree_tips
        )

        self.assertEqual(len(results), 3)
        self.assertEqual(results[0][0], "monophyletic")
        self.assertEqual(results[1][0], "insufficient_taxon_representation")
        self.assertEqual(results[2][0], "not_monophyletic")

    @patch('builtins.print')
    @patch('phykit.services.tree.hidden_paralogy_check.Phylo.read')
    def test_run_sequential_processing(self, mock_phylo_read, mock_print):
        """Test run method with sequential processing (small dataset)"""
        # Setup mock tree
        mock_master_tree = Mock()
        mock_tree = Mock()
        mock_subtree = Mock()

        # Master tree tips
        self.checker.read_tree_file = Mock(return_value=mock_master_tree)
        self.checker.get_tip_names_from_tree = Mock(
            side_effect=[
                ["taxa1", "taxa2", "taxa3"],  # Master tree tips
                ["taxa1", "taxa2"]            # Subtree tips
            ]
        )

        # Mock clades file reading
        self.checker.read_clades_file = Mock(return_value=[
            ["taxa1", "taxa2"]  # Single clade for sequential processing
        ])

        # Setup subtree
        mock_terminal1 = Mock()
        mock_terminal1.name = "taxa1"
        mock_terminal2 = Mock()
        mock_terminal2.name = "taxa2"
        mock_subtree.get_terminals.return_value = [mock_terminal1, mock_terminal2]

        mock_tree.common_ancestor.return_value = mock_subtree
        mock_tree.root_with_outgroup.return_value = None

        mock_phylo_read.return_value = mock_tree

        # Run the method
        self.checker.run()

        # Check that sequential path was taken (Phylo.read called once)
        mock_phylo_read.assert_called_once_with("test_tree.tre", "newick")

        # Check output
        mock_print.assert_called_once_with("monophyletic")

    @patch('builtins.print')
    @patch('multiprocessing.Pool')
    def test_run_parallel_processing(self, mock_pool_class, mock_print):
        """Test run method with parallel processing (large dataset)"""
        # Setup mock tree
        mock_master_tree = Mock()

        # Master tree tips
        self.checker.read_tree_file = Mock(return_value=mock_master_tree)
        self.checker.get_tip_names_from_tree = Mock(
            return_value=["taxa1", "taxa2", "taxa3", "taxa4"]
        )

        # Create many clades to trigger parallel processing (>=10)
        clades = [["taxa1", "taxa2"] for _ in range(15)]
        self.checker.read_clades_file = Mock(return_value=clades)

        # Mock multiprocessing pool
        mock_pool = MagicMock()
        mock_pool_class.return_value.__enter__.return_value = mock_pool

        # Mock pool.map results
        batch_results = [
            [["monophyletic", []] for _ in range(8)],  # First batch
            [["monophyletic", []] for _ in range(7)]   # Second batch
        ]
        mock_pool.map.return_value = batch_results

        # Run the method
        self.checker.run()

        # Check that parallel path was taken
        mock_pool_class.assert_called_once()
        mock_pool.map.assert_called_once()

        # Check that print was called 15 times
        self.assertEqual(mock_print.call_count, 15)

    @patch('builtins.print')
    @patch('phykit.services.tree.hidden_paralogy_check.Phylo.read')
    def test_run_with_insufficient_taxa(self, mock_phylo_read, mock_print):
        """Test run with insufficient taxa representation"""
        mock_master_tree = Mock()
        mock_tree = Mock()
        mock_phylo_read.return_value = mock_tree

        self.checker.read_tree_file = Mock(return_value=mock_master_tree)
        self.checker.get_tip_names_from_tree = Mock(
            return_value=["taxa1", "taxa2"]
        )

        # Clade with only one taxon in common
        self.checker.read_clades_file = Mock(return_value=[
            ["taxa1", "taxa_not_in_tree"]
        ])

        # Run the method
        self.checker.run()

        # Check output
        mock_print.assert_called_once_with("insufficient_taxon_representation")

    @patch('builtins.print')
    @patch('phykit.services.tree.hidden_paralogy_check.Phylo.read')
    def test_run_with_not_monophyletic(self, mock_phylo_read, mock_print):
        """Test run with not monophyletic result"""
        mock_master_tree = Mock()
        mock_tree = Mock()
        mock_subtree = Mock()

        self.checker.read_tree_file = Mock(return_value=mock_master_tree)
        self.checker.get_tip_names_from_tree = Mock(
            side_effect=[
                ["taxa1", "taxa2", "taxa3", "taxa4"],  # Master tree
                ["taxa1", "taxa2", "taxa3"]            # Subtree (has extra taxa3)
            ]
        )

        self.checker.read_clades_file = Mock(return_value=[
            ["taxa1", "taxa2"]  # Only asking for taxa1 and taxa2
        ])

        mock_tree.common_ancestor.return_value = mock_subtree
        mock_tree.root_with_outgroup.return_value = None

        mock_phylo_read.return_value = mock_tree

        # Run the method
        self.checker.run()

        # Check output
        mock_print.assert_called_once_with("not_monophyletic")

    @patch('multiprocessing.cpu_count')
    @patch('multiprocessing.Pool')
    def test_run_parallel_with_cpu_count(self, mock_pool_class, mock_cpu_count):
        """Test parallel processing respects CPU count limit"""
        mock_cpu_count.return_value = 16  # Many CPUs

        mock_master_tree = Mock()
        self.checker.read_tree_file = Mock(return_value=mock_master_tree)
        self.checker.get_tip_names_from_tree = Mock(return_value=["taxa1"])

        # Create many clades
        clades = [["taxa1"] for _ in range(20)]
        self.checker.read_clades_file = Mock(return_value=clades)

        mock_pool = MagicMock()
        mock_pool_class.return_value.__enter__.return_value = mock_pool
        mock_pool.map.return_value = [[["monophyletic", []]] for _ in range(20)]

        self.checker.run()

        # Should be capped at 8 workers even with 16 CPUs
        mock_pool_class.assert_called_once_with(processes=8)

    @patch('phykit.services.tree.hidden_paralogy_check.Phylo.read')
    def test_edge_cases(self, mock_phylo_read):
        """Test various edge cases"""
        # Mock tree (will be read but not used for insufficient taxa)
        mock_tree = Mock()
        mock_phylo_read.return_value = mock_tree

        # Test empty clade intersection
        clade_batch = [["taxa_not_in_tree"]]
        master_tree_tips = frozenset(["taxa1", "taxa2"])

        results = HiddenParalogyCheck._process_clade_batch(
            clade_batch, "test.tre", master_tree_tips
        )

        self.assertEqual(results[0][0], "insufficient_taxon_representation")

        # Test single taxon clade
        clade_batch = [["taxa1"]]
        master_tree_tips = frozenset(["taxa1", "taxa2"])

        results = HiddenParalogyCheck._process_clade_batch(
            clade_batch, "test.tre", master_tree_tips
        )

        self.assertEqual(results[0][0], "insufficient_taxon_representation")


if __name__ == '__main__':
    unittest.main()