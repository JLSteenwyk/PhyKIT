"""
Unit tests for PolytomyTest class
"""

import unittest
from unittest.mock import Mock, MagicMock, patch, call, mock_open
import sys
import tempfile
import os
from argparse import Namespace
from Bio import Phylo
from Bio.Phylo import Newick
from scipy.stats import _stats_py
import multiprocessing as mp

from phykit.services.tree.polytomy_test import PolytomyTest


class TestPolytomyTest(unittest.TestCase):
    """Test PolytomyTest class"""

    def setUp(self):
        """Set up test fixtures"""
        self.args = Namespace(
            trees="test_trees.txt",
            groups="test_groups.txt"
        )
        self.polytomy = PolytomyTest(self.args)

    def test_init(self):
        """Test initialization"""
        self.assertEqual(self.polytomy.trees, "test_trees.txt")
        self.assertEqual(self.polytomy.groups, "test_groups.txt")

    def test_process_args(self):
        """Test argument processing"""
        args = Namespace(
            trees="my_trees.txt",
            groups="my_groups.txt"
        )
        processed = self.polytomy.process_args(args)

        self.assertEqual(processed["trees"], "my_trees.txt")
        self.assertEqual(processed["groups"], "my_groups.txt")

    def test_read_in_groups_success(self):
        """Test reading groups file successfully"""
        test_content = """# Comment line
test1\tseq1;seq2\tseq3;seq4\tseq5;seq6\toutgroup1;outgroup2
test2\tseq7;seq8\tseq9;seq10\tseq11;seq12\toutgroup3;outgroup4
"""

        with patch("builtins.open", mock_open(read_data=test_content)):
            groups = self.polytomy.read_in_groups()

        self.assertEqual(len(groups), 2)
        self.assertEqual(groups[0][0], "test1")
        self.assertEqual(groups[0][1], ["seq1", "seq2"])
        self.assertEqual(groups[0][2], ["seq3", "seq4"])
        self.assertEqual(groups[0][3], ["seq5", "seq6"])
        self.assertEqual(groups[0][4], ["outgroup1", "outgroup2"])

    @patch('sys.exit')
    @patch('builtins.print')
    def test_read_in_groups_index_error(self, mock_print, mock_exit):
        """Test handling of incorrectly formatted groups file"""
        test_content = """test1\tincomplete_line"""

        with patch("builtins.open", mock_open(read_data=test_content)):
            self.polytomy.read_in_groups()

        mock_print.assert_any_call("test_groups.txt contains an indexing error.")
        mock_exit.assert_called_once_with(2)

    @patch('sys.exit')
    @patch('builtins.print')
    def test_read_in_groups_file_not_found(self, mock_print, mock_exit):
        """Test handling of missing groups file"""
        with patch("builtins.open", side_effect=FileNotFoundError):
            self.polytomy.read_in_groups()

        mock_print.assert_any_call("test_groups.txt corresponds to no such file.")
        mock_exit.assert_called_once_with(2)

    def test_determine_groups_of_groups(self):
        """Test determining groups of groups"""
        groups_arr = [
            ["test1", ["seq1", "seq2"], ["seq3", "seq4"], ["seq5", "seq6"], ["out1", "out2"]]
        ]

        groups_of_groups, outgroup_taxa = self.polytomy.determine_groups_of_groups(groups_arr)

        self.assertIn("test1", groups_of_groups)
        self.assertEqual(len(groups_of_groups["test1"]), 3)
        self.assertEqual(groups_of_groups["test1"][0], ["seq1", "seq2"])
        self.assertEqual(groups_of_groups["test1"][1], ["seq3", "seq4"])
        self.assertEqual(groups_of_groups["test1"][2], ["seq5", "seq6"])
        self.assertEqual(outgroup_taxa, ["out1", "out2"])

    def test_count_number_of_groups_in_triplet(self):
        """Test counting groups represented in a triplet"""
        triplet = ("seq1", "seq3", "seq5")
        groups = [["seq1", "seq2"], ["seq3", "seq4"], ["seq5", "seq6"]]

        count = self.polytomy.count_number_of_groups_in_triplet(triplet, groups)

        self.assertEqual(count, 3)

    def test_count_number_of_groups_in_triplet_partial(self):
        """Test counting groups when not all are represented"""
        triplet = ("seq1", "seq2", "seq5")
        groups = [["seq1", "seq2"], ["seq3", "seq4"], ["seq5", "seq6"]]

        count = self.polytomy.count_number_of_groups_in_triplet(triplet, groups)

        self.assertEqual(count, 2)  # Only group 0 and group 2

    def test_set_branch_lengths_in_tree_to_one(self):
        """Test setting all branch lengths to 1"""
        mock_tree = Mock(spec=Newick.Tree)
        mock_clade1 = Mock()
        mock_clade2 = Mock()
        mock_clade3 = Mock()
        mock_tree.find_clades.return_value = [mock_clade1, mock_clade2, mock_clade3]

        self.polytomy.set_branch_lengths_in_tree_to_one(mock_tree)

        self.assertEqual(mock_clade1.branch_length, 1)
        self.assertEqual(mock_clade2.branch_length, 1)
        self.assertEqual(mock_clade3.branch_length, 1)

    def test_check_if_triplet_is_a_polytomy(self):
        """Test checking if triplet is a polytomy"""
        mock_tree = Mock(spec=Newick.Tree)

        # Case 1: Polytomy (only 1 nonterminal)
        mock_tree.get_nonterminals.return_value = [Mock()]
        is_polytomy = self.polytomy.check_if_triplet_is_a_polytomy(mock_tree)
        self.assertTrue(is_polytomy)

        # Case 2: Not a polytomy (2 nonterminals)
        mock_tree.get_nonterminals.return_value = [Mock(), Mock()]
        is_polytomy = self.polytomy.check_if_triplet_is_a_polytomy(mock_tree)
        self.assertFalse(is_polytomy)

    def test_sister_relationship_counter(self):
        """Test counting sister relationships"""
        summary = {}
        tree_file = "tree1.tre"
        sisters = "0-1"

        # First occurrence
        summary = self.polytomy.sister_relationship_counter(tree_file, summary, sisters)
        self.assertEqual(summary["tree1.tre"]["0-1"], 1)

        # Second occurrence
        summary = self.polytomy.sister_relationship_counter(tree_file, summary, sisters)
        self.assertEqual(summary["tree1.tre"]["0-1"], 2)

        # Different sister relationship
        summary = self.polytomy.sister_relationship_counter(tree_file, summary, "0-2")
        self.assertEqual(summary["tree1.tre"]["0-2"], 1)

    @patch('phykit.services.tree.polytomy_test.Phylo.read')
    def test_get_triplet_tree_success(self, mock_phylo_read):
        """Test getting triplet tree successfully"""
        mock_tree = Mock()
        mock_tree.root_with_outgroup.return_value = None
        mock_phylo_read.return_value = mock_tree

        self.polytomy.prune_tree_using_taxa_list = Mock(return_value=mock_tree)

        tips = ["seq1", "seq2", "seq3", "seq4", "out1"]
        triplet = ("seq1", "seq2", "seq3")
        tree_file = "test.tre"
        outgroup_taxa = ["out1"]

        result = self.polytomy.get_triplet_tree(tips, triplet, tree_file, outgroup_taxa)

        self.assertEqual(result, mock_tree)
        mock_tree.root_with_outgroup.assert_called_once_with(["out1"])
        self.polytomy.prune_tree_using_taxa_list.assert_called_once()

    @patch('phykit.services.tree.polytomy_test.Phylo.read')
    def test_get_triplet_tree_root_error(self, mock_phylo_read):
        """Test handling of rooting error in get_triplet_tree"""
        mock_tree = Mock()
        mock_tree.root_with_outgroup.side_effect = ValueError("Cannot root")
        mock_phylo_read.return_value = mock_tree

        tips = ["seq1", "seq2", "seq3", "seq4"]
        triplet = ("seq1", "seq2", "seq3")
        tree_file = "test.tre"
        outgroup_taxa = ["out1"]

        result = self.polytomy.get_triplet_tree(tips, triplet, tree_file, outgroup_taxa)

        self.assertFalse(result)

    def test_determine_sisters_from_triplet(self):
        """Test determining which taxa are sisters"""
        groups = [["seq1", "seq2"], ["seq3", "seq4"], ["seq5", "seq6"]]

        # Test pair from groups 0 and 1
        pair = ("seq1", "seq3")
        sisters = self.polytomy.determine_sisters_from_triplet(groups, pair)
        self.assertEqual(sisters, "0-1")

        # Test pair from groups 0 and 2
        pair = ("seq2", "seq5")
        sisters = self.polytomy.determine_sisters_from_triplet(groups, pair)
        self.assertEqual(sisters, "0-2")

        # Test pair from groups 1 and 2
        pair = ("seq4", "seq6")
        sisters = self.polytomy.determine_sisters_from_triplet(groups, pair)
        self.assertEqual(sisters, "1-2")

    def test_determine_sisters_and_add_to_counter(self):
        """Test determining sisters and updating counter"""
        tip_names = ["seq1", "seq3", "seq5"]
        mock_tree = Mock()
        mock_tree.distance.return_value = 2  # Sisters
        tree_file = "test.tre"
        groups = [["seq1", "seq2"], ["seq3", "seq4"], ["seq5", "seq6"]]
        summary = {}

        self.polytomy.check_if_triplet_is_a_polytomy = Mock(return_value=False)

        summary = self.polytomy.determine_sisters_and_add_to_counter(
            tip_names, mock_tree, tree_file, groups, summary
        )

        # Should have found sister relationships
        self.assertIn("test.tre", summary)
        self.assertTrue(len(summary["test.tre"]) > 0)

    def test_determine_sisters_polytomy_case(self):
        """Test that polytomies are not counted"""
        tip_names = ["seq1", "seq3", "seq5"]
        mock_tree = Mock()
        mock_tree.distance.return_value = 2
        tree_file = "test.tre"
        groups = [["seq1", "seq2"], ["seq3", "seq4"], ["seq5", "seq6"]]
        summary = {}

        self.polytomy.check_if_triplet_is_a_polytomy = Mock(return_value=True)

        summary = self.polytomy.determine_sisters_and_add_to_counter(
            tip_names, mock_tree, tree_file, groups, summary
        )

        # Should not add anything for polytomy
        self.assertEqual(summary, {})

    def test_get_triplet_and_gene_support_freq_counts(self):
        """Test counting triplet and gene support frequencies"""
        summary = {
            "tree1": {"0-1": 5, "0-2": 3, "1-2": 2},
            "tree2": {"0-1": 2, "0-2": 8, "1-2": 1},
            "tree3": {"0-1": 1, "1-2": 4}  # Missing 0-2
        }

        triplet_counts, gene_support = self.polytomy.get_triplet_and_gene_support_freq_counts(summary)

        # Check triplet counts
        self.assertEqual(triplet_counts["g0g1_count"], 8)  # 5 + 2 + 1
        self.assertEqual(triplet_counts["g0g2_count"], 11)  # 3 + 8 + 0
        self.assertEqual(triplet_counts["g1g2_count"], 7)  # 2 + 1 + 4

        # Check gene support frequencies
        self.assertEqual(gene_support["0-1"], 1)  # tree1 max
        self.assertEqual(gene_support["0-2"], 1)  # tree2 max
        self.assertEqual(gene_support["1-2"], 1)  # tree3 max

    @patch('phykit.services.tree.polytomy_test.chisquare')
    def test_chisquare_tests(self, mock_chisquare):
        """Test chi-square statistical tests"""
        triplet_counts = {"g0g1_count": 10, "g0g2_count": 12, "g1g2_count": 8}
        gene_support = {"0-1": 5, "0-2": 6, "1-2": 4}

        mock_result1 = Mock()
        mock_result2 = Mock()
        mock_chisquare.side_effect = [mock_result1, mock_result2]

        triplet_res, gene_res = self.polytomy.chisquare_tests(triplet_counts, gene_support)

        self.assertEqual(triplet_res, mock_result1)
        self.assertEqual(gene_res, mock_result2)

        # Check the values passed to chisquare
        calls = mock_chisquare.call_args_list
        self.assertEqual(calls[0][0][0], [10, 12, 8])
        self.assertEqual(calls[1][0][0], [5, 6, 4])

    @patch('builtins.print')
    def test_print_gene_support_freq_res(self, mock_print):
        """Test printing gene support frequency results"""
        mock_result = Mock()
        mock_result.statistic = 3.4567
        mock_result.pvalue = 0.012345

        gene_support_freq = {"0-1": 5, "0-2": 6, "1-2": 4}
        trees_file_path = ["tree1.tre", "tree2.tre"]

        self.polytomy.print_gene_support_freq_res(mock_result, gene_support_freq, trees_file_path)

        # Check that results were printed
        mock_print.assert_any_call("Gene Support Frequency Results")
        mock_print.assert_any_call("==============================")
        mock_print.assert_any_call("chi-squared: 3.4567")
        mock_print.assert_any_call("p-value: 0.012345")
        mock_print.assert_any_call("total genes: 15")
        mock_print.assert_any_call("0-1: 5")
        mock_print.assert_any_call("0-2: 6")
        mock_print.assert_any_call("1-2: 4")

    @patch('builtins.print')
    def test_print_gene_support_freq_res_broken_pipe(self, mock_print):
        """Test handling of BrokenPipeError"""
        mock_print.side_effect = BrokenPipeError()

        mock_result = Mock()
        mock_result.statistic = 1.0
        mock_result.pvalue = 0.5

        gene_support_freq = {"0-1": 1, "0-2": 1, "1-2": 1}
        trees_file_path = ["tree1.tre"]

        # Should not raise exception
        try:
            self.polytomy.print_gene_support_freq_res(mock_result, gene_support_freq, trees_file_path)
        except BrokenPipeError:
            self.fail("BrokenPipeError was not caught")

    @patch('phykit.services.tree.polytomy_test.Phylo.read')
    def test_process_tree_batch(self, mock_phylo_read):
        """Test processing a batch of trees"""
        mock_tree = Mock()
        mock_phylo_read.return_value = mock_tree

        self.polytomy.get_tip_names_from_tree = Mock(return_value=["seq1", "seq2", "seq3"])
        self.polytomy.examine_all_triplets_and_sister_pairing = Mock(
            return_value={"tree1.tre": {"0-1": 5}}
        )

        tree_batch = ["tree1.tre", "tree2.tre"]
        groups_of_groups = {"test": [["seq1"], ["seq2"], ["seq3"]]}
        outgroup_taxa = ["out1"]

        result = self.polytomy._process_tree_batch(tree_batch, groups_of_groups, outgroup_taxa)

        self.assertIn("tree1.tre", result)
        self.assertEqual(result["tree1.tre"]["0-1"], 5)

    def test_process_tree_batch_exception(self):
        """Test that exceptions in tree processing are caught"""
        tree_batch = ["bad_tree.tre"]
        groups_of_groups = {"test": [["seq1"], ["seq2"], ["seq3"]]}
        outgroup_taxa = ["out1"]

        # Should not raise exception
        with patch('phykit.services.tree.polytomy_test.Phylo.read', side_effect=Exception("Bad tree")):
            result = self.polytomy._process_tree_batch(tree_batch, groups_of_groups, outgroup_taxa)
            self.assertEqual(result, {})

    @patch('phykit.services.tree.polytomy_test.read_single_column_file_to_list')
    @patch('phykit.services.tree.polytomy_test.Phylo.read')
    def test_loop_through_trees_sequential(self, mock_phylo_read, mock_read_file):
        """Test sequential processing of trees (small dataset)"""
        mock_read_file.return_value = ["tree1.tre", "tree2.tre"]  # Small dataset
        mock_tree = Mock()
        mock_phylo_read.return_value = mock_tree

        self.polytomy.get_tip_names_from_tree = Mock(return_value=["seq1", "seq2", "seq3"])

        # The method modifies and returns the summary dict, so we need to simulate that
        def examine_side_effect(tips, tree_file, summary, groups, outgroup):
            if tree_file == "tree1.tre":
                summary["tree1.tre"] = {"0-1": 3}
            elif tree_file == "tree2.tre":
                summary["tree2.tre"] = {"0-1": 2}
            return summary

        self.polytomy.examine_all_triplets_and_sister_pairing = Mock(side_effect=examine_side_effect)

        groups_of_groups = {"test": [["seq1"], ["seq2"], ["seq3"]]}
        outgroup_taxa = ["out1"]

        summary = self.polytomy.loop_through_trees_and_examine_sister_support_among_triplets(
            ["tree1.tre", "tree2.tre"], groups_of_groups, outgroup_taxa
        )

        # Should process sequentially
        self.assertEqual(len(summary), 2)
        self.assertEqual(self.polytomy.examine_all_triplets_and_sister_pairing.call_count, 2)

    @patch('multiprocessing.Pool')
    @patch('phykit.services.tree.polytomy_test.read_single_column_file_to_list')
    def test_loop_through_trees_parallel(self, mock_read_file, mock_pool_class):
        """Test parallel processing of trees (large dataset)"""
        # Create many tree files to trigger parallel processing
        trees = [f"tree{i}.tre" for i in range(20)]
        mock_read_file.return_value = trees

        mock_pool = MagicMock()
        mock_pool_class.return_value.__enter__.return_value = mock_pool

        # Mock pool.map results
        batch_results = [
            {"tree1.tre": {"0-1": 3}},
            {"tree2.tre": {"0-2": 2}}
        ]
        mock_pool.map.return_value = batch_results

        groups_of_groups = {"test": [["seq1"], ["seq2"], ["seq3"]]}
        outgroup_taxa = ["out1"]

        summary = self.polytomy.loop_through_trees_and_examine_sister_support_among_triplets(
            trees, groups_of_groups, outgroup_taxa
        )

        # Should use multiprocessing
        mock_pool_class.assert_called_once()
        mock_pool.map.assert_called_once()
        self.assertEqual(len(summary), 2)

    def test_examine_all_triplets_sequential(self):
        """Test sequential processing of triplets (small dataset)"""
        tips = ["seq1", "seq2", "seq3", "seq4", "out1"]
        tree_file = "test.tre"
        summary = {}
        groups_of_groups = {"test": [["seq1"], ["seq2"], ["seq3"]]}
        outgroup_taxa = ["out1"]

        self.polytomy.get_triplet_tree = Mock(return_value=None)

        summary = self.polytomy.examine_all_triplets_and_sister_pairing(
            tips, tree_file, summary, groups_of_groups, outgroup_taxa
        )

        # Should process sequentially for small dataset
        self.polytomy.get_triplet_tree.assert_called()

    def test_process_triplet_batch(self):
        """Test processing a batch of triplets"""
        triplet_batch = [("seq1", "seq2", "seq3")]
        tips = ["seq1", "seq2", "seq3", "seq4"]
        tree_file = "test.tre"
        groups_of_groups = {"test": [["seq1"], ["seq2"], ["seq3"]]}
        outgroup_taxa = ["out1"]

        mock_tree = Mock()
        mock_tree.get_terminals.return_value = [Mock(), Mock(), Mock()]  # 3 terminals

        self.polytomy._get_triplet_tree_cached = Mock(return_value=mock_tree)
        self.polytomy.count_number_of_groups_in_triplet = Mock(return_value=3)
        self.polytomy.get_tip_names_from_tree = Mock(return_value=["seq1", "seq2", "seq3"])
        self.polytomy.set_branch_lengths_in_tree_to_one = Mock()
        self.polytomy.determine_sisters_and_add_to_counter = Mock(
            return_value={"test.tre": {"0-1": 1}}
        )

        result = self.polytomy._process_triplet_batch(
            triplet_batch, tips, tree_file, groups_of_groups, outgroup_taxa
        )

        self.assertIn("test.tre", result)

    @patch('builtins.print')
    @patch('phykit.services.tree.polytomy_test.read_single_column_file_to_list')
    def test_run_complete_workflow(self, mock_read_file, mock_print):
        """Test complete run workflow"""
        # Mock file reading
        mock_read_file.return_value = ["tree1.tre", "tree2.tre"]

        # Mock methods
        self.polytomy.read_in_groups = Mock(return_value=[
            ["test1", ["seq1"], ["seq2"], ["seq3"], ["out1"]]
        ])
        self.polytomy.determine_groups_of_groups = Mock(return_value=(
            {"test1": [["seq1"], ["seq2"], ["seq3"]]},
            ["out1"]
        ))
        self.polytomy.loop_through_trees_and_examine_sister_support_among_triplets = Mock(return_value={
            "tree1.tre": {"0-1": 5, "0-2": 3},
            "tree2.tre": {"0-2": 4, "1-2": 2}
        })
        self.polytomy.get_triplet_and_gene_support_freq_counts = Mock(return_value=(
            {"g0g1_count": 5, "g0g2_count": 7, "g1g2_count": 2},
            {"0-1": 1, "0-2": 2, "1-2": 0}
        ))
        self.polytomy.chisquare_tests = Mock(return_value=(
            Mock(statistic=2.5, pvalue=0.1),
            Mock(statistic=1.8, pvalue=0.3)
        ))

        # Run
        self.polytomy.run()

        # Verify all methods were called
        self.polytomy.read_in_groups.assert_called_once()
        self.polytomy.determine_groups_of_groups.assert_called_once()
        self.polytomy.loop_through_trees_and_examine_sister_support_among_triplets.assert_called_once()
        self.polytomy.get_triplet_and_gene_support_freq_counts.assert_called_once()
        self.polytomy.chisquare_tests.assert_called_once()

        # Check that results were printed
        mock_print.assert_any_call("Gene Support Frequency Results")

    def test_cached_methods(self):
        """Test that cached methods work correctly"""
        # Test _count_groups_cached
        triplet = ("seq1", "seq2", "seq3")
        groups = (frozenset(["seq1"]), frozenset(["seq2"]), frozenset(["seq3"]))

        count = self.polytomy._count_groups_cached(triplet, groups)
        self.assertEqual(count, 3)

        # Call again - should use cache
        count2 = self.polytomy._count_groups_cached(triplet, groups)
        self.assertEqual(count2, 3)

    @patch('sys.exit')
    @patch('builtins.print')
    @patch('phykit.services.tree.polytomy_test.Phylo.read')
    def test_examine_all_triplets_file_not_found(self, mock_phylo_read, mock_print, mock_exit):
        """Test handling of FileNotFoundError in examine_all_triplets"""
        mock_phylo_read.side_effect = FileNotFoundError()

        tips = ["seq1", "seq2", "seq3"]
        tree_file = "missing.tre"
        summary = {}
        groups_of_groups = {"test": [["seq1"], ["seq2"], ["seq3"]]}
        outgroup_taxa = ["out1"]

        # For small dataset, it should try to read the tree and catch the error
        with self.assertRaises(FileNotFoundError):
            self.polytomy.examine_all_triplets_and_sister_pairing(
                tips, tree_file, summary, groups_of_groups, outgroup_taxa
            )


if __name__ == '__main__':
    unittest.main()