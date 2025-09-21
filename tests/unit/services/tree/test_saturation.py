"""
Unit tests for Saturation class
"""

import unittest
from unittest.mock import Mock, MagicMock, patch, call
import numpy as np
from argparse import Namespace
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import multiprocessing as mp

from phykit.services.tree.saturation import Saturation, FileFormat


class TestSaturation(unittest.TestCase):
    """Test Saturation class"""

    def setUp(self):
        """Set up test fixtures"""
        self.args = Namespace(
            tree="test_tree.tre",
            alignment="test_alignment.fasta",
            exclude_gaps=False,
            verbose=False
        )
        self.saturation = Saturation(self.args)

    def test_init(self):
        """Test initialization"""
        self.assertEqual(self.saturation.tree_file_path, "test_tree.tre")
        self.assertEqual(self.saturation.alignment_file_path, "test_alignment.fasta")
        self.assertFalse(self.saturation.exclude_gaps)
        self.assertFalse(self.saturation.verbose)

    def test_process_args(self):
        """Test argument processing"""
        args = Namespace(
            tree="my_tree.tre",
            alignment="my_alignment.fasta",
            exclude_gaps=True,
            verbose=True
        )
        processed = self.saturation.process_args(args)

        self.assertEqual(processed["tree_file_path"], "my_tree.tre")
        self.assertEqual(processed["alignment_file_path"], "my_alignment.fasta")
        self.assertTrue(processed["exclude_gaps"])
        self.assertTrue(processed["verbose"])

    def test_file_format_enum(self):
        """Test FileFormat enum"""
        self.assertEqual(FileFormat.fasta.value, "fasta")
        self.assertEqual(FileFormat.clustal.value, "clustal")
        self.assertEqual(FileFormat.phylip.value, "phylip")
        self.assertEqual(FileFormat.stockholm.value, "stockholm")

    def test_process_combo_batch_without_gaps(self):
        """Test processing combination batch without gap exclusion"""
        # Mock tree
        mock_tree = Mock()
        mock_tree.distance.return_value = 0.5

        # Create sequence arrays
        seq_arrays = {
            'seq1': np.array(['A', 'T', 'C', 'G'], dtype='U1'),
            'seq2': np.array(['A', 'T', 'G', 'G'], dtype='U1')
        }

        # Empty gap mask since exclude_gaps=False
        gap_mask = {}

        # Process batch
        combo_batch = [('seq1', 'seq2')]
        results = self.saturation._process_combo_batch(
            mock_tree, seq_arrays, gap_mask, False, combo_batch
        )

        # Check results
        self.assertEqual(len(results), 1)
        pd, ud = results[0]
        self.assertEqual(pd, 0.5)  # Patristic distance
        self.assertEqual(ud, 0.25)  # Uncorrected distance (1 mismatch out of 4)

    def test_process_combo_batch_with_gaps(self):
        """Test processing combination batch with gap exclusion"""
        # Mock tree
        mock_tree = Mock()
        mock_tree.distance.return_value = 0.3

        # Create sequence arrays with gaps
        seq_arrays = {
            'seq1': np.array(['A', 'T', '-', 'G', 'C'], dtype='U1'),
            'seq2': np.array(['A', 'G', 'C', 'G', '-'], dtype='U1')
        }

        # Gap masks
        gap_mask = {
            'seq1': np.array([False, False, True, False, False]),
            'seq2': np.array([False, False, False, False, True])
        }

        # Process batch
        combo_batch = [('seq1', 'seq2')]
        results = self.saturation._process_combo_batch(
            mock_tree, seq_arrays, gap_mask, True, combo_batch
        )

        # Check results
        self.assertEqual(len(results), 1)
        pd, ud = results[0]
        self.assertEqual(pd, 0.3)  # Patristic distance
        # Uncorrected distance: positions 0,1,3 are valid (no gaps)
        # seq1: A T G, seq2: A G G -> 1 mismatch out of 3
        self.assertAlmostEqual(ud, 1/3, places=4)

    def test_process_combo_batch_all_gaps(self):
        """Test processing when all positions have gaps"""
        mock_tree = Mock()
        mock_tree.distance.return_value = 0.2

        seq_arrays = {
            'seq1': np.array(['-', '-', '-'], dtype='U1'),
            'seq2': np.array(['A', 'T', 'G'], dtype='U1')
        }

        gap_mask = {
            'seq1': np.array([True, True, True]),
            'seq2': np.array([False, False, False])
        }

        combo_batch = [('seq1', 'seq2')]
        results = self.saturation._process_combo_batch(
            mock_tree, seq_arrays, gap_mask, True, combo_batch
        )

        # Should return NaN for uncorrected distance
        self.assertEqual(len(results), 1)
        pd, ud = results[0]
        self.assertEqual(pd, 0.2)
        self.assertTrue(np.isnan(ud))

    def test_loop_through_combos_sequential(self):
        """Test sequential processing of combinations (small dataset)"""
        # Create mock tree
        mock_tree = Mock()
        mock_tree.distance.side_effect = [0.1, 0.2, 0.3]  # 3 combinations

        # Create mock alignment
        seq1 = SeqRecord(Seq("ATCG"), id="seq1", name="seq1")
        seq2 = SeqRecord(Seq("ATGG"), id="seq2", name="seq2")
        seq3 = SeqRecord(Seq("AACG"), id="seq3", name="seq3")
        alignment = Align.MultipleSeqAlignment([seq1, seq2, seq3])

        # Mock get_gap_chars
        self.saturation.get_gap_chars = Mock(return_value={'-', 'N', '?'})

        # Small number of combinations to trigger sequential processing
        combos = [('seq1', 'seq2'), ('seq1', 'seq3'), ('seq2', 'seq3')]

        pd_list, ud_list = self.saturation.loop_through_combos_and_calculate_pds_and_pis(
            combos, alignment, mock_tree, False
        )

        # Check results
        self.assertEqual(len(pd_list), 3)
        self.assertEqual(len(ud_list), 3)
        self.assertEqual(pd_list, [0.1, 0.2, 0.3])

        # Check uncorrected distances
        # seq1 vs seq2: 1 diff out of 4 = 0.25
        # seq1 vs seq3: 1 diff out of 4 = 0.25
        # seq2 vs seq3: 2 diff out of 4 = 0.5
        self.assertAlmostEqual(ud_list[0], 0.25, places=4)
        self.assertAlmostEqual(ud_list[1], 0.25, places=4)
        self.assertAlmostEqual(ud_list[2], 0.5, places=4)

    @patch('multiprocessing.Pool')
    def test_loop_through_combos_parallel(self, mock_pool_class):
        """Test parallel processing of combinations (large dataset)"""
        # Create mock tree
        mock_tree = Mock()
        mock_tree.distance.return_value = 0.15

        # Create mock alignment
        seqs = []
        for i in range(10):
            seq = SeqRecord(Seq("ATCG"), id=f"seq{i}", name=f"seq{i}")
            seqs.append(seq)
        alignment = Align.MultipleSeqAlignment(seqs)

        # Mock get_gap_chars
        self.saturation.get_gap_chars = Mock(return_value={'-', 'N', '?'})

        # Create many combinations to trigger parallel processing
        combos = [(f"seq{i}", f"seq{j}") for i in range(10) for j in range(i+1, 10)]
        self.assertEqual(len(combos), 45)  # Should be 45 combinations

        # Add more combos to exceed threshold
        combos.extend([('seq0', 'seq1') for _ in range(10)])  # Add duplicates to reach 55

        # Mock pool
        mock_pool = MagicMock()
        mock_pool_class.return_value.__enter__.return_value = mock_pool

        # Mock pool.map to return chunk results
        mock_pool.map.return_value = [
            [(0.15, 0.0)] * 10,  # First chunk results
            [(0.15, 0.0)] * 10,  # Second chunk results
            [(0.15, 0.0)] * 10,  # Third chunk results
            [(0.15, 0.0)] * 10,  # Fourth chunk results
            [(0.15, 0.0)] * 15,  # Fifth chunk results
        ]

        pd_list, ud_list = self.saturation.loop_through_combos_and_calculate_pds_and_pis(
            combos, alignment, mock_tree, False
        )

        # Check that multiprocessing was used
        mock_pool_class.assert_called_once()
        mock_pool.map.assert_called_once()

        # Check results
        self.assertEqual(len(pd_list), 55)
        self.assertEqual(len(ud_list), 55)
        self.assertTrue(all(pd == 0.15 for pd in pd_list))
        self.assertTrue(all(ud == 0.0 for ud in ud_list))

    @patch('builtins.print')
    def test_print_res_non_verbose(self, mock_print):
        """Test result printing in non-verbose mode"""
        combos = [('seq1', 'seq2'), ('seq1', 'seq3')]
        uncorrected_distances = [0.25, 0.33]
        patristic_distances = [0.1, 0.2]
        slope = 0.8

        self.saturation.print_res(
            False, combos, uncorrected_distances, patristic_distances, slope
        )

        # In non-verbose mode, should print slope and saturation
        mock_print.assert_called_once_with("0.8\t0.2")

    @patch('builtins.print')
    def test_print_res_verbose(self, mock_print):
        """Test result printing in verbose mode"""
        combos = [('seq1', 'seq2'), ('seq1', 'seq3')]
        uncorrected_distances = [0.25, 0.33]
        patristic_distances = [0.1, 0.2]
        slope = 0.8

        self.saturation.print_res(
            True, combos, uncorrected_distances, patristic_distances, slope
        )

        # In verbose mode, should print each combination
        expected_calls = [
            call("seq1\tseq2\t0.25\t0.1"),
            call("seq1\tseq3\t0.33\t0.2")
        ]
        mock_print.assert_has_calls(expected_calls)

    @patch('builtins.print')
    def test_print_res_broken_pipe(self, mock_print):
        """Test handling of BrokenPipeError"""
        mock_print.side_effect = BrokenPipeError()

        combos = [('seq1', 'seq2')]
        uncorrected_distances = [0.25]
        patristic_distances = [0.1]
        slope = 0.8

        # Should not raise exception
        try:
            self.saturation.print_res(
                False, combos, uncorrected_distances, patristic_distances, slope
            )
        except BrokenPipeError:
            self.fail("BrokenPipeError was not caught")

    @patch('phykit.services.tree.saturation.LinearRegression')
    @patch('phykit.services.tree.saturation.get_alignment_and_format_helper')
    def test_run_complete_workflow(self, mock_get_alignment, mock_lr_class):
        """Test complete run workflow"""
        # Mock alignment
        seq1 = SeqRecord(Seq("ATCG"), id="seq1", name="seq1")
        seq2 = SeqRecord(Seq("ATGG"), id="seq2", name="seq2")
        alignment = Align.MultipleSeqAlignment([seq1, seq2])
        mock_get_alignment.return_value = (alignment, None, False)

        # Mock tree
        mock_tree = Mock()
        mock_tree.distance.return_value = 0.15
        self.saturation.read_tree_file = Mock(return_value=mock_tree)
        self.saturation.get_tip_names_from_tree = Mock(return_value=['seq1', 'seq2'])

        # Mock the distance calculation method to return patristic and uncorrected distances
        self.saturation.loop_through_combos_and_calculate_pds_and_pis = Mock(
            return_value=([0.15], [0.25])  # patristic_distances, uncorrected_distances
        )

        # Mock linear regression
        mock_model = Mock()
        mock_model.coef_ = [0.85]
        mock_lr_class.return_value = mock_model

        # Mock print_res
        self.saturation.print_res = Mock()

        # Run
        self.saturation.run()

        # Check that all methods were called
        mock_get_alignment.assert_called_once_with(self.saturation.alignment_file_path)
        self.saturation.read_tree_file.assert_called_once()
        self.saturation.get_tip_names_from_tree.assert_called_once_with(mock_tree)
        mock_model.fit.assert_called_once()
        self.saturation.print_res.assert_called_once()

        # Check print_res arguments
        call_args = self.saturation.print_res.call_args[0]
        self.assertFalse(call_args[0])  # verbose
        self.assertEqual(len(call_args[1]), 1)  # combos
        self.assertEqual(call_args[4], 0.85)  # slope

    def test_loop_through_combos_with_gaps_sequential(self):
        """Test sequential processing with gap exclusion"""
        # Create mock tree
        mock_tree = Mock()
        mock_tree.distance.side_effect = [0.1, 0.2]

        # Create alignment with gaps
        seq1 = SeqRecord(Seq("AT-G"), id="seq1", name="seq1")
        seq2 = SeqRecord(Seq("ATCG"), id="seq2", name="seq2")
        seq3 = SeqRecord(Seq("A-CG"), id="seq3", name="seq3")
        alignment = Align.MultipleSeqAlignment([seq1, seq2, seq3])

        # Mock get_gap_chars
        self.saturation.get_gap_chars = Mock(return_value={'-', 'N', '?'})

        combos = [('seq1', 'seq2'), ('seq1', 'seq3')]

        pd_list, ud_list = self.saturation.loop_through_combos_and_calculate_pds_and_pis(
            combos, alignment, mock_tree, True  # exclude_gaps=True
        )

        # Check results
        self.assertEqual(len(pd_list), 2)
        self.assertEqual(len(ud_list), 2)
        self.assertEqual(pd_list, [0.1, 0.2])

        # seq1 vs seq2: compare only positions 0,1,3 (position 2 has gap in seq1)
        # AT-G vs ATCG -> ATG vs ATG = 0 differences
        self.assertAlmostEqual(ud_list[0], 0.0, places=4)

        # seq1 vs seq3: compare only positions 0,3 (positions 1,2 have gaps)
        # AT-G vs A-CG -> AG vs AG = 0 differences
        self.assertAlmostEqual(ud_list[1], 0.0, places=4)

    @patch('multiprocessing.cpu_count')
    @patch('multiprocessing.Pool')
    def test_loop_through_combos_parallel_cpu_limit(self, mock_pool_class, mock_cpu_count):
        """Test that parallel processing respects CPU count limit"""
        mock_cpu_count.return_value = 16  # Many CPUs

        # Create mock tree
        mock_tree = Mock()
        mock_tree.distance.return_value = 0.15

        # Create mock alignment
        seqs = []
        for i in range(15):
            seq = SeqRecord(Seq("ATCG"), id=f"seq{i}", name=f"seq{i}")
            seqs.append(seq)
        alignment = Align.MultipleSeqAlignment(seqs)

        # Mock get_gap_chars
        self.saturation.get_gap_chars = Mock(return_value={'-', 'N', '?'})

        # Create many combinations
        combos = [(f"seq{i}", f"seq{j}") for i in range(15) for j in range(i+1, 15)]
        self.assertTrue(len(combos) >= 50)  # Should trigger parallel processing

        # Mock pool
        mock_pool = MagicMock()
        mock_pool_class.return_value.__enter__.return_value = mock_pool
        mock_pool.map.return_value = [[(0.15, 0.0)] * 10 for _ in range(11)]

        pd_list, ud_list = self.saturation.loop_through_combos_and_calculate_pds_and_pis(
            combos, alignment, mock_tree, False
        )

        # Check that pool was created with max 8 workers
        mock_pool_class.assert_called_once_with(processes=8)


if __name__ == '__main__':
    unittest.main()