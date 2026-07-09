"""
Unit tests for Saturation class
"""

import unittest
import builtins
import subprocess
import sys
import tempfile
from pathlib import Path
from unittest.mock import Mock, MagicMock, patch
import numpy as np
from argparse import Namespace
from Bio import Align
from Bio.Phylo.Newick import Clade, Tree
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import phykit.services.tree.saturation as saturation_module
from phykit.services.tree.saturation import Saturation


class TestSaturation(unittest.TestCase):
    """Test Saturation class"""

    def test_module_import_does_not_import_numpy_or_biopython_tree_modules(self):
        code = """
import sys
import phykit.services.tree.saturation as module
assert callable(module.print_json)
assert callable(module.get_alignment_and_format_helper)
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.files" not in sys.modules
assert "numpy" not in sys.modules
assert "multiprocessing" not in sys.modules
assert "Bio.Align" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "phykit.helpers.plot_config" not in sys.modules
"""
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_lazy_numpy_proxy_caches_resolved_attributes(self):
        lazy_np = saturation_module._LazyNumpy()

        first = lazy_np.count_nonzero
        second = lazy_np.count_nonzero

        self.assertIs(first, second)
        self.assertIs(lazy_np.__dict__["count_nonzero"], first)

    def test_lazy_multiprocessing_caches_module_and_keeps_cpu_count_patchable(self):
        import multiprocessing

        lazy_mp = saturation_module._LazyMultiprocessing()

        with patch.object(multiprocessing, "cpu_count", return_value=13):
            self.assertEqual(lazy_mp.cpu_count(), 13)

        with patch.object(multiprocessing, "cpu_count", return_value=5):
            self.assertEqual(lazy_mp.cpu_count(), 5)

        self.assertIs(lazy_mp._module, multiprocessing)

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

    def test_standard_combo_tip_derivation_does_not_scan_all_pairs(self):
        tips = [f"seq{i}" for i in range(6)]
        pairs = [(tips[i], tips[j]) for i in range(6) for j in range(i + 1, 6)]

        class StandardPairs:
            def __len__(self):
                return len(pairs)

            def __getitem__(self, index):
                if index == 0:
                    return pairs[0]
                if isinstance(index, slice):
                    return pairs[index]
                raise AssertionError("standard pair derivation should not scan pairs")

            def __iter__(self):
                raise AssertionError("standard pair derivation should not iterate")

        observed = Saturation._combo_tips_from_pairs(
            StandardPairs(),
            standard_combo_order=True,
        )

        self.assertEqual(observed, tips)

    def test_custom_combo_tip_derivation_preserves_first_seen_order(self):
        combos = [("seq2", "seq3"), ("seq1", "seq2"), ("seq4", "seq1")]

        observed = Saturation._combo_tips_from_pairs(combos)

        self.assertEqual(observed, ["seq2", "seq3", "seq1", "seq4"])

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

    def test_should_use_multiprocessing_env_controls(self):
        with patch.dict("os.environ", {"PHYKIT_DISABLE_MP": "1"}, clear=False):
            self.assertFalse(self.saturation._should_use_multiprocessing(9999))
        with patch.dict("os.environ", {"PHYKIT_FORCE_MP": "1"}, clear=False):
            self.assertTrue(self.saturation._should_use_multiprocessing(1))

    def test_plot_max_uses_array_method_for_small_plot_vectors(self):
        values = np.array([0.4, 0.1, 0.9])

        with patch(
            "phykit.services.tree.saturation.np.max",
            side_effect=AssertionError("small saturation plot maxima should use ndarray.max"),
        ):
            self.assertAlmostEqual(saturation_module._plot_max(values), 0.9)

    def test_plot_max_preserves_large_vector_numpy_path(self):
        values = np.ones(saturation_module._PLOT_DIRECT_EXTREMA_LIMIT + 1)
        calls = []
        original_max = saturation_module.np.max

        def max_spy(observed, *args, **kwargs):
            calls.append(observed.shape)
            return original_max(observed, *args, **kwargs)

        with patch("phykit.services.tree.saturation.np.max", side_effect=max_spy):
            self.assertAlmostEqual(saturation_module._plot_max(values), 1.0)

        self.assertEqual(
            calls,
            [(saturation_module._PLOT_DIRECT_EXTREMA_LIMIT + 1,)],
        )

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

    def test_sequence_to_array_uppercases_ascii_bytes(self):
        seq_arr = self.saturation._sequence_to_array("atgx-?")

        self.assertEqual(seq_arr.dtype, np.uint8)
        self.assertEqual(seq_arr.tolist(), [65, 84, 71, 88, 45, 63])

    def test_sequence_to_array_falls_back_for_non_ascii(self):
        seq_arr = self.saturation._sequence_to_array("aé")

        self.assertEqual(seq_arr.dtype.kind, "U")
        self.assertEqual(seq_arr.tolist(), ["A", "É"])

    def test_gap_values_match_byte_arrays(self):
        seq_arr = self.saturation._sequence_to_array("atgn")
        gap_values = self.saturation._gap_values_for_array(seq_arr, {"n", "?"})

        self.assertEqual(set(gap_values.tolist()), {78, 63})

    def test_gap_mask_uses_byte_lookup(self):
        seq_arr = self.saturation._sequence_to_array("AT-GN?")
        with patch("phykit.services.tree.saturation.np.isin") as mock_isin:
            gap_mask = self.saturation._gap_mask_for_array(
                seq_arr, {"-", "N", "?"}
            )

        mock_isin.assert_not_called()
        self.assertEqual(
            gap_mask.tolist(),
            [False, False, True, False, True, True],
        )

    def test_gap_mask_falls_back_for_non_ascii_arrays(self):
        seq_arr = self.saturation._sequence_to_array("AÉ-")

        gap_mask = self.saturation._gap_mask_for_array(seq_arr, {"-", "É"})

        self.assertEqual(gap_mask.tolist(), [False, True, True])

    def test_process_combo_batch_counts_byte_array_valid_positions(self):
        mock_tree = Mock()
        mock_tree.distance.return_value = 0.4
        seq_arrays = {
            "seq1": self.saturation._sequence_to_array("AT-GN"),
            "seq2": self.saturation._sequence_to_array("ATCGT"),
        }
        gap_values = self.saturation._gap_values_for_array(
            seq_arrays["seq1"], {"-", "N", "?"}
        )
        gap_mask = {
            name: np.isin(seq_arr, gap_values)
            for name, seq_arr in seq_arrays.items()
        }

        results = self.saturation._process_combo_batch(
            mock_tree, seq_arrays, gap_mask, True, [("seq1", "seq2")]
        )

        self.assertEqual(results[0][0], 0.4)
        self.assertEqual(results[0][1], 0.0)

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
    def test_loop_through_combos_fast_tree_distances(self, mock_pool_class):
        """Real trees should use cached patristic distances and preserve combo order."""
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="seq1"),
                    Clade(
                        branch_length=2.0,
                        clades=[
                            Clade(branch_length=3.0, name="seq2"),
                            Clade(branch_length=4.0, name="seq3"),
                        ],
                    ),
                ],
            )
        )
        alignment = Align.MultipleSeqAlignment([
            SeqRecord(Seq("ATCG"), id="seq1", name="seq1"),
            SeqRecord(Seq("ATGG"), id="seq2", name="seq2"),
            SeqRecord(Seq("AACG"), id="seq3", name="seq3"),
        ])
        self.saturation.MP_MIN_COMBOS = 1
        combos = [('seq2', 'seq3'), ('seq1', 'seq2'), ('seq1', 'seq3')]

        pd_list, ud_list = self.saturation.loop_through_combos_and_calculate_pds_and_pis(
            combos, alignment, tree, False
        )

        mock_pool_class.assert_not_called()
        self.assertEqual(pd_list, [tree.distance(*combo) for combo in combos])
        self.assertAlmostEqual(ud_list[0], 0.5, places=4)
        self.assertAlmostEqual(ud_list[1], 0.25, places=4)
        self.assertAlmostEqual(ud_list[2], 0.25, places=4)

    @patch('multiprocessing.Pool')
    def test_loop_through_combos_uses_ordered_fast_tree_distances(
        self, mock_pool_class
    ):
        """Ordered cached distances should be consumed directly in combo order."""
        mock_tree = Mock()
        mock_tree.distance.side_effect = AssertionError(
            "cached distances should be used"
        )
        alignment = Align.MultipleSeqAlignment([
            SeqRecord(Seq("ATCG"), id="seq1", name="seq1"),
            SeqRecord(Seq("ATGG"), id="seq2", name="seq2"),
            SeqRecord(Seq("AACG"), id="seq3", name="seq3"),
        ])
        combos = [('seq1', 'seq2'), ('seq1', 'seq3'), ('seq2', 'seq3')]
        self.saturation.calculate_pairwise_tip_distances_fast = Mock(
            return_value=(combos, [0.1, 0.2, 0.3])
        )
        self.saturation.MP_MIN_COMBOS = 1

        pd_list, ud_list = self.saturation.loop_through_combos_and_calculate_pds_and_pis(
            combos, alignment, mock_tree, False
        )

        mock_pool_class.assert_not_called()
        self.assertEqual(pd_list, [0.1, 0.2, 0.3])
        self.assertEqual(ud_list, [0.25, 0.25, 0.5])

    @patch('multiprocessing.Pool')
    def test_loop_through_standard_combos_skips_cached_combo_pairs(
        self, mock_pool_class
    ):
        """Standard-order callers should not request duplicate cached pairs."""
        mock_tree = Mock()
        mock_tree.distance.side_effect = AssertionError(
            "cached distances should be used"
        )
        alignment = Align.MultipleSeqAlignment([
            SeqRecord(Seq("ATCG"), id="seq1", name="seq1"),
            SeqRecord(Seq("ATGG"), id="seq2", name="seq2"),
            SeqRecord(Seq("AACG"), id="seq3", name="seq3"),
        ])
        combos = [('seq1', 'seq2'), ('seq1', 'seq3'), ('seq2', 'seq3')]
        self.saturation.calculate_pairwise_tip_distances_fast = Mock(
            return_value=(None, [0.1, 0.2, 0.3])
        )
        self.saturation.MP_MIN_COMBOS = 1

        pd_list, ud_list = self.saturation.loop_through_combos_and_calculate_pds_and_pis(
            combos,
            alignment,
            mock_tree,
            False,
            standard_combo_order=True,
        )

        mock_pool_class.assert_not_called()
        self.saturation.calculate_pairwise_tip_distances_fast.assert_called_once_with(
            mock_tree,
            ["seq1", "seq2", "seq3"],
            include_combos=False,
        )
        self.assertEqual(pd_list, [0.1, 0.2, 0.3])
        self.assertEqual(ud_list, [0.25, 0.25, 0.5])

    def test_calculate_uncorrected_distances_matrix_with_gaps(self):
        seq_arrays = {
            "seq1": self.saturation._sequence_to_array("AT-G"),
            "seq2": self.saturation._sequence_to_array("ATCG"),
            "seq3": self.saturation._sequence_to_array("G--G"),
            "seq4": self.saturation._sequence_to_array("----"),
            "seq5": self.saturation._sequence_to_array("NNNN"),
        }
        gap_values = self.saturation._gap_values_for_array(
            seq_arrays["seq1"], {"-", "N", "?"}
        )
        gap_mask = {
            name: np.isin(seq_arr, gap_values)
            for name, seq_arr in seq_arrays.items()
        }
        combos = [
            ("seq1", "seq2"),
            ("seq1", "seq3"),
            ("seq4", "seq5"),
        ]

        distances = self.saturation._calculate_uncorrected_distances_matrix(
            ["seq1", "seq2", "seq3", "seq4", "seq5"],
            combos,
            seq_arrays,
            gap_mask,
            True,
        )

        self.assertEqual(distances[0], 0.0)
        self.assertEqual(distances[1], 0.5)
        self.assertTrue(np.isnan(distances[2]))

    def test_standard_upper_triangle_values_preserve_pair_order(self):
        combo_tips = ["seq1", "seq2", "seq3", "seq4"]
        combos = [
            ("seq1", "seq2"),
            ("seq1", "seq3"),
            ("seq1", "seq4"),
            ("seq2", "seq3"),
            ("seq2", "seq4"),
            ("seq3", "seq4"),
        ]
        values = np.arange(16, dtype=float).reshape(4, 4)

        distances = self.saturation._standard_upper_triangle_values(
            combo_tips,
            combos,
            values,
        )

        self.assertEqual(distances, [1.0, 2.0, 3.0, 6.0, 7.0, 11.0])

    def test_standard_upper_triangle_values_reject_custom_pair_order(self):
        combo_tips = ["seq1", "seq2", "seq3"]
        combos = [("seq1", "seq3"), ("seq1", "seq2"), ("seq2", "seq3")]
        values = np.arange(9, dtype=float).reshape(3, 3)

        distances = self.saturation._standard_upper_triangle_values(
            combo_tips,
            combos,
            values,
        )

        self.assertIsNone(distances)

    def test_standard_upper_triangle_validator_avoids_tip_slices(self):
        class NoSliceList(list):
            def __getitem__(self, item):
                if isinstance(item, slice):
                    raise AssertionError("validator should stream expected pairs")
                return super().__getitem__(item)

        combo_tips = NoSliceList(["seq1", "seq2", "seq3", "seq4"])
        combos = [
            ("seq1", "seq2"),
            ("seq1", "seq3"),
            ("seq1", "seq4"),
            ("seq2", "seq3"),
            ("seq2", "seq4"),
            ("seq3", "seq4"),
        ]

        self.assertTrue(
            self.saturation._combos_are_standard_upper_triangle(
                combo_tips,
                combos,
            )
        )

    def test_standard_upper_triangle_values_trusted_order_skips_validation(self):
        combo_tips = ["seq1", "seq2", "seq3"]
        combos = [("seq1", "seq2"), ("seq1", "seq3"), ("seq2", "seq3")]
        values = np.arange(9, dtype=float).reshape(3, 3)

        with patch.object(
            self.saturation,
            "_combos_are_standard_upper_triangle",
            side_effect=AssertionError("trusted run path should not revalidate pairs"),
        ):
            distances = self.saturation._standard_upper_triangle_values(
                combo_tips,
                combos,
                values,
                standard_combo_order=True,
            )

        self.assertEqual(distances, [1.0, 2.0, 5.0])

    def test_standard_upper_triangle_gappy_distances_preserve_nan_rows(self):
        combo_tips = ["seq1", "seq2", "seq3"]
        combos = [("seq1", "seq2"), ("seq1", "seq3"), ("seq2", "seq3")]
        identity_counts = np.array(
            [
                [4.0, 2.0, 0.0],
                [2.0, 4.0, 1.0],
                [0.0, 1.0, 4.0],
            ]
        )
        adjusted_lengths = np.array(
            [
                [4.0, 4.0, 0.0],
                [4.0, 4.0, 2.0],
                [0.0, 2.0, 4.0],
            ]
        )

        distances = self.saturation._standard_upper_triangle_gappy_distances(
            combo_tips,
            combos,
            identity_counts,
            adjusted_lengths,
        )

        self.assertEqual(distances[0], 0.5)
        self.assertTrue(np.isnan(distances[1]))
        self.assertEqual(distances[2], 0.5)

    def test_standard_upper_triangle_gappy_distances_all_valid_skips_scratch_rows(self):
        combo_tips = ["seq1", "seq2", "seq3"]
        combos = [("seq1", "seq2"), ("seq1", "seq3"), ("seq2", "seq3")]
        identity_counts = np.array(
            [
                [4.0, 2.0, 1.0],
                [2.0, 4.0, 2.0],
                [1.0, 2.0, 4.0],
            ]
        )
        adjusted_lengths = np.full((3, 3), 4.0)

        with patch(
            "phykit.services.tree.saturation.np.empty",
            side_effect=AssertionError(
                "all-valid gappy distances should use direct row division"
            ),
        ):
            distances = self.saturation._standard_upper_triangle_gappy_distances(
                combo_tips,
                combos,
                identity_counts,
                adjusted_lengths,
            )

        self.assertEqual(distances, [0.5, 0.75, 0.5])

    def test_calculate_uncorrected_distances_matrix_no_gap_path_skips_valid_matrix(self):
        seq_arrays = {
            "seq1": self.saturation._sequence_to_array("ATCG"),
            "seq2": self.saturation._sequence_to_array("ATGG"),
            "seq3": self.saturation._sequence_to_array("AACG"),
        }
        combos = [
            ("seq1", "seq2"),
            ("seq1", "seq3"),
            ("seq2", "seq3"),
        ]
        with patch(
            "phykit.services.tree.saturation.np.ones",
            side_effect=AssertionError(
                "no-gap saturation matrix path should not allocate a valid matrix"
            ),
        ):
            distances = self.saturation._calculate_uncorrected_distances_matrix(
                ["seq1", "seq2", "seq3"],
                combos,
                seq_arrays,
                {},
                False,
            )

        self.assertEqual(distances, [0.25, 0.25, 0.5])

    def test_standard_no_gap_matrix_path_skips_square_distance_matrix(self):
        combo_tips = ["seq1", "seq2", "seq3"]
        combos = [
            ("seq1", "seq2"),
            ("seq1", "seq3"),
            ("seq2", "seq3"),
        ]
        seq_matrix = np.vstack([
            self.saturation._sequence_to_array("ATCG"),
            self.saturation._sequence_to_array("ATGG"),
            self.saturation._sequence_to_array("AACG"),
        ])

        with patch(
            "phykit.services.tree.saturation.np.empty",
            side_effect=AssertionError(
                "standard no-gap distances should stream the upper triangle"
            ),
        ):
            distances = self.saturation._calculate_uncorrected_distances_no_gap_matrix(
                combo_tips,
                combos,
                seq_matrix,
            )

        self.assertEqual(distances, [0.25, 0.25, 0.5])

    def test_calculate_uncorrected_distances_exclude_gaps_clean_alignment_uses_no_gap_path(self):
        seq_arrays = {
            "seq1": self.saturation._sequence_to_array("ATCG"),
            "seq2": self.saturation._sequence_to_array("ATGG"),
            "seq3": self.saturation._sequence_to_array("AACG"),
        }
        gap_mask = {
            name: np.zeros(len(seq_arr), dtype=bool)
            for name, seq_arr in seq_arrays.items()
        }
        combos = [
            ("seq1", "seq2"),
            ("seq1", "seq3"),
            ("seq2", "seq3"),
        ]
        with patch(
            "phykit.services.tree.saturation.np.vstack",
            wraps=np.vstack,
        ) as vstack_spy:
            distances = self.saturation._calculate_uncorrected_distances_matrix(
                ["seq1", "seq2", "seq3"],
                combos,
                seq_arrays,
                gap_mask,
                True,
            )

        self.assertEqual(distances, [0.25, 0.25, 0.5])
        self.assertEqual(vstack_spy.call_count, 1)

    def test_calculate_uncorrected_distances_stacks_raw_gap_masks_before_invert(self):
        seq_arrays = {
            "seq1": self.saturation._sequence_to_array("AT-G"),
            "seq2": self.saturation._sequence_to_array("ATCG"),
            "seq3": self.saturation._sequence_to_array("A-CG"),
        }
        gap_mask = {
            name: self.saturation._gap_mask_for_array(seq_arr, {"-", "N", "?"})
            for name, seq_arr in seq_arrays.items()
        }
        combos = [
            ("seq1", "seq2"),
            ("seq1", "seq3"),
            ("seq2", "seq3"),
        ]
        stacked_inputs = []

        def capture_vstack(rows, *args, **kwargs):
            stacked_inputs.append([row.copy() for row in rows])
            return np.vstack(rows, *args, **kwargs)

        with patch(
            "phykit.services.tree.saturation.np.vstack",
            side_effect=capture_vstack,
        ):
            distances = self.saturation._calculate_uncorrected_distances_matrix(
                ["seq1", "seq2", "seq3"],
                combos,
                seq_arrays,
                gap_mask,
                True,
            )

        self.assertEqual(len(stacked_inputs), 2)
        for observed, tip in zip(stacked_inputs[1], ["seq1", "seq2", "seq3"]):
            np.testing.assert_array_equal(observed, gap_mask[tip])
        self.assertEqual(distances, [0.0, 0.0, 0.0])

    @patch('multiprocessing.Pool')
    def test_loop_through_combos_uses_matrix_uncorrected_distances(
        self, mock_pool_class
    ):
        mock_tree = Mock()
        mock_tree.distance.side_effect = AssertionError(
            "cached distances should be used"
        )
        alignment = Align.MultipleSeqAlignment([
            SeqRecord(Seq("ATCG"), id="seq1", name="seq1"),
            SeqRecord(Seq("ATGG"), id="seq2", name="seq2"),
            SeqRecord(Seq("AACG"), id="seq3", name="seq3"),
        ])
        combos = [('seq1', 'seq2'), ('seq1', 'seq3'), ('seq2', 'seq3')]
        self.saturation.calculate_pairwise_tip_distances_fast = Mock(
            return_value=(combos, [0.1, 0.2, 0.3])
        )
        self.saturation._calculate_uncorrected_distances_matrix = Mock(
            return_value=[0.11, 0.22, 0.33]
        )
        self.saturation.MP_MIN_COMBOS = 1

        pd_list, ud_list = self.saturation.loop_through_combos_and_calculate_pds_and_pis(
            combos, alignment, mock_tree, False
        )

        mock_pool_class.assert_not_called()
        self.saturation._calculate_uncorrected_distances_matrix.assert_called_once()
        self.assertEqual(pd_list, [0.1, 0.2, 0.3])
        self.assertEqual(ud_list, [0.11, 0.22, 0.33])

    @patch('multiprocessing.Pool')
    def test_loop_identical_sequences_skip_uncorrected_distance_arrays(
        self, mock_pool_class
    ):
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="seq1"),
                    Clade(
                        branch_length=2.0,
                        clades=[
                            Clade(branch_length=3.0, name="seq2"),
                            Clade(branch_length=4.0, name="seq3"),
                        ],
                    ),
                ],
            )
        )
        alignment = Align.MultipleSeqAlignment([
            SeqRecord(Seq("ACGTN-?*X"), id="seq1", name="seq1"),
            SeqRecord(Seq("acgtn-?*x"), id="seq2", name="seq2"),
            SeqRecord(Seq("ACGTN-?*X"), id="seq3", name="seq3"),
        ])
        combos = [('seq1', 'seq2'), ('seq1', 'seq3'), ('seq2', 'seq3')]
        self.saturation.MP_MIN_COMBOS = 1
        with patch.object(
            self.saturation,
            "_sequence_to_array",
            side_effect=AssertionError(
                "identical sequences should skip uncorrected-distance arrays"
            ),
        ):
            pd_list, ud_list = (
                self.saturation.loop_through_combos_and_calculate_pds_and_pis(
                    combos,
                    alignment,
                    tree,
                    True,
                    is_protein=False,
                )
            )

        mock_pool_class.assert_not_called()
        self.assertEqual(pd_list, [tree.distance(*combo) for combo in combos])
        self.assertEqual(ud_list, [0.0, 0.0, 0.0])

    def test_identical_sequence_helper_does_not_slice_rows(self):
        class NoSliceList(list):
            def __getitem__(self, key):
                if isinstance(key, slice):
                    raise AssertionError("identical-sequence scan should not slice")
                return super().__getitem__(key)

        sequences = NoSliceList(["ACGT", "ACGT", "ACGT"])

        self.assertTrue(saturation_module._all_sequences_identical(sequences))

    def test_valid_site_count_for_identical_unicode_sequence_counts_gap_chars(self):
        sequence = "ACGT\u03a9N-X?"

        self.assertEqual(
            Saturation._valid_site_count_for_identical_sequence(
                sequence,
                {"-", "?", "*", "X", "N"},
            ),
            5,
        )
        self.assertEqual(
            Saturation._valid_site_count_for_identical_sequence(
                sequence,
                {"-", "?", "*", "X"},
            ),
            6,
        )

    def test_loop_identical_all_gap_sequences_return_nan_distances(self):
        tree = Tree(
            root=Clade(
                clades=[
                    Clade(branch_length=1.0, name="seq1"),
                    Clade(branch_length=1.0, name="seq2"),
                ],
            )
        )
        alignment = Align.MultipleSeqAlignment([
            SeqRecord(Seq("NN--??"), id="seq1", name="seq1"),
            SeqRecord(Seq("nn--??"), id="seq2", name="seq2"),
        ])
        combos = [('seq1', 'seq2')]
        with patch.object(
            self.saturation,
            "_sequence_to_array",
            side_effect=AssertionError(
                "identical all-gap sequences should skip arrays"
            ),
        ):
            pd_list, ud_list = (
                self.saturation.loop_through_combos_and_calculate_pds_and_pis(
                    combos,
                    alignment,
                    tree,
                    True,
                    is_protein=False,
                )
            )

        self.assertEqual(pd_list, [2.0])
        self.assertTrue(np.isnan(ud_list[0]))

    @patch('multiprocessing.Pool')
    def test_loop_through_combos_parallel(self, mock_pool_class):
        """Test parallel processing of combinations (large dataset)"""
        self.saturation.MP_MIN_COMBOS = 50
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

        # In verbose mode, pair rows should be batched into one stdout write.
        mock_print.assert_called_once_with(
            "seq1\tseq2\t0.25\t0.1\nseq1\tseq3\t0.33\t0.2"
        )

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

    @patch('phykit.services.tree.saturation.get_alignment_and_format_helper')
    def test_run_complete_workflow(self, mock_get_alignment):
        """Test complete run workflow"""
        # Mock alignment
        seq1 = SeqRecord(Seq("ATCG"), id="seq1", name="seq1")
        seq2 = SeqRecord(Seq("ATGG"), id="seq2", name="seq2")
        alignment = Align.MultipleSeqAlignment([seq1, seq2])
        mock_get_alignment.return_value = (alignment, None, False)

        # Mock tree
        mock_tree = Mock()
        mock_tree.distance.return_value = 0.15
        self.saturation.read_tree_file_unmodified = Mock(return_value=mock_tree)
        self.saturation.get_tip_names_from_tree = Mock(return_value=['seq1', 'seq2'])

        # Mock the distance calculation method to return patristic and uncorrected distances
        self.saturation.loop_through_combos_and_calculate_pds_and_pis = Mock(
            return_value=([0.15], [0.25])  # patristic_distances, uncorrected_distances
        )

        # Mock print_res
        self.saturation.print_res = Mock()

        # Run
        self.saturation.run()

        # Check that all methods were called
        mock_get_alignment.assert_called_once_with(self.saturation.alignment_file_path)
        self.saturation.read_tree_file_unmodified.assert_called_once()
        self.saturation.get_tip_names_from_tree.assert_called_once_with(mock_tree)
        self.saturation.print_res.assert_called_once()
        args = self.saturation.print_res.call_args.args
        self.assertEqual(args[0], self.saturation.verbose)
        self.assertEqual(args[1], [('seq1', 'seq2')])
        self.assertEqual(args[2], [0.25])
        self.assertEqual(args[3], [0.15])
        self.assertAlmostEqual(args[4], 0.25 / 0.15, places=6)

        # Check print_res arguments
        call_args = self.saturation.print_res.call_args[0]
        self.assertFalse(call_args[0])  # verbose
        self.assertEqual(len(call_args[1]), 1)  # combos

    @patch("phykit.services.tree.saturation.print_json")
    def test_print_res_json_non_verbose(self, mock_json):
        self.saturation.json_output = True
        self.saturation.plot = False
        self.saturation.print_res(
            False,
            [("a", "b")],
            [0.25],
            [0.5],
            0.8,
        )
        payload = mock_json.call_args.args[0]
        self.assertFalse(payload["verbose"])
        self.assertEqual(payload["summary"]["slope"], 0.8)

    @patch("phykit.services.tree.saturation.print_json")
    def test_print_res_json_verbose_with_plot(self, mock_json):
        self.saturation.json_output = True
        self.saturation.plot = True
        self.saturation.plot_output = "sat.png"
        self.saturation.print_res(
            True,
            [("a", "b")],
            [0.25],
            [0.5],
            0.8,
        )
        payload = mock_json.call_args.args[0]
        self.assertTrue(payload["verbose"])
        self.assertEqual(payload["rows"], payload["pairs"])
        self.assertEqual(payload["plot_output"], "sat.png")

    def test_plot_saturation_scatter_importerror(self):
        original_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name.startswith("matplotlib"):
                raise ImportError("no matplotlib")
            return original_import(name, globals, locals, fromlist, level)

        with patch("builtins.__import__", side_effect=fake_import):
            with self.assertRaises(SystemExit) as exc:
                self.saturation._plot_saturation_scatter(np.array([0.1]), np.array([0.2]), 0.5)
        self.assertEqual(exc.exception.code, 2)

    def test_plot_saturation_scatter_skips_redundant_tight_layout(self):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot  # noqa: F401
            from matplotlib.figure import Figure
        except ImportError:
            self.skipTest("matplotlib is not installed")

        with tempfile.TemporaryDirectory() as tmp_dir:
            output_path = Path(tmp_dir) / "saturation.png"
            self.saturation.plot_output = str(output_path)

            with patch.object(
                Figure,
                "tight_layout",
                side_effect=AssertionError("savefig bbox should handle tight cropping"),
            ):
                self.saturation._plot_saturation_scatter(
                    np.array([0.1, 0.2, 0.3]),
                    np.array([0.04, 0.08, 0.12]),
                    0.4,
                )

            self.assertTrue(output_path.exists())

    @patch('phykit.services.tree.saturation.get_alignment_and_format_helper')
    def test_run_with_plot_calls_plotter(self, mock_get_alignment):
        self.saturation.plot = True
        self.saturation._plot_saturation_scatter = Mock()

        seq1 = SeqRecord(Seq("ATCG"), id="seq1", name="seq1")
        seq2 = SeqRecord(Seq("ATGG"), id="seq2", name="seq2")
        alignment = Align.MultipleSeqAlignment([seq1, seq2])
        mock_get_alignment.return_value = (alignment, None, False)

        mock_tree = Mock()
        self.saturation.read_tree_file_unmodified = Mock(return_value=mock_tree)
        self.saturation.get_tip_names_from_tree = Mock(return_value=['seq1', 'seq2'])
        self.saturation.loop_through_combos_and_calculate_pds_and_pis = Mock(
            return_value=([0.15], [0.25])
        )
        self.saturation.print_res = Mock()

        self.saturation.run()
        self.saturation._plot_saturation_scatter.assert_called_once()

    @patch('phykit.services.tree.saturation.get_alignment_and_format_helper')
    def test_run_uses_unmodified_tree_read(self, mock_get_alignment):
        seq1 = SeqRecord(Seq("ATCG"), id="seq1", name="seq1")
        seq2 = SeqRecord(Seq("ATGG"), id="seq2", name="seq2")
        alignment = Align.MultipleSeqAlignment([seq1, seq2])
        mock_get_alignment.return_value = (alignment, None, False)
        mock_tree = Mock()
        self.saturation.read_tree_file_unmodified = Mock(return_value=mock_tree)
        self.saturation.read_tree_file = Mock(
            side_effect=AssertionError("copying tree reader should not be used")
        )
        self.saturation.get_tip_names_from_tree = Mock(return_value=['seq1', 'seq2'])
        self.saturation.loop_through_combos_and_calculate_pds_and_pis = Mock(
            return_value=([0.15], [0.25])
        )
        self.saturation.print_res = Mock()

        self.saturation.run()

        self.saturation.read_tree_file_unmodified.assert_called_once_with()
        self.saturation.loop_through_combos_and_calculate_pds_and_pis.assert_called_once()

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
        self.saturation.MP_MIN_COMBOS = 50
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
