"""
Unit tests for HiddenParalogyCheck class
"""

import multiprocessing
import unittest
from unittest.mock import Mock, MagicMock, patch, mock_open
from argparse import Namespace
from io import StringIO
import tempfile
import os
import subprocess
import sys

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
import phykit.services.tree.hidden_paralogy_check as module
from phykit.services.tree.hidden_paralogy_check import HiddenParalogyCheck


def test_module_import_does_not_import_biophylo_or_numpy():
    code = """
import sys
import phykit.services.tree.hidden_paralogy_check as module

assert hasattr(module.Phylo, "read")
assert callable(module.print_json)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
assert "multiprocessing" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_phylo_caches_resolved_reader(tmp_path):
    tree_path = tmp_path / "tree.tre"
    tree_path.write_text("((A:1,B:1):1,C:2);\n")
    lazy_phylo = module._LazyPhylo()

    first = lazy_phylo.read(str(tree_path), "newick")
    cached_module = lazy_phylo._module
    second = lazy_phylo.read(str(tree_path), "newick")

    assert cached_module is not None
    assert lazy_phylo._module is cached_module
    assert lazy_phylo.__dict__["read"] is cached_module.read
    assert [tip.name for tip in first.get_terminals()] == ["A", "B", "C"]
    assert [tip.name for tip in second.get_terminals()] == ["A", "B", "C"]


def test_lazy_multiprocessing_caches_module_and_keeps_cpu_count_patchable():
    lazy_mp = module._LazyMultiprocessing()

    with patch.object(multiprocessing, "cpu_count", return_value=11):
        assert lazy_mp.cpu_count() == 11
        assert lazy_mp._module is multiprocessing

    with patch.object(multiprocessing, "cpu_count", return_value=13) as cpu_count:
        assert lazy_mp.cpu_count() == 13

    cpu_count.assert_called_once_with()


def test_run_uses_unmodified_master_tree_read(mocker):
    args = Namespace(tree="test_tree.tre", clade="test_clades.txt", json=False)
    checker = HiddenParalogyCheck(args)
    tree = Phylo.read(StringIO("((A,B),C);"), "newick")

    read_unmodified = mocker.patch.object(
        checker, "read_tree_file_unmodified", return_value=tree
    )
    mocker.patch.object(
        checker,
        "read_tree_file",
        side_effect=AssertionError("run should use read_tree_file_unmodified"),
    )
    mocker.patch.object(checker, "get_tip_names_from_tree", return_value=["A", "B", "C"])
    build_index = mocker.patch.object(
        checker, "_build_exact_clade_index", return_value={frozenset(["A", "B"])}
    )
    mocker.patch.object(checker, "read_clades_file", return_value=[["A", "B"]])
    phylo_read = mocker.patch(
        "phykit.services.tree.hidden_paralogy_check.Phylo.read",
        side_effect=AssertionError("exact clades should not reread the tree"),
    )
    print_results = mocker.patch.object(checker, "print_results")

    checker.run()

    read_unmodified.assert_called_once_with()
    build_index.assert_called_once_with(tree)
    phylo_read.assert_not_called()
    print_results.assert_called_once_with([["monophyletic", []]])


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

    def test_read_clades_file_preserves_blank_line_entries(self):
        """Blank lines remain empty clades for compatibility."""
        test_content = "taxa1 taxa2\n\n   \ntaxa3\n"

        with patch("builtins.open", mock_open(read_data=test_content)):
            clades = self.checker.read_clades_file("test_clades.txt")

        self.assertEqual(clades, [["taxa1", "taxa2"], [], [], ["taxa3"]])

    def test_read_clades_file_iterates_without_full_file_materialization(self):
        """Large clade files should not require a full read/splitlines list."""

        class IterableOnlyHandle:
            def __enter__(self):
                return self

            def __exit__(self, *_args):
                return False

            def read(self):
                raise AssertionError("clade parser should iterate over file rows")

            def __iter__(self):
                return iter(["taxa1 taxa2\n", "\n", "taxa3 taxa4\n"])

        with patch("builtins.open", return_value=IterableOnlyHandle()):
            clades = self.checker.read_clades_file("test_clades.txt")

        self.assertEqual(clades, [["taxa1", "taxa2"], [], ["taxa3", "taxa4"]])

    @patch('sys.exit')
    @patch('builtins.print')
    def test_read_clades_file_not_found(self, mock_print, mock_exit):
        """Test handling of missing clades file"""
        with patch("builtins.open", side_effect=FileNotFoundError):
            self.checker.read_clades_file("nonexistent.txt")

        mock_print.assert_called_once_with("Clade file not found. Please check the path.")
        mock_exit.assert_called_once_with(2)

    def test_print_results(self):
        """Test printing results"""
        res_arr = [
            ["monophyletic", []],
            ["not_monophyletic", ["taxa1", "taxa2"]],
            ["insufficient_taxon_representation"]
        ]

        with patch(
            "phykit.services.tree.hidden_paralogy_check.sys.stdout.write"
        ) as write, patch("builtins.print") as mocked_print:
            self.checker.print_results(res_arr)

        mocked_print.assert_not_called()
        write.assert_called_once_with(
            "monophyletic\n"
            "not_monophyletic\n"
            "insufficient_taxon_representation\n",
        )

    def test_print_results_empty_text_output_prints_nothing(self):
        """Empty text results should not emit a blank line."""
        with patch("sys.stdout", new_callable=StringIO) as output:
            self.checker.print_results([])

        self.assertEqual(output.getvalue(), "")

    def test_print_results_json_payload(self):
        """JSON output preserves row order and unexpected taxa sorting."""
        self.checker.json_output = True
        res_arr = [
            ["monophyletic", []],
            ["not_monophyletic", ["taxa2", "taxa1"]],
            ["insufficient_taxon_representation"],
        ]

        with patch(
            "phykit.services.tree.hidden_paralogy_check.print_json"
        ) as mocked_json:
            self.checker.print_results(res_arr)

        rows = [
            {
                "clade_index": 1,
                "status": "monophyletic",
                "unexpected_taxa": [],
            },
            {
                "clade_index": 2,
                "status": "not_monophyletic",
                "unexpected_taxa": ["taxa1", "taxa2"],
            },
            {
                "clade_index": 3,
                "status": "insufficient_taxon_representation",
                "unexpected_taxa": [],
            },
        ]
        mocked_json.assert_called_once_with({"rows": rows, "clades": rows})

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

        # First read is for the optional batch exact-clade index. The mock tree
        # is nonstandard, so processing falls back to per non-exact clade reads.
        mock_master_tree = Mock()
        mock_phylo_read.side_effect = [mock_master_tree, mock_tree1, mock_tree2]

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

    def test_process_clade_batch_exact_clade_skips_reroot_and_mrca(self):
        """Exact clades in batch mode use the indexed master tree."""
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_root_with_outgroup(*_args, **_kwargs):
            raise AssertionError("root_with_outgroup should not be called")

        def fail_common_ancestor(*_args, **_kwargs):
            raise AssertionError("common_ancestor should not be called")

        tree.root_with_outgroup = fail_root_with_outgroup
        tree.common_ancestor = fail_common_ancestor

        fd, path = tempfile.mkstemp(suffix=".tre", text=True)
        os.close(fd)
        try:
            with patch(
                "phykit.services.tree.hidden_paralogy_check.Phylo.read",
                return_value=tree,
            ) as mocked_read:
                results = HiddenParalogyCheck._process_clade_batch(
                    [["A", "B"]], path, frozenset(["A", "B", "C", "D"])
                )
        finally:
            os.unlink(path)

        mocked_read.assert_called_once_with(path, "newick")
        self.assertEqual(results, [["monophyletic", []]])

    def test_process_clade_batch_nonexact_uses_direct_terminal_names(self):
        """Non-exact batch fallback avoids generic terminal materialization."""
        index_tree = Phylo.read(StringIO("((A:1,C:1):1,(B:1,D:1):1);"), "newick")
        work_tree = Phylo.read(StringIO("((A:1,C:1):1,(B:1,D:1):1);"), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("standard clades should use direct terminal names")

        original_get_terminals = TreeMixin.get_terminals
        TreeMixin.get_terminals = fail_get_terminals
        try:
            with patch(
                "phykit.services.tree.hidden_paralogy_check.Phylo.read",
                side_effect=[index_tree, work_tree],
            ):
                results = HiddenParalogyCheck._process_clade_batch(
                    [["A", "B"]], "test.tre", frozenset(["A", "B", "C", "D"])
                )
        finally:
            TreeMixin.get_terminals = original_get_terminals

        self.assertEqual(results[0][0], "not_monophyletic")
        self.assertTrue(results[0][1])

    def test_terminal_names_direct(self):
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:1);"), "newick")

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("standard clades should use direct traversal")

        original_get_terminals = TreeMixin.get_terminals
        TreeMixin.get_terminals = fail_get_terminals
        try:
            result = HiddenParalogyCheck._terminal_names_direct(tree.root)
        finally:
            TreeMixin.get_terminals = original_get_terminals

        self.assertEqual(
            result,
            {"A", "B", "C"},
        )

    def test_build_exact_clade_index_uses_direct_postorder(self):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard exact-clade index should use direct postorder")

        tree.find_clades = fail_find_clades

        exact_clades = HiddenParalogyCheck._build_exact_clade_index(tree)

        self.assertIn(frozenset({"A", "B"}), exact_clades)
        self.assertIn(frozenset({"C", "D"}), exact_clades)
        self.assertIn(frozenset({"A", "B", "C", "D"}), exact_clades)

    def test_build_exact_clade_index_handles_polytomies(self):
        tree = Phylo.read(StringIO("(A:1,B:1,C:1,D:1);"), "newick")

        exact_clades = HiddenParalogyCheck._build_exact_clade_index(tree)

        self.assertIn(frozenset({"A"}), exact_clades)
        self.assertIn(frozenset({"A", "B", "C", "D"}), exact_clades)

    @patch('phykit.services.tree.hidden_paralogy_check.Phylo.read')
    def test_process_clade_batch_uses_provided_exact_clade_index(self, mock_read):
        """Provided exact-clade indexes avoid batch tree reads."""
        mock_read.side_effect = AssertionError("exact clade should not read tree")

        results = HiddenParalogyCheck._process_clade_batch(
            [["A", "B"]],
            "unused.tre",
            frozenset(["A", "B", "C", "D"]),
            exact_clades={frozenset(["A", "B"])},
        )

        mock_read.assert_not_called()
        self.assertEqual(results, [["monophyletic", []]])

    @patch('phykit.services.tree.hidden_paralogy_check.sys.stdout.write')
    @patch('phykit.services.tree.hidden_paralogy_check.Phylo.read')
    def test_run_sequential_processing(self, mock_phylo_read, mock_write):
        """Test run method with sequential processing (small dataset)"""
        # Setup mock tree
        mock_master_tree = Mock()
        mock_tree = Mock()
        mock_subtree = Mock()

        # Master tree tips
        self.checker.read_tree_file_unmodified = Mock(return_value=mock_master_tree)
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
        mock_write.assert_called_once_with("monophyletic\n")

    @patch('phykit.services.tree.hidden_paralogy_check.sys.stdout.write')
    def test_run_sequential_exact_clade_skips_tree_reread(self, mock_write):
        """Exact clades are classified from the master tree without rerooting."""
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        self.checker.read_tree_file_unmodified = Mock(return_value=tree)
        self.checker.read_clades_file = Mock(return_value=[["A", "B"]])

        with patch(
            "phykit.services.tree.hidden_paralogy_check.Phylo.read",
            side_effect=AssertionError("exact clade should not reread tree"),
        ):
            self.checker.run()

        mock_write.assert_called_once_with("monophyletic\n")

    @patch('phykit.services.tree.hidden_paralogy_check.sys.stdout.write')
    def test_run_sequential_nonexact_uses_direct_terminal_names(self, mock_write):
        """Sequential non-exact clades avoid generic terminal materialization."""
        master_tree = Phylo.read(StringIO("((A:1,C:1):1,(B:1,D:1):1);"), "newick")
        work_tree = Phylo.read(StringIO("((A:1,C:1):1,(B:1,D:1):1);"), "newick")
        self.checker.read_tree_file_unmodified = Mock(return_value=master_tree)
        self.checker.read_clades_file = Mock(return_value=[["A", "B"]])

        def fail_get_terminals(*_args, **_kwargs):
            raise AssertionError("standard clades should use direct terminal names")

        original_get_terminals = TreeMixin.get_terminals
        TreeMixin.get_terminals = fail_get_terminals
        try:
            with patch(
                "phykit.services.tree.hidden_paralogy_check.Phylo.read",
                return_value=work_tree,
            ):
                self.checker.run()
        finally:
            TreeMixin.get_terminals = original_get_terminals

        mock_write.assert_called_once_with("not_monophyletic\n")

    @patch('phykit.services.tree.hidden_paralogy_check.sys.stdout.write')
    @patch('multiprocessing.Pool')
    def test_run_parallel_processing(self, mock_pool_class, mock_write):
        """Test run method with parallel processing (large dataset)"""
        self.checker.MP_MIN_CLADES = 10
        # Setup mock tree
        mock_master_tree = Mock()

        # Master tree tips
        self.checker.read_tree_file_unmodified = Mock(return_value=mock_master_tree)
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

        self.assertEqual(
            mock_write.call_args.args[0].splitlines(),
            ["monophyletic"] * 15,
        )

    @patch('phykit.services.tree.hidden_paralogy_check.sys.stdout.write')
    @patch('multiprocessing.Pool')
    def test_run_medium_exact_clades_skip_pool(self, mock_pool_class, mock_write):
        """Default processing keeps medium exact-clade lists sequential."""
        tree = Phylo.read(StringIO("(((A,B),(C,D)),((E,F),(G,H)));"), "newick")
        exact_clades = [["A", "B"], ["C", "D"], ["E", "F"]] * 5

        self.checker.read_tree_file_unmodified = Mock(return_value=tree)
        self.checker.get_tip_names_from_tree = Mock(
            return_value=["A", "B", "C", "D", "E", "F", "G", "H"]
        )
        self.checker.read_clades_file = Mock(return_value=exact_clades)

        self.checker.run()

        mock_pool_class.assert_not_called()
        self.assertEqual(
            mock_write.call_args.args[0].splitlines(),
            ["monophyletic"] * 15,
        )

    @patch('phykit.services.tree.hidden_paralogy_check.sys.stdout.write')
    @patch('phykit.services.tree.hidden_paralogy_check.Phylo.read')
    def test_run_with_insufficient_taxa(self, mock_phylo_read, mock_write):
        """Test run with insufficient taxa representation"""
        mock_master_tree = Mock()
        mock_tree = Mock()
        mock_phylo_read.return_value = mock_tree

        self.checker.read_tree_file_unmodified = Mock(return_value=mock_master_tree)
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
        mock_write.assert_called_once_with("insufficient_taxon_representation\n")

    @patch('phykit.services.tree.hidden_paralogy_check.sys.stdout.write')
    @patch('phykit.services.tree.hidden_paralogy_check.Phylo.read')
    def test_run_with_not_monophyletic(self, mock_phylo_read, mock_write):
        """Test run with not monophyletic result"""
        mock_master_tree = Mock()
        mock_tree = Mock()
        mock_subtree = Mock()

        self.checker.read_tree_file_unmodified = Mock(return_value=mock_master_tree)
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
        mock_write.assert_called_once_with("not_monophyletic\n")

    @patch('multiprocessing.cpu_count')
    @patch('multiprocessing.Pool')
    def test_run_parallel_with_cpu_count(self, mock_pool_class, mock_cpu_count):
        """Test parallel processing respects CPU count limit"""
        self.checker.MP_MIN_CLADES = 10
        mock_cpu_count.return_value = 16  # Many CPUs

        mock_master_tree = Mock()
        self.checker.read_tree_file_unmodified = Mock(return_value=mock_master_tree)
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

        # Duplicate and off-tree taxa should still reduce to the shared clade set.
        clade_batch = [["taxa1", "taxa1", "taxa2", "taxa_not_in_tree"]]
        master_tree_tips = frozenset(["taxa1", "taxa2", "taxa3"])
        exact_clades = {frozenset(["taxa1", "taxa2"])}

        results = HiddenParalogyCheck._process_clade_batch(
            clade_batch, "test.tre", master_tree_tips, exact_clades=exact_clades
        )

        self.assertEqual(results[0], ["monophyletic", []])


if __name__ == '__main__':
    unittest.main()
