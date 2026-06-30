"""
Unit tests for Kuhner-Felsenstein (branch score) distance.

The KF distance (Kuhner & Felsenstein 1994) is:
    KF = sqrt( sum_over_all_splits( (b1_i - b2_i)^2 ) )
where b1_i and b2_i are branch lengths for split i in each tree.
Both internal and terminal branches are included.

Expected values cross-validated against R's phangorn::KF.dist().
See tests/r_validation/validate_kf_distance.R for the R script.
"""
import copy
import math
import subprocess
import sys
import pytest
from argparse import Namespace
from io import StringIO

from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin

from phykit.services.tree.kf_distance import KuhnerFelsensteinDistance


def test_tree_distance_output_modules_defer_heavy_imports():
    lightweight_modules = [
        "phykit.services.tree.kf_distance",
        "phykit.services.tree.treeness",
        "phykit.services.tree.evolutionary_rate",
        "phykit.services.tree.collapse_branches",
        "phykit.services.tree.rf_distance",
        "phykit.services.tree.print_tree",
        "phykit.services.tree.spurious_sequence",
        "phykit.services.tree.rename_tree_tips",
    ]
    stats_modules = [
        "phykit.services.tree.patristic_distances",
        "phykit.services.tree.internal_branch_stats",
        "phykit.services.tree.terminal_branch_stats",
    ]
    code = (
        "import sys; "
        f"lightweight_modules = {lightweight_modules!r}; "
        f"stats_modules = {stats_modules!r}; "
        "modules = [__import__(module, fromlist=['*']) for module in lightweight_modules]; "
        "assert all(callable(module.print_json) for module in modules); "
        "assert 'json' not in sys.modules; "
        "assert 'phykit.helpers.json_output' not in sys.modules; "
        "assert 'Bio.Phylo' not in sys.modules; "
        "assert 'tqdm' not in sys.modules; "
        "assert 'numpy' not in sys.modules; "
        "[__import__(module) for module in stats_modules]; "
        "assert 'Bio.Phylo' not in sys.modules; "
        "assert 'tqdm' not in sys.modules"
    )

    subprocess.run([sys.executable, "-c", code], check=True)


def test_kf_distance_import_does_not_import_typing():
    code = (
        "import sys; "
        "import phykit.services.tree.kf_distance; "
        "assert 'typing' not in sys.modules; "
        "assert 'json' not in sys.modules; "
        "assert 'phykit.helpers.json_output' not in sys.modules; "
        "assert 'Bio.Phylo' not in sys.modules; "
        "assert 'numpy' not in sys.modules"
    )

    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def args():
    kwargs = dict(
        tree_zero="/some/path/to/file.tre",
        tree_one="/some/path/to/file.tre",
    )
    return Namespace(**kwargs)


class TestKuhnerFelsensteinDistance:
    def test_init_sets_tree_file_path(self, args):
        kf = KuhnerFelsensteinDistance(args)
        assert kf.tree_file_path == args.tree_zero
        assert kf.tree1_file_path == args.tree_one
        assert kf.output_file_path is None

    def test_process_args(self, args):
        kf = KuhnerFelsensteinDistance(args)
        parsed = kf.process_args(args)
        assert parsed["tree_file_path"] == args.tree_zero
        assert parsed["tree1_file_path"] == args.tree_one
        assert parsed["json_output"] is False

    def test_process_args_json(self):
        args = Namespace(
            tree_zero="/a.tre", tree_one="/b.tre", json=True
        )
        kf = KuhnerFelsensteinDistance(args)
        parsed = kf.process_args(args)
        assert parsed["json_output"] is True

    def test_calculate_kf_distance_identical_trees(self, tree_simple, args):
        kf = KuhnerFelsensteinDistance(args)
        tree_clone = copy.deepcopy(tree_simple)
        plain_kf, normalized_kf = kf.calculate_kf_distance(tree_simple, tree_clone)
        assert plain_kf == 0.0
        assert normalized_kf == 0.0

    def test_calculate_kf_distance_same_object_skips_split_extraction(
        self,
        tree_simple,
        args,
        mocker,
    ):
        kf = KuhnerFelsensteinDistance(args)
        mocker.patch.object(
            kf,
            "_get_splits_with_lengths",
            side_effect=AssertionError("same-object KF should skip split extraction"),
        )

        assert kf.calculate_kf_distance(tree_simple, tree_simple) == (0.0, 0.0)

    def test_calculate_kf_distance_different_trees(self, tree_simple, tree_simple_other, args):
        kf = KuhnerFelsensteinDistance(args)

        # Root on same taxon
        tree_simple.root_with_outgroup("dog")
        tree_simple_other.root_with_outgroup("dog")

        plain_kf, normalized_kf = kf.calculate_kf_distance(tree_simple, tree_simple_other)

        # Cross-validated against R's phangorn::KF.dist()
        assert round(plain_kf, 4) == 51.9041
        assert isinstance(normalized_kf, float)
        assert normalized_kf > 0.0

    def test_get_splits_with_lengths_includes_terminals(self, tree_simple, args):
        kf = KuhnerFelsensteinDistance(args)
        tree_simple.root_with_outgroup("dog")
        splits = kf._get_splits_with_lengths(tree_simple)

        # Should include terminal splits (single-taxon frozensets)
        terminal_names = {t.name for t in tree_simple.get_terminals()}
        for name in terminal_names:
            assert frozenset([name]) in splits

    def test_get_splits_with_lengths_includes_internals(self, tree_simple, args):
        kf = KuhnerFelsensteinDistance(args)
        tree_simple.root_with_outgroup("dog")
        splits = kf._get_splits_with_lengths(tree_simple)

        # Should include at least one multi-taxon split
        multi_taxon_splits = [s for s in splits if len(s) > 1]
        assert len(multi_taxon_splits) > 0

    def test_get_splits_with_lengths_matches_legacy_terminal_scan(
        self, tree_simple, args
    ):
        kf = KuhnerFelsensteinDistance(args)
        tree_simple.root_with_outgroup("dog")

        expected = {}
        for clade in tree_simple.find_clades(order="postorder"):
            if clade == tree_simple.root:
                continue
            tips = frozenset(t.name for t in clade.get_terminals())
            expected[tips] = (
                clade.branch_length
                if clade.branch_length is not None
                else 0.0
            )

        assert kf._get_splits_with_lengths(tree_simple) == expected

    def test_get_splits_with_lengths_avoids_repeated_terminal_scans(
        self, monkeypatch, tree_simple, args
    ):
        from Bio.Phylo.BaseTree import Clade

        kf = KuhnerFelsensteinDistance(args)
        tree_simple.root_with_outgroup("dog")
        terminal_names = {t.name for t in tree_simple.get_terminals()}

        def fail_get_terminals(self):
            raise AssertionError("get_terminals should not be called")

        monkeypatch.setattr(Clade, "get_terminals", fail_get_terminals)

        splits = kf._get_splits_with_lengths(tree_simple)

        for name in terminal_names:
            assert frozenset([name]) in splits

    def test_get_splits_with_lengths_uses_direct_postorder(self, monkeypatch, args):
        kf = KuhnerFelsensteinDistance(args)
        tree = Phylo.read(StringIO("((A:2,B:3):5,C:7);"), "newick")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError(
                "standard split extraction should use direct postorder"
            )

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        splits = kf._get_splits_with_lengths(tree)

        assert splits == {
            frozenset({"A"}): 2.0,
            frozenset({"B"}): 3.0,
            frozenset({"A", "B"}): 5.0,
            frozenset({"C"}): 7.0,
        }

    def test_get_splits_with_lengths_handles_mixed_child_counts(self, args):
        kf = KuhnerFelsensteinDistance(args)
        tree = Phylo.read(
            StringIO("(((A:1):2,(B:3,C:4,D:5):6):7,E:8);"),
            "newick",
        )

        splits = kf._get_splits_with_lengths(tree)

        assert splits == {
            frozenset({"A"}): 2.0,
            frozenset({"B"}): 3.0,
            frozenset({"C"}): 4.0,
            frozenset({"D"}): 5.0,
            frozenset({"B", "C", "D"}): 6.0,
            frozenset({"A", "B", "C", "D"}): 7.0,
            frozenset({"E"}): 8.0,
        }

    def test_get_first_tip_name_from_tree_uses_direct_tree_traversal(
        self, monkeypatch, args
    ):
        kf = KuhnerFelsensteinDistance(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:2);"), "newick")

        def fail_get_terminals(self):
            raise AssertionError("get_terminals should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_get_terminals)

        assert kf.get_first_tip_name_from_tree(tree) == "A"

    def test_get_splits_with_lengths_none_branch_length(self, args):
        """None branch lengths should be treated as 0.0."""
        from Bio import Phylo
        from io import StringIO

        kf = KuhnerFelsensteinDistance(args)
        # Tree with no branch lengths
        tree = Phylo.read(StringIO("(A,B,(C,D));"), "newick")
        splits = kf._get_splits_with_lengths(tree)

        for bl in splits.values():
            assert bl == 0.0

    def test_kf_distance_zero_length_trees(self, args):
        """Both trees with zero-length branches should return (0.0, 0.0)."""
        from Bio import Phylo
        from io import StringIO

        kf = KuhnerFelsensteinDistance(args)
        t0 = Phylo.read(StringIO("(A:0,B:0,(C:0,D:0):0);"), "newick")
        t1 = Phylo.read(StringIO("(A:0,B:0,(C:0,D:0):0);"), "newick")
        plain_kf, normalized_kf = kf.calculate_kf_distance(t0, t1)
        assert plain_kf == 0.0
        assert normalized_kf == 0.0

    def test_calculate_kf_distance_accumulates_common_and_unique_splits(
        self, mocker, args
    ):
        kf = KuhnerFelsensteinDistance(args)
        tree_zero = object()
        tree_one = object()
        split_a = frozenset({"A"})
        split_b = frozenset({"B"})
        split_ab = frozenset({"A", "B"})
        split_c = frozenset({"C"})
        mocker.patch.object(
            kf,
            "_get_splits_with_lengths",
            side_effect=[
                {split_a: 1.0, split_b: 2.0, split_ab: 5.0},
                {split_a: 1.5, split_b: 0.5, split_c: 4.0},
            ],
        )

        plain_kf, normalized_kf = kf.calculate_kf_distance(tree_zero, tree_one)

        expected_plain = math.sqrt(
            (1.0 - 1.5) ** 2
            + (2.0 - 0.5) ** 2
            + 5.0 ** 2
            + 4.0 ** 2
        )
        assert plain_kf == expected_plain
        assert normalized_kf == expected_plain / 14.0

    def test_run_text_output(self, mocker, args):
        kf = KuhnerFelsensteinDistance(args)
        tree_zero = mocker.Mock()
        tree_one = mocker.Mock()
        term = mocker.Mock()
        term.name = "A"
        tree_zero.get_terminals.return_value = [term]
        tree_one.get_terminals.return_value = [term]
        mocker.patch.object(kf, "read_tree_file", return_value=tree_zero)
        mocker.patch.object(kf, "read_tree1_file", return_value=tree_one)
        mocker.patch.object(kf, "get_tip_names_from_tree", side_effect=[["A"], ["A"]])
        mocker.patch.object(kf, "calculate_kf_distance", return_value=(5.1234, 0.3456))
        mocked_print = mocker.patch("builtins.print")

        kf.run()

        mocked_print.assert_called_once_with("5.1234\t0.3456")

    def test_run_json_output(self, mocker):
        args = Namespace(tree_zero="/a.tre", tree_one="/b.tre", json=True)
        kf = KuhnerFelsensteinDistance(args)
        tree_zero = mocker.Mock()
        tree_one = mocker.Mock()
        term = mocker.Mock()
        term.name = "A"
        tree_zero.get_terminals.return_value = [term]
        tree_one.get_terminals.return_value = [term]
        mocker.patch.object(kf, "read_tree_file", return_value=tree_zero)
        mocker.patch.object(kf, "read_tree1_file", return_value=tree_one)
        mocker.patch.object(kf, "get_tip_names_from_tree", side_effect=[["A"], ["A"]])
        mocker.patch.object(kf, "calculate_kf_distance", return_value=(5.1234, 0.3456))
        mocked_json = mocker.patch("phykit.services.tree.kf_distance.print_json")

        kf.run()

        payload = mocked_json.call_args.args[0]
        assert payload["plain_kf"] == 5.1234
        assert payload["normalized_kf"] == 0.3456

    def test_run_no_shared_taxa_exits(self, mocker, args):
        kf = KuhnerFelsensteinDistance(args)
        tree_zero = mocker.Mock()
        tree_one = mocker.Mock()
        mocker.patch.object(kf, "read_tree_file", return_value=tree_zero)
        mocker.patch.object(kf, "read_tree1_file", return_value=tree_one)
        mocker.patch.object(
            kf, "get_tip_names_from_tree",
            side_effect=[["A", "B"], ["C", "D"]],
        )

        with pytest.raises(SystemExit) as exc:
            kf.run()
        assert exc.value.code == 2
