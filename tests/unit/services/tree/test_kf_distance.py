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
import pytest
from argparse import Namespace

from phykit.services.tree.kf_distance import KuhnerFelsensteinDistance


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
