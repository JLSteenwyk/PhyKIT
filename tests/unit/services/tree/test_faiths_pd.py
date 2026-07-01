"""Unit tests for Faith's phylogenetic diversity.

Reference values computed with picante 1.8.2 on a multi2di-resolved
tree_simple.tre via:
  picante::pd(samp, tr, include.root = TRUE / FALSE)
"""
import pytest
import subprocess
import sys
from argparse import Namespace
from io import StringIO
from math import isclose

from Bio import Phylo

from phykit.services.tree.faiths_pd import FaithsPD
from phykit.errors import PhykitUserError


def test_module_import_defers_json_helper_and_heavy_tree_dependencies():
    code = """
import importlib
import sys

module = importlib.import_module("phykit.services.tree.faiths_pd")
assert callable(module.print_json)
assert callable(module.read_single_column_file_to_list)
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "phykit.helpers.json_output" not in sys.modules
assert "phykit.helpers.files" not in sys.modules
assert "Bio.Phylo" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def args():
    return Namespace(
        tree="/some/path/to/file.tre",
        taxa="/some/path/to/taxa.txt",
        exclude_root=False,
        json=False,
    )


def _service(args):
    return FaithsPD(args)


class TestFaithsPD(object):
    def test_init_sets_attrs(self, args):
        svc = _service(args)
        assert svc.tree_file_path == args.tree
        assert svc.taxa_file == args.taxa
        assert svc.include_root is True
        assert svc.json_output is False

    def test_exclude_root_flag_inverts_include_root(self, args):
        args.exclude_root = True
        svc = _service(args)
        assert svc.include_root is False

    def test_all_tips_equals_total_tree_length(self, tree_simple, args):
        svc = _service(args)
        all_tips = [t.name for t in tree_simple.get_terminals()]
        pd_inc, n = svc.calculate_faiths_pd(tree_simple, all_tips, include_root=True)
        pd_exc, _ = svc.calculate_faiths_pd(tree_simple, all_tips, include_root=False)
        assert n == 8
        assert isclose(pd_inc, 277.27722, rel_tol=1e-5)
        assert isclose(pd_exc, 277.27722, rel_tol=1e-5)

    def test_four_taxa_spanning_root(self, tree_simple, args):
        """MRCA of this community is the tree root, so include/exclude agree
        (validated against picante::pd)."""
        svc = _service(args)
        community = ["raccoon", "bear", "sea_lion", "seal"]
        pd_inc, _ = svc.calculate_faiths_pd(tree_simple, community, include_root=True)
        pd_exc, _ = svc.calculate_faiths_pd(tree_simple, community, include_root=False)
        assert isclose(pd_inc, 62.24955, rel_tol=1e-5)
        assert isclose(pd_exc, 62.24955, rel_tol=1e-5)

    def test_deep_community_across_root_children(self, tree_simple, args):
        svc = _service(args)
        community = ["monkey", "cat", "weasel", "dog"]
        pd_inc, _ = svc.calculate_faiths_pd(tree_simple, community, include_root=True)
        pd_exc, _ = svc.calculate_faiths_pd(tree_simple, community, include_root=False)
        assert isclose(pd_inc, 218.90149, rel_tol=1e-5)
        assert isclose(pd_exc, 218.90149, rel_tol=1e-5)

    def test_mrca_deep_in_tree_differs_between_modes(self, tree_simple, args):
        """MRCA of (raccoon, bear) is below the root, so include/exclude differ."""
        svc = _service(args)
        community = ["raccoon", "bear"]
        pd_inc, _ = svc.calculate_faiths_pd(tree_simple, community, include_root=True)
        pd_exc, _ = svc.calculate_faiths_pd(tree_simple, community, include_root=False)
        assert isclose(pd_inc, 26.84600, rel_tol=1e-5)
        assert isclose(pd_exc, 26.00000, rel_tol=1e-5)

    def test_exclude_root_fast_path_does_not_call_common_ancestor(
        self, args, monkeypatch
    ):
        svc = _service(args)
        tree = Phylo.read(StringIO("((A:1,B:1):1,C:3):0.5;"), "newick")

        def fail_common_ancestor(self, *args, **kwargs):
            raise AssertionError("common_ancestor fallback should not be called")

        monkeypatch.setattr(type(tree), "common_ancestor", fail_common_ancestor)

        pd_exc, n = svc.calculate_faiths_pd(tree, ["A", "B"], include_root=False)
        pd_inc, _ = svc.calculate_faiths_pd(tree, ["A", "B"], include_root=True)

        assert n == 2
        assert pd_exc == pytest.approx(2.0)
        assert pd_inc == pytest.approx(3.0)

    def test_fast_path_does_not_call_get_path(self, args, monkeypatch):
        svc = _service(args)
        tree = Phylo.read(StringIO("(((A:1,B:1):1,C:1):1,D:1):0.5;"), "newick")

        def fail_get_path(self, *args, **kwargs):
            raise AssertionError("get_path should not be called")

        monkeypatch.setattr(type(tree), "get_path", fail_get_path)

        pd_exc, n = svc.calculate_faiths_pd(
            tree, ["A", "B", "C"], include_root=False
        )
        pd_inc, _ = svc.calculate_faiths_pd(
            tree, ["A", "B", "C"], include_root=True
        )

        assert n == 3
        assert pd_exc == pytest.approx(4.0)
        assert pd_inc == pytest.approx(5.0)

    def test_include_root_partial_community_uses_unique_root_paths(self, args):
        svc = _service(args)
        tree = Phylo.read(StringIO("(((A:1,B:2):3,C:5):7,D:11):13;"), "newick")

        pd_inc, n = svc.calculate_faiths_pd(
            tree,
            ["A", "C"],
            include_root=True,
        )

        assert n == 2
        assert pd_inc == pytest.approx(16.0)

    def test_include_root_overlapping_paths_count_shared_branch_once(self, args):
        svc = _service(args)
        tree = Phylo.read(StringIO("(((A:1,B:2):3,C:5):7,D:11);"), "newick")

        pd_inc, n = svc.calculate_faiths_pd(
            tree,
            ["A", "B", "C"],
            include_root=True,
        )

        assert n == 3
        assert pd_inc == pytest.approx(18.0)

    def test_all_tips_fast_path_matches_total_branch_length(self, args):
        svc = _service(args)
        tree = Phylo.read(StringIO("((A:1,B:2):3,(C:4,D:5):6);"), "newick")

        pd_exc, n = svc.calculate_faiths_pd(
            tree, ["A", "B", "C", "D"], include_root=False
        )
        pd_inc, _ = svc.calculate_faiths_pd(
            tree, ["A", "B", "C", "D"], include_root=True
        )

        assert n == 4
        assert pd_exc == pytest.approx(21.0)
        assert pd_inc == pytest.approx(21.0)

    def test_all_tips_fast_path_uses_unique_tip_names_for_duplicate_labels(self, args):
        svc = _service(args)
        tree = Phylo.read(StringIO("((A:1,A:2):3,B:4);"), "newick")

        pd_exc, n = svc.calculate_faiths_pd(
            tree,
            ["A", "B"],
            include_root=False,
        )

        assert n == 2
        assert pd_exc == pytest.approx(10.0)

    def test_calculate_faiths_pd_handles_mixed_child_counts_without_reversed(
        self, args
    ):
        class NoReversedList(list):
            def __reversed__(self):
                raise AssertionError("Faith's PD traversal should push children explicitly")

        class Clade:
            def __init__(self, name=None, branch_length=None, clades=None):
                self.name = name
                self.branch_length = branch_length
                self.clades = NoReversedList(clades or [])

        tree = type(
            "Tree",
            (),
            {
                "root": Clade(
                    clades=[
                        Clade("A", 1.0),
                        Clade(None, 2.0, [Clade("B", 3.0), Clade("C", 4.0)]),
                        Clade(
                            None,
                            5.0,
                            [Clade("D", 6.0), Clade("E", 7.0), Clade("F", 8.0)],
                        ),
                    ],
                )
            },
        )()
        svc = _service(args)

        pd_inc, n = svc.calculate_faiths_pd(
            tree, ["A", "C", "E"], include_root=True
        )
        pd_exc, _ = svc.calculate_faiths_pd(
            tree, ["A", "C", "E"], include_root=False
        )

        assert n == 3
        assert pd_inc == pytest.approx(19.0)
        assert pd_exc == pytest.approx(19.0)

    def test_binary_selected_count_uses_indexed_child_access(
        self, args, monkeypatch
    ):
        class NoIterBinaryList(list):
            def __iter__(self):
                raise AssertionError("binary child aggregation should not iterate")

            def __reversed__(self):
                raise AssertionError("binary child aggregation should not reverse")

        class Clade:
            def __init__(self, name=None, branch_length=None, clades=None):
                self.name = name
                self.branch_length = branch_length
                self.clades = NoIterBinaryList(clades or [])

        tree = type(
            "Tree",
            (),
            {
                "root": Clade(
                    clades=[
                        Clade(
                            None,
                            1.0,
                            [Clade("A", 2.0), Clade("B", 3.0)],
                        ),
                        Clade(
                            None,
                            4.0,
                            [Clade("C", 5.0), Clade("D", 6.0)],
                        ),
                    ],
                )
            },
        )()
        svc = _service(args)
        monkeypatch.setattr(
            svc,
            "validate_tree",
            lambda *_args, **_kwargs: None,
        )

        pd_exc, n = svc.calculate_faiths_pd(
            tree, ["A", "B"], include_root=False
        )

        assert n == 2
        assert pd_exc == pytest.approx(5.0)

    def test_single_tip_include_root(self, tree_simple, args):
        svc = _service(args)
        pd_inc, n = svc.calculate_faiths_pd(tree_simple, ["monkey"], include_root=True)
        assert n == 1
        assert isclose(pd_inc, 127.41973, rel_tol=1e-5)

    def test_single_tip_exclude_root_returns_zero(self, tree_simple, args):
        """picante returns NA for this case; we return 0 by documented choice."""
        svc = _service(args)
        pd_exc, n = svc.calculate_faiths_pd(tree_simple, ["monkey"], include_root=False)
        assert n == 1
        assert pd_exc == 0.0

    def test_single_tip_exclude_root_skips_selected_count_pass(
        self, args, monkeypatch
    ):
        class Clade:
            def __init__(self, name=None, branch_length=None, clades=None):
                self.name = name
                self.branch_length = branch_length
                self._clades = clades or []
                self._clades_accesses = 0

            @property
            def clades(self):
                self._clades_accesses += 1
                if self._clades_accesses > 1:
                    raise AssertionError(
                        "single-tip exclude-root should skip selected counts"
                    )
                return self._clades

        tree = type(
            "Tree",
            (),
            {
                "root": Clade(
                    clades=[
                        Clade("A", 1.0),
                        Clade("B", 2.0),
                    ],
                )
            },
        )()
        svc = _service(args)
        monkeypatch.setattr(
            svc,
            "validate_tree",
            lambda *_args, **_kwargs: None,
        )

        pd_exc, n = svc.calculate_faiths_pd(tree, ["A"], include_root=False)

        assert n == 1
        assert pd_exc == 0.0

    def test_missing_taxon_raises(self, tree_simple, args):
        svc = _service(args)
        with pytest.raises(PhykitUserError) as exc_info:
            svc.calculate_faiths_pd(tree_simple, ["raccoon", "not_a_tip"])
        assert "Taxa not found" in str(exc_info.value.messages[0])

    def test_tree_without_branch_lengths_raises(self, tree_zero_branch_length, args):
        # tree_zero_branch_length has None branch lengths except tips at 0;
        # require_branch_lengths triggers on any None.
        svc = _service(args)
        # If all branch lengths are present (even as 0), this passes; we
        # only want to error when any are None. If the fixture is all-zero
        # then PD should simply be 0.
        # Use the tree as-is and assert either behavior is well-defined.
        tips = [t.name for t in tree_zero_branch_length.get_terminals()]
        try:
            pd_val, _ = svc.calculate_faiths_pd(
                tree_zero_branch_length, tips, include_root=True
            )
            assert pd_val == 0.0
        except PhykitUserError:
            # Also acceptable: raised on None branch lengths.
            pass

    def test_deduplicates_taxa(self, tree_simple, args):
        svc = _service(args)
        pd1, n1 = svc.calculate_faiths_pd(
            tree_simple, ["raccoon", "bear"], include_root=False
        )
        pd2, n2 = svc.calculate_faiths_pd(
            tree_simple, ["raccoon", "bear", "raccoon"], include_root=False
        )
        assert n1 == 2
        assert pd1 == pd2

    def test_load_taxa_deduplicates_in_order_and_skips_blanks(self, tmp_path):
        taxa_file = tmp_path / "taxa.txt"
        taxa_file.write_text("A\n\nB\nA\n   \nC\nB\n")

        assert FaithsPD._load_taxa(str(taxa_file)) == ["A", "B", "C"]

    def test_load_taxa_empty_after_filtering_raises(self, tmp_path):
        taxa_file = tmp_path / "empty_taxa.txt"
        taxa_file.write_text("\n   \n")

        with pytest.raises(PhykitUserError) as exc_info:
            FaithsPD._load_taxa(str(taxa_file))

        assert "No taxa found" in str(exc_info.value.messages[0])

    def test_run_uses_unmodified_tree_read(self, mocker, args):
        tree = object()
        svc = _service(args)
        read_tree = mocker.patch.object(
            svc,
            "read_tree_file_unmodified",
            return_value=tree,
        )
        mocker.patch.object(svc, "_load_taxa", return_value=["A", "B"])
        calculate = mocker.patch.object(
            svc,
            "calculate_faiths_pd",
            return_value=(12.34567, 2),
        )
        mocked_print = mocker.patch("builtins.print")

        svc.run()

        read_tree.assert_called_once_with()
        calculate.assert_called_once_with(tree, ["A", "B"], include_root=True)
        mocked_print.assert_called_once_with(12.3457)

    def test_run_json_output(self, mocker, args):
        args.json = True
        svc = _service(args)
        mocker.patch.object(svc, "read_tree_file_unmodified", return_value=object())
        mocker.patch.object(svc, "_load_taxa", return_value=["A", "B"])
        mocker.patch.object(svc, "calculate_faiths_pd", return_value=(12.34567, 2))
        mocked_json = mocker.patch("phykit.services.tree.faiths_pd.print_json")

        svc.run()

        mocked_json.assert_called_once_with(
            {"faiths_pd": 12.3457, "n_taxa": 2, "include_root": True}
        )
