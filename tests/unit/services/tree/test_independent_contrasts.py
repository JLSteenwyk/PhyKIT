"""
Unit tests for independent_contrasts (Felsenstein's PIC).

Computes phylogenetically independent contrasts for continuous traits.
Cross-validated against R's ape::pic(). Individual contrasts at
polytomy-resolved nodes may differ due to resolution order, but the
sum of squared contrasts is invariant (verified to match R within
floating-point tolerance).

See tests/r_validation/validate_pic.R for the R validation script.
"""
import pytest
import numpy as np
import subprocess
import sys
import builtins
import json
import math
from argparse import Namespace
from pathlib import Path
from Bio.Phylo.BaseTree import TreeMixin

from phykit.errors import PhykitUserError
from phykit.services.tree.independent_contrasts import IndependentContrasts
import phykit.services.tree.independent_contrasts as ic_module


def test_module_import_does_not_import_numpy():
    code = """
import sys
import phykit.services.tree.independent_contrasts as module

assert callable(module._json_dumps)
assert hasattr(module.np, "__getattr__")
assert "typing" not in sys.modules
assert "json" not in sys.modules
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_lazy_numpy_caches_resolved_attributes():
    lazy_np = ic_module._LazyNumpy()

    asarray_attr = lazy_np.asarray

    assert lazy_np.__dict__["asarray"] is asarray_attr
    assert lazy_np.asarray is asarray_attr
    assert lazy_np._module is not None


def test_small_contrast_summary_does_not_import_numpy():
    code = """
import sys
from phykit.services.tree.independent_contrasts import _contrast_summary_stats

mean_abs, variance = _contrast_summary_stats([1.0, -2.0, 3.0])
assert round(mean_abs, 6) == 2.0
assert round(variance, 6) == 6.333333
assert "numpy" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def args():
    return Namespace(
        tree="tests/sample_files/tree_simple.tre",
        trait_data="tests/sample_files/tree_simple_traits.tsv",
    )


class TestIndependentContrastsInit:
    def test_init_sets_fields(self, args):
        ic = IndependentContrasts(args)
        assert ic.tree_file_path == args.tree
        assert ic.trait_data_path == args.trait_data
        assert ic.json_output is False

    def test_process_args_json(self):
        args = Namespace(
            tree="t.tre", trait_data="d.tsv", json=True
        )
        ic = IndependentContrasts(args)
        assert ic.json_output is True


class TestPICComputation:
    def test_tree_scans_reject_legacy_and_malformed_shapes(self, args):
        ic = IndependentContrasts(args)

        assert ic._needs_default_branch_lengths(object()) is True
        assert ic._has_polytomies(object()) is None
        assert ic._postorder_clades_fast(object()) is None
        assert ic._should_use_text_summary_path(object()) is False

        class MalformedRoot:
            clades = [object()]
            branch_length = None

        class MalformedTree:
            root = MalformedRoot()

        malformed = MalformedTree()
        assert ic._needs_default_branch_lengths(malformed) is True
        assert ic._has_polytomies(malformed) is None
        assert ic._postorder_clades_fast(malformed) is None
        assert ic._should_use_text_summary_path(malformed) is False

    def test_text_summary_path_detects_deep_tree(self, args):
        ic = IndependentContrasts(args)

        class Clade:
            def __init__(self, child=None):
                self.clades = [] if child is None else [child]

        root = Clade()
        for _ in range(66):
            root = Clade(root)

        class DeepTree:
            pass

        tree = DeepTree()
        tree.root = root

        assert ic._should_use_text_summary_path(tree) is True

    def test_correct_number_of_contrasts(self, args):
        """n tips should produce n-1 contrasts."""
        ic = IndependentContrasts(args)
        tree = ic.read_tree_file()
        import copy
        tree = copy.deepcopy(tree)
        tree_tips = [t.name for t in tree.get_terminals()]
        tip_traits = ic._parse_trait_data(args.trait_data, tree_tips)
        shared = set(tip_traits.keys())
        tips_to_prune = [t for t in tree_tips if t not in shared]
        if tips_to_prune:
            tree = ic.prune_tree_using_taxa_list(tree, tips_to_prune)
        ic._resolve_polytomies(tree)
        contrasts, _ = ic._compute_pic(tree, tip_traits)
        assert len(contrasts) == len(tip_traits) - 1

    def test_sum_of_squared_contrasts_matches_r(self, args):
        """Sum of squared contrasts should match R's ape::pic() exactly.

        R value: sum(pic(trait_vec, multi2di(tree))^2) = 0.307253
        This is invariant to polytomy resolution order.
        """
        ic = IndependentContrasts(args)
        tree = ic.read_tree_file()
        import copy
        tree = copy.deepcopy(tree)
        tree_tips = [t.name for t in tree.get_terminals()]
        tip_traits = ic._parse_trait_data(args.trait_data, tree_tips)
        shared = set(tip_traits.keys())
        tips_to_prune = [t for t in tree_tips if t not in shared]
        if tips_to_prune:
            tree = ic.prune_tree_using_taxa_list(tree, tips_to_prune)
        ic._resolve_polytomies(tree)
        contrasts, _ = ic._compute_pic(tree, tip_traits)

        ss = sum(c ** 2 for c in contrasts)
        # Cross-validated against R: sum(pic(...)^2) = 0.307253
        assert ss == pytest.approx(0.307253, abs=0.001)

    def test_known_contrasts_from_resolved_clades(self, args):
        """Contrasts from fully resolved clades should match R exactly."""
        ic = IndependentContrasts(args)
        tree = ic.read_tree_file()
        import copy
        tree = copy.deepcopy(tree)
        tree_tips = [t.name for t in tree.get_terminals()]
        tip_traits = ic._parse_trait_data(args.trait_data, tree_tips)
        shared = set(tip_traits.keys())
        tips_to_prune = [t for t in tree_tips if t not in shared]
        if tips_to_prune:
            tree = ic.prune_tree_using_taxa_list(tree, tips_to_prune)
        ic._resolve_polytomies(tree)
        contrasts, node_labels = ic._compute_pic(tree, tip_traits)

        # Find contrasts by their tip composition
        contrast_by_tips = {}
        for c, tips in zip(contrasts, node_labels):
            contrast_by_tips[frozenset(tips)] = c

        # bear-raccoon: R gives -0.264757
        br = frozenset(["bear", "raccoon"])
        assert contrast_by_tips[br] == pytest.approx(-0.264757, abs=0.0001)

        # sea_lion-seal: R gives 0.085732
        ss = frozenset(["sea_lion", "seal"])
        assert contrast_by_tips[ss] == pytest.approx(0.085732, abs=0.0001)

        # monkey-cat: R gives 0.003288
        mc = frozenset(["monkey", "cat"])
        assert contrast_by_tips[mc] == pytest.approx(0.003288, abs=0.0001)

    def test_identical_traits_produce_zero_contrasts(self):
        """If all taxa have the same trait value, all contrasts should be ~0."""
        from Bio import Phylo
        from io import StringIO
        import copy

        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        args = Namespace(tree="dummy", trait_data="dummy")
        ic = IndependentContrasts(args)
        ic._resolve_polytomies(tree)
        tip_traits = {"A": 5.0, "B": 5.0, "C": 5.0, "D": 5.0}
        contrasts, _ = ic._compute_pic(tree, tip_traits)
        for c in contrasts:
            assert abs(c) < 1e-10

    def test_resolve_polytomies_uses_direct_standard_tree_traversal(self, monkeypatch):
        from Bio import Phylo
        from io import StringIO

        tree = Phylo.read(StringIO("(A:1,B:1,C:1,D:1);"), "newick")
        args = Namespace(tree="dummy", trait_data="dummy")
        ic = IndependentContrasts(args)

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic clade traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        ic._resolve_polytomies(tree)

        assert all(len(clade.clades) <= 2 for clade in ic._postorder_clades_fast(tree))

    def test_resolve_polytomies_binary_tree_does_not_import_newick(
        self, monkeypatch
    ):
        from Bio import Phylo
        from io import StringIO

        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        args = Namespace(tree="dummy", trait_data="dummy")
        ic = IndependentContrasts(args)
        original_import = builtins.__import__

        def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "Bio.Phylo.Newick" or (
                name == "Bio.Phylo" and fromlist and "Newick" in fromlist
            ):
                raise AssertionError("binary trees should not import Newick helpers")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", guarded_import)

        ic._resolve_polytomies(tree)

        assert all(len(clade.clades) <= 2 for clade in ic._postorder_clades_fast(tree))

    def test_resolve_polytomies_matches_previous_traversal_output(self):
        from Bio import Phylo
        from Bio.Phylo import Newick
        from io import StringIO

        newick = "((A:1,B:1,C:1):1,(D:1,E:1,F:1):1,G:1);"
        expected = Phylo.read(StringIO(newick), "newick")
        actual = Phylo.read(StringIO(newick), "newick")
        args = Namespace(tree="dummy", trait_data="dummy")
        ic = IndependentContrasts(args)

        stack = [expected.root]
        while stack:
            clade = stack.pop()
            children = clade.clades
            if children:
                stack.extend(reversed(children))
            while len(children) > 2:
                child1 = children.pop()
                child2 = children.pop()
                new_internal = Newick.Clade(branch_length=0.0)
                new_internal.clades = [child1, child2]
                children.append(new_internal)

        ic._resolve_polytomies(actual)

        expected_out = StringIO()
        actual_out = StringIO()
        Phylo.write(expected, expected_out, "newick")
        Phylo.write(actual, actual_out, "newick")
        assert actual_out.getvalue() == expected_out.getvalue()

    def test_boolean_tree_scans_handle_mixed_child_counts(self, monkeypatch):
        from Bio import Phylo
        from io import StringIO

        mixed = Phylo.read(
            StringIO("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);"),
            "newick",
        )
        binary = Phylo.read(
            StringIO("((A:1,B:1):1,(C:1,D:1):1);"),
            "newick",
        )

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("standard tree scan should not use find_clades")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        assert IndependentContrasts._needs_default_branch_lengths(mixed) is False
        mixed.root.clades[2].clades[1].branch_length = None
        assert IndependentContrasts._needs_default_branch_lengths(mixed) is True
        assert IndependentContrasts._has_polytomies(binary) is False
        assert IndependentContrasts._has_polytomies(mixed) is True

    def test_compute_pic_does_not_rewalk_internal_terminals(self):
        from Bio import Phylo
        from io import StringIO

        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        args = Namespace(tree="dummy", trait_data="dummy")
        ic = IndependentContrasts(args)
        tip_traits = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}

        def fail_get_terminals():
            raise AssertionError("internal get_terminals should not be called")

        for clade in tree.find_clades():
            if not clade.is_terminal():
                clade.get_terminals = fail_get_terminals

        contrasts, node_labels = ic._compute_pic(tree, tip_traits)

        assert len(contrasts) == 3
        assert frozenset(node_labels[-1]) == {"A", "B", "C", "D"}

    def test_compute_pic_uses_direct_postorder_traversal(self, monkeypatch):
        from Bio import Phylo
        from io import StringIO

        newick = "((A:1,B:1):1,(C:1,D:1):1);"
        tip_traits = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        args = Namespace(tree="dummy", trait_data="dummy")
        ic = IndependentContrasts(args)
        expected = ic._compute_pic(Phylo.read(StringIO(newick), "newick"), tip_traits)
        tree = Phylo.read(StringIO(newick), "newick")

        def fail_traversal(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "get_terminals", fail_traversal)
        monkeypatch.setattr(TreeMixin, "find_clades", fail_traversal)

        assert ic._compute_pic(tree, tip_traits) == expected

    def test_compute_pic_generic_fallback_matches_standard_tree(
        self, monkeypatch, args
    ):
        from Bio import Phylo
        from io import StringIO

        newick = "((A:1,B:2):1,(C:3,D:4):2);"
        tip_traits = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        ic = IndependentContrasts(args)
        expected = ic._compute_pic(Phylo.read(StringIO(newick), "newick"), tip_traits)
        tree = Phylo.read(StringIO(newick), "newick")
        monkeypatch.setattr(ic, "_compute_pic_standard_tree", lambda *_args: None)

        observed = ic._compute_pic(tree, tip_traits)

        np.testing.assert_allclose(observed[0], expected[0])
        assert observed[1] == expected[1]

    def test_compute_pic_generic_fallback_skips_unresolved_or_missing_nodes(
        self, monkeypatch, args
    ):
        from Bio import Phylo
        from io import StringIO

        ic = IndependentContrasts(args)
        monkeypatch.setattr(ic, "_compute_pic_standard_tree", lambda *_args: None)
        polytomy = Phylo.read(StringIO("(A:1,B:1,C:1);"), "newick")
        missing_tip = Phylo.read(StringIO("(A:1,B:1);"), "newick")

        assert ic._compute_pic(polytomy, {"A": 1.0}) == ([], [])
        assert ic._compute_pic(missing_tip, {"A": 1.0}) == ([], [])

    def test_text_summary_generic_fallback_matches_full_labels(
        self, monkeypatch, args
    ):
        from Bio import Phylo
        from io import StringIO

        tree = Phylo.read(StringIO("((D:1,B:1):1,(C:1,A:1):1);"), "newick")
        tip_traits = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
        ic = IndependentContrasts(args)
        full_contrasts, node_labels = ic._compute_pic(tree, tip_traits)
        monkeypatch.setattr(
            ic,
            "_compute_pic_text_summaries_standard_tree",
            lambda *_args: None,
        )

        contrasts, summaries = ic._compute_pic_text_summaries(tree, tip_traits)

        assert contrasts == full_contrasts
        assert summaries == [(tips[:3], len(tips)) for tips in node_labels]

    def test_legacy_polytomy_resolution_fallbacks(self, args):
        from Bio import Phylo
        from io import StringIO

        ic = IndependentContrasts(args)
        tree = Phylo.read(StringIO("(A:1,B:1,C:1,D:1);"), "newick")

        class LegacyTree:
            def find_clades(self, *call_args, **kwargs):
                return tree.find_clades(*call_args, **kwargs)

        ic._resolve_polytomies(LegacyTree())

        assert len(tree.root.clades) == 2
        assert all(len(clade.clades) <= 2 for clade in tree.find_clades())

        fallback_tree = Phylo.read(
            StringIO("(A:1,B:1,C:1,D:1);"), "newick"
        )
        ic._resolve_polytomies_fallback(fallback_tree)
        assert len(fallback_tree.root.clades) == 2
        assert all(len(clade.clades) <= 2 for clade in fallback_tree.find_clades())

    def test_malformed_polytomy_traversal_delegates_to_fallback(self, args):
        from Bio import Phylo
        from io import StringIO

        fallback_tree = Phylo.read(
            StringIO("(A:1,B:1,C:1,D:1);"), "newick"
        )

        class MalformedRoot:
            clades = [object()]

        class MalformedTree:
            root = MalformedRoot()

            def find_clades(self, *call_args, **kwargs):
                return fallback_tree.find_clades(*call_args, **kwargs)

        ic = IndependentContrasts(args)

        ic._resolve_polytomies(MalformedTree())

        assert len(fallback_tree.root.clades) == 2

    def test_standard_pic_builders_decline_legacy_tree(self, args):
        ic = IndependentContrasts(args)

        assert ic._compute_pic_standard_tree(object(), {}) is None
        assert ic._compute_pic_text_summaries_standard_tree(object(), {}) is None

    def test_standard_pic_builders_skip_unresolved_or_missing_nodes(self, args):
        from Bio import Phylo
        from io import StringIO

        ic = IndependentContrasts(args)
        polytomy = Phylo.read(StringIO("(A:1,B:1,C:1);"), "newick")
        missing_tip = Phylo.read(StringIO("(A:1,B:1);"), "newick")

        assert ic._compute_pic_standard_tree(polytomy, {"A": 1.0}) == ([], [])
        assert ic._compute_pic_standard_tree(missing_tip, {"A": 1.0}) == ([], [])
        assert ic._compute_pic_text_summaries_standard_tree(
            polytomy, {"A": 1.0}
        ) == ([], [])
        assert ic._compute_pic_text_summaries_standard_tree(
            missing_tip, {"A": 1.0}
        ) == ([], [])

    def test_direct_postorder_preserves_left_to_right_order(self):
        from Bio import Phylo
        from io import StringIO

        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        args = Namespace(tree="dummy", trait_data="dummy")
        ic = IndependentContrasts(args)

        observed = [
            clade.name if clade.name is not None else "internal"
            for clade in ic._postorder_clades_fast(tree)
        ]

        assert observed == [
            "A", "B", "internal", "C", "D", "internal", "internal",
        ]

    def test_merge_sorted_labels_preserves_sorted_descendant_order(self, args):
        ic = IndependentContrasts(args)

        assert ic._merge_sorted_labels(
            ["A", "C", "E"],
            ["B", "D", "F"],
        ) == ["A", "B", "C", "D", "E", "F"]
        assert ic._merge_sorted_labels(["A", "B"], ["C"]) == ["A", "B", "C"]
        assert ic._merge_sorted_labels(["D"], ["A", "C"]) == ["A", "C", "D"]

    def test_merge_sorted_labels_handles_empty_inputs(self, args):
        ic = IndependentContrasts(args)

        assert ic._merge_sorted_labels([], ["A", "B"]) == ["A", "B"]
        assert ic._merge_sorted_labels(["A", "B"], []) == ["A", "B"]

    def test_limited_label_merge_preserves_sorted_preview_and_count(self, args):
        ic = IndependentContrasts(args)

        assert ic._merge_limited_sorted_labels(
            ["B", "D", "F"],
            10,
            ["A", "C", "E"],
            8,
            3,
        ) == (["A", "B", "C"], 18)
        assert ic._merge_limited_sorted_labels(
            ["A", "B", "C"],
            20,
            ["D", "E", "F"],
            12,
            3,
        ) == (["A", "B", "C"], 32)

    @pytest.mark.parametrize(
        ("left", "right", "limit", "expected"),
        [
            (["A"], ["B"], 0, []),
            ([], ["A", "B"], 1, ["A"]),
            (["A", "B"], [], 1, ["A"]),
            (["A", "B"], ["C", "D"], 3, ["A", "B", "C"]),
            (["C", "D"], ["A", "B"], 3, ["A", "B", "C"]),
            (["C", "D"], ["A", "B"], 2, ["A", "B"]),
            (["A", "C"], ["B"], 3, ["A", "B", "C"]),
        ],
    )
    def test_limited_label_merge_boundary_cases(
        self, args, left, right, limit, expected
    ):
        ic = IndependentContrasts(args)

        preview, count = ic._merge_limited_sorted_labels(
            left,
            len(left),
            right,
            len(right),
            limit,
        )

        assert preview == expected
        assert count == len(left) + len(right)

    def test_compute_pic_text_summaries_match_full_label_output(self):
        from Bio import Phylo
        from io import StringIO

        tree = Phylo.read(StringIO("((D:1,B:1):1,(C:1,A:1):1);"), "newick")
        args = Namespace(tree="dummy", trait_data="dummy")
        ic = IndependentContrasts(args)
        tip_traits = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}

        full_contrasts, node_labels = ic._compute_pic(tree, tip_traits)
        summary_contrasts, node_summaries = ic._compute_pic_text_summaries(
            tree, tip_traits
        )

        assert summary_contrasts == full_contrasts
        assert node_summaries == [(tips[:3], len(tips)) for tips in node_labels]

    def test_parse_trait_data_streams_and_preserves_tree_tip_order(
        self, monkeypatch, args
    ):
        from io import StringIO

        class StreamingTraits(StringIO):
            def readlines(self, *_args, **_kwargs):
                raise AssertionError("trait parser should stream input lines")

        def open_stream(*_args, **_kwargs):
            return StreamingTraits("C\t3\nA\t1\nB\t2\nD\t4\n")

        monkeypatch.setattr("builtins.open", open_stream)
        ic = IndependentContrasts(args)

        tip_traits = ic._parse_trait_data("traits.tsv", ["B", "A", "C"])

        assert list(tip_traits) == ["B", "A", "C"]
        assert tip_traits == {"B": 2.0, "A": 1.0, "C": 3.0}

    def test_parse_trait_data_preserves_exact_tree_tip_order(
        self, monkeypatch, args
    ):
        from io import StringIO

        def open_stream(*_args, **_kwargs):
            return StringIO("A\t1\nB\t2\nC\t3\n")

        monkeypatch.setattr("builtins.open", open_stream)
        ic = IndependentContrasts(args)

        tip_traits = ic._parse_trait_data("traits.tsv", ["A", "B", "C"])

        assert list(tip_traits) == ["A", "B", "C"]
        assert tip_traits == {"A": 1.0, "B": 2.0, "C": 3.0}

    def test_parse_trait_data_skips_whitespace_prefixed_comments(
        self, monkeypatch, args
    ):
        from io import StringIO

        def open_stream(*_args, **_kwargs):
            return StringIO("   # comment\nA\t1\n\t# comment\nB\t2\nC\t3\n")

        monkeypatch.setattr("builtins.open", open_stream)
        ic = IndependentContrasts(args)

        tip_traits = ic._parse_trait_data("traits.tsv", ["A", "B", "C"])

        assert list(tip_traits) == ["A", "B", "C"]
        assert tip_traits == {"A": 1.0, "B": 2.0, "C": 3.0}

    def test_parse_trait_data_valid_rows_avoid_split(self, monkeypatch, args):
        class NoSplitLine(str):
            def strip(self):
                return self

            def split(self, *_args, **_kwargs):
                raise AssertionError("valid trait rows should not allocate splits")

        class NoSplitTraits:
            def __enter__(self):
                return iter(
                    [
                        NoSplitLine("A\t1\n"),
                        NoSplitLine("bad\t2\textra\n"),
                        NoSplitLine("B\t2\n"),
                        NoSplitLine("C\t3\n"),
                    ]
                )

            def __exit__(self, *_args):
                return False

        monkeypatch.setattr("builtins.open", lambda *_args, **_kwargs: NoSplitTraits())
        ic = IndependentContrasts(args)

        tip_traits = ic._parse_trait_data("traits.tsv", ["A", "B", "C"])

        assert list(tip_traits) == ["A", "B", "C"]
        assert tip_traits == {"A": 1.0, "B": 2.0, "C": 3.0}

    def test_parse_trait_data_skips_non_numeric_and_malformed_rows(
        self, monkeypatch, args
    ):
        from io import StringIO

        monkeypatch.setattr(
            "builtins.open",
            lambda *_args, **_kwargs: StringIO(
                "bad\tnot-a-number\nmalformed\nA\t1\nB\t2\nC\t3\n"
            ),
        )
        ic = IndependentContrasts(args)

        observed = ic._parse_trait_data("traits.tsv", ["A", "B", "C"])

        assert observed == {"A": 1.0, "B": 2.0, "C": 3.0}

    def test_parse_trait_data_reports_missing_file(self, tmp_path, args):
        ic = IndependentContrasts(args)

        with pytest.raises(PhykitUserError) as excinfo:
            ic._parse_trait_data(str(tmp_path / "missing.tsv"), ["A", "B", "C"])

        assert "not found" in excinfo.value.messages[0]

    def test_parse_trait_data_reports_insufficient_shared_taxa(
        self, monkeypatch, args
    ):
        from io import StringIO

        monkeypatch.setattr(
            "builtins.open",
            lambda *_args, **_kwargs: StringIO("A\t1\nB\t2\n"),
        )
        ic = IndependentContrasts(args)

        with pytest.raises(PhykitUserError) as excinfo:
            ic._parse_trait_data("traits.tsv", ["A", "B", "C"])

        assert "Only 2 shared taxa" in excinfo.value.messages[0]


class TestPICRun:
    @staticmethod
    def _stub_run_tail(monkeypatch, ic, tip_traits, captured):
        monkeypatch.setattr(ic, "_parse_trait_data", lambda *_args: tip_traits)
        monkeypatch.setattr(ic, "_should_use_text_summary_path", lambda _tree: False)

        def compute_pic(tree, traits):
            captured["tree"] = tree
            captured["traits"] = traits
            return [0.0], [["A", "B"]]

        monkeypatch.setattr(ic, "_compute_pic", compute_pic)
        monkeypatch.setattr(ic, "_print_text", lambda *_args, **_kwargs: None)

    def test_run_uses_fast_tip_name_helper_for_setup(self, mocker, args):
        ic = IndependentContrasts(args)
        tip_name_spy = mocker.spy(ic, "get_tip_names_from_tree")

        ic.run()

        tip_name_spy.assert_called_once()

    def test_run_text_output(self, args, capsys):
        ic = IndependentContrasts(args)
        ic.run()
        captured = capsys.readouterr()
        assert "Number of contrasts: 7" in captured.out
        assert "Mean absolute contrast:" in captured.out
        assert "Variance of contrasts:" in captured.out

    def test_run_all_shared_binary_tree_uses_read_only_tree_without_copy_or_resolve(
        self, monkeypatch, args, capsys
    ):
        from Bio import Phylo
        from io import StringIO

        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ic = IndependentContrasts(args)
        captured = {}
        tip_traits = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}

        monkeypatch.setattr(ic, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(
            ic,
            "read_tree_file",
            lambda: (_ for _ in ()).throw(
                AssertionError("mutable tree read should not be used")
            ),
        )
        monkeypatch.setattr(
            ic,
            "_fast_copy",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("binary all-shared tree should not be copied")
            ),
        )
        monkeypatch.setattr(
            ic,
            "_resolve_polytomies",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("binary all-shared tree should not be resolved")
            ),
        )
        monkeypatch.setattr(
            ic,
            "prune_tree_using_taxa_list",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                AssertionError("all-shared tree should not be pruned")
            ),
        )
        self._stub_run_tail(monkeypatch, ic, tip_traits, captured)

        ic.run()

        capsys.readouterr()
        assert captured["tree"] is tree
        assert captured["traits"] == tip_traits
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}

    def test_run_missing_branch_lengths_copies_before_validation(
        self, monkeypatch, args, capsys
    ):
        from Bio import Phylo
        from io import StringIO

        tree = Phylo.read(StringIO("((A,B):1,(C:1,D:1):1);"), "newick")
        ic = IndependentContrasts(args)
        original_fast_copy = ic._fast_copy
        copied_trees = []
        captured = {}
        tip_traits = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        monkeypatch.setattr(ic, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(ic, "_fast_copy", copy_spy)
        self._stub_run_tail(monkeypatch, ic, tip_traits, captured)

        ic.run()

        capsys.readouterr()
        assert len(copied_trees) == 1
        assert captured["tree"] is copied_trees[0]
        assert tree.root.clades[0].clades[0].branch_length is None
        assert copied_trees[0].root.clades[0].clades[0].branch_length == 1e-8

    def test_run_missing_trait_taxa_copies_before_pruning(
        self, monkeypatch, args, capsys
    ):
        class OrderedTraitValues(dict):
            def __contains__(self, key):
                raise AssertionError("ordered prune path should not scan membership")

        from Bio import Phylo
        from io import StringIO

        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        ic = IndependentContrasts(args)
        original_fast_copy = ic._fast_copy
        copied_trees = []
        captured = {}
        tip_traits = OrderedTraitValues({"A": 1.0, "B": 2.0, "C": 3.0})

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        monkeypatch.setattr(ic, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(ic_module.Tree, "_ORDERED_MAPPING_PRUNE_MIN_SIZE", 0)
        monkeypatch.setattr(ic, "_fast_copy", copy_spy)
        self._stub_run_tail(monkeypatch, ic, tip_traits, captured)

        ic.run()

        capsys.readouterr()
        assert len(copied_trees) == 1
        assert captured["tree"] is copied_trees[0]
        assert {tip.name for tip in tree.get_terminals()} == {"A", "B", "C", "D"}
        assert {tip.name for tip in copied_trees[0].get_terminals()} == {
            "A",
            "B",
            "C",
        }

    def test_run_polytomy_copies_before_resolving(
        self, monkeypatch, args, capsys
    ):
        from Bio import Phylo
        from io import StringIO

        tree = Phylo.read(StringIO("(A:1,B:1,C:1,D:1);"), "newick")
        ic = IndependentContrasts(args)
        original_fast_copy = ic._fast_copy
        copied_trees = []
        captured = {}
        tip_traits = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}

        def copy_spy(tree_to_copy):
            copied_tree = original_fast_copy(tree_to_copy)
            copied_trees.append(copied_tree)
            return copied_tree

        monkeypatch.setattr(ic, "read_tree_file_unmodified", lambda: tree)
        monkeypatch.setattr(ic, "_fast_copy", copy_spy)
        self._stub_run_tail(monkeypatch, ic, tip_traits, captured)

        ic.run()

        capsys.readouterr()
        assert len(copied_trees) == 1
        assert captured["tree"] is copied_trees[0]
        assert len(tree.root.clades) == 4
        assert all(
            len(clade.clades) <= 2
            for clade in ic._postorder_clades_fast(copied_trees[0])
        )

    def test_run_text_output_uses_summary_label_path(self, monkeypatch, args):
        ic = IndependentContrasts(args)

        def fail_compute_pic(*_args, **_kwargs):
            raise AssertionError("text output should not build full labels")

        monkeypatch.setattr(ic, "_should_use_text_summary_path", lambda _tree: True)
        monkeypatch.setattr(ic, "_compute_pic", fail_compute_pic)

        ic.run()

    def test_print_text_batches_contrast_rows(self, args, capsys):
        ic = IndependentContrasts(args)

        ic._print_text(
            [1.0, -3.0],
            [["A", "B"], ["C", "D", "E", "F"]],
        )

        captured = capsys.readouterr()
        assert captured.out == (
            "Number of contrasts: 2\n"
            "\n"
            "Node      Contrast  Tips\n"
            "------------------------------------------------------------\n"
            "1         1.000000  A, B\n"
            "2        -3.000000  C, D, E, ... (4 total)\n"
            "\n"
            "Mean absolute contrast: 2.000000\n"
            "Variance of contrasts:  8.000000\n"
        )

    def test_small_contrast_summary_stats_avoid_numpy_wrappers(self, monkeypatch):
        contrasts = [1.0, -3.0, 5.0, -7.0]
        expected_mean_abs = float(np.mean(np.abs(contrasts)))
        expected_variance = float(np.var(contrasts, ddof=1))

        def fail_wrapper(*_args, **_kwargs):
            raise AssertionError("small summary stats should use scalar reductions")

        monkeypatch.setattr(ic_module.np, "mean", fail_wrapper)
        monkeypatch.setattr(ic_module.np, "var", fail_wrapper)

        mean_abs, variance = ic_module._contrast_summary_stats(contrasts)

        assert mean_abs == pytest.approx(expected_mean_abs)
        assert variance == pytest.approx(expected_variance)

    @pytest.mark.parametrize(
        ("contrasts", "expected_mean"),
        [([], float("nan")), ([2.0], 2.0)],
    )
    def test_scalar_contrast_summary_handles_fewer_than_two_values(
        self, contrasts, expected_mean
    ):
        mean_abs, variance = ic_module._contrast_summary_stats_scalar(contrasts)

        if math.isnan(expected_mean):
            assert math.isnan(mean_abs)
        else:
            assert mean_abs == expected_mean
        assert math.isnan(variance)

    def test_print_text_summary_stats_avoid_generic_wrappers(
        self, args, monkeypatch, capsys
    ):
        ic = IndependentContrasts(args)

        def fail_wrapper(*_args, **_kwargs):
            raise AssertionError("small text summary should use scalar reductions")

        monkeypatch.setattr(ic_module.np, "mean", fail_wrapper)
        monkeypatch.setattr(ic_module.np, "var", fail_wrapper)

        ic._print_text(
            [1.0, -3.0],
            [["A", "B"], ["C", "D"]],
        )

        captured = capsys.readouterr()
        assert "Mean absolute contrast: 2.000000" in captured.out
        assert "Variance of contrasts:  8.000000" in captured.out

    def test_run_json_output(self, capsys, args):
        args.json = True
        ic = IndependentContrasts(args)
        ic.run()
        captured = capsys.readouterr()
        import json

        payload = json.loads(captured.out)
        assert payload["n_taxa"] == 8
        assert payload["n_contrasts"] == 7
        assert len(payload["contrasts"]) == 7
        # Sum of squared contrasts matches R
        ss = sum(c["contrast"] ** 2 for c in payload["contrasts"])
        assert ss == pytest.approx(0.307253, abs=0.001)

    def test_print_json_uses_scalar_summary_for_small_outputs(
        self, monkeypatch, capsys, args
    ):
        ic = IndependentContrasts(args)
        original_asarray = ic_module.np.asarray
        asarray_calls = []

        def counting_asarray(*call_args, **kwargs):
            asarray_calls.append((call_args, kwargs))
            return original_asarray(*call_args, **kwargs)

        monkeypatch.setattr(ic_module.np, "asarray", counting_asarray)

        ic._print_json(
            [1.0, -3.0],
            [["A", "B"], ["C", "D"]],
            {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0},
        )

        payload = json.loads(capsys.readouterr().out)

        assert asarray_calls == []
        assert payload["n_taxa"] == 4
        assert payload["n_contrasts"] == 2
        assert payload["mean_absolute_contrast"] == 2.0
        assert payload["variance_of_contrasts"] == 8.0
        assert payload["contrasts"] == [
            {"node": 1, "contrast": 1.0, "tips": ["A", "B"]},
            {"node": 2, "contrast": -3.0, "tips": ["C", "D"]},
        ]

    def test_print_json_ignores_broken_pipe(self, monkeypatch, args):
        ic = IndependentContrasts(args)

        def fail_dump(*_args, **_kwargs):
            raise BrokenPipeError

        monkeypatch.setattr(ic_module, "_json_dumps", fail_dump)

        ic._print_json([1.0], [["A", "B"]], {"A": 1.0, "B": 2.0})

    def test_large_contrast_summary_uses_numpy_array_path(
        self, monkeypatch, args
    ):
        ic = IndependentContrasts(args)
        original_asarray = ic_module.np.asarray
        asarray_calls = []

        def counting_asarray(*call_args, **kwargs):
            asarray_calls.append((call_args, kwargs))
            return original_asarray(*call_args, **kwargs)

        monkeypatch.setattr(ic_module.np, "asarray", counting_asarray)

        mean_abs, variance = ic_module._contrast_summary_stats([1.0] * 300)

        assert len(asarray_calls) == 1
        assert mean_abs == 1.0
        assert variance == 0.0
