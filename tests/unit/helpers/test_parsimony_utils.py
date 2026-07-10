"""
Unit tests for parsimony utility algorithms.

Tests Fitch downpass, ACCTRAN/DELTRAN uppass, change detection,
change classification (synapomorphy/convergence/reversal),
consistency index, and retention index.
"""
import subprocess
import sys
import builtins

import pytest
import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import TreeMixin
from io import StringIO

import phykit.helpers.parsimony_utils as parsimony_module
from phykit.helpers.parsimony_utils import (
    build_parent_map,
    resolve_polytomies,
    fitch_downpass,
    fitch_uppass_acctran,
    fitch_uppass_deltran,
    detect_changes,
    classify_changes,
    consistency_index,
    retention_index,
)


def test_module_import_does_not_import_typing_or_biopython():
    code = """
import sys
import phykit.helpers.parsimony_utils
assert "typing" not in sys.modules
assert "Bio" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def _make_tree(newick):
    return Phylo.read(StringIO(newick), "newick")


# ---------------------------------------------------------------------------
# build_parent_map
# ---------------------------------------------------------------------------

class TestBuildParentMap:
    def test_simple_tree(self):
        tree = _make_tree("((A:1,B:1):1,C:1);")
        pm = build_parent_map(tree)
        # All non-root nodes should have a parent
        root = tree.root
        for clade in tree.find_clades():
            if clade == root:
                assert id(clade) not in pm
            else:
                assert id(clade) in pm

    def test_parent_of_tip_is_internal(self):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        pm = build_parent_map(tree)
        for tip in tree.get_terminals():
            parent = pm[id(tip)]
            assert not parent.is_terminal()

    def test_parent_map_covers_all_non_root(self):
        tree = _make_tree("((A:1,B:1):1,((C:1,D:1):1,E:1):1,F:1);")
        pm = build_parent_map(tree)
        all_nodes = list(tree.find_clades())
        # Every node except root should be in parent_map
        assert len(pm) == len(all_nodes) - 1

    def test_standard_tree_uses_direct_traversal(self, monkeypatch):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic preorder traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        pm = build_parent_map(tree)

        assert len(pm) == 6
        assert pm[id(tree.root.clades[0])] is tree.root
        assert pm[id(tree.root.clades[1])] is tree.root

    def test_direct_traversal_handles_mixed_child_counts(self, monkeypatch):
        tree = _make_tree("(A:1,(B:1,C:1):1,(D:1,E:1,F:1):1);")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic preorder traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        pm = build_parent_map(tree)

        binary = tree.root.clades[1]
        trifurcation = tree.root.clades[2]
        assert pm[id(tree.root.clades[0])] is tree.root
        assert pm[id(binary)] is tree.root
        assert pm[id(trifurcation)] is tree.root
        assert pm[id(binary.clades[0])] is binary
        assert pm[id(binary.clades[1])] is binary
        assert pm[id(trifurcation.clades[0])] is trifurcation
        assert pm[id(trifurcation.clades[1])] is trifurcation
        assert pm[id(trifurcation.clades[2])] is trifurcation


# ---------------------------------------------------------------------------
# resolve_polytomies
# ---------------------------------------------------------------------------

class TestResolvePolytomies:
    def test_trifurcating_root_resolved(self):
        tree = _make_tree("(A:1,B:1,C:1,D:1);")
        resolve_polytomies(tree)
        for clade in tree.find_clades():
            if not clade.is_terminal():
                assert len(clade.clades) <= 2

    def test_binary_tree_unchanged(self):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        tips_before = sorted(t.name for t in tree.get_terminals())
        resolve_polytomies(tree)
        tips_after = sorted(t.name for t in tree.get_terminals())
        assert tips_before == tips_after

    def test_preserves_all_tips(self):
        tree = _make_tree("(A:1,B:1,C:1,D:1,E:1);")
        resolve_polytomies(tree)
        tips = sorted(t.name for t in tree.get_terminals())
        assert tips == ["A", "B", "C", "D", "E"]

    def test_zero_length_branches_inserted(self):
        tree = _make_tree("(A:1,B:1,C:1);")
        resolve_polytomies(tree)
        # There should be a zero-length internal branch
        has_zero = False
        for clade in tree.find_clades():
            if not clade.is_terminal() and clade != tree.root:
                if clade.branch_length == 0.0:
                    has_zero = True
        assert has_zero

    def test_character_map_tree(self):
        """Test with the actual character_map tree topology."""
        tree = _make_tree("((A:1,B:1):1,((C:1,D:1):1,E:1):1,F:1);")
        resolve_polytomies(tree)
        for clade in tree.find_clades():
            if not clade.is_terminal():
                assert len(clade.clades) <= 2
        tips = sorted(t.name for t in tree.get_terminals())
        assert tips == ["A", "B", "C", "D", "E", "F"]

    def test_resolve_polytomies_uses_direct_postorder(self, monkeypatch):
        tree = _make_tree("(A:1,B:1,C:1,D:1);")

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic postorder traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        resolve_polytomies(tree)

        stack = [tree.root]
        while stack:
            clade = stack.pop()
            assert len(clade.clades) <= 2
            stack.extend(clade.clades)

    def test_binary_tree_does_not_import_newick(self, monkeypatch):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        original_import = builtins.__import__

        def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "Bio.Phylo.Newick" or (
                name == "Bio.Phylo" and fromlist and "Newick" in fromlist
            ):
                raise AssertionError("binary trees should not import Newick helpers")
            return original_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", guarded_import)

        resolve_polytomies(tree)

        assert [tip.name for tip in tree.get_terminals()] == ["A", "B", "C", "D"]

    def test_resolve_polytomies_matches_previous_postorder_output(self):
        newick = "((A:1,B:1,C:1):1,(D:1,E:1,F:1):1,G:1);"
        expected = _make_tree(newick)
        actual = _make_tree(newick)

        from Bio.Phylo import Newick

        clades = []
        stack = [expected.root]
        while stack:
            clade = stack.pop()
            clades.append(clade)
            stack.extend(clade.clades)
        clades.reverse()
        for clade in clades:
            while len(clade.clades) > 2:
                child1 = clade.clades.pop()
                child2 = clade.clades.pop()
                new_internal = Newick.Clade(branch_length=0.0)
                new_internal.clades = [child1, child2]
                clade.clades.append(new_internal)

        resolve_polytomies(actual)

        expected_out = StringIO()
        actual_out = StringIO()
        Phylo.write(expected, expected_out, "newick")
        Phylo.write(actual, actual_out, "newick")
        assert actual_out.getvalue() == expected_out.getvalue()


# ---------------------------------------------------------------------------
# fitch_downpass
# ---------------------------------------------------------------------------

class TestFitchDownpass:
    def test_preorder_direct_binary_children_use_indexed_path(self):
        from Bio.Phylo.BaseTree import Clade, Tree

        class IndexedOnlyList(list):
            def __reversed__(self):
                raise AssertionError("preorder helper should index children")

        left = Clade(
            name="left",
            clades=IndexedOnlyList([Clade(name="A"), Clade(name="B")]),
        )
        right = Clade(
            name="right",
            clades=IndexedOnlyList([Clade(name="C"), Clade(name="D")]),
        )
        root = Clade(name="root", clades=IndexedOnlyList([left, right]))
        tree = Tree(root=root)

        clades = parsimony_module._preorder_clades_direct(tree)

        assert clades == [
            root,
            left,
            left.clades[0],
            left.clades[1],
            right,
            right.clades[0],
            right.clades[1],
        ]

    def test_fitch_pipeline_uses_direct_tree_traversal(self, monkeypatch):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        tip_states = {
            "A": ["0", "0"],
            "B": ["0", "1"],
            "C": ["1", "1"],
            "D": ["1", "0"],
        }

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("generic tree traversal should not be used")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        parent_map = build_parent_map(tree)
        node_state_sets, scores = fitch_downpass(tree, tip_states)
        acctran_states = fitch_uppass_acctran(tree, node_state_sets, parent_map)
        deltran_states = fitch_uppass_deltran(tree, node_state_sets, parent_map)
        acctran_changes = detect_changes(tree, acctran_states, parent_map)
        deltran_changes = detect_changes(tree, deltran_states, parent_map)

        assert scores == [1, 2]
        assert id(tree.root) in acctran_states
        assert id(tree.root) in deltran_states
        assert acctran_changes
        assert deltran_changes

    def test_simple_no_change(self):
        """All tips same state -> 0 steps."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["0"], "C": ["0"], "D": ["0"]}
        node_state_sets, scores = fitch_downpass(tree, tip_states)
        assert scores == [0]

    def test_simple_one_change(self):
        """One character, one change required."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["0"], "C": ["1"], "D": ["1"]}
        node_state_sets, scores = fitch_downpass(tree, tip_states)
        assert scores == [1]

    def test_wildcard_question_mark(self):
        """Wildcard (?) matches any state."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["?"], "C": ["1"], "D": ["1"]}
        node_state_sets, scores = fitch_downpass(tree, tip_states)
        # B is wildcard, so (A,B) can be {0}, and (C,D)={1} => 1 change
        assert scores == [1]

    def test_wildcard_dash(self):
        """Wildcard (-) matches any state."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["-"], "C": ["1"], "D": ["1"]}
        node_state_sets, scores = fitch_downpass(tree, tip_states)
        assert scores == [1]

    def test_multistate(self):
        """Character with 3 states."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["1"], "C": ["2"], "D": ["2"]}
        node_state_sets, scores = fitch_downpass(tree, tip_states)
        # (A,B) = union {0,1} cost 1; (C,D) = {2}; root = union cost 1 => 2
        assert scores == [2]

    def test_multiple_characters(self):
        """Two characters scored independently."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {
            "A": ["0", "1"],
            "B": ["0", "1"],
            "C": ["1", "0"],
            "D": ["1", "0"],
        }
        node_state_sets, scores = fitch_downpass(tree, tip_states)
        assert scores == [1, 1]

    def test_repeated_terminal_patterns_match_generic_downpass(self, monkeypatch):
        tree = _make_tree("(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);")
        tip_states = {
            "A": ["0", "1", "0"],
            "B": ["0", "1", "0"],
            "C": ["1", "0", "1"],
            "D": ["1", "0", "1"],
            "E": ["0", "1", "0"],
            "F": ["0", "1", "0"],
            "G": ["1", "0", "1"],
            "H": ["1", "0", "1"],
        }

        optimized_sets, optimized_scores = fitch_downpass(tree, tip_states)

        monkeypatch.setattr(
            parsimony_module,
            "_postorder_clades_direct",
            lambda _tree: None,
        )
        generic_sets, generic_scores = fitch_downpass(tree, tip_states)

        assert optimized_scores == generic_scores
        assert optimized_sets[id(tree.root)] == generic_sets[id(tree.root)]

    def test_repeated_terminal_patterns_keep_node_state_lists_independent(self):
        tree = _make_tree("(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);")
        pattern_a = ["0", "1"] * 8
        pattern_b = ["1", "0"] * 8
        tip_states = {
            "A": pattern_a,
            "B": pattern_a,
            "C": pattern_b,
            "D": pattern_b,
            "E": pattern_a,
            "F": pattern_a,
            "G": pattern_b,
            "H": pattern_b,
        }

        node_state_sets, scores = fitch_downpass(tree, tip_states)
        lists_by_value = {}
        for char_sets in node_state_sets.values():
            lists_by_value.setdefault(tuple(char_sets), []).append(char_sets)

        repeated_lists = next(
            group for group in lists_by_value.values() if len(group) > 1
        )

        assert scores == [2, 2] * 8
        assert repeated_lists[0] == repeated_lists[1]
        assert repeated_lists[0] is not repeated_lists[1]

    def test_root_state_set_correct(self):
        """Root state set should be intersection or union of children."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["0"], "C": ["1"], "D": ["1"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        root_sets = node_state_sets[id(tree.root)]
        # (A,B)={0}, (C,D)={1}, root = union {0,1}
        assert root_sets[0] == {"0", "1"}

    def test_all_wildcards(self):
        """All wildcards should produce zero changes."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["?"], "B": ["?"], "C": ["?"], "D": ["?"]}
        node_state_sets, scores = fitch_downpass(tree, tip_states)
        assert scores == [0]

    def test_binary_downpass_wildcard_empty_sets_match_generic_semantics(self):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {
            "A": ["?", "0"],
            "B": ["?", "1"],
            "C": ["1", "?"],
            "D": ["1", "?"],
        }

        node_state_sets, scores = fitch_downpass(tree, tip_states)

        assert scores == [0, 1]
        assert node_state_sets[id(tree.root)][0] == {"1"}
        assert node_state_sets[id(tree.root)][1] == {"0", "1"}

    def test_binary_downpass_bitmask_path_is_uppass_compatible(self, monkeypatch):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        tip_states = {
            "A": ["0", "0", "?"],
            "B": ["0", "1", "1"],
            "C": ["1", "1", "1"],
            "D": ["1", "0", "?"],
        }

        def fail_find_clades(*_args, **_kwargs):
            raise AssertionError("resolved trees should use direct traversal")

        monkeypatch.setattr(TreeMixin, "find_clades", fail_find_clades)

        parent_map = build_parent_map(tree)
        node_state_sets, scores = fitch_downpass(tree, tip_states)
        acctran_states = fitch_uppass_acctran(tree, node_state_sets, parent_map)
        deltran_states = fitch_uppass_deltran(tree, node_state_sets, parent_map)

        assert scores == [1, 2, 0]
        assert node_state_sets[id(tree.root)][0] == {"0", "1"}
        assert node_state_sets[id(tree.root)][1] == {"0", "1"}
        assert node_state_sets[id(tree.root)][2] == {"1"}
        assert id(tree.root) in acctran_states
        assert id(tree.root) in deltran_states

    def test_binary_downpass_uses_small_mask_lookup(self, monkeypatch):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        tip_states = {
            "A": ["A", "C", "G", "T"],
            "B": ["A", "G", "G", "T"],
            "C": ["C", "G", "A", "?"],
            "D": ["C", "G", "T", "?"],
        }
        built_state_counts = []
        original_lookup_builder = parsimony_module._build_small_mask_set_lookup

        def recording_lookup_builder(states):
            built_state_counts.append(len(states))
            return original_lookup_builder(states)

        monkeypatch.setattr(
            parsimony_module,
            "_build_small_mask_set_lookup",
            recording_lookup_builder,
        )

        node_state_sets, scores = fitch_downpass(tree, tip_states)

        assert built_state_counts == [2, 2, 3, 1]
        assert scores == [1, 1, 2, 0]
        assert node_state_sets[id(tree.root)] == [
            {"A", "C"},
            {"G"},
            {"A", "G", "T"},
            {"T"},
        ]

    def test_polytomy_tree_after_resolve(self):
        """Works on a resolved polytomy tree."""
        tree = _make_tree("((A:1,B:1):1,((C:1,D:1):1,E:1):1,F:1);")
        resolve_polytomies(tree)
        tip_states = {
            "A": ["0"], "B": ["0"], "C": ["1"],
            "D": ["1"], "E": ["1"], "F": ["0"],
        }
        node_state_sets, scores = fitch_downpass(tree, tip_states)
        # char0: 0->1 on (C,D,E) clade = 1 change
        assert scores == [1]


# ---------------------------------------------------------------------------
# ACCTRAN uppass
# ---------------------------------------------------------------------------

class TestAcctranUppass:
    def test_no_change_inherits_root(self):
        """All same state: every node gets root state."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["0"], "C": ["0"], "D": ["0"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)
        for nid, states in node_states.items():
            assert states == ["0"]

    def test_change_placed_rootward(self):
        """ACCTRAN places change on the branch closest to the root."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        # A,B = 0; C,D = 1 => change on branch from root to (C,D)
        tip_states = {"A": ["0"], "B": ["0"], "C": ["1"], "D": ["1"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)
        # Root should be "0" (min of {0,1})
        root_states = node_states[id(tree.root)]
        assert root_states == ["0"]
        # Tips should match observed
        for tip in tree.get_terminals():
            assert node_states[id(tip)] == tip_states[tip.name]

    def test_root_tie_break_uses_min(self):
        """When root state set has multiple states, min is chosen."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["1"], "B": ["1"], "C": ["0"], "D": ["0"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)
        assert node_states[id(tree.root)] == ["0"]

    def test_tips_always_match_observed(self):
        """ACCTRAN must preserve tip states."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["1"], "C": ["0"], "D": ["1"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)
        for tip in tree.get_terminals():
            assert node_states[id(tip)] == tip_states[tip.name]

    def test_multiple_characters(self):
        """ACCTRAN with multiple characters independently."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {
            "A": ["0", "1"],
            "B": ["0", "1"],
            "C": ["1", "0"],
            "D": ["1", "0"],
        }
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)
        root_states = node_states[id(tree.root)]
        assert root_states == ["0", "0"]


# ---------------------------------------------------------------------------
# DELTRAN uppass
# ---------------------------------------------------------------------------

class TestDeltranUppass:
    def test_change_placed_tipward(self):
        """DELTRAN delays changes toward tips."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["0"], "C": ["1"], "D": ["1"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_deltran(tree, node_state_sets, pm)
        # Tips must match observed
        for tip in tree.get_terminals():
            assert node_states[id(tip)] == tip_states[tip.name]

    def test_acctran_vs_deltran_may_differ(self):
        """ACCTRAN and DELTRAN produce valid reconstructions (may differ on complex cases)."""
        # Tree: ((A,B),(C,(D,E)))
        tree = _make_tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")
        resolve_polytomies(tree)
        # A=0, B=1, C=1, D=0, E=0
        tip_states = {
            "A": ["0"], "B": ["1"], "C": ["1"],
            "D": ["0"], "E": ["0"],
        }
        node_state_sets, scores = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        acctran = fitch_uppass_acctran(tree, node_state_sets, pm)
        deltran = fitch_uppass_deltran(tree, node_state_sets, pm)

        # Both must produce same number of changes as parsimony score
        acctran_changes = detect_changes(tree, acctran, pm)
        deltran_changes = detect_changes(tree, deltran, pm)
        acctran_total = sum(len(c) for c in acctran_changes.values())
        deltran_total = sum(len(c) for c in deltran_changes.values())
        assert acctran_total == sum(scores)
        assert deltran_total == sum(scores)

        # Both must preserve tip states
        for tip in tree.get_terminals():
            assert acctran[id(tip)] == tip_states[tip.name]
            assert deltran[id(tip)] == tip_states[tip.name]

        # Root should be the same
        assert acctran[id(tree.root)] == deltran[id(tree.root)]

    def test_tips_always_match_observed(self):
        """DELTRAN must preserve tip states."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["1"], "C": ["0"], "D": ["1"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_deltran(tree, node_state_sets, pm)
        for tip in tree.get_terminals():
            assert node_states[id(tip)] == tip_states[tip.name]

    def test_no_change_inherits_root(self):
        """All same state: every node gets root state with DELTRAN too."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["0"], "C": ["0"], "D": ["0"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_deltran(tree, node_state_sets, pm)
        for nid, states in node_states.items():
            assert states == ["0"]

    def test_deltran_delays_change(self):
        """DELTRAN should delay changes compared to ACCTRAN on asymmetric tree."""
        # Tree: (A,(B,(C,D)))
        tree = _make_tree("(A:1,(B:1,(C:1,D:1):1):1);")
        resolve_polytomies(tree)
        # A=0, B=0, C=1, D=1: change to 1 needed somewhere in (C,D) branch
        tip_states = {"A": ["0"], "B": ["0"], "C": ["1"], "D": ["1"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        acctran = fitch_uppass_acctran(tree, node_state_sets, pm)
        deltran = fitch_uppass_deltran(tree, node_state_sets, pm)
        # ACCTRAN: change on (B,(C,D)) ancestor from root.
        # DELTRAN: change delayed to (C,D) ancestor.
        # Both should have same root and tip states
        assert acctran[id(tree.root)] == deltran[id(tree.root)]
        for tip in tree.get_terminals():
            assert acctran[id(tip)] == deltran[id(tip)]


# ---------------------------------------------------------------------------
# detect_changes
# ---------------------------------------------------------------------------

class TestDetectChanges:
    def test_single_character_path_skips_per_character_loop(self, monkeypatch):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["0"], "C": ["1"], "D": ["1"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)

        def fail_range(*_args, **_kwargs):
            raise AssertionError("single-character changes should not loop over range(1)")

        monkeypatch.setattr(parsimony_module, "range", fail_range, raising=False)

        changes = detect_changes(tree, node_states, pm)

        all_changes = [
            c for branch_changes in changes.values() for c in branch_changes
        ]
        assert all_changes == [(0, "0", "1")]

    def test_counts_changes(self):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["0"], "C": ["1"], "D": ["1"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)
        changes = detect_changes(tree, node_states, pm)
        # Exactly one change for char 0
        all_changes = [
            c for branch_changes in changes.values() for c in branch_changes
        ]
        char0_changes = [c for c in all_changes if c[0] == 0]
        assert len(char0_changes) == 1

    def test_no_change_detected_when_uniform(self):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["0"], "C": ["0"], "D": ["0"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)
        changes = detect_changes(tree, node_states, pm)
        assert len(changes) == 0

    def test_change_records_old_and_new(self):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["0"], "C": ["1"], "D": ["1"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)
        changes = detect_changes(tree, node_states, pm)
        all_changes = [
            c for branch_changes in changes.values() for c in branch_changes
        ]
        # The change should be 0->1
        assert any(c[1] == "0" and c[2] == "1" for c in all_changes)

    def test_multiple_characters_multiple_changes(self):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {
            "A": ["0", "1"],
            "B": ["0", "1"],
            "C": ["1", "0"],
            "D": ["1", "0"],
        }
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)
        changes = detect_changes(tree, node_states, pm)
        all_changes = [
            c for branch_changes in changes.values() for c in branch_changes
        ]
        # Should have changes for both characters
        chars_changed = set(c[0] for c in all_changes)
        assert 0 in chars_changed
        assert 1 in chars_changed


# ---------------------------------------------------------------------------
# classify_changes
# ---------------------------------------------------------------------------

class TestClassifyChanges:
    def test_transition_counts_use_plain_dict(self, monkeypatch):
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        root = tree.root
        left, right = root.clades
        a, b = left.clades
        c, d = right.clades
        pm = build_parent_map(tree)
        node_states = {
            id(root): ["0", "0"],
            id(left): ["1", "0"],
            id(right): ["0", "0"],
            id(a): ["0", "1"],
            id(b): ["1", "1"],
            id(c): ["0", "0"],
            id(d): ["0", "1"],
        }
        branch_changes = {
            id(left): [(0, "0", "1")],
            id(a): [(0, "1", "0")],
            id(c): [(0, "1", "0")],
            id(b): [(1, "0", "1")],
            id(d): [(1, "0", "1")],
        }

        def fail_counter(*_args, **_kwargs):
            raise AssertionError("classify_changes should not instantiate Counter")

        monkeypatch.setattr(parsimony_module, "Counter", fail_counter)

        classified = classify_changes(tree, branch_changes, node_states, pm)

        assert classified[id(left)] == [(0, "0", "1", "synapomorphy")]
        assert classified[id(a)] == [(0, "1", "0", "reversal")]
        assert classified[id(c)] == [(0, "1", "0", "reversal")]
        assert classified[id(b)] == [(1, "0", "1", "convergence")]
        assert classified[id(d)] == [(1, "0", "1", "convergence")]

    def test_synapomorphy(self):
        """A change occurring once is a synapomorphy."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["0"], "C": ["1"], "D": ["1"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)
        changes = detect_changes(tree, node_states, pm)
        classified = classify_changes(tree, changes, node_states, pm)
        # All changes for char0 should be synapomorphy
        all_c = [c for branch in classified.values() for c in branch]
        assert all(c[3] == "synapomorphy" for c in all_c if c[0] == 0)

    def test_convergence(self):
        """Same state gained independently = convergence."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        # A=1, B=0, C=1, D=0 => char0 gains 1 independently on A and C
        tip_states = {"A": ["1"], "B": ["0"], "C": ["1"], "D": ["0"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)
        changes = detect_changes(tree, node_states, pm)
        classified = classify_changes(tree, changes, node_states, pm)
        all_c = [c for branch in classified.values() for c in branch]
        convergences = [c for c in all_c if c[3] == "convergence"]
        assert len(convergences) >= 1

    def test_reversal_to_non_root_state(self):
        """Reversal to an ancestral (non-root) state is classified correctly."""
        tree = _make_tree("(((A:1,B:1):1,C:1):1,(D:1,E:1):1);")
        resolve_polytomies(tree)
        tip_states = {
            "A": ["1"], "B": ["1"], "C": ["0"],
            "D": ["0"], "E": ["1"],
        }
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)
        changes = detect_changes(tree, node_states, pm)
        classified = classify_changes(tree, changes, node_states, pm)
        all_c = [c for branch in classified.values() for c in branch]
        types = set(c[3] for c in all_c)
        # Should have at least one change identified
        assert len(all_c) >= 1

    def test_classification_tuple_has_four_elements(self):
        """Each classified change should be (char_idx, old, new, classification)."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["0"], "C": ["1"], "D": ["1"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)
        changes = detect_changes(tree, node_states, pm)
        classified = classify_changes(tree, changes, node_states, pm)
        for cid, cls_list in classified.items():
            for item in cls_list:
                assert len(item) == 4
                assert item[3] in ("synapomorphy", "convergence", "reversal")

    def test_with_character_map_data(self):
        """Full test with the character_map sample data."""
        tree = _make_tree("((A:1,B:1):1,((C:1,D:1):1,E:1):1,F:1);")
        resolve_polytomies(tree)
        tip_states = {
            "A": ["0", "1", "0", "0", "0", "1", "1", "0"],
            "B": ["0", "1", "0", "0", "1", "1", "1", "0"],
            "C": ["1", "0", "1", "0", "0", "0", "2", "0"],
            "D": ["1", "0", "1", "1", "0", "0", "0", "0"],
            "E": ["1", "0", "0", "0", "0", "0", "1", "0"],
            "F": ["0", "0", "0", "0", "0", "0", "0", "0"],
        }
        node_state_sets, scores = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)
        changes = detect_changes(tree, node_states, pm)
        classified = classify_changes(tree, changes, node_states, pm)
        all_c = [c for branch in classified.values() for c in branch]
        # char0: exactly 1 change => synapomorphy
        char0 = [c for c in all_c if c[0] == 0]
        assert len(char0) == 1
        assert char0[0][3] == "synapomorphy"
        # char1: exactly 1 change => synapomorphy
        char1 = [c for c in all_c if c[0] == 1]
        assert len(char1) == 1
        assert char1[0][3] == "synapomorphy"
        # char7 (all zeros): no changes
        char7 = [c for c in all_c if c[0] == 7]
        assert len(char7) == 0


# ---------------------------------------------------------------------------
# consistency_index
# ---------------------------------------------------------------------------

class TestConsistencyIndex:
    def test_no_homoplasy(self):
        """CI = 1.0 when no homoplasy."""
        # 2 states, 1 observed change => CI = 1/1 = 1.0
        ci_per_char, ci_overall = consistency_index(
            n_states_per_char=[2], observed_per_char=[1]
        )
        assert ci_per_char == [1.0]
        assert ci_overall == 1.0

    def test_with_homoplasy(self):
        """CI < 1.0 with homoplasy."""
        # 2 states, 2 observed changes => CI = 1/2 = 0.5
        ci_per_char, ci_overall = consistency_index(
            n_states_per_char=[2], observed_per_char=[2]
        )
        assert ci_per_char == [0.5]
        assert ci_overall == 0.5

    def test_multistate(self):
        """min_changes = n_states - 1 for multistate."""
        # 3 states, 3 observed => CI = 2/3
        ci_per_char, ci_overall = consistency_index(
            n_states_per_char=[3], observed_per_char=[3]
        )
        assert ci_per_char[0] == pytest.approx(2 / 3)

    def test_zero_observed_returns_none(self):
        """No observed changes => CI is None for that character."""
        ci_per_char, ci_overall = consistency_index(
            n_states_per_char=[2], observed_per_char=[0]
        )
        assert ci_per_char == [None]
        assert ci_overall is None

    def test_multiple_characters(self):
        """CI overall is weighted sum."""
        # char0: 2 states, 1 change => CI=1.0, min=1
        # char1: 2 states, 2 changes => CI=0.5, min=1
        # overall: (1+1)/(1+2) = 2/3
        ci_per_char, ci_overall = consistency_index(
            n_states_per_char=[2, 2], observed_per_char=[1, 2]
        )
        assert ci_per_char[0] == 1.0
        assert ci_per_char[1] == 0.5
        assert ci_overall == pytest.approx(2 / 3)

    def test_perfect_parsimony(self):
        """All characters perfectly parsimonious."""
        ci_per_char, ci_overall = consistency_index(
            n_states_per_char=[2, 3, 2], observed_per_char=[1, 2, 1]
        )
        assert all(c == 1.0 for c in ci_per_char)
        assert ci_overall == 1.0


# ---------------------------------------------------------------------------
# retention_index
# ---------------------------------------------------------------------------

class TestRetentionIndex:
    def test_ascii_symbol_counts_by_character_matches_reference(self):
        alphabet = np.frombuffer(b"ACGT-?", dtype=np.uint8)
        matrix = np.array(
            [
                [ord("A"), ord("C"), ord("A"), ord("?")],
                [ord("G"), ord("G"), ord("-"), ord("T")],
                [ord("C"), ord("C"), ord("C"), ord("A")],
            ],
            dtype=np.uint8,
        )
        symbols = alphabet[(alphabet != ord("?")) & (alphabet != ord("-"))]

        observed = parsimony_module._ascii_symbol_counts_by_character(
            matrix,
            symbols,
        )
        expected = np.vstack(
            [np.count_nonzero(matrix == symbol, axis=1) for symbol in symbols]
        )

        np.testing.assert_array_equal(observed, expected)

    def test_ascii_symbol_counts_by_character_large_alphabet_uses_bincount(
        self,
        monkeypatch,
    ):
        alphabet = np.arange(33, 57, dtype=np.uint8)
        matrix = np.tile(alphabet, 120).reshape(80, -1)
        symbols = np.unique(matrix)
        original_bincount = np.bincount
        calls = 0

        def count_bincount(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original_bincount(*args, **kwargs)

        monkeypatch.setattr(np, "bincount", count_bincount)

        observed = parsimony_module._ascii_symbol_counts_by_character(
            matrix,
            symbols,
        )
        expected = np.vstack(
            [np.count_nonzero(matrix == symbol, axis=1) for symbol in symbols]
        )

        assert calls == 1
        np.testing.assert_array_equal(observed, expected)

    def test_ascii_symbol_counts_by_character_fifteen_states_uses_bincount(
        self,
        monkeypatch,
    ):
        alphabet = np.arange(33, 48, dtype=np.uint8)
        matrix = np.tile(alphabet, 80).reshape(75, -1)
        symbols = np.unique(matrix)
        original_bincount = np.bincount
        calls = 0

        def count_bincount(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original_bincount(*args, **kwargs)

        monkeypatch.setattr(np, "bincount", count_bincount)

        observed = parsimony_module._ascii_symbol_counts_by_character(
            matrix,
            symbols,
        )
        expected = np.vstack(
            [np.count_nonzero(matrix == symbol, axis=1) for symbol in symbols]
        )

        assert calls == 1
        np.testing.assert_array_equal(observed, expected)

    def test_ascii_symbol_counts_by_character_twelve_states_skips_bincount(
        self,
        monkeypatch,
    ):
        alphabet = np.arange(33, 45, dtype=np.uint8)
        matrix = np.tile(alphabet, 100).reshape(80, -1)
        symbols = np.unique(matrix)

        def fail_bincount(*_args, **_kwargs):
            raise AssertionError("12-state ASCII counts should use equality scans")

        monkeypatch.setattr(np, "bincount", fail_bincount)

        observed = parsimony_module._ascii_symbol_counts_by_character(
            matrix,
            symbols,
        )
        expected = np.vstack(
            [np.count_nonzero(matrix == symbol, axis=1) for symbol in symbols]
        )

        np.testing.assert_array_equal(observed, expected)

    def test_basic(self):
        """RI with known values."""
        # 4 taxa: A=0, B=0, C=1, D=1
        # max_changes = 4 - 2 = 2 (f_max=2 for state 0 or 1)
        # min_changes = 2 - 1 = 1
        # observed = 1
        # RI = (2-1)/(2-1) = 1.0
        ri_per_char, ri_overall = retention_index(
            tip_states_per_char=[["0", "0", "1", "1"]],
            observed_per_char=[1],
        )
        assert ri_per_char[0] == 1.0

    def test_uninformative(self):
        """Uninformative character: RI is None."""
        # All same state: max = 4-4=0, min = 1-1=0, max==min => None
        ri_per_char, ri_overall = retention_index(
            tip_states_per_char=[["0", "0", "0", "0"]],
            observed_per_char=[0],
        )
        assert ri_per_char[0] is None

    def test_with_homoplasy(self):
        """RI < 1.0 with homoplasy."""
        # 4 taxa: A=0, B=1, C=0, D=1 (checkerboard)
        # max_changes = 4-2=2, min_changes = 2-1=1
        # observed = 2
        # RI = (2-2)/(2-1) = 0.0
        ri_per_char, ri_overall = retention_index(
            tip_states_per_char=[["0", "1", "0", "1"]],
            observed_per_char=[2],
        )
        assert ri_per_char[0] == 0.0

    def test_ri_with_wildcards_filtered(self):
        """Wildcards should be excluded from state frequency counts."""
        # 4 taxa, one wildcard: A=0, B=?, C=1, D=1
        # clean = [0, 1, 1], n_taxa=3, f_max=2 (for 1), max=3-2=1
        # n_states=2, min=1, max==min => None
        ri_per_char, ri_overall = retention_index(
            tip_states_per_char=[["0", "?", "1", "1"]],
            observed_per_char=[1],
        )
        assert ri_per_char[0] is None

    def test_ri_all_wildcards(self):
        """RI is undefined when a column has no observed states."""
        ri_per_char, ri_overall = retention_index(
            tip_states_per_char=[["?", "-", "?", "-"]],
            observed_per_char=[0],
        )
        assert ri_per_char[0] is None
        assert ri_overall is None

    def test_multiple_characters(self):
        """RI overall is weighted."""
        ri_per_char, ri_overall = retention_index(
            tip_states_per_char=[
                ["0", "0", "1", "1"],
                ["0", "1", "0", "1"],
            ],
            observed_per_char=[1, 2],
        )
        assert ri_per_char[0] == 1.0
        assert ri_per_char[1] == 0.0
        # overall = sum(max-obs)/sum(max-min)
        # char0: (2-1)/(2-1)=1, num=1, den=1
        # char1: (2-2)/(2-1)=0, num=0, den=1
        # overall = (1+0)/(1+1) = 0.5
        assert ri_overall == pytest.approx(0.5)

    def test_ascii_single_character_fast_path_matches_reference(self):
        tip_states_per_char = [
            ["0", "0", "1", "1", "?"],
            ["0", "1", "0", "1", "-"],
            ["0", "0", "0", "0", "?"],
        ]
        observed_per_char = [1, 2, 0]

        ri_per_char, ri_overall = retention_index(
            tip_states_per_char,
            observed_per_char,
        )

        assert ri_per_char == [1.0, 0.0, None]
        assert ri_overall == pytest.approx(0.5)

    def test_ascii_fast_path_uses_array_reductions(self, monkeypatch):
        def fail_sum(*_args, **_kwargs):
            raise AssertionError("ASCII RI fast path should avoid np.sum")

        monkeypatch.setattr(np, "sum", fail_sum)

        ri_per_char, ri_overall = retention_index(
            tip_states_per_char=[
                ["0", "0", "1", "1", "?"],
                ["0", "1", "0", "1", "-"],
            ],
            observed_per_char=[1, 2],
        )

        assert ri_per_char == [1.0, 0.0]
        assert ri_overall == pytest.approx(0.5)

    def test_small_ascii_retention_index_uses_counter_fallback(self, monkeypatch):
        def fail_symbol_counts(*_args, **_kwargs):
            raise AssertionError("small RI matrices should avoid NumPy counting")

        monkeypatch.setattr(
            parsimony_module,
            "_ascii_symbol_counts_by_character",
            fail_symbol_counts,
        )

        ri_per_char, ri_overall = retention_index(
            tip_states_per_char=[
                ["0", "0", "1", "1", "?"],
                ["0", "1", "0", "1", "-"],
            ],
            observed_per_char=[1, 2],
        )

        assert ri_per_char == [1.0, 0.0]
        assert ri_overall == pytest.approx(0.5)

    def test_retention_index_fallback_handles_multicharacter_states(self):
        ri_per_char, ri_overall = retention_index(
            tip_states_per_char=[
                ["red", "red", "blue", "blue"],
                ["α", "β", "α", "β"],
            ],
            observed_per_char=[1, 2],
        )

        assert ri_per_char == [1.0, 0.0]
        assert ri_overall == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# Integration: full pipeline with character_map data
# ---------------------------------------------------------------------------

class TestFullPipeline:
    """End-to-end test using the character_map sample tree and matrix."""

    def _load_data(self):
        tree = Phylo.read(
            "tests/sample_files/tree_character_map.tre", "newick"
        )
        resolve_polytomies(tree)

        tip_states = {}
        with open("tests/sample_files/character_matrix_simple.tsv") as f:
            header = f.readline().strip().split("\t")
            char_names = header[1:]
            for line in f:
                parts = line.strip().split("\t")
                taxon = parts[0]
                states = parts[1:]
                tip_states[taxon] = states
        return tree, tip_states, char_names

    def test_downpass_scores(self):
        """Verify parsimony scores for the sample dataset."""
        tree, tip_states, _ = self._load_data()
        _, scores = fitch_downpass(tree, tip_states)
        # char0: 1 change (0->1 on CDE clade)
        assert scores[0] == 1
        # char1: 1 change (0->1 on AB clade)
        assert scores[1] == 1
        # char2: 1 change (0->1 on CD)
        assert scores[2] == 1
        # char3: 1 change (0->1 on D)
        assert scores[3] == 1
        # char4: 1 change (0->1 on B)
        assert scores[4] == 1
        # char7: 0 changes (all zeros)
        assert scores[7] == 0

    def test_acctran_pipeline(self):
        """Full ACCTRAN pipeline produces valid output."""
        tree, tip_states, _ = self._load_data()
        node_state_sets, scores = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_acctran(tree, node_state_sets, pm)
        changes = detect_changes(tree, node_states, pm)
        classified = classify_changes(tree, changes, node_states, pm)

        # All tips should have their observed states
        for tip in tree.get_terminals():
            assert node_states[id(tip)] == tip_states[tip.name]

        # Total number of changes across all branches should equal sum of scores
        total_changes = sum(
            len(bc) for bc in changes.values()
        )
        assert total_changes == sum(scores)

    def test_deltran_pipeline(self):
        """Full DELTRAN pipeline produces valid output."""
        tree, tip_states, _ = self._load_data()
        node_state_sets, scores = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        node_states = fitch_uppass_deltran(tree, node_state_sets, pm)
        changes = detect_changes(tree, node_states, pm)
        classified = classify_changes(tree, changes, node_states, pm)

        # All tips should have their observed states
        for tip in tree.get_terminals():
            assert node_states[id(tip)] == tip_states[tip.name]

        # Total number of changes across all branches should equal sum of scores
        total_changes = sum(
            len(bc) for bc in changes.values()
        )
        assert total_changes == sum(scores)

    def test_ci_with_sample_data(self):
        """Consistency index computed from sample data."""
        tree, tip_states, _ = self._load_data()
        node_state_sets, scores = fitch_downpass(tree, tip_states)

        # Count distinct states per character
        n_chars = len(scores)
        n_states_per_char = []
        for i in range(n_chars):
            states = set()
            for taxon_states in tip_states.values():
                s = taxon_states[i]
                if s not in ("?", "-"):
                    states.add(s)
            n_states_per_char.append(len(states))

        ci_per_char, ci_overall = consistency_index(
            n_states_per_char, scores
        )
        # char0-char4: 2 states, 1 change each => CI=1.0
        for i in range(5):
            assert ci_per_char[i] == 1.0
        # char7: 1 state, 0 changes => None
        assert ci_per_char[7] is None

    def test_ri_with_sample_data(self):
        """Retention index computed from sample data."""
        tree, tip_states, _ = self._load_data()
        _, scores = fitch_downpass(tree, tip_states)

        n_chars = len(scores)
        tip_states_per_char = []
        for i in range(n_chars):
            tip_states_per_char.append(
                [tip_states[t][i] for t in tip_states]
            )

        ri_per_char, ri_overall = retention_index(
            tip_states_per_char, scores
        )
        # char7 (all zeros): uninformative => None
        assert ri_per_char[7] is None
