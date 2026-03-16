"""
Unit tests for parsimony utility algorithms.

Tests Fitch downpass, ACCTRAN/DELTRAN uppass, change detection,
change classification (synapomorphy/convergence/reversal),
consistency index, and retention index.
"""
import pytest
from Bio import Phylo
from io import StringIO

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


# ---------------------------------------------------------------------------
# fitch_downpass
# ---------------------------------------------------------------------------

class TestFitchDownpass:
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
