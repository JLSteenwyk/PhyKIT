# Character Map Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `character_map` command that maps synapomorphies and homoplasies onto a phylogeny using Fitch parsimony with ACCTRAN/DELTRAN optimization.

**Architecture:** Reusable parsimony algorithms live in `phykit/helpers/parsimony_utils.py` (Fitch downpass, ACCTRAN/DELTRAN uppass, change detection, classification, CI/RI). The service class `phykit/services/tree/character_map.py` handles I/O, orchestration, and cladogram/phylogram plotting. Cross-validated against R's phangorn.

**Tech Stack:** Python, BioPython (tree parsing), matplotlib (plotting), phangorn (R, validation only)

**Spec:** `docs/superpowers/specs/2026-03-16-character-map-design.md`

---

## Chunk 1: Test Data and Parsimony Utilities (Core Algorithms)

### Task 1: Create test dataset

**Files:**
- Create: `tests/sample_files/character_matrix_simple.tsv`
- Create: `tests/sample_files/tree_character_map.tre`

This dataset is used by all subsequent tests. 6 taxa, 8 characters with known properties:
- Characters with unique changes (synapomorphies)
- Characters with convergent gains (homoplasies)
- Characters with reversals
- A multistate character (3 states)
- A character with missing data (?)
- An uninformative character (same state everywhere)

The tree: `((A:1,B:1):1,((C:1,D:1):1,E:1):1,F:1);` (trifurcating root, 6 taxa)

- [ ] **Step 1: Create the tree file**

```
# tests/sample_files/tree_character_map.tre
((A:1,B:1):1,((C:1,D:1):1,E:1):1,F:1);
```

- [ ] **Step 2: Create the character matrix**

```tsv
taxon	char0	char1	char2	char3	char4	char5	char6	char7
A	0	1	0	0	0	1	1	0
B	0	1	0	0	1	1	1	0
C	1	0	1	0	0	0	2	0
D	1	0	1	1	0	0	0	0
E	1	0	0	0	0	0	1	0
F	0	0	0	0	0	0	0	0
```

Character design:
- `char0`: 0→1 on branch to (C,D,E) clade = synapomorphy
- `char1`: 0→1 on branch to (A,B) clade = synapomorphy
- `char2`: 0→1 on branch to (C,D) = synapomorphy
- `char3`: 0→1 on branch to D only = autapomorphy (synapomorphy at terminal)
- `char4`: 0→1 on branch to B only = autapomorphy
- `char5`: A,B=1, C,D,E=0, F=0 — ACCTRAN places 0→1 on (A,B) stem, then
  1→0 on (C,D,E) stem = reversal to ancestral state 0
- `char6`: multistate — 3 states (0,1,2), multiple changes
- `char7`: uninformative — all zeros

- [ ] **Step 3: Commit test data**

```bash
git add tests/sample_files/character_matrix_simple.tsv tests/sample_files/tree_character_map.tre
git commit -m "test: add sample data for character_map command"
```

---

### Task 2: Fitch downpass and helper utilities

**Files:**
- Create: `phykit/helpers/parsimony_utils.py`
- Create: `tests/unit/helpers/test_parsimony_utils.py`

- [ ] **Step 1: Write failing tests for `build_parent_map` and `resolve_polytomies`**

In `tests/unit/helpers/test_parsimony_utils.py`:

```python
import pytest
from Bio import Phylo
from io import StringIO

from phykit.helpers.parsimony_utils import (
    build_parent_map,
    resolve_polytomies,
)


def _make_tree(newick):
    return Phylo.read(StringIO(newick), "newick")


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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/helpers/test_parsimony_utils.py -v`
Expected: FAIL (module not found)

- [ ] **Step 3: Implement `build_parent_map` and `resolve_polytomies`**

Create `phykit/helpers/parsimony_utils.py`:

```python
"""
Reusable parsimony algorithms for discrete character analysis.

Provides generalized Fitch downpass, ACCTRAN/DELTRAN uppass,
character change detection, synapomorphy/homoplasy classification,
and consistency/retention index computation.
"""
from typing import Dict, List, Optional, Set, Tuple


def build_parent_map(tree) -> Dict[int, object]:
    """Build dict mapping node id -> parent clade."""
    parent_map = {}
    for clade in tree.find_clades(order="preorder"):
        for child in clade.clades:
            parent_map[id(child)] = clade
    return parent_map


def resolve_polytomies(tree) -> None:
    """Resolve multifurcations by inserting zero-length branches."""
    from Bio.Phylo import Newick
    for clade in tree.find_clades(order="postorder"):
        while len(clade.clades) > 2:
            child1 = clade.clades.pop()
            child2 = clade.clades.pop()
            new_internal = Newick.Clade(branch_length=0.0)
            new_internal.clades = [child1, child2]
            clade.clades.append(new_internal)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/helpers/test_parsimony_utils.py -v`
Expected: PASS

- [ ] **Step 5: Write failing tests for `fitch_downpass`**

Append to `tests/unit/helpers/test_parsimony_utils.py`:

```python
from phykit.helpers.parsimony_utils import fitch_downpass


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

    def test_wildcard_handled(self):
        """Wildcard (?) matches any state."""
        tree = _make_tree("((A:1,B:1):1,(C:1,D:1):1);")
        resolve_polytomies(tree)
        tip_states = {"A": ["0"], "B": ["?"], "C": ["1"], "D": ["1"]}
        node_state_sets, scores = fitch_downpass(tree, tip_states)
        # B is wildcard, so (A,B) can be {0}, and (C,D)={1} => 1 change
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
```

- [ ] **Step 6: Implement `fitch_downpass`**

Add to `phykit/helpers/parsimony_utils.py`:

```python
def fitch_downpass(
    tree, tip_states: Dict[str, List[str]]
) -> Tuple[Dict[int, List[Set[str]]], List[int]]:
    """Generalized Fitch downpass for discrete characters.

    Args:
        tree: BioPython tree (must be bifurcating).
        tip_states: Dict mapping taxon name -> list of state strings,
            one per character. Use "?" or "-" for wildcards.

    Returns:
        node_state_sets: Dict[node_id -> List[Set[str]]] per-character
            state sets at each node.
        scores: List[int] parsimony score per character.
    """
    # Determine number of characters from first taxon
    n_chars = len(next(iter(tip_states.values())))
    wildcard = {"?", "-"}

    # Collect all observed states per character (for wildcard expansion)
    all_states_per_char: List[Set[str]] = [set() for _ in range(n_chars)]
    for states in tip_states.values():
        for i, s in enumerate(states):
            if s not in wildcard:
                all_states_per_char[i].add(s)

    node_state_sets: Dict[int, List[Set[str]]] = {}
    scores = [0] * n_chars

    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            char_sets = []
            for i, s in enumerate(tip_states[clade.name]):
                if s in wildcard:
                    char_sets.append(set(all_states_per_char[i]))
                else:
                    char_sets.append({s})
            node_state_sets[id(clade)] = char_sets
        else:
            child_state_lists = [
                node_state_sets[id(c)] for c in clade.clades
            ]
            char_sets = []
            for i in range(n_chars):
                sets = [cs[i] for cs in child_state_lists]
                intersection = sets[0]
                for s in sets[1:]:
                    intersection = intersection & s
                if intersection:
                    char_sets.append(intersection)
                else:
                    union = sets[0]
                    for s in sets[1:]:
                        union = union | s
                    char_sets.append(union)
                    scores[i] += 1
            node_state_sets[id(clade)] = char_sets

    return node_state_sets, scores
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `pytest tests/unit/helpers/test_parsimony_utils.py -v`
Expected: PASS

- [ ] **Step 8: Commit**

```bash
git add phykit/helpers/parsimony_utils.py tests/unit/helpers/test_parsimony_utils.py
git commit -m "feat: add parsimony_utils with Fitch downpass, build_parent_map, resolve_polytomies"
```

---

### Task 3: ACCTRAN and DELTRAN uppass

**Files:**
- Modify: `phykit/helpers/parsimony_utils.py`
- Modify: `tests/unit/helpers/test_parsimony_utils.py`

- [ ] **Step 1: Write failing tests for ACCTRAN uppass**

Append to test file:

```python
from phykit.helpers.parsimony_utils import fitch_uppass_acctran


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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/helpers/test_parsimony_utils.py::TestAcctranUppass -v`
Expected: FAIL

- [ ] **Step 3: Implement `fitch_uppass_acctran`**

Add to `phykit/helpers/parsimony_utils.py`:

```python
def fitch_uppass_acctran(
    tree,
    node_state_sets: Dict[int, List[Set[str]]],
    parent_map: Dict[int, object],
) -> Dict[int, List[str]]:
    """ACCTRAN uppass: assign final states pushing changes rootward.

    Args:
        tree: bifurcating BioPython tree.
        node_state_sets: from fitch_downpass.
        parent_map: from build_parent_map.

    Returns:
        Dict[node_id -> List[str]] final state per character per node.
    """
    n_chars = len(next(iter(node_state_sets.values())))
    node_states: Dict[int, List[str]] = {}

    for clade in tree.find_clades(order="preorder"):
        cid = id(clade)
        state_sets = node_state_sets[cid]

        if clade == tree.root:
            node_states[cid] = [min(ss) for ss in state_sets]
        else:
            parent = parent_map[cid]
            parent_finals = node_states[id(parent)]
            finals = []
            for i in range(n_chars):
                if parent_finals[i] in state_sets[i]:
                    finals.append(parent_finals[i])
                else:
                    finals.append(min(state_sets[i]))
            node_states[cid] = finals

    return node_states
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/helpers/test_parsimony_utils.py::TestAcctranUppass -v`
Expected: PASS

- [ ] **Step 5: Write failing tests for DELTRAN uppass**

```python
from phykit.helpers.parsimony_utils import fitch_uppass_deltran


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

    def test_acctran_vs_deltran_differ(self):
        """ACCTRAN and DELTRAN produce different ancestral states."""
        # Tree: ((A,B),(C,(D,E)))
        tree = _make_tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")
        resolve_polytomies(tree)
        # A=0, B=1, C=1, D=0, E=0
        tip_states = {"A": ["0"], "B": ["1"], "C": ["1"], "D": ["0"], "E": ["0"]}
        node_state_sets, _ = fitch_downpass(tree, tip_states)
        pm = build_parent_map(tree)
        acctran = fitch_uppass_acctran(tree, node_state_sets, pm)
        deltran = fitch_uppass_deltran(tree, node_state_sets, pm)
        # They should differ at some internal node
        internal_ids = [
            id(c) for c in tree.find_clades() if not c.is_terminal()
        ]
        differs = any(
            acctran[nid] != deltran[nid] for nid in internal_ids
        )
        assert differs
```

- [ ] **Step 6: Implement `fitch_uppass_deltran`**

Add to `phykit/helpers/parsimony_utils.py`:

```python
def fitch_uppass_deltran(
    tree,
    node_state_sets: Dict[int, List[Set[str]]],
    parent_map: Dict[int, object],
) -> Dict[int, List[str]]:
    """DELTRAN uppass: assign final states pushing changes tipward.

    Follows Swofford & Maddison (1987). Uses the parent's downpass
    state set (not just final state) to resolve conflicts.
    """
    n_chars = len(next(iter(node_state_sets.values())))
    node_states: Dict[int, List[str]] = {}

    for clade in tree.find_clades(order="preorder"):
        cid = id(clade)
        state_sets = node_state_sets[cid]

        if clade == tree.root:
            node_states[cid] = [min(ss) for ss in state_sets]
        else:
            parent = parent_map[cid]
            pid = id(parent)
            parent_finals = node_states[pid]
            parent_downpass = node_state_sets[pid]
            finals = []
            for i in range(n_chars):
                if parent_finals[i] in state_sets[i]:
                    # Parent's state is compatible — inherit (delay change)
                    finals.append(parent_finals[i])
                else:
                    # Must change. Prefer state shared with parent's
                    # downpass set for continuity.
                    shared = state_sets[i] & parent_downpass[i]
                    if shared:
                        finals.append(min(shared))
                    else:
                        finals.append(min(state_sets[i]))
            node_states[cid] = finals

    return node_states
```

- [ ] **Step 7: Run all parsimony utils tests**

Run: `pytest tests/unit/helpers/test_parsimony_utils.py -v`
Expected: PASS

- [ ] **Step 8: Commit**

```bash
git add phykit/helpers/parsimony_utils.py tests/unit/helpers/test_parsimony_utils.py
git commit -m "feat: add ACCTRAN and DELTRAN uppass to parsimony_utils"
```

---

### Task 4: Change detection, classification, CI/RI

**Files:**
- Modify: `phykit/helpers/parsimony_utils.py`
- Modify: `tests/unit/helpers/test_parsimony_utils.py`

- [ ] **Step 1: Write failing tests for change detection and classification**

```python
from phykit.helpers.parsimony_utils import (
    detect_changes,
    classify_changes,
    consistency_index,
    retention_index,
)


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
        all_changes = [c for branch_changes in changes.values() for c in branch_changes]
        char0_changes = [c for c in all_changes if c[0] == 0]
        assert len(char0_changes) == 1


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
        # Tree: ((A,B),(C,(D,E)))
        tree = _make_tree("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")
        resolve_polytomies(tree)
        # A=1, B=1, C=0, D=0, E=0
        # ACCTRAN: root=0, (A,B)=1 (change), (C,(D,E))=0
        # Then within (A,B), no further change.
        # Now use a second scenario: root=0, gain 1 at subtree, then revert
        # A=0, B=0, C=1, D=0, E=0
        # root=0, (C,(D,E)) ancestor gets {0,1}, ACCTRAN picks 0
        # C gets 1 (change), D,E stay 0
        # That's a synapomorphy. Let's build a reversal:
        # A=1, B=1, C=1, D=0, E=0
        # root={0,1}, ACCTRAN root=0
        # (A,B) ancestor: {1}, change 0->1
        # (C,(D,E)) ancestor: {0,1}, inherits 0
        # C: {1}, change 0->1 (convergence with A,B!)
        # (D,E): {0}, inherits 0
        # Still no reversal. Need deeper tree.
        tree2 = _make_tree("(((A:1,B:1):1,C:1):1,(D:1,E:1):1);")
        resolve_polytomies(tree2)
        # A=1, B=1, C=0, D=0, E=0
        # root=0, ((A,B),C) ancestor={0,1} ACCTRAN->0
        # (A,B) ancestor={1}, change 0->1 (synapomorphy)
        # C={0}, inherits 0
        # Now: A=1, B=0, C=0, D=0, E=0
        # (A,B) ancestor={0,1}, ACCTRAN inherits from parent (0? no, union)
        # Actually simpler: just check that the function finds reversals
        # when a descendant returns to a state seen at an ancestor
        tip_states = {"A": ["1"], "B": ["1"], "C": ["0"], "D": ["0"], "E": ["1"]}
        node_state_sets, _ = fitch_downpass(tree2, tip_states)
        pm = build_parent_map(tree2)
        node_states = fitch_uppass_acctran(tree2, node_state_sets, pm)
        changes = detect_changes(tree2, node_states, pm)
        classified = classify_changes(tree2, changes, node_states, pm)
        all_c = [c for branch in classified.values() for c in branch]
        types = set(c[3] for c in all_c)
        # Should have at least synapomorphy or convergence/reversal
        assert len(all_c) >= 1


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
        # All same state: max = min = 0
        ri_per_char, ri_overall = retention_index(
            tip_states_per_char=[["0", "0", "0", "0"]],
            observed_per_char=[0],
        )
        assert ri_per_char[0] is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/helpers/test_parsimony_utils.py -k "Detect or Classify or Consistency or Retention" -v`
Expected: FAIL

- [ ] **Step 3: Implement `detect_changes`, `classify_changes`, `consistency_index`, `retention_index`**

Add to `phykit/helpers/parsimony_utils.py`:

```python
def detect_changes(
    tree,
    node_states: Dict[int, List[str]],
    parent_map: Dict[int, object],
) -> Dict[int, List[Tuple[int, str, str]]]:
    """Detect character state changes on each branch.

    Returns:
        Dict[child_node_id -> List[(char_idx, old_state, new_state)]]
    """
    n_chars = len(next(iter(node_states.values())))
    changes: Dict[int, List[Tuple[int, str, str]]] = {}

    for clade in tree.find_clades(order="preorder"):
        if clade == tree.root:
            continue
        cid = id(clade)
        pid = id(parent_map[cid])
        branch_changes = []
        for i in range(n_chars):
            old = node_states[pid][i]
            new = node_states[cid][i]
            if old != new:
                branch_changes.append((i, old, new))
        if branch_changes:
            changes[cid] = branch_changes

    return changes


def classify_changes(
    tree,
    branch_changes: Dict[int, List[Tuple[int, str, str]]],
    node_states: Dict[int, List[str]],
    parent_map: Dict[int, object],
) -> Dict[int, List[Tuple[int, str, str, str]]]:
    """Classify each change as synapomorphy, convergence, or reversal.

    Returns:
        Dict[child_node_id -> List[(char_idx, old, new, classification)]]
    """
    from collections import Counter

    # Count how many times each (char, new_state) appears
    transition_counts: Counter = Counter()
    for changes in branch_changes.values():
        for char_idx, old, new in changes:
            transition_counts[(char_idx, new)] += 1

    result: Dict[int, List[Tuple[int, str, str, str]]] = {}
    for cid, changes in branch_changes.items():
        classified = []
        for char_idx, old, new in changes:
            if transition_counts[(char_idx, new)] == 1:
                classification = "synapomorphy"
            else:
                # Walk ancestor chain to check if new_state is ancestral
                is_reversal = False
                node = parent_map.get(cid)
                while node is not None:
                    if node_states[id(node)][char_idx] == new:
                        is_reversal = True
                        break
                    node = parent_map.get(id(node))
                classification = "reversal" if is_reversal else "convergence"
            classified.append((char_idx, old, new, classification))
        result[cid] = classified

    return result


def consistency_index(
    n_states_per_char: List[int],
    observed_per_char: List[int],
) -> Tuple[List[Optional[float]], float]:
    """Compute per-character and overall consistency index.

    CI_i = (n_states_i - 1) / observed_i
    CI_overall = sum(min_i) / sum(observed_i)
    """
    ci_per_char = []
    sum_min = 0
    sum_obs = 0
    for n_states, observed in zip(n_states_per_char, observed_per_char):
        min_changes = n_states - 1
        if observed == 0:
            ci_per_char.append(None)
        else:
            ci_per_char.append(min_changes / observed)
            sum_min += min_changes
            sum_obs += observed

    ci_overall = sum_min / sum_obs if sum_obs > 0 else None
    return ci_per_char, ci_overall


def retention_index(
    tip_states_per_char: List[List[str]],
    observed_per_char: List[int],
) -> Tuple[List[Optional[float]], Optional[float]]:
    """Compute per-character and overall retention index (Farris 1989).

    max_changes_i = n_taxa - f_max_i (count of most frequent state)
    RI_i = (max_i - observed_i) / (max_i - min_i)
    """
    from collections import Counter

    ri_per_char = []
    sum_num = 0  # sum(max - observed)
    sum_den = 0  # sum(max - min)

    for states, observed in zip(tip_states_per_char, observed_per_char):
        # Filter out wildcards
        clean = [s for s in states if s not in ("?", "-")]
        if not clean:
            ri_per_char.append(None)
            continue

        counts = Counter(clean)
        n_taxa = len(clean)
        n_states = len(counts)
        f_max = max(counts.values())
        max_changes = n_taxa - f_max
        min_changes = n_states - 1

        if max_changes == min_changes:
            ri_per_char.append(None)
        else:
            ri = (max_changes - observed) / (max_changes - min_changes)
            ri_per_char.append(ri)
            sum_num += max_changes - observed
            sum_den += max_changes - min_changes

    ri_overall = sum_num / sum_den if sum_den > 0 else None
    return ri_per_char, ri_overall
```

- [ ] **Step 4: Run all tests**

Run: `pytest tests/unit/helpers/test_parsimony_utils.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add phykit/helpers/parsimony_utils.py tests/unit/helpers/test_parsimony_utils.py
git commit -m "feat: add change detection, classification, CI/RI to parsimony_utils"
```

---

## Chunk 2: Service Class, CLI Registration, Plotting

### Task 5: Character matrix parsing and service skeleton

**Files:**
- Create: `phykit/services/tree/character_map.py`
- Create: `tests/unit/services/tree/test_character_map.py`

- [ ] **Step 1: Write failing test for matrix parsing**

```python
import pytest
from argparse import Namespace
from pathlib import Path

from phykit.services.tree.character_map import CharacterMap


class TestParseCharacterMatrix:
    def test_parses_simple_matrix(self, tmp_path):
        matrix_file = tmp_path / "matrix.tsv"
        matrix_file.write_text(
            "taxon\tc0\tc1\n"
            "A\t0\t1\n"
            "B\t1\t0\n"
            "C\t0\t1\n"
        )
        args = Namespace(
            tree="tests/sample_files/tree_character_map.tre",
            data=str(matrix_file),
            output=str(tmp_path / "out.png"),
        )
        svc = CharacterMap(args)
        char_names, tip_states = svc._parse_character_matrix(str(matrix_file))
        assert char_names == ["c0", "c1"]
        assert tip_states["A"] == ["0", "1"]
        assert tip_states["B"] == ["1", "0"]

    def test_missing_data_preserved(self, tmp_path):
        matrix_file = tmp_path / "matrix.tsv"
        matrix_file.write_text(
            "taxon\tc0\n"
            "A\t0\n"
            "B\t?\n"
        )
        args = Namespace(
            tree="tests/sample_files/tree_character_map.tre",
            data=str(matrix_file),
            output=str(tmp_path / "out.png"),
        )
        svc = CharacterMap(args)
        _, tip_states = svc._parse_character_matrix(str(matrix_file))
        assert tip_states["B"] == ["?"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/services/tree/test_character_map.py -v`
Expected: FAIL (module not found)

- [ ] **Step 3: Implement service skeleton with matrix parser**

Create `phykit/services/tree/character_map.py`:

```python
"""
Character map: synapomorphy and homoplasy visualization on a phylogeny.

Maps character state changes onto branches using Fitch parsimony with
ACCTRAN or DELTRAN optimization. Produces a cladogram (default) or
phylogram with color-coded circles: blue = synapomorphy, red = convergence,
gray = reversal.
"""
from typing import Dict, List, Tuple

from .base import Tree
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig
from ...helpers.parsimony_utils import (
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
from ...errors import PhykitUserError


class CharacterMap(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.data_path = parsed["data_path"]
        self.output_path = parsed["output_path"]
        self.optimization = parsed["optimization"]
        self.phylogram = parsed["phylogram"]
        self.characters_filter = parsed["characters_filter"]
        self.verbose = parsed["verbose"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> Dict:
        chars_str = getattr(args, "characters", None)
        chars_filter = None
        if chars_str:
            chars_filter = [int(c.strip()) for c in chars_str.split(",")]

        return dict(
            tree_file_path=args.tree,
            data_path=args.data,
            output_path=args.output,
            optimization=getattr(args, "optimization", "acctran"),
            phylogram=getattr(args, "phylogram", False),
            characters_filter=chars_filter,
            verbose=getattr(args, "verbose", False),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    def _parse_character_matrix(
        self, path: str
    ) -> Tuple[List[str], Dict[str, List[str]]]:
        """Parse TSV character matrix.

        Returns (character_names, {taxon: [state0, state1, ...]}).
        """
        try:
            with open(path) as f:
                lines = [
                    l.strip() for l in f if l.strip() and not l.startswith("#")
                ]
        except FileNotFoundError:
            raise PhykitUserError(
                [f"{path} corresponds to no such file or directory."],
                code=2,
            )

        if len(lines) < 2:
            raise PhykitUserError(
                ["Character matrix must have a header row and at least one data row."],
                code=2,
            )

        header = lines[0].split("\t")
        char_names = header[1:]  # first column is taxon name

        tip_states: Dict[str, List[str]] = {}
        for line in lines[1:]:
            parts = line.split("\t")
            taxon = parts[0]
            states = parts[1:]
            if len(states) != len(char_names):
                raise PhykitUserError(
                    [
                        f"Taxon '{taxon}' has {len(states)} states "
                        f"but header defines {len(char_names)} characters."
                    ],
                    code=2,
                )
            tip_states[taxon] = states

        return char_names, tip_states

    def run(self) -> None:
        # Placeholder — implemented in Task 6
        pass
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/unit/services/tree/test_character_map.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add phykit/services/tree/character_map.py tests/unit/services/tree/test_character_map.py
git commit -m "feat: add CharacterMap service skeleton with TSV matrix parser"
```

---

### Task 6: Service `run()` method — orchestration and text/JSON output

**Files:**
- Modify: `phykit/services/tree/character_map.py`
- Modify: `tests/unit/services/tree/test_character_map.py`

- [ ] **Step 1: Write failing test for JSON output with ground truth**

Append to test file:

```python
class TestCharacterMapRun:
    def test_json_output(self, mocker, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_character_map.tre",
            data="tests/sample_files/character_matrix_simple.tsv",
            output=str(tmp_path / "out.png"),
            optimization="acctran",
            json=True,
        )
        svc = CharacterMap(args)
        mocked_json = mocker.patch(
            "phykit.services.tree.character_map.print_json"
        )
        svc.run()

        payload = mocked_json.call_args.args[0]
        assert payload["optimization"] == "acctran"
        assert payload["n_characters"] == 8
        assert "ci" in payload
        assert "ri" in payload
        assert "tree_length" in payload
        assert len(payload["characters"]) == 8

    def test_creates_png(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_character_map.tre",
            data="tests/sample_files/character_matrix_simple.tsv",
            output=str(tmp_path / "out.png"),
        )
        svc = CharacterMap(args)
        svc.run()
        assert Path(tmp_path / "out.png").exists()

    def test_too_few_shared_taxa(self, tmp_path):
        tree_file = tmp_path / "small.tre"
        tree_file.write_text("(X:1,Y:1);")
        matrix_file = tmp_path / "matrix.tsv"
        matrix_file.write_text("taxon\tc0\nX\t0\nY\t1\n")
        args = Namespace(
            tree=str(tree_file),
            data=str(matrix_file),
            output=str(tmp_path / "out.png"),
        )
        svc = CharacterMap(args)
        with pytest.raises(SystemExit):
            svc.run()
```

- [ ] **Step 2: Implement `run()` method**

Replace the placeholder `run()` in `character_map.py`:

```python
    def run(self) -> None:
        import copy

        tree = self.read_tree_file()
        char_names, tip_states_raw = self._parse_character_matrix(self.data_path)

        # Prune to shared taxa
        tree_tips = set(t.name for t in tree.get_terminals())
        data_taxa = set(tip_states_raw.keys())
        shared = tree_tips & data_taxa

        if len(shared) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa between tree and data.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        tips_to_prune = [t for t in tree_tips if t not in shared]
        if tips_to_prune:
            tree = self.prune_tree_using_taxa_list(tree, tips_to_prune)

        tip_states = {k: v for k, v in tip_states_raw.items() if k in shared}

        # Resolve polytomies
        resolve_polytomies(tree)

        # Ladderize if requested
        if self.plot_config.ladderize:
            tree.ladderize()

        # Run Fitch downpass
        node_state_sets, scores = fitch_downpass(tree, tip_states)

        # Uppass
        parent_map = build_parent_map(tree)
        if self.optimization == "deltran":
            node_states = fitch_uppass_deltran(tree, node_state_sets, parent_map)
        else:
            node_states = fitch_uppass_acctran(tree, node_state_sets, parent_map)

        # Detect and classify changes
        changes = detect_changes(tree, node_states, parent_map)
        classified = classify_changes(tree, changes, node_states, parent_map)

        # Compute CI/RI
        n_chars = len(char_names)
        n_states_per_char = []
        tip_states_per_char = []
        for i in range(n_chars):
            states_for_char = [
                tip_states[taxon][i]
                for taxon in tip_states
                if tip_states[taxon][i] not in ("?", "-")
            ]
            n_states_per_char.append(len(set(states_for_char)))
            tip_states_per_char.append(states_for_char)

        ci_per_char, ci_overall = consistency_index(n_states_per_char, scores)
        ri_per_char, ri_overall = retention_index(tip_states_per_char, scores)

        # Count parsimony-informative characters
        from collections import Counter
        n_informative = 0
        for i in range(n_chars):
            counts = Counter(tip_states_per_char[i])
            if sum(1 for c in counts.values() if c >= 2) >= 2:
                n_informative += 1

        tree_length = sum(scores)

        # Plot
        self._plot_character_map(tree, classified, char_names, parent_map)

        # Output
        if self.json_output:
            self._print_json(
                char_names, scores, ci_per_char, ci_overall,
                ri_per_char, ri_overall, classified, tree, parent_map,
                n_informative, tree_length,
            )
        elif self.verbose:
            self._print_verbose(
                char_names, scores, ci_per_char, ri_per_char,
                classified, tree, parent_map, n_informative,
                tree_length, ci_overall, ri_overall,
            )
        else:
            self._print_summary(
                n_chars, n_informative, tree_length,
                ci_overall, ri_overall,
            )
```

Also add the output helper methods (`_print_summary`, `_print_verbose`, `_print_json`) and a placeholder `_plot_character_map`. These are straightforward text formatting — implement them following the spec's output format. The plot method is implemented in Task 7.

- [ ] **Step 3: Run tests**

Run: `pytest tests/unit/services/tree/test_character_map.py -v`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add phykit/services/tree/character_map.py tests/unit/services/tree/test_character_map.py
git commit -m "feat: add character_map run() with full parsimony pipeline and output"
```

---

### Task 7: Cladogram/phylogram plotting

**Files:**
- Modify: `phykit/services/tree/character_map.py`
- Modify: `tests/unit/services/tree/test_character_map.py`

- [ ] **Step 1: Write failing test for cladogram and phylogram modes**

```python
class TestPlotModes:
    def test_cladogram_default(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_character_map.tre",
            data="tests/sample_files/character_matrix_simple.tsv",
            output=str(tmp_path / "clado.png"),
        )
        svc = CharacterMap(args)
        svc.run()
        assert Path(tmp_path / "clado.png").exists()
        assert Path(tmp_path / "clado.png").stat().st_size > 0

    def test_phylogram_mode(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_character_map.tre",
            data="tests/sample_files/character_matrix_simple.tsv",
            output=str(tmp_path / "phylo.png"),
            phylogram=True,
        )
        svc = CharacterMap(args)
        svc.run()
        assert Path(tmp_path / "phylo.png").exists()

    def test_character_filter(self, tmp_path):
        args = Namespace(
            tree="tests/sample_files/tree_character_map.tre",
            data="tests/sample_files/character_matrix_simple.tsv",
            output=str(tmp_path / "filtered.png"),
            characters="0,1",
        )
        svc = CharacterMap(args)
        svc.run()
        assert Path(tmp_path / "filtered.png").exists()
```

- [ ] **Step 2: Implement `_plot_character_map`**

Replace the placeholder. Follow the phylogram drawing pattern from `quartet_pie.py` (node_x/node_y computation, horizontal + vertical branch lines). For cladogram mode, compute x from topological depth. Place colored circles along branches at evenly-spaced positions.

Key implementation notes:
- Use `matplotlib.use("Agg")` for headless rendering
- Compute `node_x` based on depth (cladogram) or branch length (phylogram)
- Compute `node_y` from tip indices + internal node averaging
- For each branch with classified changes, place circles at even intervals
- Color: `#2b8cbe` (synapomorphy), `#d62728` (convergence), `#969696` (reversal)
- Annotate: character index above, state transition below
- If `--characters` filter is set, only draw changes for those characters
- Add legend with `matplotlib.patches.Patch`

The method is ~120 lines following the established PhyKIT pattern. Use `config.resolve(n_rows=len(tips))` for auto-scaling.

- [ ] **Step 3: Run tests**

Run: `pytest tests/unit/services/tree/test_character_map.py -v`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add phykit/services/tree/character_map.py tests/unit/services/tree/test_character_map.py
git commit -m "feat: add cladogram/phylogram plotting to character_map"
```

---

### Task 8: CLI registration and integration tests

**Files:**
- Modify: `phykit/cli_registry.py`
- Modify: `phykit/service_factories.py`
- Modify: `phykit/services/tree/__init__.py`
- Modify: `phykit/phykit.py`
- Modify: `setup.py` (console_scripts entry points)
- Create: `tests/integration/tree/test_character_map_integration.py`

- [ ] **Step 1: Register the command**

In `phykit/cli_registry.py`, add:
```python
"charmap": "character_map",
"synapomorphy_map": "character_map",
```

In `phykit/service_factories.py`, add:
```python
CharacterMap = _LazyServiceFactory("phykit.services.tree.character_map", "CharacterMap")
```

In `phykit/services/tree/__init__.py`, add to the `_EXPORTS` dict:
```python
"CharacterMap": "character_map",
```

In `setup.py`, add to the `console_scripts` list in `entry_points`:
```python
"pk_character_map = phykit.phykit:character_map",
"pk_charmap = phykit.phykit:character_map",
"pk_synapomorphy_map = phykit.phykit:character_map",
```

- [ ] **Step 2: Add command handler in `phykit/phykit.py`**

Add a `def character_map(argv):` static method following the pattern from `parsimony_score`. Include full help text describing the command, its arguments, and the polytomy/ACCTRAN/DELTRAN behavior. Add to the menu listing. Add `[--ladderize]` to the usage line.

Also add a module-level function at the bottom of `phykit/phykit.py` (required by console_scripts):
```python
def character_map():
    Phykit.character_map(sys.argv[1:])
```

- [ ] **Step 3: Write integration tests**

Create `tests/integration/tree/test_character_map_integration.py`:

```python
from mock import patch
from pathlib import Path
import json
import sys

import pytest

from phykit.phykit import Phykit

here = Path(__file__)


@pytest.mark.integration
class TestCharacterMap:
    @patch("builtins.print")
    def test_basic_invocation(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap.png")
        testargs = [
            "phykit", "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_alias_charmap(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap.png")
        testargs = [
            "phykit", "charmap",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_alias_synapomorphy_map(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap.png")
        testargs = [
            "phykit", "synapomorphy_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_json_ground_truth(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap.png")
        testargs = [
            "phykit", "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output, "--json",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["n_characters"] == 8
        assert payload["optimization"] == "acctran"
        assert "ci" in payload
        assert "ri" in payload

    @patch("builtins.print")
    def test_deltran_option(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap.png")
        testargs = [
            "phykit", "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output, "--optimization", "deltran",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_ladderize_flag(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap.png")
        testargs = [
            "phykit", "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output, "--ladderize",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_character_filter(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap.png")
        testargs = [
            "phykit", "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output, "--characters", "0,1",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_verbose_output(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap.png")
        testargs = [
            "phykit", "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output, "--verbose",
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        # Verbose output should contain per-character detail
        output_text = " ".join(
            str(c) for call in mocked_print.call_args_list for c in call.args
        )
        assert "Character" in output_text or "char" in output_text

    @patch("builtins.print")
    def test_pdf_output(self, mocked_print, tmp_path):
        output = str(tmp_path / "charmap.pdf")
        testargs = [
            "phykit", "character_map",
            "-t", f"{here.parent.parent.parent}/sample_files/tree_character_map.tre",
            "-d", f"{here.parent.parent.parent}/sample_files/character_matrix_simple.tsv",
            "-o", output,
        ]
        with patch.object(sys, "argv", testargs):
            Phykit()
        assert Path(output).exists()

    @patch("builtins.print")
    def test_character_filter_preserves_ci_ri(self, mocked_print, tmp_path):
        """--characters filter should not change CI/RI values."""
        base = f"{here.parent.parent.parent}/sample_files"
        # Run without filter
        out1 = str(tmp_path / "full.png")
        with patch.object(sys, "argv", [
            "phykit", "character_map",
            "-t", f"{base}/tree_character_map.tre",
            "-d", f"{base}/character_matrix_simple.tsv",
            "-o", out1, "--json",
        ]):
            Phykit()
        full = json.loads(mocked_print.call_args.args[0])

        mocked_print.reset_mock()

        # Run with filter
        out2 = str(tmp_path / "filtered.png")
        with patch.object(sys, "argv", [
            "phykit", "character_map",
            "-t", f"{base}/tree_character_map.tre",
            "-d", f"{base}/character_matrix_simple.tsv",
            "-o", out2, "--json", "--characters", "0,1",
        ]):
            Phykit()
        filtered = json.loads(mocked_print.call_args.args[0])

        assert full["ci"] == filtered["ci"]
        assert full["ri"] == filtered["ri"]
```

- [ ] **Step 4: Run all tests**

Run: `pytest tests/unit/helpers/test_parsimony_utils.py tests/unit/services/tree/test_character_map.py tests/integration/tree/test_character_map_integration.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add phykit/cli_registry.py phykit/service_factories.py phykit/services/tree/__init__.py phykit/phykit.py tests/integration/tree/test_character_map_integration.py
git commit -m "feat: register character_map command with CLI, add integration tests"
```

---

## Chunk 3: Validation, Documentation, Release

### Task 9: R cross-validation

**Files:**
- Create: `tests/r_validation/validate_character_map.R`

- [ ] **Step 1: Write R validation script**

```r
#!/usr/bin/env Rscript
# validate_character_map.R
#
# Cross-validate character_map CI/RI against phangorn.
#
# Requires: ape, phangorn
#
# Usage:
#   cd tests/r_validation
#   Rscript validate_character_map.R

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})

tree <- read.tree("../sample_files/tree_character_map.tre")
tree <- multi2di(tree)

# Read character matrix as phyDat
# Convert TSV to a matrix phangorn can read
data_raw <- read.table(
  "../sample_files/character_matrix_simple.tsv",
  header = TRUE, sep = "\t", row.names = 1,
  colClasses = "character"
)

# Convert to phyDat
data_matrix <- as.matrix(data_raw)
levels <- sort(unique(as.vector(data_matrix)))
levels <- levels[levels != "?" & levels != "-"]
phydat <- phyDat(data_matrix, type = "USER", levels = levels, ambiguity = "?")

# Compute parsimony
fitch_score <- parsimony(tree, phydat, method = "fitch")

# Compute CI and RI
ci <- CI(tree, phydat)
ri <- RI(tree, phydat)

cat("=== KEY VALUES FOR PYTHON TESTS ===\n")
cat(sprintf("Fitch parsimony score (tree length): %d\n", fitch_score))
cat(sprintf("Consistency index (CI): %.6f\n", ci))
cat(sprintf("Retention index (RI): %.6f\n", ri))

# Ancestral reconstruction (ACCTRAN)
anc <- ancestral.pars(tree, phydat, type = "ACCTRAN")
cat("\n=== ANCESTRAL STATES (ACCTRAN) ===\n")
for (i in seq_along(tree$node.label)) {
  node_idx <- length(tree$tip.label) + i
  cat(sprintf("Node %d:", node_idx))
  for (ch in seq_len(ncol(data_matrix))) {
    state_probs <- anc[[node_idx]][ch, ]
    state <- levels[which.max(state_probs)]
    cat(sprintf(" %s", state))
  }
  cat("\n")
}
```

- [ ] **Step 2: Commit**

```bash
git add tests/r_validation/validate_character_map.R
git commit -m "test: add R validation script for character_map CI/RI"
```

---

### Task 10: Documentation and version bump

**Files:**
- Modify: `docs/usage/index.rst`
- Modify: `docs/change_log/index.rst`
- Modify: `phykit/version.py`

- [ ] **Step 1: Add command documentation to `docs/usage/index.rst`**

Add a new section following the existing pattern (see `quartet_pie` section for reference):

```rst
.. _cmd-character_map:

Character map (synapomorphy/homoplasy mapping)
##############################################
Function names: character_map; charmap; synapomorphy_map |br|
Command line interface: pk_character_map; pk_charmap; pk_synapomorphy_map

Map character state changes onto a phylogeny using Fitch parsimony with
ACCTRAN (default) or DELTRAN optimization. Produces a cladogram (default)
or phylogram with color-coded circles on each branch showing synapomorphies
(blue), convergences (red), and reversals (gray).

Input: a Newick tree file and a TSV character matrix (header row with
character names, one row per taxon with discrete states). Missing data
(``?`` or ``-``) is treated as a wildcard.

Reports the consistency index (CI), retention index (RI), and tree length.

[usage block, options, shared plot options]
```

- [ ] **Step 2: Add changelog entry**

Add to `docs/change_log/index.rst` under the 2.1.49 section (or bump to 2.1.50):

```rst
* Added character map command (``character_map`` / ``charmap`` /
  ``synapomorphy_map``): maps synapomorphies and homoplasies onto a
  phylogeny using Fitch parsimony with ACCTRAN or DELTRAN optimization

  - Color-coded circles on branches: blue (synapomorphy), red
    (convergence), gray (reversal)
  - Supports cladogram (default) and phylogram layouts
  - Reports consistency index (CI) and retention index (RI)
  - Optional ``--characters`` filter to display specific characters
  - Cross-validated against R's phangorn package
```

- [ ] **Step 3: Bump version**

Update `phykit/version.py` to the next version.

- [ ] **Step 4: Run full test suite**

Run: `pytest tests/unit/helpers/test_parsimony_utils.py tests/unit/services/tree/test_character_map.py tests/integration/tree/test_character_map_integration.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add docs/usage/index.rst docs/change_log/index.rst phykit/version.py
git commit -m "docs: add character_map documentation and changelog entry"
```
