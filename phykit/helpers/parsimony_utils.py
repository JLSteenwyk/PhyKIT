"""
Reusable parsimony algorithms for discrete character analysis.

Provides generalized Fitch downpass, ACCTRAN/DELTRAN uppass,
character change detection, synapomorphy/homoplasy classification,
and consistency/retention index computation.
"""
from collections import Counter
from typing import Dict, List, Optional, Set, Tuple


def build_parent_map(tree) -> Dict[int, object]:
    """Build dict mapping node id -> parent clade."""
    parent_map = {}
    for clade in tree.find_clades(order="preorder"):
        for child in clade.clades:
            parent_map[id(child)] = clade
    return parent_map


def resolve_polytomies(tree) -> None:
    """Resolve multifurcations by inserting zero-length branches.

    Mutates the tree in place, converting any node with >2 children
    into a series of bifurcating nodes connected by zero-length branches.
    """
    from Bio.Phylo import Newick

    for clade in tree.find_clades(order="postorder"):
        while len(clade.clades) > 2:
            child1 = clade.clades.pop()
            child2 = clade.clades.pop()
            new_internal = Newick.Clade(branch_length=0.0)
            new_internal.clades = [child1, child2]
            clade.clades.append(new_internal)


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
                # If all children have empty sets (all wildcards, no
                # observed states), treat as universal match: no change.
                if all(len(s) == 0 for s in sets):
                    char_sets.append(set())
                    continue
                # If some children are empty (wildcard) and some not,
                # the non-empty ones define the intersection.
                non_empty = [s for s in sets if len(s) > 0]
                intersection = non_empty[0]
                for s in non_empty[1:]:
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


def fitch_uppass_acctran(
    tree,
    node_state_sets: Dict[int, List[Set[str]]],
    parent_map: Dict[int, object],
) -> Dict[int, List[str]]:
    """ACCTRAN uppass: assign final states pushing changes rootward.

    Root gets min(root_state_set) for tie-breaking.
    If parent state is in child state set, inherit; else use min(child set).

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


def fitch_uppass_deltran(
    tree,
    node_state_sets: Dict[int, List[Set[str]]],
    parent_map: Dict[int, object],
) -> Dict[int, List[str]]:
    """DELTRAN uppass: assign final states pushing changes tipward.

    Follows Swofford & Maddison (1987). Uses the parent's downpass
    state set (not just final state) to resolve conflicts.

    If parent final state is in child state set, inherit (delay change).
    Else: prefer min(child_set & parent_downpass_set), fallback min(child_set).

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
            pid = id(parent)
            parent_finals = node_states[pid]
            parent_downpass = node_state_sets[pid]
            finals = []
            for i in range(n_chars):
                if parent_finals[i] in state_sets[i]:
                    # Parent's state is compatible -- inherit (delay change)
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


def detect_changes(
    tree,
    node_states: Dict[int, List[str]],
    parent_map: Dict[int, object],
) -> Dict[int, List[Tuple[int, str, str]]]:
    """Detect character state changes on each branch.

    Compares each node's final states to its parent's final states.

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

    A change to state X for character i is:
    - synapomorphy: if (i, X) appears only once across all branches
    - reversal: if (i, X) appears multiple times AND X is found at
      an ancestor node
    - convergence: if (i, X) appears multiple times AND X is NOT found
      at any ancestor node

    Returns:
        Dict[child_node_id -> List[(char_idx, old, new, classification)]]
    """
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
) -> Tuple[List[Optional[float]], Optional[float]]:
    """Compute per-character and overall consistency index.

    CI_i = (n_states_i - 1) / observed_i
    CI_overall = sum(min_i) / sum(observed_i)

    Returns None for characters with 0 observed changes.
    """
    ci_per_char: List[Optional[float]] = []
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
    min_changes_i = n_states_i - 1
    RI_i = (max_i - observed_i) / (max_i - min_i)

    Returns None for uninformative characters where max == min.
    """
    ri_per_char: List[Optional[float]] = []
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
