"""
Reusable parsimony algorithms for discrete character analysis.

Provides generalized Fitch downpass, ACCTRAN/DELTRAN uppass,
character change detection, synapomorphy/homoplasy classification,
and consistency/retention index computation.
"""
from collections import Counter


_RETENTION_INDEX_NUMPY_MIN_CELLS = 4096


def build_parent_map(tree) -> dict[int, object]:
    """Build dict mapping node id -> parent clade."""
    direct_result = _build_parent_map_direct(tree)
    if direct_result is not None:
        return direct_result

    parent_map = {}
    for clade in tree.find_clades(order="preorder"):
        for child in clade.clades:
            parent_map[id(child)] = clade
    return parent_map


def _build_parent_map_direct(tree):
    try:
        root = tree.root
        root.clades
    except AttributeError:
        return None

    parent_map = {}
    stack = [root]
    try:
        pop = stack.pop
        extend = stack.extend
        while stack:
            clade = pop()
            children = clade.clades
            for child in children:
                parent_map[id(child)] = clade
            if children:
                extend(children)
    except AttributeError:
        return None
    return parent_map


def _preorder_clades_direct(tree):
    try:
        root = tree.root
        root.clades
    except AttributeError:
        return None

    clades = []
    stack = [root]
    try:
        pop = stack.pop
        append = stack.append
        append_clade = clades.append
        while stack:
            clade = pop()
            append_clade(clade)
            children = clade.clades
            if children:
                child_count = len(children)
                if child_count == 2:
                    append(children[1])
                    append(children[0])
                else:
                    for index in range(child_count - 1, -1, -1):
                        append(children[index])
    except AttributeError:
        return None
    return clades


def _postorder_clades_direct(tree):
    try:
        root = tree.root
        root.clades
    except AttributeError:
        return None

    clades = []
    stack = [root]
    try:
        pop = stack.pop
        extend = stack.extend
        append = clades.append
        while stack:
            clade = pop()
            append(clade)
            extend(clade.clades)
    except AttributeError:
        return None
    clades.reverse()
    return clades


def resolve_polytomies(tree) -> None:
    """Resolve multifurcations by inserting zero-length branches.

    Mutates the tree in place, converting any node with >2 children
    into a series of bifurcating nodes connected by zero-length branches.
    """
    try:
        root = tree.root
        root.clades
    except AttributeError:
        from Bio.Phylo import Newick

        clades = tree.find_clades(order="postorder")
        for clade in clades:
            while len(clade.clades) > 2:
                child1 = clade.clades.pop()
                child2 = clade.clades.pop()
                new_internal = Newick.Clade(branch_length=0.0)
                new_internal.clades = [child1, child2]
                clade.clades.append(new_internal)
        return

    stack = [root]
    newick_clade = None
    try:
        pop = stack.pop
        extend = stack.extend
        while stack:
            clade = pop()
            children = clade.clades
            if children:
                extend(children)
            while len(children) > 2:
                if newick_clade is None:
                    from Bio.Phylo import Newick

                    newick_clade = Newick.Clade
                child1 = children.pop()
                child2 = children.pop()
                new_internal = newick_clade(branch_length=0.0)
                new_internal.clades = [child1, child2]
                children.append(new_internal)
    except AttributeError:
        from Bio.Phylo import Newick

        clades = tree.find_clades(order="postorder")
        for clade in clades:
            while len(clade.clades) > 2:
                child1 = clade.clades.pop()
                child2 = clade.clades.pop()
                new_internal = Newick.Clade(branch_length=0.0)
                new_internal.clades = [child1, child2]
                clade.clades.append(new_internal)


def fitch_downpass(
    tree, tip_states: dict[str, list[str]]
) -> tuple[dict[int, list[set[str]]], list[int]]:
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
    all_states_per_char: list[set[str]] = [set() for _ in range(n_chars)]
    for states in tip_states.values():
        for i, s in enumerate(states):
            if s not in wildcard:
                all_states_per_char[i].add(s)

    clades = _postorder_clades_direct(tree)
    if clades is None:
        clades = tree.find_clades(order="postorder")
    else:
        bitmask_result = _fitch_downpass_bitmask(
            clades,
            tip_states,
            all_states_per_char,
            n_chars,
            wildcard,
        )
        if bitmask_result is not None:
            return bitmask_result

    node_state_sets: dict[int, list[set[str]]] = {}
    scores = [0] * n_chars

    for clade in clades:
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
            if len(child_state_lists) == 2:
                left_states, right_states = child_state_lists
                for i in range(n_chars):
                    left = left_states[i]
                    right = right_states[i]
                    if not left and not right:
                        char_sets.append(set())
                        continue
                    if not left:
                        char_sets.append(right)
                        continue
                    if not right:
                        char_sets.append(left)
                        continue

                    intersection = left & right
                    if intersection:
                        char_sets.append(intersection)
                    else:
                        char_sets.append(left | right)
                        scores[i] += 1
            else:
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


def _fitch_downpass_bitmask(
    clades,
    tip_states: dict[str, list[str]],
    all_states_per_char: list[set[str]],
    n_chars: int,
    wildcard: set[str],
):
    states_by_char = [tuple(sorted(states)) for states in all_states_per_char]
    state_to_mask = [
        {state: 1 << index for index, state in enumerate(states)}
        for states in states_by_char
    ]
    all_masks = [
        (1 << len(states)) - 1 if states else 0
        for states in states_by_char
    ]
    tip_mask_cache = {}

    node_masks: dict[int, list[int]] = {}
    scores = [0] * n_chars

    for clade in clades:
        children = clade.clades
        cid = id(clade)
        if not children:
            states = tip_states[clade.name]
            states_key = tuple(states)
            masks = tip_mask_cache.get(states_key)
            if masks is None:
                masks = [
                    all_masks[i] if state in wildcard else state_to_mask[i][state]
                    for i, state in enumerate(states)
                ]
                tip_mask_cache[states_key] = masks
            node_masks[cid] = masks
            continue

        if len(children) > 2:
            return None
        if len(children) == 1:
            node_masks[cid] = node_masks[id(children[0])]
            continue

        left_masks = node_masks[id(children[0])]
        right_masks = node_masks[id(children[1])]
        char_masks = [0] * n_chars
        for i in range(n_chars):
            left = left_masks[i]
            right = right_masks[i]
            if not left and not right:
                continue
            if not left:
                char_masks[i] = right
                continue
            if not right:
                char_masks[i] = left
                continue

            intersection = left & right
            if intersection:
                char_masks[i] = intersection
            else:
                char_masks[i] = left | right
                scores[i] += 1
        node_masks[cid] = char_masks

    node_state_sets = {}
    use_mask_tuple_cache = (
        n_chars >= 16
        and len(tip_mask_cache) * 4 <= len(tip_states)
    )
    small_state_chars = [len(states) <= 6 for states in states_by_char]

    if all(small_state_chars):
        small_mask_set_lookups = [
            _build_small_mask_set_lookup(states) for states in states_by_char
        ]
        if not use_mask_tuple_cache:
            for cid, masks in node_masks.items():
                node_state_sets[cid] = [
                    small_mask_set_lookups[i][mask]
                    for i, mask in enumerate(masks)
                ]
            return node_state_sets, scores

        mask_tuple_cache = {}
        for cid, masks in node_masks.items():
            key = tuple(masks)
            cached_sets = mask_tuple_cache.get(key)
            if cached_sets is None:
                cached_sets = tuple(
                    small_mask_set_lookups[i][mask]
                    for i, mask in enumerate(masks)
                )
                mask_tuple_cache[key] = cached_sets
            node_state_sets[cid] = list(cached_sets)
        return node_state_sets, scores

    if not any(small_state_chars):
        mask_set_caches = [{} for _ in range(n_chars)]
        if not use_mask_tuple_cache:
            for cid, masks in node_masks.items():
                char_sets = []
                for i, mask in enumerate(masks):
                    cache = mask_set_caches[i]
                    state_set = cache.get(mask)
                    if state_set is None:
                        states = states_by_char[i]
                        state_set = frozenset(
                            states[index]
                            for index in range(len(states))
                            if mask & (1 << index)
                        )
                        cache[mask] = state_set
                    char_sets.append(state_set)
                node_state_sets[cid] = char_sets
            return node_state_sets, scores

        mask_tuple_cache = {}
        for cid, masks in node_masks.items():
            key = tuple(masks)
            cached_sets = mask_tuple_cache.get(key)
            if cached_sets is None:
                char_sets = []
                for i, mask in enumerate(masks):
                    cache = mask_set_caches[i]
                    state_set = cache.get(mask)
                    if state_set is None:
                        states = states_by_char[i]
                        state_set = frozenset(
                            states[index]
                            for index in range(len(states))
                            if mask & (1 << index)
                        )
                        cache[mask] = state_set
                    char_sets.append(state_set)
                cached_sets = tuple(char_sets)
                mask_tuple_cache[key] = cached_sets
            char_sets = list(cached_sets)
            node_state_sets[cid] = char_sets
        return node_state_sets, scores

    small_mask_set_lookups = [
        _build_small_mask_set_lookup(states) if is_small else None
        for states, is_small in zip(states_by_char, small_state_chars)
    ]
    mask_set_caches = [
        {} if lookup is None else None for lookup in small_mask_set_lookups
    ]
    if not use_mask_tuple_cache:
        for cid, masks in node_masks.items():
            char_sets = []
            for i, mask in enumerate(masks):
                lookup = small_mask_set_lookups[i]
                if lookup is not None:
                    state_set = lookup[mask]
                else:
                    cache = mask_set_caches[i]
                    state_set = cache.get(mask)
                    if state_set is None:
                        states = states_by_char[i]
                        state_set = frozenset(
                            states[index]
                            for index in range(len(states))
                            if mask & (1 << index)
                        )
                        cache[mask] = state_set
                char_sets.append(state_set)
            node_state_sets[cid] = char_sets
        return node_state_sets, scores

    mask_tuple_cache = {}
    for cid, masks in node_masks.items():
        key = tuple(masks)
        cached_sets = mask_tuple_cache.get(key)
        if cached_sets is None:
            char_sets = []
            for i, mask in enumerate(masks):
                lookup = small_mask_set_lookups[i]
                if lookup is not None:
                    state_set = lookup[mask]
                else:
                    cache = mask_set_caches[i]
                    state_set = cache.get(mask)
                    if state_set is None:
                        states = states_by_char[i]
                        state_set = frozenset(
                            states[index]
                            for index in range(len(states))
                            if mask & (1 << index)
                        )
                        cache[mask] = state_set
                char_sets.append(state_set)
            cached_sets = tuple(char_sets)
            mask_tuple_cache[key] = cached_sets
        char_sets = list(cached_sets)
        node_state_sets[cid] = char_sets

    return node_state_sets, scores


def _build_small_mask_set_lookup(states: tuple[str, ...]):
    if len(states) > 6:
        return None

    lookup = []
    for mask in range(1 << len(states)):
        lookup.append(_mask_to_state_set(states, mask))
    return lookup


def _mask_to_state_set(states: tuple[str, ...], mask: int):
    return frozenset(
        states[index]
        for index in range(len(states))
        if mask & (1 << index)
    )


def fitch_uppass_acctran(
    tree,
    node_state_sets: dict[int, list[set[str]]],
    parent_map: dict[int, object],
) -> dict[int, list[str]]:
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
    node_states: dict[int, list[str]] = {}

    clades = _preorder_clades_direct(tree)
    if clades is None:
        clades = tree.find_clades(order="preorder")

    for clade in clades:
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
    node_state_sets: dict[int, list[set[str]]],
    parent_map: dict[int, object],
) -> dict[int, list[str]]:
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
    node_states: dict[int, list[str]] = {}

    clades = _preorder_clades_direct(tree)
    if clades is None:
        clades = tree.find_clades(order="preorder")

    for clade in clades:
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
    node_states: dict[int, list[str]],
    parent_map: dict[int, object],
) -> dict[int, list[tuple[int, str, str]]]:
    """Detect character state changes on each branch.

    Compares each node's final states to its parent's final states.

    Returns:
        Dict[child_node_id -> List[(char_idx, old_state, new_state)]]
    """
    n_chars = len(next(iter(node_states.values())))
    changes: dict[int, list[tuple[int, str, str]]] = {}

    clades = _preorder_clades_direct(tree)
    if clades is None:
        clades = tree.find_clades(order="preorder")

    if n_chars == 1:
        root = tree.root
        for clade in clades:
            if clade == root:
                continue
            cid = id(clade)
            old = node_states[id(parent_map[cid])][0]
            new = node_states[cid][0]
            if old != new:
                changes[cid] = [(0, old, new)]
        return changes

    for clade in clades:
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
    branch_changes: dict[int, list[tuple[int, str, str]]],
    node_states: dict[int, list[str]],
    parent_map: dict[int, object],
) -> dict[int, list[tuple[int, str, str, str]]]:
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
    transition_counts: dict[tuple[int, str], int] = {}
    for changes in branch_changes.values():
        for char_idx, old, new in changes:
            key = (char_idx, new)
            transition_counts[key] = transition_counts.get(key, 0) + 1

    result: dict[int, list[tuple[int, str, str, str]]] = {}
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
    n_states_per_char: list[int],
    observed_per_char: list[int],
) -> tuple[list[float | None], float | None]:
    """Compute per-character and overall consistency index.

    CI_i = (n_states_i - 1) / observed_i
    CI_overall = sum(min_i) / sum(observed_i)

    Returns None for characters with 0 observed changes.
    """
    ci_per_char: list[float | None] = []
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
    tip_states_per_char: list[list[str]],
    observed_per_char: list[int],
) -> tuple[list[float | None], float | None]:
    """Compute per-character and overall retention index (Farris 1989).

    max_changes_i = n_taxa - f_max_i (count of most frequent state)
    min_changes_i = n_states_i - 1
    RI_i = (max_i - observed_i) / (max_i - min_i)

    Returns None for uninformative characters where max == min.
    """
    direct_result = _retention_index_ascii_single_char(
        tip_states_per_char,
        observed_per_char,
    )
    if direct_result is not None:
        return direct_result

    ri_per_char: list[float | None] = []
    sum_num = 0  # sum(max - observed)
    sum_den = 0  # sum(max - min)

    for states, observed in zip(tip_states_per_char, observed_per_char):
        counts = Counter(states)
        wildcard_count = counts.pop("?", 0) + counts.pop("-", 0)
        n_taxa = len(states) - wildcard_count
        if n_taxa == 0:
            ri_per_char.append(None)
            continue

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


def _retention_index_ascii_single_char(
    tip_states_per_char: list[list[str]],
    observed_per_char: list[int],
) -> tuple[list[float | None], float | None] | None:
    n_chars = min(len(tip_states_per_char), len(observed_per_char))
    if n_chars == 0:
        return [], None

    columns = tip_states_per_char[:n_chars]
    n_taxa = len(columns[0])
    if any(len(column) != n_taxa for column in columns):
        return None
    if n_chars * n_taxa < _RETENTION_INDEX_NUMPY_MIN_CELLS:
        return None

    try:
        data = "".join("".join(column) for column in columns).encode("ascii")
    except UnicodeEncodeError:
        return None
    if len(data) != n_chars * n_taxa:
        return None

    import numpy as np

    matrix = np.frombuffer(data, dtype=np.uint8).reshape(n_chars, n_taxa)
    symbols = np.unique(matrix)
    symbols = symbols[(symbols != ord("?")) & (symbols != ord("-"))]
    if symbols.size == 0:
        return [None] * n_chars, None

    symbol_counts = _ascii_symbol_counts_by_character(matrix, symbols)
    observed_taxa = symbol_counts.sum(axis=0)
    n_states = np.count_nonzero(symbol_counts, axis=0)
    f_max = np.max(symbol_counts, axis=0)
    max_changes = observed_taxa - f_max
    min_changes = n_states - 1
    observed = np.asarray(observed_per_char[:n_chars], dtype=np.int64)

    valid = (observed_taxa > 0) & (max_changes != min_changes)
    ri_per_char: list[float | None] = [None] * n_chars
    if valid.any():
        numerators = max_changes[valid] - observed[valid]
        denominators = max_changes[valid] - min_changes[valid]
        for index, value in zip(np.flatnonzero(valid), numerators / denominators):
            ri_per_char[int(index)] = float(value)

        sum_den = int(denominators.sum())
        ri_overall = (
            float(numerators.sum() / sum_den) if sum_den > 0 else None
        )
    else:
        ri_overall = None

    return ri_per_char, ri_overall


def _ascii_symbol_counts_by_character(matrix, symbols):
    import numpy as np

    if symbols.size >= 15:
        n_chars = matrix.shape[0]
        max_code = int(matrix.max()) + 1
        encoded = matrix.astype(np.int64)
        encoded += (np.arange(n_chars, dtype=np.int64) * max_code)[:, None]
        counts = np.bincount(
            encoded.ravel(),
            minlength=n_chars * max_code,
        ).reshape(n_chars, max_code)
        return counts[:, symbols].T

    return np.vstack(
        [np.count_nonzero(matrix == symbol, axis=1) for symbol in symbols]
    )
