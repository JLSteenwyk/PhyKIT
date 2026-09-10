"""Strict, unrooted four-group classification for genomic concordance maps."""

from collections import defaultdict
import math

from Bio import Phylo

from ..errors import PhykitUserError
from .quartet_utils import _collect_clade_tip_sets, canonical_split


CLASSES = (
    "topology_1", "topology_2", "topology_3", "unresolved",
    "insufficient_sampling", "incompatible_groups",
)


def fail(message):
    raise PhykitUserError([message], code=2)


def validate_groups(groups):
    if len(groups) != 4 or any(not taxa for taxa in groups.values()):
        fail("Exactly four named, nonempty taxon groups are required.")
    seen = set()
    for name, taxa in groups.items():
        if not isinstance(name, str) or not name.strip():
            fail("Group names must be nonempty strings.")
        if not isinstance(taxa, (list, tuple, set, frozenset)):
            fail(f"Group {name} must contain a list of taxon names.")
        if any(not isinstance(t, str) or not t.strip() for t in taxa):
            fail(f"Group {name} contains an invalid taxon name.")
        if len(set(taxa)) != len(taxa) or seen.intersection(taxa):
            fail(f"Taxa must be unique and groups disjoint: {name}.")
        seen.update(taxa)
    return {name: frozenset(taxa) for name, taxa in groups.items()}


def read_tree(path):
    try:
        tree = Phylo.read(path, "newick")
    except (OSError, ValueError, IndexError) as exc:
        fail(f"Cannot read exactly one Newick tree from {path}: {exc}")
    names = [tip.name for tip in tree.get_terminals()]
    if any(not name for name in names) or len(set(names)) != len(names):
        fail(f"Tree {path} must have unique, nonempty tip names.")
    return tree


def classify_tree(tree, groups, min_support=None, support_scale=100,
                  missing_support="collapse"):
    """Require each sampled group to form an unrooted clan, then resolve it.

    Repeated representations of an induced edge (degree-two roots or paths
    through excluded taxa) use the minimum available numeric support. When
    filtering, any absent label on that path follows missing_support.
    Branch lengths never enter classification.
    """
    groups = validate_groups(groups)
    if support_scale not in (1, 100):
        fail("Support scale must be 1 or 100.")
    if min_support is not None and (
        not math.isfinite(min_support) or not 0 <= min_support <= support_scale
    ):
        fail("Minimum support must be finite and within the selected scale.")
    if missing_support not in ("collapse", "keep", "error"):
        fail("Missing-support policy must be collapse, keep, or error.")
    names = [tip.name for tip in tree.get_terminals()]
    if any(not name for name in names) or len(set(names)) != len(names):
        fail("Gene trees must have unique, nonempty tip names.")
    wanted = frozenset().union(*groups.values())
    sampled = [taxa.intersection(names) for taxa in groups.values()]
    counts = dict(zip(groups, map(len, sampled)))
    result = {
        "classification": "insufficient_sampling", "sampled_counts": counts,
        "sampling_status": "complete" if all(sampled) else "missing_groups",
        "ignored_taxa": sorted(set(names) - wanted), "resolution_support": None,
        "resolution_support_missing": False,
    }
    if not all(sampled):
        return result
    all_taxa = frozenset().union(*sampled)
    clade_tips = _collect_clade_tip_sets(tree, all_taxa)
    evidence = defaultdict(list)
    for clade in tree.find_clades():
        if clade is tree.root or clade.is_terminal():
            continue
        side = clade_tips[id(clade)]
        if not side or side == all_taxa:
            continue
        value = clade.confidence
        if value is not None and (
            not math.isfinite(value) or not 0 <= value <= support_scale
        ):
            fail(f"Support {value} is outside the selected 0-{support_scale} scale.")
        evidence[canonical_split(side, all_taxa)].append(value)
    retained = {}
    for split, values in evidence.items():
        known = [v for v in values if v is not None]
        support = min(known) if known else None
        absent = len(known) != len(values)
        if min_support is not None:
            if absent and missing_support == "error":
                fail("Missing branch support encountered with --missing-support error.")
            if absent and missing_support == "collapse":
                continue
            if support is not None and support < min_support:
                continue
        retained[split] = (support, absent)

    # Missing a clan-defining edge may mean a polytomy, not contradictory data.
    uncertain_group = False
    for group in sampled:
        if len(group) == 1 or canonical_split(group, all_taxa) in retained:
            continue
        outside = all_taxa - group
        for split in retained:
            left, right = tuple(split)
            if left & group and right & group and left & outside and right & outside:
                result["classification"] = "incompatible_groups"
                return result
        uncertain_group = True
    result["classification"] = "unresolved"
    if uncertain_group:
        return result
    for index, partner in enumerate((1, 2, 3), 1):
        split = canonical_split(sampled[0] | sampled[partner], all_taxa)
        if split in retained:
            result["classification"] = f"topology_{index}"
            result["resolution_support"], result["resolution_support_missing"] = retained[split]
            break
    return result


def groups_from_branch(tree, branch_taxa):
    """Derive four groups around an edge identified by either complete side.

    Suppress artificial degree-two vertices so selection is root invariant.
    Group A/B lie on the selected side, C/D on its complement; components
    within a side are sorted lexicographically by their complete tip lists.
    """
    graph = defaultdict(set)
    terminals = {c: c.name for c in tree.get_terminals()}
    if len(set(terminals.values())) != len(terminals) or any(
        not name for name in terminals.values()
    ):
        fail("Reference tips must have unique, nonempty names.")
    for clade in tree.find_clades():
        for child in clade.clades:
            graph[clade].add(child)
            graph[child].add(clade)
    for node in list(graph):
        if node in terminals:
            continue
        neighbors = graph[node]
        if len(neighbors) == 1:
            neighbor = next(iter(neighbors))
            graph[neighbor].remove(node)
            del graph[node]
        elif len(neighbors) == 2:
            left, right = neighbors
            graph[left].remove(node)
            graph[right].remove(node)
            graph[left].add(right)
            graph[right].add(left)
            del graph[node]

    def component(start, blocked):
        found = set()
        stack = [(start, blocked)]
        while stack:
            node, parent = stack.pop()
            if node in terminals:
                found.add(terminals[node])
            stack.extend((n, node) for n in graph[node] if n is not parent)
        return frozenset(found)

    target = frozenset(branch_taxa)
    for left, neighbors in graph.items():
        for right in neighbors:
            if component(left, right) != target:
                continue
            if len(graph[left]) != 3 or len(graph[right]) != 3:
                fail("Selected branch must have two degree-three internal endpoints.")
            sides = []
            for node, other in ((left, right), (right, left)):
                parts = [component(n, node) for n in graph[node] if n is not other]
                sides.extend(sorted(parts, key=lambda x: tuple(sorted(x))))
            return validate_groups(dict(zip(("A", "B", "C", "D"), sides)))
    fail("Branch taxa do not identify a complete side of a reference-tree edge.")
