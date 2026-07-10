"""
Shared utilities for quartet concordance factor computation and ASTRAL parsing.

Provides gene concordance factor (gCF/gDF1/gDF2) computation via the
four-group bipartition decomposition, and parsing of ASTRAL -t 2 or
wASTRAL --support 3 q1/q2/q3 annotations from Newick node labels.
"""
from collections import Counter

from ..errors import PhykitUserError


def canonical_split(tips: frozenset, all_taxa: frozenset) -> frozenset:
    """Normalize a bipartition to a canonical frozenset-of-frozensets."""
    complement = all_taxa - tips
    return frozenset((tips, complement))


def _collect_clade_tip_sets(tree, allowed_taxa: frozenset = None) -> dict[int, frozenset]:
    direct_clade_tips = _collect_clade_tip_sets_direct(tree, allowed_taxa)
    if direct_clade_tips is not None:
        return direct_clade_tips

    clade_tips: dict[int, frozenset] = {}

    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            if allowed_taxa is None or clade.name in allowed_taxa:
                clade_tips[id(clade)] = frozenset({clade.name})
            else:
                clade_tips[id(clade)] = frozenset()
        else:
            tips = frozenset()
            for child in clade.clades:
                tips = tips | clade_tips.get(id(child), frozenset())
            clade_tips[id(clade)] = tips

    return clade_tips


def _collect_clade_tip_sets_direct(
    tree, allowed_taxa: frozenset = None
) -> dict[int, frozenset] | None:
    try:
        root = tree.root
        root.clades
    except AttributeError:
        return None

    preorder = []
    stack = [root]
    append = preorder.append
    pop = stack.pop
    extend = stack.extend
    try:
        while stack:
            clade = pop()
            append(clade)
            children = clade.clades
            if children:
                extend(children)
    except AttributeError:
        return None

    clade_tips: dict[int, frozenset] = {}
    try:
        for clade in reversed(preorder):
            children = clade.clades
            if children:
                child_count = len(children)
                if child_count == 2:
                    tips = (
                        clade_tips[id(children[0])]
                        | clade_tips[id(children[1])]
                    )
                elif child_count == 1:
                    tips = clade_tips[id(children[0])]
                else:
                    tips = frozenset().union(
                        *(clade_tips[id(child)] for child in children)
                    )
                clade_tips[id(clade)] = tips
            elif allowed_taxa is None or clade.name in allowed_taxa:
                clade_tips[id(clade)] = frozenset({clade.name})
            else:
                clade_tips[id(clade)] = frozenset()
    except (AttributeError, KeyError, TypeError):
        return None
    return clade_tips


def _preorder_clades(tree):
    direct_clades = _preorder_clades_direct(tree)
    if direct_clades is not None:
        return direct_clades
    return list(tree.find_clades(order="preorder"))


def _preorder_clades_direct(tree) -> list | None:
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
        clades_append = clades.append
        while stack:
            clade = pop()
            clades_append(clade)
            children = clade.clades
            child_count = len(children)
            if child_count == 2:
                append(children[1])
                append(children[0])
            elif child_count:
                for idx in range(child_count - 1, -1, -1):
                    append(children[idx])
    except AttributeError:
        return None
    return clades


def compute_gcf_per_node(
    species_tree, gene_trees: list
) -> dict[int, tuple[float, float, float, int, int, int]]:
    """Compute (gCF, gDF1, gDF2, concordant, disc1, disc2) per internal node.

    Uses the four-group decomposition (C1, C2, S, D) around each internal
    branch to identify the concordant bipartition and the two NNI alternatives.

    Returns dict mapping clade id -> (gCF, gDF1, gDF2, n_conc, n_d1, n_d2).
    """
    species_clade_tips = _collect_clade_tip_sets(species_tree)
    try:
        all_taxa = species_clade_tips[id(species_tree.root)]
    except (AttributeError, KeyError):
        all_taxa = frozenset(t.name for t in species_tree.get_terminals())
    species_preorder = _preorder_clades(species_tree)

    # Extract bipartitions from all gene trees and count how many trees support
    # each split. Each gene tree contributes at most once per split.
    gt_split_counts = Counter()
    for gt in gene_trees:
        gt_clade_tips = _collect_clade_tip_sets(gt, all_taxa)
        shared = gt_clade_tips.get(id(gt.root), frozenset())
        if len(shared) < 4:
            continue
        splits = set()
        for clade in _preorder_clades(gt):
            if not clade.clades:
                continue
            tips = gt_clade_tips.get(id(clade), frozenset())
            if len(tips) <= 1 or tips == shared:
                continue
            splits.add(canonical_split(tips, shared))
        gt_split_counts.update(splits)

    result = {}
    for clade in species_preorder:
        if not clade.clades or clade == species_tree.root:
            continue

        children = clade.clades
        if len(children) < 2:
            continue

        C1 = species_clade_tips.get(id(children[0]), frozenset())
        C2 = species_clade_tips.get(id(children[1]), frozenset())
        # Handle polytomies: merge extra children into C2
        for c in children[2:]:
            C2 = C2 | species_clade_tips.get(id(c), frozenset())

        remaining = all_taxa - C1 - C2
        if not remaining:
            continue

        concordant_bp = canonical_split(C1 | C2, all_taxa)
        nni1_bp = canonical_split(remaining | C2, all_taxa)
        nni2_bp = canonical_split(C1 | remaining, all_taxa)

        conc = gt_split_counts[concordant_bp]
        d1 = gt_split_counts[nni1_bp]
        d2 = gt_split_counts[nni2_bp]

        total = conc + d1 + d2
        if total > 0:
            gcf = conc / total
            gdf1 = d1 / total
            gdf2 = d2 / total
        else:
            gcf, gdf1, gdf2 = 1.0, 0.0, 0.0

        result[id(clade)] = (gcf, gdf1, gdf2, conc, d1, d2)

    return result


def parse_astral_annotations(
    tree,
) -> dict[int, tuple[float, float, float]]:
    """Parse q1/q2/q3 annotations from ASTRAL/wASTRAL Newick node labels.

    Supports ASTRAL -t 2 and wASTRAL --support 3 output formats:
      '[q1=0.5;q2=0.3;q3=0.2;f1=50;...]'
      '[CULength=2.6;f1=108;...;q1=0.96;q2=0.03;q3=0.01]'
      'q1=0.5;q2=0.3;q3=0.2'

    Returns dict mapping clade id -> (q1, q2, q3).
    Only nodes with valid q1/q2/q3 are included.
    """
    result = {}
    for clade in _preorder_clades(tree):
        if clade.is_terminal():
            continue
        label = clade.name or clade.comment or ""
        qs = _parse_qs_from_label(str(label))
        if qs is not None:
            result[id(clade)] = qs
    return result


def parse_astral_branch_info(
    tree,
) -> dict[int, dict[str, float]]:
    """Parse f1 (concordant count) and pp1 (LPP) from ASTRAL/wASTRAL labels.

    Returns dict mapping clade id -> {"f1": ..., "pp1": ...}.
    Keys are only present if the value was found in the annotation.
    """
    result = {}
    for clade in _preorder_clades(tree):
        if clade.is_terminal():
            continue
        label = clade.name or clade.comment or ""
        info = _parse_branch_info_from_label(str(label))
        if info:
            result[id(clade)] = info
    return result


def _parse_branch_info_from_label(label: str) -> dict[str, float]:
    """Extract f1, pp1 (and optionally f2, f3, pp2, pp3) from an ASTRAL label."""
    s = label.strip("'\"").strip("[]")
    info = {}
    for part in s.split(";"):
        if "=" not in part:
            continue
        key, val = part.split("=", 1)
        key = key.strip()
        try:
            if key in ("f1", "f2", "f3", "pp1", "pp2", "pp3"):
                info[key] = float(val)
        except ValueError:
            continue
    return info


def _parse_qs_from_label(label: str) -> tuple[float, float, float] | None:
    """Extract q1, q2, q3 from an ASTRAL node label string."""
    s = label.strip("'\"").strip("[]")
    q1 = q2 = q3 = None
    for part in s.split(";"):
        if "=" not in part:
            continue
        key, val = part.split("=", 1)
        key = key.strip()
        try:
            if key == "q1":
                q1 = float(val)
            elif key == "q2":
                q2 = float(val)
            elif key == "q3":
                q3 = float(val)
        except ValueError:
            continue
    if q1 is not None and q2 is not None and q3 is not None:
        return (q1, q2, q3)
    return None
