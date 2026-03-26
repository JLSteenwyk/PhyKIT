"""
Shared utilities for quartet concordance factor computation and ASTRAL parsing.

Provides gene concordance factor (gCF/gDF1/gDF2) computation via the
four-group bipartition decomposition, and parsing of ASTRAL -t 2 or
wASTRAL --support 3 q1/q2/q3 annotations from Newick node labels.
"""
from typing import Dict, List, Optional, Tuple

from ..errors import PhykitUserError


def canonical_split(tips: frozenset, all_taxa: frozenset) -> frozenset:
    """Normalize a bipartition to a canonical frozenset-of-frozensets."""
    complement = all_taxa - tips
    return frozenset([tips, complement])


def compute_gcf_per_node(
    species_tree, gene_trees: list
) -> Dict[int, Tuple[float, float, float, int, int, int]]:
    """Compute (gCF, gDF1, gDF2, concordant, disc1, disc2) per internal node.

    Uses the four-group decomposition (C1, C2, S, D) around each internal
    branch to identify the concordant bipartition and the two NNI alternatives.

    Returns dict mapping clade id -> (gCF, gDF1, gDF2, n_conc, n_d1, n_d2).
    """
    all_taxa = frozenset(t.name for t in species_tree.get_terminals())

    # Build parent map
    parent_map = {}
    for clade in species_tree.find_clades(order="preorder"):
        for child in clade.clades:
            parent_map[id(child)] = clade

    # Extract bipartitions from all gene trees
    gt_splits = []
    for gt in gene_trees:
        gt_taxa = frozenset(t.name for t in gt.get_terminals())
        shared = gt_taxa & all_taxa
        if len(shared) < 4:
            gt_splits.append(set())
            continue
        splits = set()
        for clade in gt.get_nonterminals():
            tips = frozenset(t.name for t in clade.get_terminals()) & shared
            if len(tips) <= 1 or tips == shared:
                continue
            splits.add(canonical_split(tips, shared))
        gt_splits.append(splits)

    result = {}
    for clade in species_tree.find_clades(order="preorder"):
        if clade.is_terminal() or clade == species_tree.root:
            continue

        children = clade.clades
        if len(children) < 2:
            continue

        C1 = frozenset(t.name for t in children[0].get_terminals())
        C2 = frozenset(t.name for t in children[1].get_terminals())
        # Handle polytomies: merge extra children into C2
        for c in children[2:]:
            C2 = C2 | frozenset(t.name for t in c.get_terminals())

        remaining = all_taxa - C1 - C2
        if not remaining:
            continue

        concordant_bp = canonical_split(C1 | C2, all_taxa)
        nni1_bp = canonical_split(remaining | C2, all_taxa)
        nni2_bp = canonical_split(C1 | remaining, all_taxa)

        conc = sum(1 for s in gt_splits if concordant_bp in s)
        d1 = sum(1 for s in gt_splits if nni1_bp in s)
        d2 = sum(1 for s in gt_splits if nni2_bp in s)

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
) -> Dict[int, Tuple[float, float, float]]:
    """Parse q1/q2/q3 annotations from ASTRAL/wASTRAL Newick node labels.

    Supports ASTRAL -t 2 and wASTRAL --support 3 output formats:
      '[q1=0.5;q2=0.3;q3=0.2;f1=50;...]'
      '[CULength=2.6;f1=108;...;q1=0.96;q2=0.03;q3=0.01]'
      'q1=0.5;q2=0.3;q3=0.2'

    Returns dict mapping clade id -> (q1, q2, q3).
    Only nodes with valid q1/q2/q3 are included.
    """
    result = {}
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            continue
        label = clade.name or clade.comment or ""
        qs = _parse_qs_from_label(str(label))
        if qs is not None:
            result[id(clade)] = qs
    return result


def parse_astral_branch_info(
    tree,
) -> Dict[int, Dict[str, float]]:
    """Parse f1 (concordant count) and pp1 (LPP) from ASTRAL/wASTRAL labels.

    Returns dict mapping clade id -> {"f1": ..., "pp1": ...}.
    Keys are only present if the value was found in the annotation.
    """
    result = {}
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            continue
        label = clade.name or clade.comment or ""
        info = _parse_branch_info_from_label(str(label))
        if info:
            result[id(clade)] = info
    return result


def _parse_branch_info_from_label(label: str) -> Dict[str, float]:
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


def _parse_qs_from_label(label: str) -> Optional[Tuple[float, float, float]]:
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
