"""Shared VCV matrix utilities for phylogenetic comparative methods."""

import pickle
import sys
from io import StringIO
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from Bio import Phylo

from ...errors import PhykitUserError


def build_vcv_matrix(tree, ordered_names: List[str]) -> np.ndarray:
    """Build variance-covariance matrix from a single phylogenetic tree.

    VCV[i,i] = root_to_tip_distance(i)
    VCV[i,j] = (root_to_tip(i) + root_to_tip(j) - pairwise_dist(i,j)) / 2
    """
    n = len(ordered_names)
    vcv = np.zeros((n, n))

    root_to_tip = {}
    for name in ordered_names:
        root_to_tip[name] = tree.distance(tree.root, name)

    for i in range(n):
        for j in range(i, n):
            if i == j:
                vcv[i, j] = root_to_tip[ordered_names[i]]
            else:
                d_ij = tree.distance(ordered_names[i], ordered_names[j])
                shared_path = (
                    root_to_tip[ordered_names[i]]
                    + root_to_tip[ordered_names[j]]
                    - d_ij
                ) / 2.0
                vcv[i, j] = shared_path
                vcv[j, i] = shared_path

    return vcv


def parse_gene_trees(gene_trees_path: str) -> List:
    """Parse a multi-Newick file, returning a list of Bio.Phylo trees.

    Supports two formats:
    - Inline Newick strings (lines starting with '(')
    - File paths (one per line, relative to the gene trees file's parent dir)

    Lines starting with '#' are treated as comments.
    """
    try:
        lines = Path(gene_trees_path).read_text().splitlines()
    except FileNotFoundError:
        raise PhykitUserError(
            [
                f"{gene_trees_path} corresponds to no such file or directory.",
                "Please check filename and pathing",
            ],
            code=2,
        )

    cleaned = [l.strip() for l in lines if l.strip() and not l.strip().startswith("#")]

    if not cleaned:
        raise PhykitUserError(
            [
                "Gene trees file is empty or contains only comments.",
                "Please provide at least one gene tree.",
            ],
            code=2,
        )

    trees = []
    for line in cleaned:
        if line.startswith("("):
            trees.append(Phylo.read(StringIO(line), "newick"))
        else:
            tree_path = Path(gene_trees_path).parent / line
            trees.append(Phylo.read(str(tree_path), "newick"))

    return trees


def build_discordance_vcv(
    species_tree, gene_trees: List, ordered_names: List[str]
) -> Tuple[np.ndarray, dict]:
    """Build discordance-aware VCV by averaging per-gene-tree VCVs.

    1. Find shared taxa (intersection of species tree tips, gene tree tips,
       and ordered_names)
    2. Prune gene trees to shared taxa
    3. Build VCV from each gene tree
    4. Average and correct to nearest PSD

    Returns (vcv_matrix, metadata_dict) where metadata includes:
      n_gene_trees, n_shared_taxa, psd_corrected, min_eigenvalue_pre_correction
    """
    # Get species tree tip names
    species_tips = set(t.name for t in species_tree.get_terminals())
    ordered_set = set(ordered_names)

    # Find taxa shared across all gene trees
    shared_taxa = species_tips & ordered_set
    for gt in gene_trees:
        gt_tips = set(t.name for t in gt.get_terminals())
        shared_taxa = shared_taxa & gt_tips

    if len(shared_taxa) < 3:
        raise PhykitUserError(
            [
                f"Only {len(shared_taxa)} taxa shared across species tree "
                f"and all gene trees.",
                "At least 3 shared taxa are required.",
            ],
            code=2,
        )

    # Use sorted shared taxa for consistent ordering
    shared_ordered = sorted(shared_taxa)

    # Build VCV from each gene tree
    n = len(shared_ordered)
    vcv_sum = np.zeros((n, n))
    n_used = 0

    for gt in gene_trees:
        # Validate that gene tree has branch lengths
        for clade in gt.find_clades():
            if clade.branch_length is None and clade != gt.root:
                raise PhykitUserError(
                    [
                        "Gene tree contains branches without lengths.",
                        "All gene trees must have branch lengths for "
                        "discordance-aware VCV computation.",
                    ],
                    code=2,
                )

        # Prune gene tree to shared taxa
        gt_copy = pickle.loads(pickle.dumps(gt, protocol=pickle.HIGHEST_PROTOCOL))
        gt_tips = set(t.name for t in gt_copy.get_terminals())
        tips_to_remove = gt_tips - shared_taxa
        if tips_to_remove:
            for tip_name in tips_to_remove:
                gt_copy.prune(tip_name)

        vcv_g = build_vcv_matrix(gt_copy, shared_ordered)
        vcv_sum += vcv_g
        n_used += 1

    # Average
    vcv_avg = vcv_sum / n_used

    # Ensure PSD
    vcv_psd, was_corrected, min_eval = _nearest_psd(vcv_avg)

    # If ordered_names differs from shared_ordered, we need to return a VCV
    # that matches ordered_names (which may be a superset). But per the plan,
    # we only use shared taxa — the caller must handle subsetting.
    # Re-map to the original ordered_names ordering if they match shared_ordered.
    # If ordered_names has taxa not in shared_taxa, we need to subset ordered_names.
    # The caller is responsible for using the returned shared_ordered.

    metadata = {
        "n_gene_trees": n_used,
        "n_shared_taxa": len(shared_ordered),
        "shared_taxa": shared_ordered,
        "psd_corrected": was_corrected,
        "min_eigenvalue_pre_correction": float(min_eval),
    }

    return vcv_psd, metadata


def _nearest_psd(matrix: np.ndarray) -> Tuple[np.ndarray, bool, float]:
    """Clip negative eigenvalues to ensure positive semi-definiteness.

    Returns (corrected_matrix, was_corrected, min_eigenvalue_before_correction).
    """
    eigvals, eigvecs = np.linalg.eigh(matrix)
    min_eval = float(np.min(eigvals))

    if min_eval >= 0:
        return matrix, False, min_eval

    eigvals_clipped = np.maximum(eigvals, 0)
    corrected = eigvecs @ np.diag(eigvals_clipped) @ eigvecs.T
    # Ensure symmetry (numerical precision)
    corrected = (corrected + corrected.T) / 2.0

    return corrected, True, min_eval
