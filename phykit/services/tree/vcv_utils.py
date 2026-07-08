"""Shared VCV matrix utilities for phylogenetic comparative methods."""

from __future__ import annotations

import os

from ...errors import PhykitUserError
from ...helpers.tree_paths import build_object_parent_map, path_from_root
from .base import Tree

_path_isabs = os.path.isabs


class _LazyPickle:
    _module = None

    def _load(self):
        module = self._module
        if module is None:
            import pickle as module

            self._module = module
        return module

    def __getattr__(self, name):
        value = getattr(self._load(), name)
        setattr(self, name, value)
        return value

    def dumps(self, *args, **kwargs):
        module = self._load()
        dumps = module.dumps
        self.dumps = dumps
        if "loads" not in self.__dict__:
            self.loads = module.loads

        return dumps(*args, **kwargs)

    def loads(self, *args, **kwargs):
        module = self._load()
        loads = module.loads
        self.loads = loads
        if "dumps" not in self.__dict__:
            self.dumps = module.dumps

        return loads(*args, **kwargs)


class _LazyNumpy:
    _module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = _np
            self._module = module

        value = getattr(module, name)
        setattr(self, name, value)
        return value


pickle = _LazyPickle()
np = _LazyNumpy()
_EIGVALSH = None


class _LazyPhylo:
    def read(self, *args, **kwargs):
        from Bio import Phylo as _Phylo

        self.read = _Phylo.read
        return self.read(*args, **kwargs)


Phylo = _LazyPhylo()


def eigvalsh(*args, **kwargs):
    global _EIGVALSH

    if _EIGVALSH is None:
        from scipy.linalg import eigvalsh as _eigvalsh

        _EIGVALSH = _eigvalsh

    return _EIGVALSH(*args, **kwargs)


def _get_tip_names(tree) -> list[str]:
    names = Tree.calculate_terminal_names_fast(tree)
    if names is not None:
        return names
    return [tip.name for tip in tree.get_terminals()]


def _build_vcv_matrix_fallback(tree, ordered_names: list[str]) -> np.ndarray:
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


def _build_vcv_matrix_fast(tree, ordered_names: list[str]) -> np.ndarray | None:
    try:
        root = tree.root
    except AttributeError:
        return None

    fast_vcv = _build_vcv_from_descendant_indices(
        ordered_names,
        root,
        lambda clade: clade.branch_length or 0.0,
    )
    if fast_vcv is not None:
        return fast_vcv

    try:
        terminals = tree.get_terminals()
        parent_map = build_object_parent_map(tree)
    except AttributeError:
        return None

    try:
        terminal_by_name = {terminal.name: terminal for terminal in terminals}
    except TypeError:
        return None
    if any(name not in terminal_by_name for name in ordered_names):
        return None

    return _build_vcv_from_branch_accumulation(
        ordered_names,
        terminal_by_name,
        root,
        parent_map,
        lambda clade: clade.branch_length or 0.0,
    )


def _build_vcv_from_branch_accumulation(
    ordered_names: list[str],
    terminal_by_name: dict[str, object],
    root,
    parent_map: dict[object, object],
    branch_length_getter,
) -> np.ndarray | None:
    n = len(ordered_names)
    vcv = np.zeros((n, n))
    clade_to_indices = {}
    asarray = np.asarray
    intp_dtype = np.intp

    for idx, name in enumerate(ordered_names):
        path = path_from_root(terminal_by_name[name], root, parent_map)
        if path is None:
            return None
        for clade in path[1:]:
            clade_to_indices.setdefault(clade, []).append(idx)

    for clade, indices in clade_to_indices.items():
        branch_length = branch_length_getter(clade)
        if branch_length == 0.0:
            continue
        if len(indices) == 1:
            vcv[indices[0], indices[0]] += branch_length
            continue
        idx = asarray(indices, dtype=intp_dtype)
        vcv[idx[:, None], idx] += branch_length

    return vcv


def _build_vcv_from_descendant_indices(
    ordered_names: list[str],
    root,
    branch_length_getter,
) -> np.ndarray | None:
    n = len(ordered_names)
    name_to_index = {name: idx for idx, name in enumerate(ordered_names)}
    if len(name_to_index) != n:
        return None

    try:
        root.clades
    except AttributeError:
        return None

    vcv = np.zeros((n, n))
    preorder = []
    stack = [root]
    try:
        pop = stack.pop
        append = stack.append
        extend = stack.extend
        append_preorder = preorder.append
        while stack:
            clade = pop()
            append_preorder(clade)
            children = clade.clades
            if children:
                if len(children) == 2:
                    append(children[1])
                    append(children[0])
                else:
                    extend(reversed(children))
    except AttributeError:
        return None

    descendant_indices = {}
    seen_ordered_tips = set()
    asarray = np.asarray
    intp_dtype = np.intp
    for clade in reversed(preorder):
        children = clade.clades
        if children:
            indices = []
            for child in children:
                child_indices = descendant_indices.get(child)
                if child_indices:
                    indices.extend(child_indices)
        else:
            idx = name_to_index.get(clade.name)
            if idx is None:
                indices = []
            elif idx in seen_ordered_tips:
                return None
            else:
                seen_ordered_tips.add(idx)
                indices = [idx]

        descendant_indices[clade] = indices
        if clade is root or not indices:
            continue
        branch_length = branch_length_getter(clade)
        if branch_length == 0.0:
            continue
        if len(indices) == 1:
            vcv[indices[0], indices[0]] += branch_length
            continue
        idx = asarray(indices, dtype=intp_dtype)
        vcv[idx[:, None], idx] += branch_length

    if len(seen_ordered_tips) != n:
        return None
    return vcv


def build_vcv_matrix(tree, ordered_names: list[str]) -> np.ndarray:
    fast_vcv = _build_vcv_matrix_fast(tree, ordered_names)
    if fast_vcv is not None:
        return fast_vcv
    return _build_vcv_matrix_fallback(tree, ordered_names)


def build_transformed_vcv_matrix(
    tree,
    ordered_names: list[str],
    branch_length_transform,
) -> np.ndarray:
    try:
        root = tree.root
    except AttributeError:
        return _build_transformed_vcv_matrix_fallback(
            tree, ordered_names, branch_length_transform
        )

    vcv = _build_vcv_from_descendant_indices(
        ordered_names,
        root,
        lambda clade: branch_length_transform(clade.branch_length or 0.0),
    )
    if vcv is not None:
        return vcv

    try:
        terminals = tree.get_terminals()
        parent_map = build_object_parent_map(tree)
    except AttributeError:
        return _build_transformed_vcv_matrix_fallback(
            tree, ordered_names, branch_length_transform
        )

    try:
        terminal_by_name = {terminal.name: terminal for terminal in terminals}
    except TypeError:
        return _build_transformed_vcv_matrix_fallback(
            tree, ordered_names, branch_length_transform
        )
    if any(name not in terminal_by_name for name in ordered_names):
        return _build_transformed_vcv_matrix_fallback(
            tree, ordered_names, branch_length_transform
        )

    vcv = _build_vcv_from_branch_accumulation(
        ordered_names,
        terminal_by_name,
        root,
        parent_map,
        lambda clade: branch_length_transform(clade.branch_length or 0.0),
    )
    if vcv is None:
        return _build_transformed_vcv_matrix_fallback(
            tree, ordered_names, branch_length_transform
        )
    return vcv


def _build_transformed_vcv_matrix_fallback(
    tree,
    ordered_names: list[str],
    branch_length_transform,
) -> np.ndarray:
    transformed = pickle.loads(pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL))
    for clade in transformed.find_clades():
        if clade.branch_length is not None and clade.branch_length > 0:
            clade.branch_length = branch_length_transform(clade.branch_length)
    return build_vcv_matrix(transformed, ordered_names)


def parse_gene_trees(gene_trees_path: str) -> list:
    """Parse a multi-Newick file, returning a list of Bio.Phylo trees.

    Supports two formats:
    - Inline Newick strings (lines starting with '(')
    - File paths (one per line, relative to the gene trees file's parent dir)

    Lines starting with '#' are treated as comments.
    """
    from io import StringIO
    from pathlib import Path

    source = Path(gene_trees_path)
    try:
        with source.open() as handle:
            cleaned = [
                stripped
                for line in handle
                if (stripped := line.strip())
                and stripped[0] != "#"
            ]
    except FileNotFoundError:
        raise PhykitUserError(
            [
                f"{gene_trees_path} corresponds to no such file or directory.",
                "Please check filename and pathing",
            ],
            code=2,
        )

    if not cleaned:
        raise PhykitUserError(
            [
                "Gene trees file is empty or contains only comments.",
                "Please provide at least one gene tree.",
            ],
            code=2,
        )

    trees = []
    parent_str = str(source.parent)
    parent_prefix = "" if parent_str == "." else parent_str + os.sep
    for line in cleaned:
        if line.startswith("("):
            trees.append(Phylo.read(StringIO(line), "newick"))
        else:
            tree_path = line if _path_isabs(line) else parent_prefix + line
            trees.append(Phylo.read(tree_path, "newick"))

    return trees


def _copy_prune_gene_tree_to_shared_taxa(gene_tree, shared_taxa):
    try:
        root = gene_tree.root
        root.clades
    except AttributeError:
        root = None

    if root is not None:
        stack = [root]
        try:
            pop = stack.pop
            extend = stack.extend
            contains = shared_taxa.__contains__
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    extend(children)
                elif not contains(clade.name):
                    break
            else:
                return gene_tree
        except AttributeError:
            pass
    else:
        try:
            if all(tip.name in shared_taxa for tip in gene_tree.get_terminals()):
                return gene_tree
        except AttributeError:
            pass

    gt_copy = pickle.loads(pickle.dumps(gene_tree, protocol=pickle.HIGHEST_PROTOCOL))

    try:
        root = gt_copy.root
        root.clades
    except AttributeError:
        root = None

    if root is not None:
        tips_to_remove = []
        target_ids = set()
        terminal_count = 0
        stack = [root]
        try:
            pop = stack.pop
            append = stack.append
            append_remove = tips_to_remove.append
            add_target = target_ids.add
            contains = shared_taxa.__contains__
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    child_count = len(children)
                    if child_count == 2:
                        append(children[1])
                        append(children[0])
                    else:
                        for index in range(child_count - 1, -1, -1):
                            append(children[index])
                else:
                    terminal_count += 1
                    if not contains(clade.name):
                        append_remove(clade)
                        add_target(id(clade))
        except AttributeError:
            tips_to_remove = None

        if tips_to_remove is not None:
            if (
                tips_to_remove
                and len(tips_to_remove) < terminal_count
                and Tree._prune_terminal_objects_batch_standard_tree(
                    gt_copy,
                    target_ids,
                )
            ):
                return gt_copy
            for tip in tips_to_remove:
                gt_copy.prune(tip)
            return gt_copy

    tips_to_remove = [
        tip for tip in gt_copy.get_terminals() if tip.name not in shared_taxa
    ]
    for tip in tips_to_remove:
        gt_copy.prune(tip)
    return gt_copy


def _build_pruned_subset_vcv_matrix(gene_tree, ordered_names: list[str]):
    """Build the VCV equivalent of pruning a gene tree to ordered_names.

    Branches ancestral to every retained taxon are skipped because pruning
    extra taxa can collapse those stems above the retained subtree root.
    """
    try:
        root = gene_tree.root
        root.clades
    except AttributeError:
        return None

    n = len(ordered_names)
    name_to_index = {name: idx for idx, name in enumerate(ordered_names)}
    if len(name_to_index) != n:
        return None

    preorder = []
    stack = [root]
    try:
        pop = stack.pop
        append = stack.append
        extend = stack.extend
        append_preorder = preorder.append
        while stack:
            clade = pop()
            append_preorder(clade)
            children = clade.clades
            if children:
                if len(children) == 2:
                    append(children[1])
                    append(children[0])
                else:
                    extend(reversed(children))
    except AttributeError:
        return None

    vcv = np.zeros((n, n))
    descendant_indices = {}
    seen_ordered_tips = set()
    asarray = np.asarray
    intp_dtype = np.intp
    for clade in reversed(preorder):
        children = clade.clades
        if children:
            indices = []
            for child in children:
                child_indices = descendant_indices.get(child)
                if child_indices:
                    indices.extend(child_indices)
        else:
            idx = name_to_index.get(clade.name)
            if idx is None:
                indices = []
            elif idx in seen_ordered_tips:
                return None
            else:
                seen_ordered_tips.add(idx)
                indices = [idx]

        descendant_indices[clade] = indices
        if clade is root or not indices or len(indices) == n:
            continue

        branch_length = clade.branch_length or 0.0
        if branch_length == 0.0:
            continue
        if len(indices) == 1:
            vcv[indices[0], indices[0]] += branch_length
            continue
        idx = asarray(indices, dtype=intp_dtype)
        vcv[idx[:, None], idx] += branch_length

    if len(seen_ordered_tips) != n:
        return None
    return vcv


def _gene_tree_has_missing_branch_lengths(gene_tree) -> bool:
    try:
        root = gene_tree.root
        root.clades
    except AttributeError:
        root = None

    if root is not None:
        stack = [root]
        pop = stack.pop
        extend = stack.extend
        try:
            while stack:
                clade = pop()
                if clade is not root and clade.branch_length is None:
                    return True
                children = clade.clades
                if children:
                    extend(children)
            return False
        except AttributeError:
            pass

    for clade in gene_tree.find_clades():
        if clade.branch_length is None and clade != gene_tree.root:
            return True
    return False


def _get_tip_names_and_missing_branch_lengths(gene_tree) -> tuple[list[str], bool]:
    try:
        root = gene_tree.root
        root.clades
    except AttributeError:
        root = None

    if root is not None:
        tips = []
        missing_branch_length = False
        stack = [root]
        try:
            pop = stack.pop
            append = stack.append
            tips_append = tips.append
            while stack:
                clade = pop()
                if clade is not root and clade.branch_length is None:
                    missing_branch_length = True
                children = clade.clades
                child_count = len(children)
                if child_count == 2:
                    append(children[1])
                    append(children[0])
                elif child_count:
                    for idx in range(child_count - 1, -1, -1):
                        append(children[idx])
                else:
                    tips_append(clade.name)
            return tips, missing_branch_length
        except AttributeError:
            pass

    return _get_tip_names(gene_tree), _gene_tree_has_missing_branch_lengths(gene_tree)


def build_discordance_vcv(
    species_tree, gene_trees: list, ordered_names: list[str]
) -> tuple[np.ndarray, dict]:
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
    species_tips = set(_get_tip_names(species_tree))
    ordered_set = set(ordered_names)

    # Find taxa shared across all gene trees
    shared_taxa = species_tips & ordered_set
    gene_tree_tip_sets = []
    gene_tree_missing_branch_lengths = []
    for gt in gene_trees:
        gt_tip_names, has_missing_branch_lengths = (
            _get_tip_names_and_missing_branch_lengths(gt)
        )
        gt_tips = set(gt_tip_names)
        gene_tree_tip_sets.append(gt_tips)
        gene_tree_missing_branch_lengths.append(has_missing_branch_lengths)
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

    for gt, gt_tips, has_missing_branch_lengths in zip(
        gene_trees,
        gene_tree_tip_sets,
        gene_tree_missing_branch_lengths,
    ):
        if has_missing_branch_lengths:
            raise PhykitUserError(
                [
                    "Gene tree contains branches without lengths.",
                    "All gene trees must have branch lengths for "
                    "discordance-aware VCV computation.",
                ],
                code=2,
            )

        if gt_tips == shared_taxa:
            vcv_g = build_vcv_matrix(gt, shared_ordered)
        else:
            vcv_g = _build_pruned_subset_vcv_matrix(gt, shared_ordered)
            if vcv_g is None:
                gt_copy = _copy_prune_gene_tree_to_shared_taxa(gt, shared_taxa)
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


def _nearest_psd(matrix: np.ndarray) -> tuple[np.ndarray, bool, float]:
    """Clip negative eigenvalues to ensure positive semi-definiteness.

    Returns (corrected_matrix, was_corrected, min_eigenvalue_before_correction).
    """
    try:
        min_eval = float(
            eigvalsh(
                matrix,
                subset_by_index=(0, 0),
                check_finite=False,
            )[0]
        )
    except (TypeError, ValueError, np.linalg.LinAlgError):
        eigvals, eigvecs = np.linalg.eigh(matrix)
        return _nearest_psd_from_eigendecomposition(matrix, eigvals, eigvecs)

    if min_eval >= 0:
        return matrix, False, min_eval

    eigvals, eigvecs = np.linalg.eigh(matrix)
    return _nearest_psd_from_eigendecomposition(matrix, eigvals, eigvecs)


def _nearest_psd_from_eigendecomposition(
    matrix: np.ndarray, eigvals: np.ndarray, eigvecs: np.ndarray
) -> tuple[np.ndarray, bool, float]:
    min_eval = float(eigvals.min())

    if min_eval >= 0:
        return matrix, False, min_eval

    eigvals_clipped = np.maximum(eigvals, 0)
    corrected = (eigvecs * eigvals_clipped) @ eigvecs.T
    # Ensure symmetry (numerical precision)
    corrected = (corrected + corrected.T) / 2.0

    return corrected, True, min_eval
