# Spectral Discordance Decomposition Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a `spectral_discordance` command that decomposes gene-tree space via PCA on a bipartition matrix, with spectral clustering and eigengap auto-K detection.

**Architecture:** Each gene tree is encoded as a binary/weighted vector over the union of all observed bipartitions. PCA via SVD gives ordination scores plus directly interpretable bipartition loading vectors. Spectral clustering on the normalized graph Laplacian with eigengap heuristic auto-detects groups of genes sharing alternative topologies.

**Tech Stack:** numpy (SVD, eigendecomposition), scipy (cdist), matplotlib (scatter + eigengap plots), Bio.Phylo (tree I/O)

---

## Reference Files (read-only)

- `phykit/services/tree/discordance_asymmetry.py:318-412` — `_parse_gene_trees`, `_canonical_split`, `_extract_splits` patterns to reuse
- `phykit/services/tree/vcv_utils.py:45-84` — `parse_gene_trees` utility
- `phykit/services/tree/phylogenetic_ordination.py` — existing PCA pattern (`np.linalg.eigh`, `np.linalg.svd`)
- `phykit/services/tree/rf_distance.py` — RF distance reference
- `tests/sample_files/gene_trees_simple.nwk` — 10 gene trees (7 concordant, 3 alternative topologies)
- `tests/sample_files/tree_simple.tre` — species tree for these gene trees

## R Validation Script

Run this R script to generate reference values for test assertions:

```r
library(phangorn)
library(ape)

# Read trees
gene_files <- readLines("tests/sample_files/gene_trees_simple.nwk")
gene_files <- gene_files[nchar(trimws(gene_files)) > 0]
gene_trees <- lapply(gene_files, function(x) read.tree(text=x))
species_tree <- read.tree("tests/sample_files/tree_simple.tre")

# Prune all to shared taxa
shared <- Reduce(intersect, lapply(gene_trees, function(t) t$tip.label))
gene_trees <- lapply(gene_trees, function(t) drop.tip(t, setdiff(t$tip.label, shared)))
species_tree <- drop.tip(species_tree, setdiff(species_tree$tip.label, shared))

# Root all trees consistently
gene_trees <- lapply(gene_trees, function(t) root(t, shared[1], resolve.root=TRUE))
species_tree <- root(species_tree, shared[1], resolve.root=TRUE)

# Compute RF distance matrix
n <- length(gene_trees)
rf_mat <- matrix(0, n, n)
for(i in 1:(n-1)) {
  for(j in (i+1):n) {
    rf_mat[i,j] <- RF.dist(gene_trees[[i]], gene_trees[[j]], normalize=TRUE)
    rf_mat[j,i] <- rf_mat[i,j]
  }
}
cat("RF distance matrix (normalized):\n")
print(round(rf_mat, 6))

# Build bipartition presence/absence matrix
all_bips <- list()
for(i in 1:n) {
  pp <- prop.part(gene_trees[[i]])
  for(j in seq_along(pp)) {
    bip <- sort(attr(pp, "labels")[pp[[j]]])
    if(length(bip) > 1 && length(bip) < length(shared)) {
      all_bips[[length(all_bips)+1]] <- paste(bip, collapse=",")
    }
  }
}
unique_bips <- sort(unique(unlist(all_bips)))
cat("\nNumber of unique bipartitions:", length(unique_bips), "\n")

# Binary matrix
X <- matrix(0, n, length(unique_bips))
for(i in 1:n) {
  pp <- prop.part(gene_trees[[i]])
  for(j in seq_along(pp)) {
    bip <- paste(sort(attr(pp, "labels")[pp[[j]]]), collapse=",")
    idx <- match(bip, unique_bips)
    if(!is.na(idx)) X[i,idx] <- 1
  }
}

# PCA
pca <- prcomp(X, center=TRUE, scale.=FALSE)
cat("\nVariance explained (first 5 PCs):\n")
ve <- pca$sdev^2 / sum(pca$sdev^2)
print(round(ve[1:min(5,length(ve))], 6))
cat("\nPC scores (first 3 PCs):\n")
print(round(pca$x[,1:min(3,ncol(pca$x))], 6))
```

Run this and record the values for test assertions. Key targets:
- Number of unique bipartitions across the 10 gene trees
- Variance explained by first 3 PCs
- Sign-invariant check: PC scores should match up to sign flip per column

---

### Task 1: Create test gene tree sample file validation and test scaffold

**Files:**
- Verify: `tests/sample_files/gene_trees_simple.nwk` (exists, 10 trees)
- Verify: `tests/sample_files/tree_simple.tre` (exists)
- Create: `tests/unit/services/tree/test_spectral_discordance.py`

**Step 1: Create test file with imports and fixtures**

```python
import json
import os
import sys
import tempfile

import pytest
import numpy as np
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from phykit.services.tree.spectral_discordance import SpectralDiscordance
from phykit.errors import PhykitUserError

here = Path(__file__)
SAMPLE_FILES = here.parent.parent.parent.parent / "sample_files"
TREE_SIMPLE = str(SAMPLE_FILES / "tree_simple.tre")
GENE_TREES = str(SAMPLE_FILES / "gene_trees_simple.nwk")


@pytest.fixture
def default_args():
    return Namespace(
        tree=TREE_SIMPLE,
        gene_trees=GENE_TREES,
        metric="nrf",
        clusters=None,
        n_pcs=None,
        top_loadings=5,
        plot=None,
        json=False,
    )


@pytest.fixture
def no_species_tree_args():
    return Namespace(
        tree=None,
        gene_trees=GENE_TREES,
        metric="nrf",
        clusters=None,
        n_pcs=None,
        top_loadings=5,
        plot=None,
        json=False,
    )


@pytest.fixture
def wrf_args():
    return Namespace(
        tree=TREE_SIMPLE,
        gene_trees=GENE_TREES,
        metric="wrf",
        clusters=None,
        n_pcs=None,
        top_loadings=5,
        plot=None,
        json=False,
    )
```

**Step 2: Verify the import fails (class doesn't exist yet)**

Run: `pytest tests/unit/services/tree/test_spectral_discordance.py --collect-only 2>&1 | head -5`
Expected: ImportError for `spectral_discordance`

---

### Task 2: Create SpectralDiscordance service skeleton

**Files:**
- Create: `phykit/services/tree/spectral_discordance.py`

**Step 1: Write minimal class skeleton**

```python
import sys
from io import StringIO
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from Bio import Phylo

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class SpectralDiscordance(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        tree_file = parsed["tree_file_path"]
        if tree_file is not None:
            super().__init__(tree_file_path=tree_file)
        else:
            # No species tree — skip Tree.__init__ tree loading
            self.tree_file_path = None
        self.gene_trees_path = parsed["gene_trees_path"]
        self.metric = parsed["metric"]
        self.n_clusters = parsed["n_clusters"]
        self.n_pcs = parsed["n_pcs"]
        self.top_loadings = parsed["top_loadings"]
        self.plot_output = parsed["plot_output"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=getattr(args, "tree", None),
            gene_trees_path=args.gene_trees,
            metric=getattr(args, "metric", "nrf"),
            n_clusters=getattr(args, "clusters", None),
            n_pcs=getattr(args, "n_pcs", None),
            top_loadings=getattr(args, "top_loadings", 5),
            plot_output=getattr(args, "plot", None),
            json_output=getattr(args, "json", False),
        )

    def run(self) -> None:
        raise NotImplementedError("spectral_discordance not yet implemented")
```

**Step 2: Verify test file can now import the class**

Run: `pytest tests/unit/services/tree/test_spectral_discordance.py --collect-only`
Expected: collected 0 items (no tests yet, but no import error)

**Step 3: Commit**

```
feat(spectral_discordance): add service skeleton and test scaffold
```

---

### Task 3: Implement gene tree loading and bipartition extraction

**Files:**
- Modify: `phykit/services/tree/spectral_discordance.py`
- Modify: `tests/unit/services/tree/test_spectral_discordance.py`

**Step 1: Write failing tests for gene tree loading and bipartition extraction**

Add to test file:

```python
class TestGeneTreeLoading:
    def test_loads_all_gene_trees(self, default_args):
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        assert len(trees) == 10

    def test_file_not_found(self, default_args):
        svc = SpectralDiscordance(default_args)
        with pytest.raises(PhykitUserError):
            svc._parse_gene_trees("no_such_file.nwk")

    def test_shared_taxa(self, default_args):
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        shared = svc._get_shared_taxa(trees, species_tree=None)
        assert len(shared) == 8  # all 8 taxa shared across all gene trees


class TestBipartitionExtraction:
    def test_canonical_split_smaller_side(self, default_args):
        svc = SpectralDiscordance(default_args)
        all_taxa = frozenset({"A", "B", "C", "D", "E"})
        result = svc._canonical_split(frozenset({"A", "B"}), all_taxa)
        assert result == frozenset({"A", "B"})

    def test_canonical_split_larger_side_flips(self, default_args):
        svc = SpectralDiscordance(default_args)
        all_taxa = frozenset({"A", "B", "C", "D", "E"})
        result = svc._canonical_split(frozenset({"A", "B", "C", "D"}), all_taxa)
        assert result == frozenset({"E"})

    def test_extract_splits_returns_frozensets(self, default_args):
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        shared = svc._get_shared_taxa(trees, species_tree=None)
        all_taxa_fs = frozenset(shared)
        splits = svc._extract_splits(trees[0], all_taxa_fs)
        assert isinstance(splits, set)
        for s in splits:
            assert isinstance(s, frozenset)

    def test_extract_splits_no_trivial(self, default_args):
        """No single-taxon or all-taxa splits."""
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        shared = svc._get_shared_taxa(trees, species_tree=None)
        all_taxa_fs = frozenset(shared)
        splits = svc._extract_splits(trees[0], all_taxa_fs)
        for s in splits:
            assert len(s) > 1 or len(all_taxa_fs) - len(s) > 1
            assert s != all_taxa_fs

    def test_bipartition_matrix_shape(self, default_args):
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        shared = svc._get_shared_taxa(trees, species_tree=None)
        X, bip_index, _ = svc._build_bipartition_matrix(trees, shared, metric="nrf")
        assert X.shape[0] == 10  # 10 gene trees
        assert X.shape[1] == len(bip_index)
        assert X.shape[1] > 0

    def test_bipartition_matrix_binary_for_nrf(self, default_args):
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        shared = svc._get_shared_taxa(trees, species_tree=None)
        X, _, _ = svc._build_bipartition_matrix(trees, shared, metric="nrf")
        assert set(np.unique(X)).issubset({0.0, 1.0})

    def test_bipartition_matrix_wrf_has_branch_lengths(self, default_args):
        svc = SpectralDiscordance(default_args)
        trees = svc._parse_gene_trees(GENE_TREES)
        shared = svc._get_shared_taxa(trees, species_tree=None)
        X, _, _ = svc._build_bipartition_matrix(trees, shared, metric="wrf")
        # wrf matrix should have some entries > 1 (branch lengths)
        assert np.max(X) > 1.0

    def test_species_tree_flags(self, default_args):
        svc = SpectralDiscordance(default_args)
        species_tree = svc.read_tree_file()
        trees = svc._parse_gene_trees(GENE_TREES)
        shared = svc._get_shared_taxa(trees, species_tree=species_tree)
        _, bip_index, sp_flags = svc._build_bipartition_matrix(
            trees, shared, metric="nrf", species_tree=species_tree
        )
        # At least some bipartitions should be in the species tree
        assert any(sp_flags.values())
        # And some should not (alternative topologies)
        assert not all(sp_flags.values())
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/services/tree/test_spectral_discordance.py -v 2>&1 | tail -20`
Expected: FAILED (methods don't exist)

**Step 3: Implement the methods**

Add to `SpectralDiscordance` class in `spectral_discordance.py`:

```python
    def _parse_gene_trees(self, path: str) -> list:
        try:
            lines = Path(path).read_text().splitlines()
        except FileNotFoundError:
            raise PhykitUserError(
                [
                    f"{path} corresponds to no such file or directory.",
                    "Please check filename and pathing",
                ],
                code=2,
            )

        cleaned = [l.strip() for l in lines if l.strip() and not l.strip().startswith("#")]
        trees = []
        for line in cleaned:
            if line.startswith("("):
                trees.append(Phylo.read(StringIO(line), "newick"))
            else:
                tree_path = Path(path).parent / line
                trees.append(Phylo.read(str(tree_path), "newick"))

        if len(trees) < 2:
            raise PhykitUserError(
                ["At least 2 gene trees are required."],
                code=2,
            )
        return trees

    def _get_shared_taxa(self, gene_trees, species_tree=None) -> set:
        taxa_sets = []
        for tree in gene_trees:
            taxa_sets.append(frozenset(
                t.name for t in tree.get_terminals()
            ))
        shared = set.intersection(*[set(ts) for ts in taxa_sets])
        if species_tree is not None:
            sp_taxa = set(t.name for t in species_tree.get_terminals())
            shared = shared & sp_taxa
        if len(shared) < 4:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa across gene trees.",
                    "At least 4 shared taxa are required.",
                ],
                code=2,
            )
        return shared

    @staticmethod
    def _canonical_split(taxa_side, all_taxa):
        complement = all_taxa - taxa_side
        if len(taxa_side) < len(complement):
            return frozenset(taxa_side)
        elif len(taxa_side) > len(complement):
            return frozenset(complement)
        else:
            return min(frozenset(taxa_side), frozenset(complement),
                       key=lambda s: sorted(s))

    def _extract_splits(self, tree, all_taxa_fs) -> set:
        splits = set()
        tip_sets = {}
        stack = [(tree.root, False)]
        while stack:
            node, children_done = stack[-1]
            if not node.clades:
                stack.pop()
                name = node.name
                tip_sets[id(node)] = frozenset((name,)) if name in all_taxa_fs else frozenset()
            elif children_done:
                stack.pop()
                merged = frozenset().union(*(tip_sets[id(c)] for c in node.clades))
                tip_sets[id(node)] = merged
                if len(merged) > 1 and merged != all_taxa_fs:
                    splits.add(self._canonical_split(merged, all_taxa_fs))
            else:
                stack[-1] = (node, True)
                for child in reversed(node.clades):
                    stack.append((child, False))
        return splits

    def _extract_splits_with_lengths(self, tree, all_taxa_fs) -> Dict[frozenset, float]:
        """Extract bipartitions with associated branch lengths for wrf."""
        split_lengths = {}
        tip_sets = {}
        stack = [(tree.root, False)]
        while stack:
            node, children_done = stack[-1]
            if not node.clades:
                stack.pop()
                name = node.name
                tip_sets[id(node)] = frozenset((name,)) if name in all_taxa_fs else frozenset()
            elif children_done:
                stack.pop()
                merged = frozenset().union(*(tip_sets[id(c)] for c in node.clades))
                tip_sets[id(node)] = merged
                if len(merged) > 1 and merged != all_taxa_fs:
                    canon = self._canonical_split(merged, all_taxa_fs)
                    bl = node.branch_length if node.branch_length else 0.0
                    split_lengths[canon] = bl
            else:
                stack[-1] = (node, True)
                for child in reversed(node.clades):
                    stack.append((child, False))
        return split_lengths

    def _build_bipartition_matrix(
        self, gene_trees, shared_taxa, metric="nrf", species_tree=None
    ) -> Tuple[np.ndarray, List[frozenset], Dict[frozenset, bool]]:
        all_taxa_fs = frozenset(shared_taxa)

        # Prune gene trees to shared taxa
        pruned_trees = []
        for gt in gene_trees:
            import copy
            gt_copy = copy.deepcopy(gt)
            tips_to_remove = [
                t.name for t in gt_copy.get_terminals()
                if t.name not in shared_taxa
            ]
            if tips_to_remove:
                for tip_name in tips_to_remove:
                    gt_copy.prune(tip_name)
            pruned_trees.append(gt_copy)

        # Extract bipartitions
        if metric == "wrf":
            per_tree_splits = []
            all_bips = set()
            for gt in pruned_trees:
                sl = self._extract_splits_with_lengths(gt, all_taxa_fs)
                per_tree_splits.append(sl)
                all_bips.update(sl.keys())
        else:
            per_tree_splits = []
            all_bips = set()
            for gt in pruned_trees:
                splits = self._extract_splits(gt, all_taxa_fs)
                per_tree_splits.append(splits)
                all_bips.update(splits)

        # Stable sorted index
        bip_index = sorted(all_bips, key=lambda s: (len(s), sorted(s)))

        # Build matrix
        G = len(gene_trees)
        B = len(bip_index)
        X = np.zeros((G, B))

        bip_to_col = {bip: i for i, bip in enumerate(bip_index)}
        for g, splits in enumerate(per_tree_splits):
            if metric == "wrf":
                for bip, length in splits.items():
                    X[g, bip_to_col[bip]] = length
            else:
                for bip in splits:
                    X[g, bip_to_col[bip]] = 1.0

        # Species tree flags
        sp_flags = {bip: False for bip in bip_index}
        if species_tree is not None:
            import copy
            sp_copy = copy.deepcopy(species_tree)
            sp_tips_to_remove = [
                t.name for t in sp_copy.get_terminals()
                if t.name not in shared_taxa
            ]
            if sp_tips_to_remove:
                for tip_name in sp_tips_to_remove:
                    sp_copy.prune(tip_name)
            sp_splits = self._extract_splits(sp_copy, all_taxa_fs)
            for bip in bip_index:
                if bip in sp_splits:
                    sp_flags[bip] = True

        return X, bip_index, sp_flags
```

**Step 4: Run tests**

Run: `pytest tests/unit/services/tree/test_spectral_discordance.py -v`
Expected: All tests in TestGeneTreeLoading and TestBipartitionExtraction PASS

**Step 5: Commit**

```
feat(spectral_discordance): implement gene tree loading and bipartition matrix
```

---

### Task 4: Implement PCA via SVD

**Files:**
- Modify: `phykit/services/tree/spectral_discordance.py`
- Modify: `tests/unit/services/tree/test_spectral_discordance.py`

**Step 1: Write failing tests for PCA**

```python
class TestPCA:
    def _get_pca(self, args):
        svc = SpectralDiscordance(args)
        trees = svc._parse_gene_trees(GENE_TREES)
        species_tree = svc.read_tree_file() if svc.tree_file_path else None
        shared = svc._get_shared_taxa(trees, species_tree)
        X, bip_index, sp_flags = svc._build_bipartition_matrix(
            trees, shared, metric=args.metric, species_tree=species_tree
        )
        scores, var_explained, loadings = svc._run_pca(X)
        return scores, var_explained, loadings, bip_index

    def test_scores_shape(self, default_args):
        scores, _, _, _ = self._get_pca(default_args)
        assert scores.shape[0] == 10  # 10 gene trees
        assert scores.shape[1] <= 10  # at most G components

    def test_variance_explained_sums_to_one(self, default_args):
        _, ve, _, _ = self._get_pca(default_args)
        assert pytest.approx(np.sum(ve), abs=1e-6) == 1.0

    def test_variance_explained_descending(self, default_args):
        _, ve, _, _ = self._get_pca(default_args)
        for i in range(len(ve) - 1):
            assert ve[i] >= ve[i + 1] - 1e-10

    def test_loadings_shape(self, default_args):
        scores, _, loadings, bip_index = self._get_pca(default_args)
        # loadings: n_components x n_bipartitions
        assert loadings.shape[1] == len(bip_index)
        assert loadings.shape[0] == scores.shape[1]

    def test_first_pc_captures_most_variance(self, default_args):
        _, ve, _, _ = self._get_pca(default_args)
        # With 7 concordant + 3 discordant trees, PC1 should capture >30% of variance
        assert ve[0] > 0.3

    def test_pca_wrf(self, wrf_args):
        """Weighted RF PCA should also work."""
        scores, ve, loadings, _ = self._get_pca(wrf_args)
        assert scores.shape[0] == 10
        assert pytest.approx(np.sum(ve), abs=1e-6) == 1.0
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/services/tree/test_spectral_discordance.py::TestPCA -v`
Expected: FAILED (method `_run_pca` doesn't exist)

**Step 3: Implement PCA**

Add to `SpectralDiscordance`:

```python
    def _run_pca(
        self, X: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """PCA via SVD on centered bipartition matrix.

        Returns (scores, variance_explained, loadings).
        scores: G x n_components
        variance_explained: n_components (sums to 1.0)
        loadings: n_components x B
        """
        # Center columns
        X_centered = X - X.mean(axis=0)

        # Handle degenerate case: all trees identical
        if np.allclose(X_centered, 0):
            n_comp = min(X.shape)
            return (
                np.zeros((X.shape[0], n_comp)),
                np.zeros(n_comp),
                np.zeros((n_comp, X.shape[1])),
            )

        U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
        scores = U * S  # G x min(G, B)
        total_var = np.sum(S ** 2)
        var_explained = S ** 2 / total_var if total_var > 0 else np.zeros_like(S)
        loadings = Vt  # n_components x B

        return scores, var_explained, loadings
```

**Step 4: Run tests**

Run: `pytest tests/unit/services/tree/test_spectral_discordance.py::TestPCA -v`
Expected: All PASS

**Step 5: Commit**

```
feat(spectral_discordance): implement PCA via SVD on bipartition matrix
```

---

### Task 5: Implement spectral clustering with eigengap

**Files:**
- Modify: `phykit/services/tree/spectral_discordance.py`
- Modify: `tests/unit/services/tree/test_spectral_discordance.py`

**Step 1: Write failing tests for spectral clustering**

```python
class TestSpectralClustering:
    def _get_clusters(self, args, n_clusters=None):
        svc = SpectralDiscordance(args)
        trees = svc._parse_gene_trees(GENE_TREES)
        species_tree = svc.read_tree_file() if svc.tree_file_path else None
        shared = svc._get_shared_taxa(trees, species_tree)
        X, _, _ = svc._build_bipartition_matrix(trees, shared, metric=args.metric)
        X_centered = X - X.mean(axis=0)
        labels, K, eigengaps = svc._spectral_cluster(X_centered, n_clusters=n_clusters)
        return labels, K, eigengaps

    def test_labels_shape(self, default_args):
        labels, _, _ = self._get_clusters(default_args)
        assert len(labels) == 10

    def test_labels_valid_range(self, default_args):
        labels, K, _ = self._get_clusters(default_args)
        assert all(0 <= l < K for l in labels)

    def test_k_at_least_2(self, default_args):
        _, K, _ = self._get_clusters(default_args)
        assert K >= 2

    def test_eigengaps_positive(self, default_args):
        _, _, eigengaps = self._get_clusters(default_args)
        assert all(g >= -1e-10 for g in eigengaps)

    def test_override_k(self, default_args):
        labels, K, _ = self._get_clusters(default_args, n_clusters=3)
        assert K == 3
        assert len(set(labels)) <= 3

    def test_concordant_trees_cluster_together(self, default_args):
        """Trees 0-5 are concordant and 6-8 are discordant.
        They should mostly end up in different clusters."""
        labels, K, _ = self._get_clusters(default_args)
        # At least 2 clusters
        assert K >= 2
        # The concordant trees (0-5) should mostly share a cluster
        concordant_labels = [labels[i] for i in range(6)]
        from collections import Counter
        most_common = Counter(concordant_labels).most_common(1)[0][1]
        assert most_common >= 4  # at least 4 of 6 concordant in same cluster
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/services/tree/test_spectral_discordance.py::TestSpectralClustering -v`
Expected: FAILED

**Step 3: Implement spectral clustering**

Add to `SpectralDiscordance`:

```python
    def _spectral_cluster(
        self, X_centered: np.ndarray, n_clusters: int = None,
    ) -> Tuple[np.ndarray, int, np.ndarray]:
        """Spectral clustering with eigengap auto-K detection.

        Returns (labels, K, eigengaps).
        """
        from scipy.spatial.distance import pdist, squareform

        G = X_centered.shape[0]

        # Pairwise Euclidean distances
        dists = squareform(pdist(X_centered, metric="euclidean"))

        # Gaussian kernel bandwidth (median heuristic)
        upper_tri = dists[np.triu_indices(G, k=1)]
        sigma = np.median(upper_tri)
        if sigma == 0:
            nonzero = upper_tri[upper_tri > 0]
            sigma = np.mean(nonzero) if len(nonzero) > 0 else 1.0

        # Similarity matrix
        W = np.exp(-dists ** 2 / (2 * sigma ** 2))
        np.fill_diagonal(W, W.diagonal() + 1e-10)  # prevent isolated nodes

        # Degree matrix and normalized Laplacian
        d = np.sum(W, axis=1)
        d_inv_sqrt = 1.0 / np.sqrt(d)
        D_inv_sqrt = np.diag(d_inv_sqrt)
        L_norm = np.eye(G) - D_inv_sqrt @ W @ D_inv_sqrt

        # Eigendecompose (symmetric, so use eigh for stability)
        eigenvalues, eigenvectors = np.linalg.eigh(L_norm)

        # Eigengap heuristic
        max_k = min(20, G // 2)
        if max_k < 2:
            max_k = 2
        eigengaps = np.diff(eigenvalues[:max_k + 1])

        if n_clusters is not None:
            K = n_clusters
        else:
            # argmax of gaps starting from index 1 (skip trivial λ_0≈0)
            # K = i means use first i eigenvectors
            K = int(np.argmax(eigengaps[1:max_k]) + 2)  # +2: skip index 0 gap, 1-indexed
            K = max(K, 2)

        # K-means on first K eigenvectors
        vecs = eigenvectors[:, :K]
        # Normalize rows
        row_norms = np.linalg.norm(vecs, axis=1, keepdims=True)
        row_norms[row_norms == 0] = 1.0
        vecs = vecs / row_norms

        labels = self._kmeans(vecs, K)

        return labels, K, eigengaps

    @staticmethod
    def _kmeans(X: np.ndarray, K: int, max_iter: int = 300) -> np.ndarray:
        """Simple K-means clustering (avoids sklearn dependency)."""
        G = X.shape[0]
        rng = np.random.RandomState(42)

        # K-means++ initialization
        centers = [X[rng.randint(G)]]
        for _ in range(1, K):
            dists = np.array([
                min(np.sum((x - c) ** 2) for c in centers)
                for x in X
            ])
            probs = dists / dists.sum() if dists.sum() > 0 else np.ones(G) / G
            centers.append(X[rng.choice(G, p=probs)])
        centers = np.array(centers)

        for _ in range(max_iter):
            # Assign
            dists_to_centers = np.array([
                np.sum((X - c) ** 2, axis=1) for c in centers
            ]).T
            labels = np.argmin(dists_to_centers, axis=1)

            # Update
            new_centers = np.array([
                X[labels == k].mean(axis=0) if np.any(labels == k) else centers[k]
                for k in range(K)
            ])
            if np.allclose(new_centers, centers):
                break
            centers = new_centers

        return labels
```

**Step 4: Run tests**

Run: `pytest tests/unit/services/tree/test_spectral_discordance.py::TestSpectralClustering -v`
Expected: All PASS

**Step 5: Commit**

```
feat(spectral_discordance): implement spectral clustering with eigengap auto-K
```

---

### Task 6: Implement run(), text output, and JSON output

**Files:**
- Modify: `phykit/services/tree/spectral_discordance.py`
- Modify: `tests/unit/services/tree/test_spectral_discordance.py`

**Step 1: Write failing tests for run()**

```python
class TestRun:
    @patch("builtins.print")
    def test_text_output(self, mocked_print, default_args):
        svc = SpectralDiscordance(default_args)
        svc.run()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Spectral Discordance" in all_output
        assert "Variance explained" in all_output
        assert "cluster" in all_output.lower()

    @patch("builtins.print")
    def test_json_output(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            metric="nrf",
            clusters=None,
            n_pcs=None,
            top_loadings=5,
            plot=None,
            json=True,
        )
        svc = SpectralDiscordance(args)
        svc.run()

        payload = json.loads(mocked_print.call_args.args[0])
        assert "scores" in payload
        assert "variance_explained" in payload
        assert "cluster_assignments" in payload
        assert "top_loadings" in payload
        assert "eigengap_values" in payload
        assert "n_clusters" in payload
        assert "metric" in payload

    @patch("builtins.print")
    def test_no_species_tree(self, mocked_print, no_species_tree_args):
        svc = SpectralDiscordance(no_species_tree_args)
        svc.run()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "Spectral Discordance" in all_output

    @patch("builtins.print")
    def test_wrf_metric(self, mocked_print, wrf_args):
        svc = SpectralDiscordance(wrf_args)
        svc.run()

        all_output = " ".join(
            str(call.args[0]) for call in mocked_print.call_args_list if call.args
        )
        assert "wrf" in all_output.lower() or "weighted" in all_output.lower()

    @patch("builtins.print")
    def test_json_loadings_have_species_tree_flag(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            metric="nrf",
            clusters=None,
            n_pcs=None,
            top_loadings=5,
            plot=None,
            json=True,
        )
        svc = SpectralDiscordance(args)
        svc.run()

        payload = json.loads(mocked_print.call_args.args[0])
        pc1_loadings = payload["top_loadings"]["PC1"]
        assert len(pc1_loadings) > 0
        assert "in_species_tree" in pc1_loadings[0]
        assert "bipartition" in pc1_loadings[0]
        assert "loading" in pc1_loadings[0]
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/unit/services/tree/test_spectral_discordance.py::TestRun -v`
Expected: FAILED (run() raises NotImplementedError)

**Step 3: Implement run(), output formatting**

Replace the `run()` method and add output helpers:

```python
    def run(self) -> None:
        # Load gene trees
        gene_trees = self._parse_gene_trees(self.gene_trees_path)

        # Load species tree if provided
        species_tree = None
        if self.tree_file_path is not None:
            species_tree = self.read_tree_file()

        # Shared taxa
        shared = self._get_shared_taxa(gene_trees, species_tree)

        # Build bipartition matrix
        X, bip_index, sp_flags = self._build_bipartition_matrix(
            gene_trees, shared, metric=self.metric, species_tree=species_tree
        )

        G = X.shape[0]
        if G < 5:
            print(
                f"Warning: only {G} gene trees; clustering may be unreliable.",
                file=sys.stderr,
            )

        # PCA
        scores, var_explained, loadings = self._run_pca(X)
        n_pcs = self.n_pcs if self.n_pcs else min(10, scores.shape[1])
        n_pcs = min(n_pcs, scores.shape[1])

        # Spectral clustering
        X_centered = X - X.mean(axis=0)
        labels, K, eigengaps = self._spectral_cluster(
            X_centered, n_clusters=self.n_clusters
        )

        # Top loadings per PC
        top_loadings = self._get_top_loadings(
            loadings, bip_index, sp_flags, n_pcs, self.top_loadings
        )

        # Plot
        if self.plot_output:
            self._plot_scatter(scores, labels, K, self.plot_output + "_scatter.png")
            self._plot_eigengap(eigengaps, K, self.plot_output + "_eigengap.png")

        # Output
        if self.json_output:
            result = self._format_json(
                scores, var_explained, top_loadings, labels, K,
                eigengaps, n_pcs, bip_index, sp_flags
            )
            print_json(result)
        else:
            self._print_text(
                scores, var_explained, top_loadings, labels, K,
                eigengaps, n_pcs, G, len(shared), len(bip_index)
            )

    def _bipartition_to_str(self, bip: frozenset) -> str:
        return "{" + ",".join(sorted(bip)) + "}"

    def _get_top_loadings(
        self, loadings, bip_index, sp_flags, n_pcs, top_n,
    ) -> Dict[str, list]:
        result = {}
        for pc in range(min(n_pcs, loadings.shape[0])):
            pc_loadings = loadings[pc]
            top_idx = np.argsort(np.abs(pc_loadings))[::-1][:top_n]
            entries = []
            for idx in top_idx:
                bip = bip_index[idx]
                entries.append({
                    "bipartition": self._bipartition_to_str(bip),
                    "loading": float(pc_loadings[idx]),
                    "in_species_tree": sp_flags.get(bip, False),
                })
            result[f"PC{pc + 1}"] = entries
        return result

    def _format_json(
        self, scores, var_explained, top_loadings, labels, K,
        eigengaps, n_pcs, bip_index, sp_flags,
    ) -> Dict:
        G = scores.shape[0]
        score_dict = {}
        for g in range(G):
            row = {}
            for pc in range(min(n_pcs, scores.shape[1])):
                row[f"PC{pc + 1}"] = float(scores[g, pc])
            row["cluster"] = int(labels[g])
            score_dict[f"gene_tree_{g + 1}"] = row

        return {
            "metric": self.metric,
            "n_gene_trees": G,
            "n_bipartitions": len(bip_index),
            "n_clusters": int(K),
            "scores": score_dict,
            "variance_explained": {
                f"PC{i + 1}": float(var_explained[i])
                for i in range(min(n_pcs, len(var_explained)))
            },
            "top_loadings": top_loadings,
            "cluster_assignments": [int(l) for l in labels],
            "eigengap_values": [float(g) for g in eigengaps],
        }

    def _print_text(
        self, scores, var_explained, top_loadings, labels, K,
        eigengaps, n_pcs, n_trees, n_taxa, n_bips,
    ) -> None:
        print("Spectral Discordance Decomposition")
        print(f"\nMetric: {self.metric}")
        print(f"Gene trees: {n_trees}")
        print(f"Shared taxa: {n_taxa}")
        print(f"Unique bipartitions: {n_bips}")
        print(f"Clusters detected: {K}")

        # Variance explained
        print("\nVariance explained:")
        for i in range(min(n_pcs, len(var_explained))):
            bar = "#" * int(var_explained[i] * 50)
            print(f"  PC{i + 1:<3d} {var_explained[i]:>8.4f}  {bar}")

        # Top loadings
        for pc_name, entries in top_loadings.items():
            print(f"\nTop loadings for {pc_name}:")
            for e in entries:
                flag = " *" if e["in_species_tree"] else ""
                print(f"  {e['loading']:>8.4f}  {e['bipartition']}{flag}")

        # Scores table
        show_pcs = min(n_pcs, scores.shape[1], 5)
        header = f"  {'Gene tree':<14s}"
        for pc in range(show_pcs):
            header += f"{'PC' + str(pc + 1):>10s}"
        header += f"{'Cluster':>10s}"
        print(f"\nPC scores:")
        print(header)
        for g in range(scores.shape[0]):
            row = f"  {'gene_tree_' + str(g + 1):<14s}"
            for pc in range(show_pcs):
                row += f"{scores[g, pc]:>10.4f}"
            row += f"{labels[g]:>10d}"
            print(row)
```

**Step 4: Run tests**

Run: `pytest tests/unit/services/tree/test_spectral_discordance.py::TestRun -v`
Expected: All PASS

**Step 5: Commit**

```
feat(spectral_discordance): implement run() with text and JSON output
```

---

### Task 7: Implement plots (scatter + eigengap)

**Files:**
- Modify: `phykit/services/tree/spectral_discordance.py`
- Modify: `tests/unit/services/tree/test_spectral_discordance.py`

**Step 1: Write failing test**

```python
class TestPlot:
    @patch("builtins.print")
    def test_plots_created(self, mocked_print):
        try:
            import matplotlib
        except ImportError:
            pytest.skip("matplotlib not installed")

        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = os.path.join(tmpdir, "test")
            args = Namespace(
                tree=TREE_SIMPLE,
                gene_trees=GENE_TREES,
                metric="nrf",
                clusters=None,
                n_pcs=None,
                top_loadings=5,
                plot=prefix,
                json=False,
            )
            svc = SpectralDiscordance(args)
            svc.run()

            scatter = prefix + "_scatter.png"
            eigengap = prefix + "_eigengap.png"
            assert os.path.exists(scatter), f"Missing {scatter}"
            assert os.path.exists(eigengap), f"Missing {eigengap}"
            assert os.path.getsize(scatter) > 0
            assert os.path.getsize(eigengap) > 0
```

**Step 2: Run to verify it fails**

**Step 3: Implement plot methods**

```python
    def _plot_scatter(
        self, scores, labels, K, output_path,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print(
                "matplotlib is required for plotting. Install matplotlib and retry."
            )
            raise SystemExit(2)

        fig, ax = plt.subplots(figsize=(8, 6))

        cmap = plt.get_cmap("tab10")
        for k in range(K):
            mask = labels == k
            ax.scatter(
                scores[mask, 0], scores[mask, 1],
                c=[cmap(k)], label=f"Cluster {k + 1}",
                s=60, edgecolors="black", linewidths=0.5, alpha=0.8,
            )
            # Label each point
            idxs = np.where(mask)[0]
            for idx in idxs:
                ax.annotate(
                    str(idx + 1), (scores[idx, 0], scores[idx, 1]),
                    fontsize=7, ha="center", va="bottom",
                    xytext=(0, 4), textcoords="offset points",
                )

        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.set_title("Gene Tree Space — Spectral Discordance")
        ax.legend()
        fig.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved scatter plot: {output_path}")

    def _plot_eigengap(
        self, eigengaps, K, output_path,
    ) -> None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print(
                "matplotlib is required for plotting. Install matplotlib and retry."
            )
            raise SystemExit(2)

        fig, ax = plt.subplots(figsize=(8, 4))

        n = len(eigengaps)
        x = np.arange(1, n + 1)
        colors = ["tab:red" if i == K - 1 else "tab:blue" for i in range(n)]
        ax.bar(x, eigengaps, color=colors, edgecolor="black", linewidth=0.5)

        ax.set_xlabel("Eigenvalue index (i)")
        ax.set_ylabel(r"Eigengap ($\lambda_{i+1} - \lambda_i$)")
        ax.set_title(f"Eigengap Heuristic — K = {K} clusters")
        ax.axvline(K, color="tab:red", linestyle="--", alpha=0.7, label=f"K = {K}")
        ax.legend()
        fig.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved eigengap plot: {output_path}")
```

**Step 4: Run tests**

Run: `pytest tests/unit/services/tree/test_spectral_discordance.py::TestPlot -v`
Expected: PASS

**Step 5: Commit**

```
feat(spectral_discordance): add scatter and eigengap plots
```

---

### Task 8: Register the command (service factory, CLI, setup.py)

**Files:**
- Modify: `phykit/service_factories.py:101` — add LazyServiceFactory
- Modify: `phykit/cli_registry.py:159` — add aliases
- Modify: `phykit/phykit.py:277` — help banner
- Modify: `phykit/phykit.py:5001` — @staticmethod handler after discordance_asymmetry
- Modify: `phykit/phykit.py:5436` — module-level wrapper before treeness
- Modify: `setup.py:254` — console_scripts entry points

**Step 1: Add to service_factories.py**

Insert before `SERVICE_FACTORIES` dict (after line 100):

```python
SpectralDiscordance = _LazyServiceFactory("phykit.services.tree.spectral_discordance", "SpectralDiscordance")
```

**Step 2: Add aliases to cli_registry.py**

Insert before `# Helper aliases` (after line 159):

```python
    "spec_disc": "spectral_discordance",
    "sd": "spectral_discordance",
```

**Step 3: Add help banner entry to phykit.py**

Insert before the `treeness` entry (around line 277):

```
                spectral_discordance (alias: spec_disc; sd)
                    - PCA + spectral clustering of gene tree space via
                      bipartition decomposition
```

**Step 4: Add @staticmethod handler to Phykit class**

Insert after the `discordance_asymmetry` handler (after line 5001):

```python
    @staticmethod
    def spectral_discordance(argv):
        parser = _new_parser(
            description=textwrap.dedent(
                f"""\
                {help_header}

                Spectral discordance decomposition — decompose gene tree
                space via PCA on a bipartition presence/absence (or
                branch-length) matrix, with spectral clustering and
                automatic cluster detection via the eigengap heuristic.

                Each gene tree is encoded as a vector over the union of
                all bipartitions observed across gene trees. PCA reveals
                the axes of topological variation, with loading vectors
                identifying which bipartitions drive each PC. Spectral
                clustering groups genes sharing similar topologies.

                Two metrics are available:
                - nrf (default): binary presence/absence (normalized RF)
                - wrf: branch-length weighted

                Aliases:
                  spectral_discordance, spec_disc, sd
                Command line interfaces:
                  pk_spectral_discordance, pk_spec_disc, pk_sd

                Usage:
                phykit spectral_discordance -g <gene_trees> [-t <tree>] [--metric nrf|wrf] [--clusters K] [--n-pcs N] [--top-loadings N] [--plot <prefix>] [--json]

                Options
                =====================================================
                -g/--gene-trees             file of gene trees (one
                                            Newick per line, or file
                                            of filenames)

                -t/--tree                   species tree (optional; flags
                                            species-tree bipartitions in
                                            loading output)

                --metric                    distance metric: nrf or wrf
                                            (default: nrf)

                --clusters                  override auto-detected K

                --n-pcs                     number of PCs to report
                                            (default: min(10, G-1))

                --top-loadings              top bipartitions per PC
                                            (default: 5)

                --plot                      output prefix for plots
                                            (generates _scatter.png and
                                            _eigengap.png)

                --json                      output results as JSON
                """
            ),
        )
        parser.add_argument(
            "-g", "--gene-trees", type=str, required=True, help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "-t", "--tree", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--metric", type=str, required=False, default="nrf",
            choices=["nrf", "wrf"], help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--clusters", type=int, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--n-pcs", type=int, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--top-loadings", type=int, required=False, default=5,
            help=SUPPRESS, metavar=""
        )
        parser.add_argument(
            "--plot", type=str, required=False, default=None,
            help=SUPPRESS, metavar=""
        )
        _add_json_argument(parser)
        _run_service(parser, argv, SpectralDiscordance)
```

**Step 5: Add module-level wrapper**

Insert before `# Alignment- and tree-based functions` (around line 5440):

```python
def spectral_discordance(argv=None):
    Phykit.spectral_discordance(sys.argv[1:])
```

**Step 6: Add console_scripts to setup.py**

Insert before `"pk_saturation` (around line 254):

```python
            "pk_spectral_discordance = phykit.phykit:spectral_discordance",
            "pk_spec_disc = phykit.phykit:spectral_discordance",
            "pk_sd = phykit.phykit:spectral_discordance",
```

**Step 7: Smoke test**

Run: `python -m phykit spectral_discordance -g tests/sample_files/gene_trees_simple.nwk -t tests/sample_files/tree_simple.tre`
Expected: text output with scores, clusters, loadings

**Step 8: Commit**

```
feat(spectral_discordance): register CLI command with aliases and entry points
```

---

### Task 9: Add R validation tests

**Files:**
- Modify: `tests/unit/services/tree/test_spectral_discordance.py`

**Step 1: Add R validation test class**

The key mathematical property to validate: PCA on the binary bipartition matrix should produce the same variance-explained proportions and the same scores (up to sign flips) as R's `prcomp()` on the same matrix. We validate structure rather than exact values since the bipartition indexing may differ.

```python
class TestRValidation:
    """Validate mathematical properties that must match R's prcomp().

    Properties validated:
    - Euclidean distances between rows of the centered bipartition matrix
      equal sqrt(RF) for the binary (nrf) case
    - Variance explained sums to 1 and is in descending order
    - Reconstruction: scores @ loadings ≈ X_centered
    - Orthogonality: loadings are orthonormal rows
    """

    def _setup(self):
        args = Namespace(
            tree=TREE_SIMPLE,
            gene_trees=GENE_TREES,
            metric="nrf",
            clusters=None,
            n_pcs=None,
            top_loadings=5,
            plot=None,
            json=False,
        )
        svc = SpectralDiscordance(args)
        trees = svc._parse_gene_trees(GENE_TREES)
        species_tree = svc.read_tree_file()
        shared = svc._get_shared_taxa(trees, species_tree)
        X, bip_index, sp_flags = svc._build_bipartition_matrix(
            trees, shared, metric="nrf", species_tree=species_tree
        )
        scores, ve, loadings = svc._run_pca(X)
        return svc, X, scores, ve, loadings, bip_index, trees, shared

    def test_euclidean_equals_sqrt_rf(self):
        """Euclidean distance in binary bipartition space = sqrt(RF)."""
        svc, X, _, _, _, bip_index, _, _ = self._setup()
        from scipy.spatial.distance import pdist, squareform
        euclidean = squareform(pdist(X, metric="euclidean"))
        manhattan = squareform(pdist(X, metric="cityblock"))  # L1 = RF
        # Euclidean = sqrt(L1) for binary data
        for i in range(X.shape[0]):
            for j in range(i + 1, X.shape[0]):
                expected = np.sqrt(manhattan[i, j])
                assert euclidean[i, j] == pytest.approx(expected, abs=1e-10)

    def test_reconstruction(self):
        """scores @ loadings ≈ X_centered (full rank)."""
        _, X, scores, _, loadings, _, _, _ = self._setup()
        X_centered = X - X.mean(axis=0)
        reconstructed = scores @ loadings
        assert np.allclose(reconstructed, X_centered, atol=1e-10)

    def test_loadings_orthonormal(self):
        """Loading vectors should be orthonormal."""
        _, _, _, _, loadings, _, _, _ = self._setup()
        # Vt @ V = I (rows of loadings are orthonormal)
        product = loadings @ loadings.T
        assert np.allclose(product, np.eye(loadings.shape[0]), atol=1e-10)

    def test_scores_uncorrelated(self):
        """PC scores should be uncorrelated (covariance is diagonal)."""
        _, _, scores, _, _, _, _, _ = self._setup()
        cov = np.cov(scores.T)
        # Off-diagonal should be ~0
        mask = ~np.eye(cov.shape[0], dtype=bool)
        assert np.allclose(cov[mask], 0, atol=1e-10)
```

**Step 2: Run tests**

Run: `pytest tests/unit/services/tree/test_spectral_discordance.py::TestRValidation -v`
Expected: All PASS

**Step 3: Commit**

```
test(spectral_discordance): add R-equivalent mathematical validation tests
```

---

### Task 10: Update docs and version bump

**Files:**
- Modify: `docs/usage/index.rst:1127` — add command reference section before `.. _cmd-concordance_asr:`
- Modify: `docs/tutorials/index.rst:1841` — add tutorial section before tutorial 12
- Modify: `phykit/version.py` — bump version

**Step 1: Add command reference to docs/usage/index.rst**

Insert before `.. _cmd-concordance_asr:` (line 1127):

```rst
.. _cmd-spectral_discordance:

Spectral discordance decomposition
###################################
Function names: spectral_discordance; spec_disc; sd |br|
Command line interface: pk_spectral_discordance; pk_spec_disc; pk_sd

Decompose gene tree space via PCA on a bipartition presence/absence (or
branch-length) matrix, with spectral clustering and automatic cluster
detection via the eigengap heuristic.

Each gene tree is encoded as a vector over the union of all observed
bipartitions. PCA reveals the axes of topological variation, with loading
vectors directly identifying which bipartitions drive each PC. Spectral
clustering groups genes sharing similar topologies; the number of clusters
is auto-detected from the eigengap of the graph Laplacian.

Two metrics are available:

- **nrf** (default): binary presence/absence, consistent with normalized Robinson-Foulds distance
- **wrf**: branch-length weighted

.. code-block:: shell

   phykit spectral_discordance -g <gene_trees> [-t <tree>] [--metric nrf|wrf] [--clusters K] [--n-pcs N] [--top-loadings N] [--plot <prefix>] [--json]

Options: |br|
*-g/\\-\\-gene-trees*: file of gene trees (one Newick per line, or file of filenames) |br|
*-t/\\-\\-tree*: species tree (optional; flags species-tree bipartitions in loadings) |br|
*--metric*: distance metric: ``nrf`` or ``wrf`` (default: ``nrf``) |br|
*--clusters*: override auto-detected number of clusters |br|
*--n-pcs*: number of PCs to report (default: min(10, G-1)) |br|
*--top-loadings*: top bipartitions per PC to display (default: 5) |br|
*--plot*: output prefix for plots (generates ``_scatter.png`` and ``_eigengap.png``) |br|
*--json*: output results as JSON

|
```

**Step 2: Add tutorial section to docs/tutorials/index.rst**

Insert before `12. Testing for rate heterogeneity` (line 1841). This becomes a new section, and renumber 12 onward.

Keep the tutorial concise — 3 steps: run basic analysis, interpret loadings, generate plots.

**Step 3: Bump version**

```python
__version__ = "2.1.36"
```

**Step 4: Run full test suite**

Run: `pytest tests/unit/services/tree/test_spectral_discordance.py -v`
Run: `pytest tests/unit/services/tree/test_ancestral_reconstruction.py -v`  (regression)
Expected: All pass

**Step 5: CLI smoke test**

```bash
python -m phykit spectral_discordance \
    -g tests/sample_files/gene_trees_simple.nwk \
    -t tests/sample_files/tree_simple.tre \
    --metric nrf
```

**Step 6: Commit**

```
feat(spectral_discordance): add docs, tutorial, and bump version to 2.1.36
```
