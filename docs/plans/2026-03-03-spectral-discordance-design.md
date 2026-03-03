# Spectral Discordance Decomposition — Design

## Problem

Gene tree discordance is commonly summarized per-branch (gCF/gDF) but this loses the global structure of tree-space variation. Users want to ask: which groups of genes share alternative topologies, and what genomic features predict tree topology variation?

## Solution

A new `spectral_discordance` command that:
1. Encodes each gene tree as a binary/weighted vector over the union of all bipartitions
2. Runs PCA on the bipartition matrix to ordinate gene trees in low-dimensional space
3. Performs spectral clustering with automatic K detection via eigengap heuristic
4. Reports PC scores, variance explained, bipartition loadings, and cluster assignments

## Command Interface

```
phykit spectral_discordance \
    -g <gene_trees>              # required
    [-t <species_tree>]          # optional: flag species-tree bipartitions in loadings
    [--metric nrf|wrf]           # default: nrf (normalized RF = binary; wrf = branch lengths)
    [--clusters K]               # override auto-detected K
    [--n-pcs N]                  # PCs to report (default: min(10, G-1))
    [--top-loadings N]           # top bipartitions per PC (default: 5)
    [--plot <prefix>]            # <prefix>_scatter.png + <prefix>_eigengap.png
    [--json]
```

Aliases: `spec_disc`, `sd`
Entry points: `pk_spectral_discordance`, `pk_spec_disc`, `pk_sd`

## Algorithm

### Step 1 — Load and normalize

- Load gene trees via `vcv_utils.parse_gene_trees()`; optionally load species tree
- Compute shared taxa intersection across all gene trees (and species tree if provided)
- Prune all trees to shared taxa; error if < 4 shared taxa

### Step 2 — Build bipartition matrix X (G × B)

- Extract canonical bipartitions per tree via postorder traversal (reuse `_canonical_split` pattern from `discordance_asymmetry.py`)
- Collect union of all unique bipartitions; sort into stable index
- `nrf` mode: X[g,b] = 1 if bipartition b present in gene tree g, else 0
- `wrf` mode: X[g,b] = branch_length of bipartition b in gene tree g, else 0
- If species tree provided, record which bipartitions appear in it

### Step 3 — PCA via SVD

- Center X (subtract column means); centering preserves pairwise Euclidean distances
- `U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)`
- Scores: `U * S` (shape G × min(G,B))
- Variance explained: `S² / sum(S²)`
- Loadings: rows of Vt — each entry maps a bipartition to a PC
- Note: Euclidean distance in binary bipartition space = sqrt(RF). PCA preserves Euclidean distances, so the ordination is monotonically consistent with RF.

### Step 4 — Spectral clustering with eigengap auto-K

- Pairwise Euclidean distances from centered X
- Gaussian kernel: `W = exp(-D² / (2σ²))`, σ = median of pairwise distances
  - If σ = 0 (>50% identical trees): fall back to mean of nonzero distances
  - Add epsilon (1e-10) to W diagonal to prevent isolated nodes
- Degree matrix: D_deg_ii = Σ_j W_ij
- Normalized Laplacian: `L_norm = I - D_deg^{-1/2} W D_deg^{-1/2}`
- Eigendecompose L_norm; eigenvalues ascending: λ_1 ≤ λ_2 ≤ ...
- Eigengap: K = argmax_i(λ_{i+1} - λ_i) for i = 1..min(20, G/2), minimum K = 2
- `--clusters` overrides auto K
- K-means on first K eigenvectors of L_norm
- Warn if G < 5

### Step 5 — Output

**Text (default):** Summary header, then TSV table of gene_tree / PC1 / PC2 / ... / cluster. Then top-N bipartition loadings per PC (bipartition as taxon set, `*` if in species tree).

**JSON:** Structured output with: scores, variance_explained, top_loadings (per PC: bipartition string + value + in_species_tree flag), cluster_assignments, eigengap_values, K, metric used.

**Plots:** (1) PC1 vs PC2 scatter colored by cluster with legend; (2) eigengap bar chart with chosen K highlighted.

## Mathematical Notes

- For `nrf`: Euclidean distance = sqrt(RF), NOT RF itself. The relationship is monotonic, so PCA ordination is consistent with RF ranking.
- For `wrf`: Euclidean and weighted-RF (L1) have no closed-form link. PCA preserves Euclidean distances in branch-length space, which is a useful but distinct metric from wRF.
- Bipartitions use canonical form: frozenset of the smaller taxon side, with lexicographic tie-breaking (matching `discordance_asymmetry._canonical_split`).

## R Validation Targets

- **RF distances:** `phangorn::RF.dist()` for pairwise normalized RF matrix
- **PCA of tree space:** `phangorn::treedist()` + `stats::prcomp()` on the bipartition matrix
- **Spectral clustering eigenvalues:** `kernlab::specc()` or manual Laplacian eigendecomposition

Specific validation: for a set of gene trees, the variance explained per PC and the PC scores should match `prcomp()` on the same binary bipartition matrix to machine precision.

## Edge Cases

- All gene trees identical: zero variance; report gracefully with warning
- G < 5: warn that clustering may be unreliable
- σ = 0 in kernel: fall back to mean of nonzero distances
- Disconnected similarity graph: epsilon on diagonal prevents this
- Missing taxa across gene trees: prune to shared set first

## Files to Create/Modify

1. **Create** `phykit/services/tree/spectral_discordance.py` — core implementation
2. **Modify** `phykit/service_factories.py` — add LazyServiceFactory entry
3. **Modify** `phykit/cli_registry.py` — add aliases
4. **Modify** `phykit/phykit.py` — add @staticmethod handler, help banner entry, module-level wrapper
5. **Modify** `setup.py` — add console_scripts entry points
6. **Create** `tests/unit/services/tree/test_spectral_discordance.py` — tests
7. **Create** `tests/sample_files/gene_trees_simple.txt` — test gene tree file
8. **Modify** `docs/usage/index.rst` — command reference
9. **Modify** `docs/tutorials/index.rst` — tutorial section
