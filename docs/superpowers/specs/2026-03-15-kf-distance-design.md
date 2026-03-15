# Kuhner-Felsenstein (Branch Score) Distance

**Date:** 2026-03-15
**Status:** Approved
**Scope:** New `kf_distance` command

## Problem

PhyKIT has Robinson-Foulds distance for comparing tree topologies, but RF ignores branch lengths entirely. The Kuhner-Felsenstein (KF) distance (Kuhner & Felsenstein 1994), also called the branch score distance, measures topological *and* branch-length differences between two trees. It is widely used alongside RF and is available in R via `phangorn::KF.dist()`.

## Design

### Command interface

```
phykit kf_distance <tree_zero> <tree_one> [--json]
```

**Canonical method name:** `kf_distance` (the static method on the `Phykit` class).

**Aliases in `cli_registry.py`:** `kuhner_felsenstein_distance`, `kf_dist`, `kf` — all mapping to `"kf_distance"`. Note: `kf_distance` itself is NOT in the alias registry (it resolves directly via `hasattr`).

**Output (text):** `<plain_kf>\t<normalized_kf>`

**Output (--json):**
```json
{"plain_kf": 42.1234, "normalized_kf": 0.3456}
```

(JSON keys `plain_kf` / `normalized_kf` match the RF distance convention of `plain_rf` / `normalized_rf`.)

### Algorithm

The KF distance is defined as:

```
KF = sqrt( sum_over_all_splits( (b1_i - b2_i)^2 ) )
```

where `b1_i` and `b2_i` are the branch lengths for split `i` in tree 1 and tree 2 respectively. If a split is absent from one tree, its branch length is 0 in that tree.

**Splits include both internal and terminal (external) branches.** Each terminal branch induces a trivial split: `{taxon}` vs. the rest. This matches `phangorn::KF.dist()`, which sums over all branches.

**Branch length handling:** `None` branch lengths (from Newick trees without branch lengths) are treated as 0.0.

**Normalization:** `normalized_KF = KF / (total_branch_length_tree1 + total_branch_length_tree2)`. This is a convenience measure — it is *not* formally bounded in [0, 1] in all cases (unlike RF normalization which has a clean upper bound). If the denominator is 0 (both trees have zero-length branches), `normalized_KF` is reported as 0.0.

**Preprocessing** (same as RF distance):
1. Parse both trees
2. Prune to shared taxa (error if no shared taxa: "Trees share no common taxa.")
3. Root on the same taxon (first shared terminal). Rooting ensures bipartitions are computed in a consistent orientation. The root branch itself has no parent edge, so only the child branches of the root contribute splits.

**Bipartition extraction:**
For each branch in the tree (both internal and terminal, excluding the root node itself), extract:
- The split: a frozenset of tip names descending from that branch
- The branch length (or 0.0 if `None`)

Represent each tree's splits as a dict: `{frozenset_of_tips: branch_length}`.

Compute KF over the union of all splits from both trees.

### Implementation

**New file:** `phykit/services/tree/kf_distance.py`

```python
class KuhnerFelsensteinDistance(Tree):
    def __init__(self, args):
        # Same pattern as RobinsonFouldsDistance

    def run(self):
        # Parse, prune, root (reuse same logic as RF)
        # Error if no shared taxa
        # Extract splits with branch lengths
        # Compute KF distance

    def process_args(self, args):
        # tree_zero, tree_one, json_output

    def _get_splits_with_lengths(self, tree):
        """Return dict mapping frozenset(tip_names) -> branch_length.

        Includes both internal and terminal branches.
        None branch lengths are treated as 0.0.
        """

    def calculate_kf_distance(self, tree_zero, tree_one):
        """Return (plain_kf, normalized_kf) tuple."""
```

### Registration (5 files to update + module-level function)

1. **`phykit/service_factories.py`** — add `KuhnerFelsensteinDistance` factory
2. **`phykit/services/tree/__init__.py`** — add to `_EXPORTS`
3. **`phykit/cli_registry.py`** — add aliases: `kuhner_felsenstein_distance`, `kf_dist`, `kf` → `"kf_distance"`
4. **`phykit/phykit.py`** — add `kf_distance` static method with parser and help text; add entry to the main help text banner listing all tree commands; add module-level function `def kf_distance(argv=None): Phykit.kf_distance(sys.argv[1:])` at the bottom of the file
5. **`setup.py`** — add `pk_kf_distance`, `pk_kuhner_felsenstein_distance`, `pk_kf_dist`, `pk_kf` entry points referencing `phykit.phykit:kf_distance`

### R validation

**New file:** `tests/r_validation/validate_kf_distance.R`

Uses the same sample trees as the RF distance tests. Computes KF distance with `phangorn::KF.dist()` and prints expected values. The unit/integration tests assert PhyKIT output matches R output.

Script structure (following existing pattern — uses `suppressPackageStartupMessages`, relative paths from `tests/r_validation/`, prints `KEY VALUES FOR PYTHON TESTS` section):

```r
#!/usr/bin/env Rscript
# validate_kf_distance.R
#
# Compute KF (branch score) distance using phangorn::KF.dist()
# and compare against PhyKIT's kf_distance command.
#
# The KF distance (Kuhner & Felsenstein 1994) is:
#   KF = sqrt( sum_over_all_splits( (b1_i - b2_i)^2 ) )
# where b1_i and b2_i are the branch lengths for split i in
# each tree, and 0 is used for splits absent from a tree.
#
# Requires: ape, phangorn
#
# Usage:
#   cd tests/r_validation
#   Rscript validate_kf_distance.R

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})

t0 <- read.tree("../sample_files/tree_simple.tre")
t1 <- read.tree("../sample_files/tree_simple_other_topology.tre")

# Prune to shared taxa
shared <- intersect(t0$tip.label, t1$tip.label)
t0 <- keep.tip(t0, shared)
t1 <- keep.tip(t1, shared)

kf <- KF.dist(t0, t1)
total_bl <- sum(t0$edge.length) + sum(t1$edge.length)
normalized_kf <- as.numeric(kf) / total_bl

cat("=== KEY VALUES FOR PYTHON TESTS ===\n")
cat(sprintf("Number of shared taxa: %d\n", length(shared)))
cat(sprintf("Total branch length tree 0: %.6f\n", sum(t0$edge.length)))
cat(sprintf("Total branch length tree 1: %.6f\n", sum(t1$edge.length)))
cat(sprintf("KF distance (plain): %.4f\n", kf))
cat(sprintf("KF distance (normalized): %.4f\n", normalized_kf))
```

### Changelog entry

Following the established pattern (cross-validation against R package):

```
Added Kuhner-Felsenstein (branch score) distance command
(``kf_distance`` / ``kf``):

* Computes the KF distance between two phylogenies, incorporating
  both topology and branch length differences (Kuhner & Felsenstein
  1994)
* Reports plain and normalized KF distance
* Includes both internal and terminal branch lengths in the
  computation, matching the standard definition
* Supports ``--json`` output
* Cross-validated against R's phangorn::KF.dist(); R validation
  script provided in ``tests/r_validation/validate_kf_distance.R``
```

### Testing

**Unit tests** (`tests/unit/services/tree/test_kf_distance.py`):
- Test init and process_args
- Test `_get_splits_with_lengths` returns correct splits and lengths (both internal and terminal)
- Test `_get_splits_with_lengths` treats `None` branch lengths as 0.0
- Test `calculate_kf_distance` on identical trees (KF=0)
- Test `calculate_kf_distance` on different trees (match R output)
- Test normalized_kf division-by-zero guard (both trees zero-length)
- Test run with text and JSON output

**Integration tests** (`tests/integration/tree/test_kf_distance_integration.py`):
- End-to-end CLI test with sample trees (expected values from R validation)
- Test aliases (`kuhner_felsenstein_distance`, `kf_dist`, `kf`)
- Test JSON output
- Test trees with different taxa (pruning)

### Files changed

| File | Change |
|------|--------|
| `phykit/services/tree/kf_distance.py` | **New** — KF distance service |
| `phykit/service_factories.py` | Add factory |
| `phykit/services/tree/__init__.py` | Add to exports |
| `phykit/cli_registry.py` | Add aliases |
| `phykit/phykit.py` | Add CLI method, help text, main banner entry, module-level function |
| `setup.py` | Add entry points |
| `docs/usage/index.rst` | Add command documentation |
| `docs/change_log/index.rst` | Add changelog entry |
| `tests/unit/services/tree/test_kf_distance.py` | **New** — unit tests |
| `tests/integration/tree/test_kf_distance_integration.py` | **New** — integration tests |
| `tests/r_validation/validate_kf_distance.R` | **New** — R validation script |
