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

**Aliases:** `kf_distance`, `kuhner_felsenstein_distance`, `kf_dist`, `kf`

**Output (text):** `<kf_distance>\t<normalized_kf_distance>`

**Output (--json):**
```json
{"kf_distance": 42.1234, "normalized_kf_distance": 0.3456}
```

### Algorithm

The KF distance is defined as:

```
KF = sqrt( sum_over_all_splits( (b1_i - b2_i)^2 ) )
```

where `b1_i` and `b2_i` are the branch lengths for split `i` in tree 1 and tree 2 respectively. If a split is absent from one tree, its branch length is 0 in that tree.

**Normalization:** `normalized_KF = KF / (sum_of_all_branch_lengths_tree1 + sum_of_all_branch_lengths_tree2)`. This gives a 0-1 measure where 0 means identical trees.

**Preprocessing** (same as RF distance):
1. Parse both trees
2. Prune to shared taxa
3. Root on the same taxon (first shared terminal)

**Bipartition extraction:**
For each internal branch (non-terminal, non-root), extract:
- The split: a frozenset of tip names on one side
- The branch length

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
        # Extract splits with branch lengths
        # Compute KF distance

    def process_args(self, args):
        # tree_zero, tree_one, json_output

    def _get_splits_with_lengths(self, tree):
        """Return dict mapping frozenset(tip_names) -> branch_length for each internal branch."""

    def calculate_kf_distance(self, tree_zero, tree_one):
        """Return (kf, normalized_kf) tuple."""
```

### Registration (5 files to update)

1. **`phykit/service_factories.py`** — add `KuhnerFelsensteinDistance` factory
2. **`phykit/services/tree/__init__.py`** — add to `_EXPORTS`
3. **`phykit/cli_registry.py`** — add aliases
4. **`phykit/phykit.py`** — add static method with parser and help text
5. **`setup.py`** — add `pk_kf_distance`, `pk_kf_dist`, `pk_kf` entry points

### R validation

**New file:** `tests/r_validation/validate_kf_distance.R`

Uses the same sample trees as the RF distance tests. Computes KF distance with `phangorn::KF.dist()` and prints expected values. The unit/integration tests assert PhyKIT output matches R output.

Script structure (following existing pattern):
```r
#!/usr/bin/env Rscript
# validate_kf_distance.R
#
# Compute KF (branch score) distance using phangorn::KF.dist()
# and compare against PhyKIT's kf_distance command.
#
# Requires: ape, phangorn
#
# Usage:
#   Rscript tests/r_validation/validate_kf_distance.R

library(ape)
library(phangorn)

t0 <- read.tree("tests/sample_files/tree_simple.tre")
t1 <- read.tree("tests/sample_files/tree_simple_other_topology.tre")

# Prune to shared taxa
shared <- intersect(t0$tip.label, t1$tip.label)
t0 <- keep.tip(t0, shared)
t1 <- keep.tip(t1, shared)

kf <- KF.dist(t0, t1)
cat(sprintf("KF distance: %.4f\n", kf))
```

### Changelog entry

Following the established pattern (cross-validation against R package):

```
**2.1.40**:
Added Kuhner-Felsenstein (branch score) distance command
(``kf_distance`` / ``kf``):

* Computes the KF distance between two phylogenies, incorporating
  both topology and branch length differences
* Reports plain and normalized KF distance
* Supports ``--json`` output
* Cross-validated against R's phangorn::KF.dist(); R validation
  script provided in ``tests/r_validation/validate_kf_distance.R``
```

### Testing

**Unit tests** (`tests/unit/services/tree/test_kf_distance.py`):
- Test init and process_args
- Test `_get_splits_with_lengths` returns correct splits and lengths
- Test `calculate_kf_distance` on identical trees (KF=0)
- Test `calculate_kf_distance` on different trees (match R output)
- Test run with text and JSON output

**Integration tests** (`tests/integration/tree/test_kf_distance_integration.py`):
- End-to-end CLI test with sample trees
- Test aliases
- Test JSON output
- Expected values validated against `phangorn::KF.dist()`

### Files changed

| File | Change |
|------|--------|
| `phykit/services/tree/kf_distance.py` | **New** — KF distance service |
| `phykit/service_factories.py` | Add factory |
| `phykit/services/tree/__init__.py` | Add to exports |
| `phykit/cli_registry.py` | Add aliases |
| `phykit/phykit.py` | Add CLI method with help text |
| `setup.py` | Add entry points |
| `docs/usage/index.rst` | Add command documentation |
| `docs/change_log/index.rst` | Add changelog entry |
| `tests/unit/services/tree/test_kf_distance.py` | **New** — unit tests |
| `tests/integration/tree/test_kf_distance_integration.py` | **New** — integration tests |
| `tests/r_validation/validate_kf_distance.R` | **New** — R validation script |
