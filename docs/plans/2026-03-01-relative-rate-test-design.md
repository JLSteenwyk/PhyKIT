# Design: Relative Rate Test (`relative_rate_test`)

## Summary

Add Tajima's relative rate test (1993) to PhyKIT. The test determines whether
two lineages have evolved at equal rates since diverging, using a simple
model-free chi-squared test on site pattern counts. Supports single-alignment
and batch (multi-gene) modes with multiple testing correction.

## Motivation

Relative rate tests are simple but surprisingly hard to run at scale with
existing tools. MEGA requires manual 3-taxon selection through a GUI.
R's `pegas::rr.test()` handles one triplet at a time. No tool automatically
tests all pairwise ingroup combinations across hundreds of genes.

## Algorithm

### Tajima's Relative Rate Test (Equation 4)

Given three aligned sequences — two ingroup taxa (x, y) and one outgroup (o):

1. At each alignment column, classify the site (skip gaps/ambiguous characters):
   - x differs from o, but y matches o → **m1** (unique change on lineage x)
   - y differs from o, but x matches o → **m2** (unique change on lineage y)
   - Both differ from o, or neither differs → uninformative, skip

2. Test statistic: `chi2 = (m1 - m2)^2 / (m1 + m2)`, df = 1

3. P-value: `1 - chi2_cdf(chi2, df=1)`

4. Edge case: if m1 + m2 = 0, report chi2 = 0, p = 1.0

### Outgroup Determination

The outgroup is inferred from a rooted tree as the earliest-diverging taxon
(the taxon on the other side of the root from the largest clade). All
remaining taxa are treated as ingroup; all C(n-1, 2) pairwise combinations
are tested.

### Multiple Testing Correction

- **Bonferroni**: p_adj = min(p * n_tests, 1.0)
- **Benjamini-Hochberg FDR**: rank-based adjustment

Both are reported for each pairwise test.

## CLI Interface

```
# Single alignment
phykit relative_rate_test -a <alignment> -t <rooted_tree> [--json] [-v]

# Batch mode (multiple alignments, one shared tree)
phykit relative_rate_test -l <alignment_list> -t <rooted_tree> [--json] [-v]
```

**Aliases**: `relative_rate_test`, `rrt`, `tajima_rrt`

### Arguments

| Flag | Description |
|------|-------------|
| `-a/--alignment` | Single alignment file (FASTA, any BioPython-supported format) |
| `-l/--alignment-list` | File with one alignment path per line (batch mode) |
| `-t/--tree` | Rooted tree file (Newick) |
| `-v/--verbose` | Single: per-site detail. Batch: per-gene results |
| `--json` | JSON output |

`-a` and `-l` are mutually exclusive.

## Output

### Single mode (text)

```
Outgroup: O
Number of ingroup taxa: 5
Number of pairwise tests: 10
---
taxon1  taxon2  m1  m2  chi2     p_value   p_bonf   p_fdr    significant
A       B       45  12  19.105   1.24e-05  1.24e-04 6.20e-05 *
A       C       30  28  0.069    0.793     1.000    0.881
```

### Batch mode (text)

```
Number of alignments: 100
Outgroup: O
Number of ingroup taxa: 5
Number of pairwise tests per alignment: 10
---
taxon1  taxon2  n_reject  n_total  pct_reject  median_chi2
A       B       45        100      45.0        8.32
A       C       2         100      2.0         0.41
```

### JSON

Structured dict with all fields for programmatic consumption.

## Validation

### Primary: R's `pegas::rr.test()`

```r
library(pegas)
library(ape)
aln <- read.FASTA("test.fa")
result <- rr.test(aln[["A"]], aln[["B"]], aln[["outgroup"]])
cat("Chi-squared:", result$Chi, "\n")
cat("P-value:", result$Pval, "\n")
```

Compare m1, m2, chi2, and p-value to machine precision.

### Secondary: MEGA

Manual validation for a few test cases.

### Test cases

1. Equal rates (m1 ~ m2): p >> 0.05
2. Highly unequal rates: p << 0.05
3. All identical sequences: m1 = m2 = 0, p = 1.0
4. Gapped sites: verify consistent gap handling
5. Batch mode: verify cross-gene aggregation

## Files

### Create

- `phykit/services/tree/relative_rate_test.py` (~200-250 lines)
- `tests/unit/services/tree/test_relative_rate_test.py`
- `tests/integration/tree/test_relative_rate_test_integration.py`
- Test alignment + tree files in `tests/sample_files/`

### Modify

- `phykit/phykit.py` — handler + wrapper + help text
- `phykit/cli_registry.py` — aliases
- `phykit/service_factories.py` — lazy factory
- `phykit/services/tree/__init__.py` — exports
- `setup.py` — entry points
- `docs/usage/index.rst` — usage documentation
- `docs/change_log/index.rst` — changelog entry
- `phykit/version.py` — version bump

## References

- Tajima, F. (1993). Simple methods for testing the molecular evolutionary
  clock hypothesis. *Genetics*, 135, 599-607.
- R package `pegas`: `rr.test()` function
  (https://search.r-project.org/CRAN/refmans/pegas/html/rr.test.html)
