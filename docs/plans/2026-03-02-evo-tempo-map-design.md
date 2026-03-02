# Evolutionary Tempo Mapping — Design

## Goal

Detect rate-topology associations in phylogenomic datasets by comparing branch length distributions between concordant and discordant gene trees at each species tree branch.

## Motivation

Under the multispecies coalescent, discordant gene trees should have shorter internal branches near the discordant node (because coalescence happened deeper, in the ancestral population). Deviations from this expectation — e.g., anomalously long branches in discordant gene trees — suggest substitution rate heterogeneity correlated with topology, which could indicate adaptive evolution, different selective pressures in hybridizing lineages, or systematic error from model misspecification.

## Command

```
phykit evo_tempo_map -t <species_tree> -g <gene_trees> [--plot <output>] [-v] [--json]
```

Aliases: `evo_tempo_map`, `etm`

| Flag | Description |
|------|-------------|
| `-t/--tree` | Species tree (Newick) |
| `-g/--gene-trees` | Multi-Newick file of gene trees with branch lengths |
| `--plot` | Optional: output path for box/strip plot (PNG) |
| `-v/--verbose` | Print per-gene-tree classification details |
| `--json` | JSON output |

## Algorithm

### Per species tree branch

1. Extract the bipartition B_s defined by the branch (e.g., {cat,monkey}|{rest})
2. Compute the two NNI alternative bipartitions B_nni1, B_nni2
3. For each gene tree:
   - Extract all bipartitions (pruned to shared taxa)
   - If B_s present → concordant; extract that branch's length
   - If B_nni1 or B_nni2 present → discordant; extract that branch's length
   - Otherwise → uninformative for this branch (skip)
4. Collect concordant lengths L_c and discordant lengths L_d
5. Summary stats: mean, median, std for each group
6. Mann-Whitney U test on L_c vs L_d
7. Permutation test (1000 permutations): shuffle labels, recompute difference in medians

### Global summary

- Treeness (internal/total branch length) for concordant vs discordant gene trees
- Benjamini-Hochberg FDR correction across branches

### Concordance determination

Bipartition matching (same as gCF in concordance_asr): a gene tree is concordant at a species tree branch if it contains the same bipartition. Uses canonical split normalization.

### Branch mapping

For the branch length comparison, extract the branch in each gene tree that defines the matching bipartition (concordant) or the corresponding NNI alternative bipartition (discordant). This is the homologous internal branch.

## Output

### Text (default)

```
branch                  n_conc  n_disc  median_conc  median_disc  U_pval    perm_pval
{cat,monkey}            8       2       20.50        15.00        0.0312    0.0280
{sea_lion,seal}         7       3       7.50         6.50         0.4521    0.4600
...
---
Global treeness: concordant=0.145, discordant=0.098 (U p=0.023)
Branches tested: 5, significant (FDR<0.05): 1
```

### JSON

```json
{
  "branches": [
    {
      "split": ["cat", "monkey"],
      "n_concordant": 8,
      "n_discordant": 2,
      "concordant_lengths": {"mean": 20.5, "median": 20.5, "std": 0.35},
      "discordant_lengths": {"mean": 15.0, "median": 15.0, "std": 0.7},
      "mann_whitney_U": 2.0,
      "mann_whitney_p": 0.0312,
      "permutation_p": 0.028,
      "fdr_p": 0.156
    }
  ],
  "global": {
    "treeness_concordant": {"mean": 0.145, "median": 0.142},
    "treeness_discordant": {"mean": 0.098, "median": 0.095},
    "treeness_U_p": 0.023,
    "n_gene_trees": 10,
    "n_branches_tested": 5,
    "n_significant_fdr05": 1
  }
}
```

### Plot

Grouped box/strip plot: one panel per species tree branch (x-axis), concordant vs discordant branch lengths on y-axis, colored by group. Significant branches (FDR < 0.05) marked with asterisk.

## Edge Cases

- Gene tree missing taxa: prune to shared taxa; skip gene tree for a branch if < 2 taxa on either side after pruning
- No discordant gene trees for a branch: report concordant stats only, p-value = NA
- < 3 gene trees per group: report warning, still compute but note low power
- Gene trees without branch lengths: raise PhykitUserError
- Multifurcating gene trees: bipartition matching still works (bipartition either exists or doesn't)

## Verification Strategy

- Reuse and cross-check bipartition matching against `concordance_asr` gCF computation
- Reuse and cross-check treeness calculation against `treeness` command
- Construct hand-verified test cases with known concordant/discordant gene trees
- Test with identical gene trees → all concordant, no discordant, tests should be NA
- Test with known branch length differences → tests should detect them
