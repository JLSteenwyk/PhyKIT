# Quartet Pie Chart Visualization

**Date:** 2026-03-15
**Status:** Approved
**Scope:** New `quartet_pie` command

## Problem

Users want to visualize incomplete lineage sorting (ILS) and introgression signals on phylogenies using pie charts at internal nodes showing quartet support proportions — similar to what phytop does with ASTRAL output. PhyKIT already computes gene concordance factors (gCF/gDF1/gDF2) in `concordance_asr` and can draw pie charts in `ancestral_reconstruction`, but there is no dedicated command that combines these into a single visualization.

## Design

### Command interface

```
phykit quartet_pie -t <species_tree> [-g <gene_trees>] -o <output>
    [--annotate] [--json] [plot options...]
```

**Aliases:** `quartet_pie`, `qpie`, `quartet_pie_chart`

**Input modes:**

1. **Native mode** (`-g` provided): PhyKIT computes gCF/gDF1/gDF2 from species tree + gene trees via bipartition matching. Uses the four-group decomposition (C1, C2, S, D) around each internal branch, same algorithm as `concordance_asr._compute_gcf_per_node()`.

2. **ASTRAL mode** (no `-g`): Parses q1/q2/q3 annotations from ASTRAL's `-t 2` output. ASTRAL annotates internal nodes with `[q1=0.5;q2=0.3;q3=0.2;...]` or `'[q1=0.5;q2=0.3;q3=0.2]'` in the Newick string. The parser strips quotes and brackets and extracts `q1`, `q2`, `q3` values. If no ASTRAL annotations are found and no `-g` is provided, exit with a clear error.

**Output:** Figure file (format from extension: `.png`, `.pdf`, `.svg`, `.jpg`). Optional `--json` output with per-node concordance values.

### Visualization

- Horizontal phylogram with branch lengths to scale
- Pie chart at each internal node (excluding root) showing three proportions:
  - **Concordant / q1** — default blue (`#2b8cbe`)
  - **Discordant alt 1 / q2** — default red (`#d62728`)
  - **Discordant alt 2 / q3** — default gray (`#969696`)
- Colors overridable via `--colors "#2b8cbe,#d62728,#969696"`
- Tip labels on the right
- Legend showing what each color represents
- Pie radius scales with figure width (`max_x * 0.02`)
- `--annotate` flag adds numeric gCF/gDF1/gDF2 text near each pie
- All shared `PlotConfig` options available

### Algorithm (native mode)

For each internal branch in the species tree:
1. Identify the four groups: C1, C2 (children), S (sister), D (rest)
2. Construct the three possible bipartitions:
   - Concordant: C1∪C2 | S∪D
   - NNI alt 1: S∪C2 | C1∪D
   - NNI alt 2: C1∪S | C2∪D
3. Count how many gene trees contain each bipartition
4. Compute proportions: gCF, gDF1, gDF2 (sum to 1.0)

This is the same algorithm used in `concordance_asr` and `discordance_asymmetry`. The implementation will be extracted into a standalone function in `phykit/helpers/quartet_utils.py` to avoid importing the full `ConcordanceAsr` class.

### Algorithm (ASTRAL mode)

For each internal node in the species tree:
1. Parse the node label/comment for `q1=`, `q2=`, `q3=` annotations
2. Strip surrounding quotes and brackets
3. Extract float values
4. If any node lacks annotations, skip it (no pie drawn)

### Validation

**Ground truth:** The gCF/gDF computation will be validated against:
1. **Internal consistency:** The existing `concordance_asr` command computes the same values — the `quartet_pie` native mode must produce identical results on the same input.
2. **ASTRAL cross-validation:** When ASTRAL is available, running `astral -i gene_trees.nwk -t 2` produces q1/q2/q3 annotations. These should agree with PhyKIT's gCF/gDF1/gDF2 within sampling tolerance.

**Sample data ground truth** (from `tree_simple.tre` + `gene_trees_simple.nwk`):

| Node | gCF | gDF1 | gDF2 | Concordant | Disc1 | Disc2 |
|------|-----|------|------|------------|-------|-------|
| bear+raccoon (2 taxa) | 1.000 | 0.000 | 0.000 | 7 | 0 | 0 |
| sea_lion+seal (2 taxa) | 1.000 | 0.000 | 0.000 | 9 | 0 | 0 |
| monkey+cat (2 taxa) | 1.000 | 0.000 | 0.000 | 10 | 0 | 0 |
| monkey+cat+weasel (3 taxa) | 0.474 | 0.526 | 0.000 | 9 | 10 | 0 |
| 5-taxon clade | 0.250 | 0.375 | 0.375 | 6 | 9 | 9 |

### Implementation

**New file:** `phykit/helpers/quartet_utils.py`
- `compute_gcf_per_node(species_tree, gene_trees)` — returns `Dict[int, Tuple[float, float, float]]` mapping clade id to (gCF, gDF1, gDF2)
- `parse_astral_annotations(tree)` — returns same dict from ASTRAL node labels
- `canonical_split(tips, all_taxa)` — normalize bipartition representation

**New file:** `phykit/services/tree/quartet_pie.py`
- `QuartetPie(Tree)` — service class
- `_plot_quartet_pie(tree, proportions, output_path)` — phylogram + pie charts
- Uses `PlotConfig` for all visual parameters
- Tree layout reuses the preorder traversal + node positioning approach from `discordance_asymmetry._plot()`
- Pie rendering uses `matplotlib.patches.Wedge` (same as `ancestral_reconstruction._plot_discrete_asr()`)

### Registration

1. **`phykit/service_factories.py`** — add `QuartetPie` factory
2. **`phykit/services/tree/__init__.py`** — add to `_EXPORTS`
3. **`phykit/cli_registry.py`** — add aliases: `qpie`, `quartet_pie_chart` → `"quartet_pie"`
4. **`phykit/phykit.py`** — add `quartet_pie` static method, help banner, module-level function
5. **`setup.py`** — add `pk_quartet_pie`, `pk_qpie` entry points

### JSON output

```json
{
  "n_taxa": 8,
  "n_gene_trees": 10,
  "input_mode": "native",
  "nodes": [
    {
      "node_tips": ["bear", "raccoon"],
      "gCF": 1.0,
      "gDF1": 0.0,
      "gDF2": 0.0,
      "concordant_count": 7,
      "disc1_count": 0,
      "disc2_count": 0
    },
    ...
  ],
  "output_file": "quartet_pie.png"
}
```

### Changelog entry

```
Added quartet pie chart visualization command (``quartet_pie`` / ``qpie``):

* Draws a phylogram with pie charts at internal nodes showing gene
  concordance (gCF) and discordance (gDF1, gDF2) proportions
* Native mode: computes quartet proportions from species tree +
  gene trees via bipartition matching
* ASTRAL mode: parses q1/q2/q3 annotations from ASTRAL -t 2 output
* Optional ``--annotate`` flag adds numeric values near each pie
* Default colors: blue (concordant), red (discordant alt 1),
  gray (discordant alt 2); overridable via ``--colors``
* Supports all shared plot options and ``--json`` output
* Validated against concordance_asr gCF values
```

### Testing

**Unit tests** (`tests/unit/helpers/test_quartet_utils.py`):
- Test `compute_gcf_per_node` against ground truth table above
- Test `parse_astral_annotations` with multiple ASTRAL label formats
- Test `canonical_split` correctness

**Unit tests** (`tests/unit/services/tree/test_quartet_pie.py`):
- Test init and process_args
- Test ASTRAL mode detection (no `-g` flag)
- Test native mode (with `-g` flag)
- Test plot file creation
- Test JSON output structure
- Test `--annotate` flag

**Integration tests** (`tests/integration/tree/test_quartet_pie_integration.py`):
- End-to-end CLI with sample trees
- Test aliases
- Test JSON output values match ground truth table

### Files changed

| File | Change |
|------|--------|
| `phykit/helpers/quartet_utils.py` | **New** — gCF computation and ASTRAL parsing |
| `phykit/services/tree/quartet_pie.py` | **New** — QuartetPie service with pie chart plotting |
| `phykit/service_factories.py` | Add factory |
| `phykit/services/tree/__init__.py` | Add to exports |
| `phykit/cli_registry.py` | Add aliases |
| `phykit/phykit.py` | Add CLI method, help text, help banner, module-level function |
| `setup.py` | Add entry points |
| `docs/usage/index.rst` | Add command documentation |
| `docs/change_log/index.rst` | Add changelog entry |
| `tests/unit/helpers/test_quartet_utils.py` | **New** — unit tests for shared module |
| `tests/unit/services/tree/test_quartet_pie.py` | **New** — unit tests for service |
| `tests/integration/tree/test_quartet_pie_integration.py` | **New** — integration tests |
