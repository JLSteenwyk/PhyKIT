# Character Map Design Spec

## Overview

A new PhyKIT command (`character_map` / `charmap` / `synapomorphy_map`) that maps
synapomorphies and homoplasies onto a phylogeny using Fitch parsimony with
ACCTRAN/DELTRAN optimization. Produces a cladogram (or phylogram) with color-coded
circles on each branch showing where character state changes occur.

Note: `cmap` is already taken by `cont_map` and must not be used as an alias.

Motivated by a user request from an entomologist/systematist who needs a modern
replacement for Winclada-style character mapping on phylogenies.

## CLI Interface

```
phykit character_map -t <tree> -d <data> -o <output>
    [--optimization acctran|deltran]
    [--phylogram]
    [--characters 1,3,7,12]
    [--verbose]
    [--json]
    [shared plot options: --fig-width, --fig-height, --dpi, --colors,
     --ladderize, --no-title, --title, --legend-position,
     --ylabel-fontsize, --title-fontsize, --axis-fontsize]
```

**Aliases:** `character_map`, `charmap`, `synapomorphy_map`

### Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `-t/--tree` | Yes | — | Newick tree file |
| `-d/--data` | Yes | — | TSV character matrix (header + rows). Named `-d/--data` rather than `--trait-data` because this is a multi-character matrix, not a single-trait file. |
| `-o/--output` | Yes | — | Output figure path (.png, .pdf, .svg) |
| `--optimization` | No | `acctran` | ACCTRAN or DELTRAN ancestral state optimization |
| `--phylogram` | No | off | Draw phylogram instead of cladogram |
| `--characters` | No | all | Comma-separated character indices to display (0-based, matching the TSV header) |
| `--verbose` | No | off | Print per-character detail |
| `--json` | No | off | Output structured JSON |

### Input Format

**Character matrix (TSV):**

```
taxon	char0	char1	char2	char3	char4
Taxon_A	0	1	0	1	2
Taxon_B	0	1	1	1	0
Taxon_C	1	0	0	0	2
Taxon_D	1	0	1	0	1
```

- First row: header with character names (or indices)
- Subsequent rows: `taxon_name<TAB>state0<TAB>state1<TAB>...`
- States: any discrete values (0, 1, 2, etc.)
- Missing data: `?` or `-` treated as wildcards (match any state)

## Validation & Preprocessing

Before running the algorithm:

1. **Taxon matching:** Compute the intersection of tree tips and matrix taxa.
   Prune the tree to shared taxa and filter the matrix to shared taxa.
   Require at least 3 shared taxa (raise `PhykitUserError` otherwise).
2. **Polytomy resolution:** Resolve any multifurcations by inserting
   zero-length branches (matching the existing pattern in `parsimony_score.py`
   `_resolve_polytomies()`). This ensures the Fitch algorithm and change
   detection operate on a bifurcating tree.
3. **Ladderization:** If `--ladderize` is set, ladderize the tree after
   polytomy resolution and before plotting. Ladderization does not affect
   parsimony computation.
4. **Parsimony-informative characters:** A character is parsimony-informative
   if at least 2 states each occur in at least 2 taxa. Uninformative
   characters are still included in the analysis but reported separately
   in the summary.

## Architecture

### New Files

| File | Purpose |
|------|---------|
| `phykit/helpers/parsimony_utils.py` | Generalized Fitch downpass, ACCTRAN/DELTRAN uppass, change detection, synapomorphy classification, CI/RI |
| `phykit/services/tree/character_map.py` | Service class: orchestration and plotting |
| `tests/unit/helpers/test_parsimony_utils.py` | Unit tests for parsimony utilities |
| `tests/unit/services/tree/test_character_map.py` | Unit tests for service |
| `tests/integration/tree/test_character_map_integration.py` | CLI integration tests |
| `tests/r_validation/validate_character_map.R` | Cross-validation against phangorn |
| `tests/sample_files/character_matrix_simple.tsv` | Test dataset |

### Modified Files

| File | Change |
|------|--------|
| `phykit/phykit.py` | Add command handler + help text |
| `phykit/cli_registry.py` | Add aliases |
| `phykit/service_factories.py` | Add lazy factory |
| `phykit/services/tree/__init__.py` | Export CharacterMap |
| `docs/usage/index.rst` | Documentation |
| `docs/change_log/index.rst` | Changelog entry |

### Component Separation

**`parsimony_utils.py`** contains reusable algorithms with no plotting or I/O:

- `fitch_downpass(tree, tip_states)` → per-node state sets, parsimony score
- `fitch_uppass_acctran(tree, node_state_sets, parent_map)` → per-node final states
- `fitch_uppass_deltran(tree, node_state_sets, parent_map)` → per-node final states
- `detect_changes(tree, node_states, parent_map)` → per-branch list of (char_idx, old, new)
- `classify_changes(tree, branch_changes, node_states)` → each change tagged as synapomorphy/convergence/reversal
- `consistency_index(n_states_per_char, observed_changes)` → per-char and overall CI
- `retention_index(tip_states, observed_changes)` → per-char and overall RI
- `resolve_polytomies(tree)` → insert zero-length branches at multifurcations
- `build_parent_map(tree)` → dict mapping node id → parent clade

**`character_map.py`** handles I/O, argument parsing, orchestration, and plotting.

## Algorithm Detail

### Step 1: Generalized Fitch Downpass

For each character independently, postorder traversal:

```
for each node (leaves first):
    if terminal:
        if state is wildcard (? or -):
            state_set = set of all observed states for this character
        else:
            state_set = {observed_state}
    else:
        child_sets = [state_set of each child]
        intersection = intersection of all child_sets
        if intersection is non-empty:
            state_set = intersection  (no cost)
        else:
            state_set = union of all child_sets  (cost += 1)
```

Returns: `Dict[node_id → List[Set[str]]]` (per-character state sets at each node)
and `List[int]` (parsimony score per character).

### Step 2: ACCTRAN Uppass

Preorder traversal (root first):

```
root: final_state = min(root_state_set)  # deterministic tie-breaking

for each node in preorder (skip root):
    parent_state = final_state[parent]
    if parent_state in node_state_set:
        final_state[node] = parent_state  # inherit, no change
    else:
        final_state[node] = min(node_state_set)  # change happens here
```

ACCTRAN pushes changes rootward: if the parent's state is not compatible
with the child, the change is on the branch to that child. This favors
reversals over convergence.

### Step 3: DELTRAN Uppass

Single preorder pass following Swofford & Maddison (1987):

```
root: final_state = min(root_state_set)

for each node in preorder (skip root):
    parent_final = final_state[parent]
    if parent_final in node_state_set:
        # Parent's state is compatible — inherit it (delay the change)
        final_state[node] = parent_final
    else:
        # Parent's state is NOT compatible — change must occur here.
        # Prefer a state from the intersection of node's downpass set
        # and the parent's downpass set (to maintain continuity with
        # the parent's lineage). If that intersection is empty,
        # pick from the node's own downpass set.
        shared = node_state_set & parent_state_set
        if shared:
            final_state[node] = min(shared)
        else:
            final_state[node] = min(node_state_set)
```

DELTRAN pushes changes tipward: changes are delayed as long as possible.
This favors convergence over reversals. The key difference from ACCTRAN
is that DELTRAN uses the parent's *downpass state set* (not just its
final state) when resolving conflicts, which allows it to choose states
that defer the change.

### Step 4: Change Detection

For each branch (parent → child), for each character:

```
if final_state[parent][char] != final_state[child][char]:
    record (char_index, old_state, new_state) on this branch
```

### Step 5: Synapomorphy vs Homoplasy Classification

For each (character, new_state) pair across all branches:

```
count = number of branches where character transitions TO new_state
if count == 1:
    classification = "synapomorphy"
else:
    # Determine if convergence or reversal by checking whether
    # new_state was the ancestral state at any ancestor of this branch.
    # Walk from the parent up to the root; if any ancestor's final state
    # for this character equals new_state, it is a reversal.
    is_reversal = False
    ancestor = parent
    while ancestor is not None:
        if final_state[ancestor][character] == new_state:
            is_reversal = True
            break
        ancestor = parent_map[ancestor]  # walk toward root

    if is_reversal:
        classification = "reversal"
    else:
        classification = "convergence"
```

This correctly handles cases where a reversal returns to a non-root
ancestral state (e.g., 0→1→2→1 where the reversal to 1 is not the
root state 0).

### Step 6: CI and RI

**Per-character consistency index:**
```
ci_i = min_changes_i / observed_changes_i
```
Where `min_changes_i = number_of_states_i - 1` (minimum possible on any tree).

**Per-character retention index (Farris 1989):**
```
ri_i = (max_changes_i - observed_i) / (max_changes_i - min_changes_i)
```
Where `max_changes_i = n_taxa - f_max_i` and `f_max_i` is the number of
tips with the most frequent state for character `i`. This is the standard
formula: the maximum number of steps any character can require equals the
number of taxa minus the count of the most common state.

Handle division by zero: if `max_changes_i == min_changes_i`, RI is
undefined for that character (set to `None`).

**Overall (ensemble):**
```
CI = sum(min_i) / sum(observed_i)
RI = sum(max_i - observed_i) / sum(max_i - min_i)
```

Uninformative characters (where `max == min`) are excluded from the
ensemble RI denominator.

## Visualization

### Layout

- **Default: cladogram** — branches stretched so all tips align at the right edge
- **Optional: phylogram** (`--phylogram`) — branch length proportional to actual values
- Standard PhyKIT phylogram drawing code (node_x, node_y computation, horizontal + vertical lines)

### Cladogram Branch Length Computation

All tips are placed at the same x-coordinate (`max_x`). Internal nodes
are placed based on their depth from the root:

```
max_depth = max depth (in edges) from root to any tip
step = max_x / max_depth

for each node in preorder:
    if root:
        x = 0
    elif terminal:
        x = max_x  # all tips aligned at right edge
    else:
        # internal node: placed at its topological depth
        x = depth_from_root * step
```

### Character Change Circles

For each branch with character changes:

- Circles placed at evenly-spaced positions along the branch
- Circle radius: auto-scaled based on number of tips and figure size
- Fill color by classification:
  - Synapomorphy: blue (`#2b8cbe`, default)
  - Convergence: red (`#d62728`, default)
  - Reversal: gray (`#969696`, default)
  - Customizable via `--colors "blue,red,gray"`
- Character number rendered above the circle (small font)
- State transition rendered below the circle (e.g., "0→1")

### Character Filtering

`--characters 1,3,7,12` displays only the specified character indices.
All characters are still used for parsimony computation and CI/RI;
the filter only affects which changes are drawn.

### Legend

Located at the upper right (default), showing:
- Blue circle: Synapomorphy
- Red circle: Convergence
- Gray circle: Reversal

Customizable position via `--legend-position`.

## Text/JSON Output

### Default Text Output

```
Character Map (ACCTRAN)

Characters: 10
Parsimony-informative: 7
Tree length (total steps): 14
Consistency index (CI): 0.833
Retention index (RI): 0.750

Figure saved: output.png
```

### Verbose Output (`--verbose`)

Adds per-character detail:

```
Character 0: steps=2, CI=0.500, RI=0.333
  Branch (root → node1): 0 → 1 [synapomorphy]
  Branch (node2 → Taxon_C): 0 → 1 [convergence]

Character 1: steps=1, CI=1.000, RI=1.000
  Branch (root → node1): 0 → 1 [synapomorphy]
...
```

### JSON Output (`--json`)

```json
{
  "n_characters": 10,
  "n_informative": 7,
  "tree_length": 14,
  "ci": 0.833,
  "ri": 0.750,
  "optimization": "acctran",
  "output_file": "output.png",
  "characters": [
    {
      "index": 0,
      "name": "char0",
      "steps": 2,
      "ci": 0.500,
      "ri": 0.333,
      "changes": [
        {"branch": "root → node1", "from": "0", "to": "1", "type": "synapomorphy"},
        {"branch": "node2 → Taxon_C", "from": "0", "to": "1", "type": "convergence"}
      ]
    }
  ]
}
```

## Validation

Cross-validate against R's phangorn package:

1. Create a small test dataset (6 taxa, 10 characters with known properties)
2. R script (`tests/r_validation/validate_character_map.R`) computes:
   - ACCTRAN ancestral states via `ancestral.pars(tree, data, type="ACCTRAN")`
   - Per-branch state changes
   - CI via `CI(tree, data)`
   - RI via `RI(tree, data)`
3. PhyKIT JSON output compared against R output
4. Must match exactly on: ancestral states, change list, CI, RI

## Test Plan

### Unit Tests (`test_parsimony_utils.py`)

- `test_fitch_downpass_simple` — 4-taxon tree, known state sets
- `test_fitch_downpass_with_wildcards` — missing data handled correctly
- `test_fitch_downpass_multistate` — characters with 3+ states
- `test_acctran_uppass_basic` — verify rootward change placement
- `test_deltran_uppass_basic` — verify tipward change placement
- `test_acctran_vs_deltran_differ` — same data, different change placement
- `test_detect_changes_counts` — correct number of changes per branch
- `test_classify_synapomorphy` — unique change classified correctly
- `test_classify_convergence` — repeated gain classified correctly
- `test_classify_reversal` — return to ancestral (non-root) state classified correctly
- `test_ci_perfect` — no homoplasy → CI = 1.0
- `test_ci_with_homoplasy` — CI < 1.0
- `test_ci_multistate` — min_changes = n_states - 1
- `test_ri_computation` — known RI values using Farris formula
- `test_ri_uninformative_character` — handle division by zero
- `test_resolve_polytomies` — multifurcation correctly resolved

### Unit Tests (`test_character_map.py`)

- `test_init_sets_fields` — argument parsing
- `test_parse_character_matrix` — TSV parsing, missing data
- `test_creates_png` — plot file created
- `test_creates_pdf` — PDF output works
- `test_cladogram_mode` — default layout
- `test_phylogram_mode` — `--phylogram` flag
- `test_character_filter` — `--characters` flag
- `test_json_output` — JSON structure and values
- `test_too_few_tips` — validation error
- `test_taxon_mismatch_prunes` — tree/data mismatch handled by pruning

### Integration Tests (`test_character_map_integration.py`)

- `test_basic_invocation` — end-to-end with sample data
- `test_alias_charmap` — alias works
- `test_alias_synapomorphy_map` — alias works
- `test_json_ground_truth` — JSON values match expected
- `test_deltran_option` — `--optimization deltran`
- `test_ladderize_flag` — `--ladderize` works
- `test_character_filter` — `--characters 1,3` shows subset

### R Validation (`validate_character_map.R`)

- Load same tree + matrix used in unit tests
- Compare CI, RI, ancestral states, change counts
- Assert exact match
