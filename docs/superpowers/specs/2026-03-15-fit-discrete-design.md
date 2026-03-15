# fitDiscrete Command + Shared Discrete Model Utilities

**Date:** 2026-03-15
**Status:** Approved
**Scope:** Extract shared discrete trait code, add `fit_discrete` command

## Problem

1. PhyKIT has `fit_continuous` for comparing continuous trait evolution models but no equivalent for discrete traits. R's `geiger::fitDiscrete()` fills this role and is widely used.

2. The Q-matrix fitting code (`_build_q_matrix`, `_felsenstein_pruning`, `_fit_q_matrix`, `_matrix_exp`) is duplicated between `stochastic_character_map.py` and `ancestral_reconstruction.py`. Adding a third copy in `fit_discrete.py` would worsen the duplication.

## Design

### Part 1: Extract shared discrete model utilities

**New file:** `phykit/helpers/discrete_models.py`

Extract these functions from `stochastic_character_map.py` (they are identical in `ancestral_reconstruction.py`):

| Function | Source (SCM lines) | Description |
|----------|-------------------|-------------|
| `build_q_matrix(params, k, model)` | 216-241 | Build Q-matrix from parameters for ER/SYM/ARD |
| `matrix_exp(Q, t)` | 243-244 | Wrapper around `scipy.linalg.expm(Q * t)` |
| `felsenstein_pruning(tree, tip_states, Q, pi, states)` | 246-276 | Postorder traversal computing conditional likelihoods + log-likelihood |
| `fit_q_matrix(tree, tip_states, states, model)` | 278-338 | Multi-start ML optimization of Q-matrix parameters |

These become **module-level functions** (not methods) since they don't need `self`. They only use `numpy`, `scipy.linalg.expm`, and `scipy.optimize.minimize`.

**Update `stochastic_character_map.py`:** Replace the duplicated methods with imports from `discrete_models.py`. The methods become thin wrappers calling the shared functions. For example:

```python
from ...helpers.discrete_models import fit_q_matrix

# In _fit_q_matrix method:
def _fit_q_matrix(self, tree, tip_states, states, model):
    return fit_q_matrix(tree, tip_states, states, model)
```

Or simply call the functions directly in `run()`.

**Update `ancestral_reconstruction.py`:** Same treatment — replace duplicated methods with imports.

**Backward compatibility:** The public API (CLI commands) doesn't change. Only internal method calls are redirected.

### Part 2: `fit_discrete` command

**Command interface:**

```
phykit fit_discrete -t <tree> -d <trait_data> -c <trait_column>
    [--models ER,SYM,ARD] [--json]
```

**Aliases:** `fit_discrete`, `fitdiscrete`, `fd`

**Output (text):** Model comparison table (following `fit_continuous` format):

```
Model   lnL        AIC      ΔAIC    AIC_w   BIC      ΔBIC    n_params  n_states
ER      -12.3456   26.6912  0.0000  0.7234  27.3844  0.0000  1         3
SYM     -11.8901   29.7802  3.0890  0.1533  31.8598  4.4754  3         3
ARD     -10.5678   33.1356  6.4444  0.0293  37.2948  9.9104  6         3
```

Where:
- `lnL`: log-likelihood
- `AIC`: Akaike Information Criterion = -2*lnL + 2*n_params
- `ΔAIC`: difference from best model's AIC
- `AIC_w`: Akaike weights = exp(-ΔAIC/2) / sum(exp(-ΔAIC_j/2))
- `BIC`: Bayesian Information Criterion = -2*lnL + n_params*ln(n_obs)
- `ΔBIC`: difference from best model's BIC
- `n_params`: number of free rate parameters in the Q-matrix
- `n_states`: number of discrete states

`n_obs` for BIC is the number of tips that have trait data AND are present in the tree (the intersection count), matching `fit_continuous` convention.

Parameter counts assume a fixed equal root prior (pi = 1/k). Only Q-matrix rate parameters are counted:
- ER: 1 parameter
- SYM: k(k-1)/2 parameters
- ARD: k(k-1) parameters

**Output (--json):**
```json
{
  "n_taxa": 8,
  "n_states": 3,
  "states": ["carnivore", "herbivore", "omnivore"],
  "models": [
    {
      "model": "ER",
      "lnL": -12.3456,
      "aic": 26.6912,
      "delta_aic": 0.0,
      "aic_weight": 0.7234,
      "bic": 27.3844,
      "delta_bic": 0.0,
      "n_params": 1,
      "q_matrix": [[...], [...], [...]],
      "rates": {"all": 0.1234}
    },
    ...
  ]
}
```

**Default models:** `ER,SYM,ARD` (all three). User can select a subset with `--models ER,ARD`.

**Model validation:** Invalid model names in `--models` produce a clear error at argument parsing time (not deep in optimization).

### Algorithm

For each model (ER, SYM, ARD):
1. Build initial Q-matrix from parameters
2. Optimize parameters via ML (multi-start, `scipy.optimize.minimize`)
3. Compute log-likelihood via Felsenstein pruning
4. Compute AIC = -2*lnL + 2*k_params
5. Compute BIC = -2*lnL + k_params*ln(n_taxa) (n_taxa = number of tips with data)

After fitting all models:
1. Compute ΔAIC = AIC_i - min(AIC)
2. Compute Akaike weights: w_i = exp(-ΔAIC_i/2) / sum(exp(-ΔAIC_j/2))
3. Sort by AIC (best first)

### Implementation

**New file:** `phykit/services/tree/fit_discrete.py`

```python
class FitDiscrete(Tree):
    def __init__(self, args):
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.trait_data_path = parsed["trait_data_path"]
        self.trait_column = parsed["trait_column"]
        self.selected_models = parsed["models"]
        self.json_output = parsed["json_output"]

    def run(self):
        tree = self.read_tree_file()
        tip_states = self._parse_trait_data(...)
        states = sorted(set(tip_states.values()))

        results = []
        for model_name in self.selected_models:
            result = self._fit_model(tree, tip_states, states, model_name)
            results.append(result)

        results = self._compute_model_comparison(results, len(tip_states))
        # Output text or JSON

    def _fit_model(self, tree, tip_states, states, model_name):
        Q, lnL = fit_q_matrix(tree, tip_states, states, model_name)
        n_params = _count_params(len(states), model_name)
        aic = -2 * lnL + 2 * n_params
        return {"model": model_name, "lnL": lnL, "aic": aic, ...}

    def _compute_model_comparison(self, results, n_obs):
        # ΔAIC, AIC weights, BIC, sort by AIC
```

### Trait parsing

There are three discrete trait parsers across the codebase:
1. `StochasticCharacterMap._parse_discrete_trait_file(path, column, tree_tips)` — multi-column TSV with header
2. `AncestralReconstruction._parse_discrete_trait_data_single(path, tree_tips)` — two-column TSV, no header
3. `AncestralReconstruction._parse_discrete_trait_data_multi(path, tree_tips, trait_column)` — multi-column TSV with header

Extract a single unified function into `discrete_models.py`:
```python
def parse_discrete_traits(path, tree_tips, trait_column=None):
    """Parse discrete trait data from a TSV file.

    If trait_column is None: expects 2-column format (taxon, state).
    If trait_column is given: expects multi-column with header, extracts named column.
    Returns {taxon: state} dict for shared taxa (intersection with tree_tips).
    """
```

`StochasticCharacterMap`, `AncestralReconstruction`, and `DensityMap` will all import from this shared function. `DensityMap` currently reaches into SCM internals (`scm._parse_discrete_trait_file`), so it must also be updated.

### Tree validation

`FitDiscrete.run()` must validate the tree before fitting: at least 3 tips, all branches must have lengths. This follows the pattern in `StochasticCharacterMap._validate_tree` and `FitContinuous._validate_tree`.

### Registration

1. **`phykit/service_factories.py`** — add `FitDiscrete` factory
2. **`phykit/services/tree/__init__.py`** — add to `_EXPORTS`
3. **`phykit/cli_registry.py`** — add aliases: `fitdiscrete`, `fd` → `"fit_discrete"`
4. **`phykit/phykit.py`** — add `fit_discrete` static method, help banner entry, module-level function
5. **`setup.py`** — add `pk_fit_discrete`, `pk_fitdiscrete`, `pk_fd` entry points

### R validation

**New file:** `tests/r_validation/validate_fit_discrete.R`

Uses `geiger::fitDiscrete()` on the same sample tree and traits. Fits ER, SYM, ARD models and prints lnL, AIC, and AIC weights for cross-validation.

```r
suppressPackageStartupMessages({
  library(ape)
  library(geiger)
})

tree <- read.tree("../sample_files/tree_simple.tre")
traits <- read.delim("../sample_files/tree_simple_discrete_traits.tsv")
trait_vec <- setNames(traits$diet, traits$taxon)

# Use equal root prior to match PhyKIT's fixed flat prior
cat("=== KEY VALUES FOR PYTHON TESTS ===\n")
for (model in c("ER", "SYM", "ARD")) {
  fit <- fitDiscrete(tree, trait_vec, model = model, type = "equal")
  n <- length(trait_vec)
  bic <- -2 * fit$opt$lnL + fit$opt$k * log(n)
  cat(sprintf("Model: %s  lnL: %.4f  AIC: %.4f  BIC: %.4f  k: %d\n",
              model, fit$opt$lnL, fit$opt$aic, bic, fit$opt$k))
}
```

**Tolerance:** Log-likelihood values should agree within 0.1 of R's values, since multi-start optimization in PhyKIT and geiger may converge to slightly different local optima. AIC/BIC are deterministic given lnL and n_params.

### Changelog entry

```
Added discrete trait model comparison command
(``fit_discrete`` / ``fd``):

* Compares ER (Equal Rates), SYM (Symmetric), and ARD (All Rates
  Different) Mk models of discrete character evolution
* Reports log-likelihood, AIC, ΔAIC, Akaike weights, and number
  of parameters for each model
* Extracts shared Q-matrix fitting code from stochastic_character_map
  and ancestral_reconstruction into phykit/helpers/discrete_models.py,
  eliminating code duplication
* Supports ``--json`` output with full Q-matrix and rate parameters
* Cross-validated against R's geiger::fitDiscrete(); R validation
  script provided in ``tests/r_validation/validate_fit_discrete.R``
```

### Testing

**Unit tests** (`tests/unit/helpers/test_discrete_models.py`):
- Test `build_q_matrix` for ER, SYM, ARD (correct shape, rows sum to 0)
- Test `felsenstein_pruning` returns valid log-likelihood
- Test `fit_q_matrix` recovers known parameters on simple trees

**Unit tests** (`tests/unit/services/tree/test_fit_discrete.py`):
- Test init and process_args
- Test `_fit_model` returns correct structure
- Test `_compute_model_comparison` produces correct ΔAIC and weights
- Test run with text and JSON output

**Integration tests** (`tests/integration/tree/test_fit_discrete_integration.py`):
- End-to-end CLI test with sample data (values validated against R)
- Test aliases
- Test JSON output
- Test `--models` flag for subset selection

### Files changed

| File | Change |
|------|--------|
| `phykit/helpers/discrete_models.py` | **New** — shared Q-matrix, Felsenstein pruning, trait parsing |
| `phykit/services/tree/fit_discrete.py` | **New** — FitDiscrete service |
| `phykit/services/tree/stochastic_character_map.py` | Replace duplicated methods with imports |
| `phykit/services/tree/ancestral_reconstruction.py` | Replace duplicated methods with imports |
| `phykit/services/tree/density_map.py` | Update to import trait parser from shared module |
| `phykit/service_factories.py` | Add factory |
| `phykit/services/tree/__init__.py` | Add to exports (only if other newer services like FitContinuous are also listed; otherwise skip for consistency) |
| `phykit/cli_registry.py` | Add aliases |
| `phykit/phykit.py` | Add CLI method, help text, help banner, module-level function |
| `setup.py` | Add entry points |
| `docs/usage/index.rst` | Add command documentation |
| `docs/change_log/index.rst` | Add changelog entry |
| `tests/unit/helpers/test_discrete_models.py` | **New** — unit tests for shared module |
| `tests/unit/services/tree/test_fit_discrete.py` | **New** — unit tests for service |
| `tests/integration/tree/test_fit_discrete_integration.py` | **New** — integration tests |
| `tests/r_validation/validate_fit_discrete.R` | **New** — R validation script |
