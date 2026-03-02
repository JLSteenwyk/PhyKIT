# Design: Phylogenetic Effect Size (R² Variance Decomposition)

## Problem

Phylogenetic comparative methods report p-values but lack standardized effect size measures. Lambda is a model parameter, not a variance decomposition. A proper R² analog — decomposing total trait variance into phylogenetic, predictor, and residual components — fills a decades-old gap.

## Approach

**Approach A (GLS Sums-of-Squares)** for continuous-response commands, supplemented by **Approach B (Likelihood-Based pseudo-R²)** for GLMs where SS decomposition doesn't apply.

Effect sizes are added directly to existing command outputs (text and JSON). No new standalone command.

## Metrics Per Command

### 1. `phylogenetic_signal` — R²_phylo

Fraction of trait variance explained by phylogeny.

```
σ²_BM  = (e' C⁻¹ e) / n          # MLE under Brownian motion (already computed)
σ²_WN  = Σ(x - x̄)² / n           # MLE under white noise (no phylogeny)

R²_phylo = 1 - (σ²_BM / σ²_WN)
```

- R²_phylo → 1: phylogeny explains all variance
- R²_phylo → 0: phylogeny explains nothing
- R²_phylo < 0: possible (phylogeny actively misleads), valid

**Output additions:**
- Text: `R²_phylo: 0.xxx`
- JSON: `"r_squared_phylo": float`

**R validation:** `phytools::phylosig()` returns `$sig2` (= σ²_BM). σ²_WN = `var(x)*(n-1)/n`. Verify ratio.

### 2. `phylogenetic_regression` (PGLS) — Three-Way Decomposition

Already computed: SS_resid, SS_total (GLS-weighted), R²_pred.

New computation:
```
σ²_ols = var(y) * (n-1)/n                    # MLE under no phylogeny
σ²_gls = sigma2                               # from existing GLS fit

R²_total = 1 - σ²_gls / σ²_ols              # total explained (phylo + predictors)
R²_pred  = existing r_squared                 # predictors given phylogeny
R²_phylo = R²_total - R²_pred                # phylogeny's unique contribution
```

**Output additions:**
- Text: `R²_phylo`, `R²_total` alongside existing R²
- JSON: `"r_squared_phylo": float`, `"r_squared_total": float`
- Keep existing `"r_squared"` key (= R²_pred) for backward compatibility

**R validation:** `caper::pgls()` reports `$r.squared`. Existing R²_pred must match. Full decomposition validated by comparing `caper::pgls()` σ² vs `lm()` σ².

### 3. `phylogenetic_glm` — McFadden's Pseudo-R²

SS decomposition doesn't apply to binary/count data. Use likelihood-based R².

```
R²_McFadden = 1 - (LL_full / LL_null)
```

LL_null = intercept-only model with same phylogenetic correlation structure.

- Poisson GEE: re-run `_fit_poisson_gee` with intercept-only design matrix
- Logistic MPLE: re-run `_fit_logistic_mple` with intercept-only design matrix

**Output additions:**
- Text: `Pseudo-R² (McFadden): 0.xxx`
- JSON: `"pseudo_r_squared_mcfadden": float`, `"ll_null": float`

**R validation:** `phylolm::phyloglm()` reports log-likelihoods. Compute McFadden's R² manually from `logLik(full)` / `logLik(null)`.

### 4. `fit_continuous` — Variance Ratio per Model

White Noise model already fitted as one of the 7 models.

```
R²_model = 1 - (σ²_model / σ²_White)
```

If White isn't in the selected models subset, fit it silently as baseline.

**Output additions:**
- Text: R² column in model comparison table
- JSON: `"r_squared": float` per model entry

**R validation:** `geiger::fitContinuous()` reports `$opt$sigsq`. Get σ²_White from `fitContinuous(model="white")`. Verify ratios.

### 5. `rate_heterogeneity` — R²_regime

Both single-rate and multi-rate σ² already computed.

```
σ²_multi_weighted = Σ(n_r/n * σ²_r)          # weighted by tips per regime
R²_regime = 1 - (σ²_multi_weighted / σ²_single)
```

**Output additions:**
- Text: `R²_regime: 0.xxx`
- JSON: `"r_squared_regime": float`

**R validation:** `phytools::brownie.lite()` reports `$sig2.single` and `$sig2.multiple`. Compute ratio.

### 6. `ouwie` — R²_OU per Model

BM1 as baseline. If BM1 isn't in selected models, fit it silently.

```
R²_OU = 1 - (σ²_model / σ²_BM1)
```

For multi-regime models, σ²_model is the weighted average of per-regime σ² values.

**Output additions:**
- Text: R² column in model comparison table
- JSON: `"r_squared": float` per model entry

**R validation:** `OUwie::OUwie()` reports `$sigma.sq`. Compare vs BM1.

## Files to Modify

| File | Change |
|------|--------|
| `phykit/services/tree/phylogenetic_signal.py` | Add R²_phylo computation and output |
| `phykit/services/tree/phylogenetic_regression.py` | Add three-way decomposition and output |
| `phykit/services/tree/phylogenetic_glm.py` | Add null model fit and McFadden's pseudo-R² |
| `phykit/services/tree/fit_continuous.py` | Add per-model R² vs White Noise |
| `phykit/services/tree/rate_heterogeneity.py` | Add R²_regime computation |
| `phykit/services/tree/ouwie.py` | Add per-model R² vs BM1 |
| `tests/unit/services/tree/test_phylogenetic_signal.py` | Effect size tests |
| `tests/unit/services/tree/test_phylogenetic_regression.py` | Effect size tests |
| `tests/unit/services/tree/test_phylogenetic_glm.py` | Effect size tests |
| `tests/unit/services/tree/test_fit_continuous.py` | Effect size tests |
| `tests/unit/services/tree/test_rate_heterogeneity.py` | Effect size tests |
| `tests/unit/services/tree/test_ouwie.py` | Effect size tests |
| `tests/r_validation/validate_signal_r2.R` | R reference values |
| `tests/r_validation/validate_pgls_r2.R` | R reference values |
| `tests/r_validation/validate_glm_pseudo_r2.R` | R reference values |
| `tests/r_validation/validate_fit_continuous_r2.R` | R reference values |
| `tests/r_validation/validate_brownie_r2.R` | R reference values |
| `tests/r_validation/validate_ouwie_r2.R` | R reference values |
| `docs/usage/index.rst` | Document effect size outputs |
| `docs/change_log/index.rst` | Changelog entry |

## R Validation Strategy

Each command gets an R script in `tests/r_validation/` that:
1. Loads the same sample data used in Python tests
2. Fits the equivalent R model
3. Extracts the component values (σ², LL, etc.)
4. Computes the R² metric from components
5. Prints expected values as comments

Python unit tests include hardcoded expected values from these R scripts, with tolerances appropriate for each metric (typically abs=1e-4 for variance ratios).

## Edge Cases

- **R²_phylo < 0**: Valid for phylogenetic signal — means phylogeny is misleading. Report as-is.
- **σ²_WN = 0**: Constant trait — R²_phylo undefined. Report as NaN with warning.
- **White Noise not in selected models** (fit_continuous): Fit silently as baseline, don't include in model comparison table.
- **BM1 not in selected models** (ouwie): Fit silently as baseline.
- **Intercept-only GLM fails to converge**: Report pseudo-R² as NaN with warning.
- **Single-predictor PGLS**: R²_pred = R²_total is possible when phylogeny adds nothing.
- **Lambda = 0 PGLS**: R²_phylo should be ~0, R²_pred should match OLS R².
