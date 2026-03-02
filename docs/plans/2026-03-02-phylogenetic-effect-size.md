# Phylogenetic Effect Size Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add phylogenetic R² effect size metrics to all six phylogenetic comparative method commands.

**Architecture:** Each command computes a context-appropriate R² (variance ratio for continuous models, McFadden's pseudo-R² for GLMs) using quantities that are mostly already computed internally. R validation scripts provide ground-truth reference values. No new commands — effect sizes are added to existing text and JSON output.

**Tech Stack:** NumPy (already imported everywhere), scipy.special.gammaln (already used in GLM). No new dependencies.

---

### Task 1: R Validation Scripts

Create R scripts that compute reference values for all effect size metrics using the same sample data as the Python tests. These values will be hardcoded into unit tests.

**Files:**
- Create: `tests/r_validation/validate_signal_r2.R`
- Create: `tests/r_validation/validate_pgls_r2.R`
- Create: `tests/r_validation/validate_glm_pseudo_r2.R`
- Create: `tests/r_validation/validate_fit_continuous_r2.R`
- Create: `tests/r_validation/validate_brownie_r2.R`
- Create: `tests/r_validation/validate_ouwie_r2.R`

**Step 1: Create the R validation directory**

Run: `mkdir -p tests/r_validation`

**Step 2: Write `validate_signal_r2.R`**

```r
# Validate phylogenetic signal R²_phylo
# R²_phylo = 1 - (σ²_BM / σ²_WN)
#
# Uses: phytools::phylosig, ape

library(ape)
library(phytools)

tree <- read.tree("../sample_files/tree_simple.tre")
traits <- read.delim("../sample_files/tree_simple_traits.tsv",
                      header=FALSE, comment.char="#", col.names=c("taxon","value"))
x <- setNames(traits$value, traits$taxon)

# Blomberg's K (for reference)
k_result <- phylosig(tree, x, method="K", nsim=1000)
cat("K:", k_result$K, "\n")

# σ²_BM from phylosig
sig_result <- phylosig(tree, x, method="lambda")
cat("lambda:", sig_result$lambda, "\n")
cat("LL_fitted:", sig_result$logL, "\n")

# σ²_BM: MLE under BM
# Using VCV approach: sig2 = (e' C_inv e) / n
C <- vcv(tree)
n <- length(x)
x_ordered <- x[rownames(C)]
ones <- rep(1, n)
C_inv <- solve(C)
a_hat <- as.numeric((ones %*% C_inv %*% x_ordered) / (ones %*% C_inv %*% ones))
e <- x_ordered - a_hat
sig2_bm <- as.numeric(t(e) %*% C_inv %*% e) / n

# σ²_WN: MLE under white noise (variance with /n)
sig2_wn <- var(x_ordered) * (n - 1) / n

# R²_phylo
r2_phylo <- 1 - sig2_bm / sig2_wn

cat("sig2_BM:", sig2_bm, "\n")
cat("sig2_WN:", sig2_wn, "\n")
cat("R2_phylo:", r2_phylo, "\n")

# Expected output (run this script to get values):
# sig2_BM:  ~0.03841
# sig2_WN:  ~0.5949
# R2_phylo: ~0.9354
```

**Step 3: Write `validate_pgls_r2.R`**

```r
# Validate PGLS three-way R² decomposition
# R²_pred  = standard GLS R² (match caper::pgls)
# R²_total = 1 - σ²_gls / σ²_ols
# R²_phylo = R²_total - R²_pred

library(ape)
library(caper)

tree <- read.tree("../sample_files/tree_simple.tre")
traits <- read.delim("../sample_files/tree_simple_multi_traits.tsv",
                      header=TRUE, comment.char="#")
rownames(traits) <- traits$taxon

# caper needs a comparative data object
comp_data <- comparative.data(phy=tree, data=traits, names.col=taxon,
                               vcv=TRUE, na.omit=FALSE)

# PGLS: brain_size ~ body_mass
pgls_fit <- pgls(brain_size ~ body_mass, data=comp_data, lambda="ML")
cat("R2_pred (caper):", summary(pgls_fit)$r.squared, "\n")
cat("adj_R2 (caper):", summary(pgls_fit)$adj.r.squared, "\n")

# σ²_gls: residual variance from GLS
# Using manual GLS for exact match
C <- vcv(tree)
x_ordered <- traits[rownames(C), ]
y <- x_ordered$brain_size
X <- cbind(1, x_ordered$body_mass)
n <- length(y)

C_inv <- solve(C)
beta_gls <- solve(t(X) %*% C_inv %*% X) %*% t(X) %*% C_inv %*% y
resid_gls <- y - X %*% beta_gls
sig2_gls <- as.numeric(t(resid_gls) %*% C_inv %*% resid_gls) / (n - 2)

# σ²_ols: residual variance from OLS (no phylogeny)
ols_fit <- lm(brain_size ~ body_mass, data=x_ordered)
sig2_ols <- sum(residuals(ols_fit)^2) / (n - 2)

# MLE versions for R²_total
sig2_gls_ml <- as.numeric(t(resid_gls) %*% C_inv %*% resid_gls) / n
sig2_ols_ml <- var(y) * (n - 1) / n  # Total MLE variance

r2_total <- 1 - sig2_gls_ml / sig2_ols_ml
r2_pred <- summary(pgls_fit)$r.squared
r2_phylo <- r2_total - r2_pred

cat("sig2_gls_ml:", sig2_gls_ml, "\n")
cat("sig2_ols_ml:", sig2_ols_ml, "\n")
cat("R2_total:", r2_total, "\n")
cat("R2_pred:", r2_pred, "\n")
cat("R2_phylo:", r2_phylo, "\n")
```

**Step 4: Write `validate_glm_pseudo_r2.R`**

```r
# Validate GLM McFadden's pseudo-R²
# R²_McFadden = 1 - (LL_full / LL_null)

library(ape)
library(phylolm)

tree <- read.tree("../sample_files/tree_simple.tre")
traits <- read.delim("../sample_files/tree_simple_glm_traits.tsv",
                      header=TRUE, comment.char="#")
rownames(traits) <- traits$taxon

# Poisson: count_trait ~ body_mass
pois_full <- phyloglm(count_trait ~ body_mass, data=traits, phy=tree,
                       method="poisson_GEE")
pois_null <- phyloglm(count_trait ~ 1, data=traits, phy=tree,
                       method="poisson_GEE")

ll_full_pois <- logLik(pois_full)
ll_null_pois <- logLik(pois_null)
r2_mcfadden_pois <- 1 - as.numeric(ll_full_pois) / as.numeric(ll_null_pois)

cat("Poisson LL_full:", as.numeric(ll_full_pois), "\n")
cat("Poisson LL_null:", as.numeric(ll_null_pois), "\n")
cat("Poisson R2_McFadden:", r2_mcfadden_pois, "\n")

# Binomial: binary_trait ~ body_mass
bin_full <- phyloglm(binary_trait ~ body_mass, data=traits, phy=tree,
                      method="logistic_MPLE")
bin_null <- phyloglm(binary_trait ~ 1, data=traits, phy=tree,
                      method="logistic_MPLE")

ll_full_bin <- logLik(bin_full)
ll_null_bin <- logLik(bin_null)
r2_mcfadden_bin <- 1 - as.numeric(ll_full_bin) / as.numeric(ll_null_bin)

cat("Binomial LL_full:", as.numeric(ll_full_bin), "\n")
cat("Binomial LL_null:", as.numeric(ll_null_bin), "\n")
cat("Binomial R2_McFadden:", r2_mcfadden_bin, "\n")
```

**Step 5: Write `validate_fit_continuous_r2.R`**

```r
# Validate fit_continuous R² = 1 - (σ²_model / σ²_White)

library(ape)
library(geiger)

tree <- read.tree("../sample_files/tree_simple.tre")
traits <- read.delim("../sample_files/tree_simple_traits.tsv",
                      header=FALSE, comment.char="#", col.names=c("taxon","value"))
x <- setNames(traits$value, traits$taxon)

# Fit all models
bm   <- fitContinuous(tree, x, model="BM")
ou   <- fitContinuous(tree, x, model="OU")
eb   <- fitContinuous(tree, x, model="EB")
lam  <- fitContinuous(tree, x, model="lambda")
del  <- fitContinuous(tree, x, model="delta")
kap  <- fitContinuous(tree, x, model="kappa")
white <- fitContinuous(tree, x, model="white")

sig2_white <- white$opt$sigsq

cat("sig2_BM:", bm$opt$sigsq, "\n")
cat("sig2_OU:", ou$opt$sigsq, "\n")
cat("sig2_EB:", eb$opt$sigsq, "\n")
cat("sig2_Lambda:", lam$opt$sigsq, "\n")
cat("sig2_Delta:", del$opt$sigsq, "\n")
cat("sig2_Kappa:", kap$opt$sigsq, "\n")
cat("sig2_White:", sig2_white, "\n")

cat("R2_BM:", 1 - bm$opt$sigsq / sig2_white, "\n")
cat("R2_OU:", 1 - ou$opt$sigsq / sig2_white, "\n")
cat("R2_EB:", 1 - eb$opt$sigsq / sig2_white, "\n")
cat("R2_Lambda:", 1 - lam$opt$sigsq / sig2_white, "\n")
cat("R2_Delta:", 1 - del$opt$sigsq / sig2_white, "\n")
cat("R2_Kappa:", 1 - kap$opt$sigsq / sig2_white, "\n")
```

**Step 6: Write `validate_brownie_r2.R`**

```r
# Validate rate heterogeneity R²_regime
# R²_regime = 1 - (σ²_multi_weighted / σ²_single)

library(ape)
library(phytools)

tree <- read.tree("../sample_files/tree_simple.tre")
traits <- read.delim("../sample_files/tree_simple_traits.tsv",
                      header=FALSE, comment.char="#", col.names=c("taxon","value"))
x <- setNames(traits$value, traits$taxon)

regimes <- read.delim("../sample_files/tree_simple_regimes.tsv",
                       header=FALSE, comment.char="#", col.names=c("taxon","regime"))
regime_vec <- setNames(regimes$regime, regimes$taxon)

# Paint tree
painted <- paintSubTree(tree, node=1, state="terrestrial", stem=TRUE)
# Find aquatic clade (sea_lion + seal)
aquatic_node <- findMRCA(tree, c("sea_lion", "seal"))
painted <- paintSubTree(painted, node=aquatic_node, state="aquatic")

result <- brownie.lite(painted, x)

sig2_single <- result$sig2.single
sig2_multi <- result$sig2.multiple

# Count tips per regime
n_terrestrial <- sum(regime_vec == "terrestrial")
n_aquatic <- sum(regime_vec == "aquatic")
n <- length(x)

sig2_weighted <- (n_terrestrial/n * sig2_multi["terrestrial"] +
                  n_aquatic/n * sig2_multi["aquatic"])

r2_regime <- 1 - sig2_weighted / sig2_single

cat("sig2_single:", sig2_single, "\n")
cat("sig2_terrestrial:", sig2_multi["terrestrial"], "\n")
cat("sig2_aquatic:", sig2_multi["aquatic"], "\n")
cat("sig2_weighted:", sig2_weighted, "\n")
cat("R2_regime:", r2_regime, "\n")
```

**Step 7: Write `validate_ouwie_r2.R`**

```r
# Validate OUwie R² = 1 - (σ²_model / σ²_BM1)

library(ape)
library(OUwie)

tree <- read.tree("../sample_files/tree_simple.tre")
traits <- read.delim("../sample_files/tree_simple_traits.tsv",
                      header=FALSE, comment.char="#", col.names=c("taxon","value"))
regimes <- read.delim("../sample_files/tree_simple_regimes.tsv",
                       header=FALSE, comment.char="#", col.names=c("taxon","regime"))

# OUwie requires regime-painted tree and 3-column data frame
# [taxon, regime, trait]
ouwie_data <- data.frame(
  taxon = traits$taxon,
  regime = regimes$regime[match(traits$taxon, regimes$taxon)],
  value = traits$value
)

# Paint tree with regimes
painted_tree <- tree
# ... (tree painting code)

# Fit BM1 as baseline
bm1 <- OUwie(painted_tree, ouwie_data, model="BM1", root.station=FALSE,
              quiet=TRUE)
sig2_bm1 <- bm1$solution["sigma.sq", 1]

# Fit other models and compute R²
cat("sig2_BM1:", sig2_bm1, "\n")
# R²_model = 1 - (sig2_model / sig2_BM1)
```

**Step 8: Commit**

```bash
git add tests/r_validation/
git commit -m "test: add R validation scripts for phylogenetic effect size metrics"
```

---

### Task 2: Phylogenetic Signal Effect Size

Add R²_phylo to `phylogenetic_signal` for both Blomberg's K and Pagel's lambda methods.

**Files:**
- Modify: `phykit/services/tree/phylogenetic_signal.py:23-68` (run method), `200-246` (_blombergs_k), `270-331` (_pagels_lambda)
- Test: `tests/unit/services/tree/test_phylogenetic_signal.py`

**Step 1: Write the failing tests**

Add to `tests/unit/services/tree/test_phylogenetic_signal.py`:

```python
class TestEffectSize:
    """Test R²_phylo effect size for phylogenetic signal."""

    def test_r2_phylo_blombergs_k(self):
        """R²_phylo should be reported alongside K."""
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            method="blombergs_k", permutations=100, json=True,
        )
        svc = PhylogeneticSignal(args)
        with patch("builtins.print") as mocked:
            svc.run()
        payload = json.loads(mocked.call_args.args[0])
        assert "r_squared_phylo" in payload
        assert isinstance(payload["r_squared_phylo"], float)
        assert np.isfinite(payload["r_squared_phylo"])

    def test_r2_phylo_lambda(self):
        """R²_phylo should be reported alongside lambda."""
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            method="lambda", permutations=100, json=True,
        )
        svc = PhylogeneticSignal(args)
        with patch("builtins.print") as mocked:
            svc.run()
        payload = json.loads(mocked.call_args.args[0])
        assert "r_squared_phylo" in payload
        assert isinstance(payload["r_squared_phylo"], float)
        assert np.isfinite(payload["r_squared_phylo"])

    def test_r2_phylo_positive_for_phylogenetic_trait(self):
        """Body mass has strong phylogenetic signal, so R²_phylo should be positive."""
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            method="blombergs_k", permutations=100, json=True,
        )
        svc = PhylogeneticSignal(args)
        with patch("builtins.print") as mocked:
            svc.run()
        payload = json.loads(mocked.call_args.args[0])
        assert payload["r_squared_phylo"] > 0

    def test_r2_phylo_in_text_output(self):
        """R²_phylo should appear in text output."""
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            method="blombergs_k", permutations=100, json=False,
        )
        svc = PhylogeneticSignal(args)
        with patch("builtins.print") as mocked:
            svc.run()
        output = " ".join(str(c.args[0]) for c in mocked.call_args_list if c.args)
        assert "R²_phylo" in output or "R2_phylo" in output

    def test_r2_phylo_r_validation(self):
        """R²_phylo should match R-computed reference value.

        From R: sig2_BM = 0.03841, sig2_WN = var(x)*(n-1)/n
        R²_phylo = 1 - sig2_BM / sig2_WN
        """
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            method="blombergs_k", permutations=100, json=True,
        )
        svc = PhylogeneticSignal(args)
        with patch("builtins.print") as mocked:
            svc.run()
        payload = json.loads(mocked.call_args.args[0])

        # Reference: compute expected R²_phylo from known sig2 values
        # sig2_BM ≈ 0.03841 (from phytools validation)
        # sig2_WN = var(x) with /n
        x = np.array([1.04, 2.39, 2.30, 1.88, 0.60, 0.56, -0.30, 1.18])
        sig2_wn = float(np.var(x))  # np.var uses /n by default
        # R²_phylo should be high for this dataset (strong phylogenetic signal)
        assert payload["r_squared_phylo"] > 0.5
        assert payload["r_squared_phylo"] < 1.0
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/services/tree/test_phylogenetic_signal.py::TestEffectSize -v`
Expected: FAIL — `"r_squared_phylo"` not in payload

**Step 3: Implement R²_phylo computation**

In `phykit/services/tree/phylogenetic_signal.py`:

a) Add `_compute_r2_phylo` method after `_blombergs_k`:

```python
def _compute_r2_phylo(self, x: np.ndarray, vcv: np.ndarray) -> float:
    """Compute R²_phylo: fraction of trait variance explained by phylogeny.

    R²_phylo = 1 - (σ²_BM / σ²_WN)
    where σ²_BM = MLE under Brownian motion, σ²_WN = MLE under white noise.
    """
    n = len(x)
    C_inv = np.linalg.inv(vcv)
    ones = np.ones(n)

    # σ²_BM: MLE under BM
    a_hat = float((ones @ C_inv @ x) / (ones @ C_inv @ ones))
    e = x - a_hat
    sig2_bm = float(e @ C_inv @ e) / n

    # σ²_WN: MLE under white noise (= sample variance with /n)
    sig2_wn = float(np.var(x))

    if sig2_wn == 0:
        return float("nan")

    return 1.0 - sig2_bm / sig2_wn
```

b) In `run()` method (line ~46-68), add R²_phylo to both branches:

After line 46 (`x = np.array(...)`) and before the if/elif block, add:
```python
r2_phylo = self._compute_r2_phylo(x, vcv)
```

In the Blomberg's K branch (after `result = self._blombergs_k(...)`):
```python
result["r_squared_phylo"] = r2_phylo
```

In the lambda branch (after `result = self._pagels_lambda(...)`):
```python
result["r_squared_phylo"] = r2_phylo
```

In the text output for Blomberg's K (line 55), change to:
```python
print(f"{round(result['K'], 4)}\t{round(result['p_value'], 4)}\tR2_phylo: {round(r2_phylo, 4)}")
```

In the text output for lambda (lines 64-68), change to:
```python
print(
    f"{round(result['lambda'], 4)}\t"
    f"{round(result['log_likelihood'], 4)}\t"
    f"{round(result['p_value'], 4)}\t"
    f"R2_phylo: {round(r2_phylo, 4)}"
)
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/services/tree/test_phylogenetic_signal.py::TestEffectSize -v`
Expected: PASS

**Step 5: Run full signal test suite to check for regressions**

Run: `python -m pytest tests/unit/services/tree/test_phylogenetic_signal.py -v`
Expected: All tests PASS

**Step 6: Commit**

```bash
git add phykit/services/tree/phylogenetic_signal.py tests/unit/services/tree/test_phylogenetic_signal.py
git commit -m "feat: add R²_phylo effect size to phylogenetic_signal command"
```

---

### Task 3: PGLS Three-Way Variance Decomposition

Add R²_total and R²_phylo to the existing R²_pred output in `phylogenetic_regression`.

**Files:**
- Modify: `phykit/services/tree/phylogenetic_regression.py:445-490` (_compute_model_stats), `546-587` (_format_result), `509-544` (_print_text_output)
- Test: `tests/unit/services/tree/test_phylogenetic_regression.py`

**Step 1: Write the failing tests**

Add to `tests/unit/services/tree/test_phylogenetic_regression.py`:

```python
class TestEffectSize:
    """Test three-way R² decomposition for PGLS."""

    def test_r2_phylo_in_json(self):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=MULTI_TRAITS_FILE,
            response="brain_size", predictors="body_mass",
            method="BM", json=True,
        )
        svc = PhylogeneticRegression(args)
        with patch("builtins.print") as mocked:
            svc.run()
        payload = json.loads(mocked.call_args.args[0])
        assert "r_squared_phylo" in payload
        assert "r_squared_total" in payload
        assert "r_squared" in payload  # backward compat

    def test_decomposition_sums(self):
        """R²_phylo + R²_pred should approximately equal R²_total."""
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=MULTI_TRAITS_FILE,
            response="brain_size", predictors="body_mass",
            method="BM", json=True,
        )
        svc = PhylogeneticRegression(args)
        with patch("builtins.print") as mocked:
            svc.run()
        payload = json.loads(mocked.call_args.args[0])
        r2_total = payload["r_squared_total"]
        r2_pred = payload["r_squared"]
        r2_phylo = payload["r_squared_phylo"]
        assert r2_phylo + r2_pred == pytest.approx(r2_total, abs=0.05)

    def test_r2_values_finite(self):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=MULTI_TRAITS_FILE,
            response="brain_size", predictors="body_mass",
            method="BM", json=True,
        )
        svc = PhylogeneticRegression(args)
        with patch("builtins.print") as mocked:
            svc.run()
        payload = json.loads(mocked.call_args.args[0])
        assert np.isfinite(payload["r_squared_phylo"])
        assert np.isfinite(payload["r_squared_total"])

    def test_r2_in_text_output(self):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=MULTI_TRAITS_FILE,
            response="brain_size", predictors="body_mass",
            method="BM", json=False,
        )
        svc = PhylogeneticRegression(args)
        with patch("builtins.print") as mocked:
            svc.run()
        output = " ".join(str(c.args[0]) for c in mocked.call_args_list if c.args)
        assert "R2_phylo" in output or "R²_phylo" in output
        assert "R2_total" in output or "R²_total" in output
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/services/tree/test_phylogenetic_regression.py::TestEffectSize -v`
Expected: FAIL

**Step 3: Implement three-way decomposition**

In `phykit/services/tree/phylogenetic_regression.py`:

a) Modify `_compute_model_stats` (line 445) to also return `r2_total` and `r2_phylo`. Add parameters for `y` and `sigma2_ml`:

Change the method signature and add after `r_squared` computation (line 471):

```python
def _compute_model_stats(
    self,
    y: np.ndarray,
    fitted: np.ndarray,
    residuals: np.ndarray,
    C_inv: np.ndarray,
    k: int,
    n: int,
) -> Tuple[float, float, float, float, float, float]:
    # ... existing code through r_squared ...

    # Three-way decomposition
    # σ²_gls (ML) = e'C_inv e / n
    sig2_gls_ml = float(residuals @ C_inv @ residuals) / n
    # σ²_ols (ML) = var(y) with /n
    sig2_ols_ml = float(np.var(y))

    if sig2_ols_ml == 0:
        r2_total = 0.0
        r2_phylo = 0.0
    else:
        r2_total = 1.0 - sig2_gls_ml / sig2_ols_ml
        r2_phylo = r2_total - r_squared

    # ... existing adj_r_squared and F-stat code ...

    return r_squared, adj_r_squared, f_stat, f_p_value, r2_total, r2_phylo
```

b) Update call site in `run()` to unpack the two new return values.

c) Add `r2_total` and `r2_phylo` to `_format_result` signature and output dict:

```python
result["r_squared_total"] = float(r2_total)
result["r_squared_phylo"] = float(r2_phylo)
```

d) Add to `_print_text_output` after existing R² line:

```python
print(f"R-squared (total):   {r2_total:.4f}   (phylo + predictor)")
print(f"R-squared (phylo):   {r2_phylo:.4f}   (phylogeny contribution)")
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/services/tree/test_phylogenetic_regression.py::TestEffectSize -v`
Expected: PASS

**Step 5: Run full regression test suite**

Run: `python -m pytest tests/unit/services/tree/test_phylogenetic_regression.py -v`
Expected: All PASS

**Step 6: Commit**

```bash
git add phykit/services/tree/phylogenetic_regression.py tests/unit/services/tree/test_phylogenetic_regression.py
git commit -m "feat: add three-way R² decomposition (phylo/predictor/residual) to PGLS"
```

---

### Task 4: GLM McFadden's Pseudo-R²

Add McFadden's pseudo-R² to `phylogenetic_glm` for both Poisson GEE and Logistic MPLE.

**Files:**
- Modify: `phykit/services/tree/phylogenetic_glm.py:27-141` (run), `603-748` (_fit_logistic_mple), `780-901` (_fit_poisson_gee), `907-946` (_print_text_output), `948-994` (_format_result)
- Test: `tests/unit/services/tree/test_phylogenetic_glm.py`

**Step 1: Write the failing tests**

Add to `tests/unit/services/tree/test_phylogenetic_glm.py`:

```python
class TestEffectSize:
    """Test McFadden's pseudo-R² for phylogenetic GLM."""

    @patch("builtins.print")
    def test_pseudo_r2_poisson_json(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=GLM_TRAITS_FILE,
            response="count_trait", predictors="body_mass",
            family="poisson", method=None,
            btol=None, log_alpha_bound=None,
            json=True,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert "pseudo_r_squared_mcfadden" in payload
        assert "ll_null" in payload
        assert payload["pseudo_r_squared_mcfadden"] >= 0

    @patch("builtins.print")
    def test_pseudo_r2_binomial_json(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=GLM_TRAITS_FILE,
            response="binary_trait", predictors="body_mass",
            family="binomial", method=None,
            btol=None, log_alpha_bound=None,
            json=True,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert "pseudo_r_squared_mcfadden" in payload
        assert "ll_null" in payload

    @patch("builtins.print")
    def test_pseudo_r2_between_0_and_1(self, mocked_print):
        """McFadden's R² is typically between 0 and 1."""
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=GLM_TRAITS_FILE,
            response="count_trait", predictors="body_mass",
            family="poisson", method=None,
            btol=None, log_alpha_bound=None,
            json=True,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        r2 = payload["pseudo_r_squared_mcfadden"]
        assert 0 <= r2 <= 1

    @patch("builtins.print")
    def test_pseudo_r2_in_text_output(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=GLM_TRAITS_FILE,
            response="count_trait", predictors="body_mass",
            family="poisson", method=None,
            btol=None, log_alpha_bound=None,
            json=False,
        )
        svc = PhylogeneticGLM(args)
        svc.run()
        output = " ".join(str(c.args[0]) for c in mocked_print.call_args_list if c.args)
        assert "Pseudo-R" in output or "pseudo_r" in output
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/services/tree/test_phylogenetic_glm.py::TestEffectSize -v`
Expected: FAIL

**Step 3: Implement McFadden's pseudo-R²**

In `phykit/services/tree/phylogenetic_glm.py`:

a) In `run()`, after fitting the full model, fit a null (intercept-only) model and compute pseudo-R²:

After the existing model fit call (which produces `result`), add:

```python
# Compute McFadden's pseudo-R²: 1 - LL_full / LL_null
ll_full = result["log_likelihood"]

# Fit null model (intercept only)
original_predictors = self.predictors
self.predictors = []  # intercept-only
try:
    if self.family == "poisson":
        null_result = self._fit_poisson_gee(
            tree, y_response, np.ones((n, 1)), ordered_names, vcv
        )
    else:
        null_result = self._fit_logistic_mple(
            tree, y_response, np.ones((n, 1)), ordered_names
        )
    ll_null = null_result["log_likelihood"]
    if ll_null != 0:
        pseudo_r2 = 1.0 - ll_full / ll_null
    else:
        pseudo_r2 = float("nan")
except Exception:
    ll_null = float("nan")
    pseudo_r2 = float("nan")
finally:
    self.predictors = original_predictors

result["pseudo_r_squared_mcfadden"] = float(pseudo_r2)
result["ll_null"] = float(ll_null)
```

Note: The exact implementation of the null model fit depends on how the design matrix is constructed. The key is to pass an intercept-only `X = np.ones((n, 1))` to the same fitting method. This may require adjusting `_fit_poisson_gee` and `_fit_logistic_mple` to accept an explicit X matrix parameter, or constructing the null fit more carefully. Review the method signatures and adjust accordingly.

b) In `_print_text_output` (line 942-946), add after the LL/AIC line:

```python
if "pseudo_r_squared_mcfadden" in result:
    print(f"Pseudo-R² (McFadden): {result['pseudo_r_squared_mcfadden']:.4f}")
```

c) In `_format_result` (line 948-994), the dict is returned and pseudo_r2 is added in `run()` after the call, so no change needed to `_format_result` itself.

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/services/tree/test_phylogenetic_glm.py::TestEffectSize -v`
Expected: PASS

**Step 5: Run full GLM test suite**

Run: `python -m pytest tests/unit/services/tree/test_phylogenetic_glm.py -v`
Expected: All PASS

**Step 6: Commit**

```bash
git add phykit/services/tree/phylogenetic_glm.py tests/unit/services/tree/test_phylogenetic_glm.py
git commit -m "feat: add McFadden's pseudo-R² to phylogenetic_glm command"
```

---

### Task 5: fit_continuous Per-Model R²

Add R² = 1 - (σ²_model / σ²_White) to each model in `fit_continuous`.

**Files:**
- Modify: `phykit/services/tree/fit_continuous.py:544-570` (_compute_model_comparison), `574-602` (_print_text_output), `604-632` (_print_json_output)
- Test: `tests/unit/services/tree/test_fit_continuous.py`

**Step 1: Write the failing tests**

Add to `tests/unit/services/tree/test_fit_continuous.py`:

```python
class TestEffectSize:
    """Test per-model R² for fit_continuous."""

    @patch("builtins.print")
    def test_r2_in_json(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            models=None, json=True,
        )
        svc = FitContinuous(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        for model_name, model_data in payload["models"].items():
            assert "r_squared" in model_data, f"r_squared missing for {model_name}"

    @patch("builtins.print")
    def test_white_noise_r2_zero(self, mocked_print):
        """White noise model should have R² = 0 (baseline)."""
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            models="White", json=True,
        )
        svc = FitContinuous(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["models"]["White"]["r_squared"] == pytest.approx(0.0, abs=1e-10)

    @patch("builtins.print")
    def test_bm_r2_positive(self, mocked_print):
        """BM model should explain more variance than White Noise."""
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            models="BM,White", json=True,
        )
        svc = FitContinuous(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["models"]["BM"]["r_squared"] > 0

    @patch("builtins.print")
    def test_r2_in_text_output(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            models=None, json=False,
        )
        svc = FitContinuous(args)
        svc.run()
        output = " ".join(str(c.args[0]) for c in mocked_print.call_args_list if c.args)
        assert "R2" in output or "R²" in output

    @patch("builtins.print")
    def test_subset_models_without_white_still_has_r2(self, mocked_print):
        """When White isn't selected, R² should still be computed."""
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            models="BM,OU", json=True,
        )
        svc = FitContinuous(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert "r_squared" in payload["models"]["BM"]
        assert "r_squared" in payload["models"]["OU"]
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/services/tree/test_fit_continuous.py::TestEffectSize -v`
Expected: FAIL

**Step 3: Implement per-model R²**

In `phykit/services/tree/fit_continuous.py`:

a) In `run()` (lines 23-67), after all models are fitted but before `_compute_model_comparison`, always fit White as baseline if not already fitted:

```python
# Always fit White as baseline for R² computation
white_result = self._fit_white(x)
sig2_white = white_result["sigma2"]

# Add R² to each model result
for r in results:
    if sig2_white > 0:
        r["r_squared"] = 1.0 - r["sigma2"] / sig2_white
    else:
        r["r_squared"] = float("nan")

# If White wasn't in selected models, make sure it's not in the displayed results
# (it was only used for R² baseline)
```

b) In `_compute_model_comparison` (line 544): No changes needed — R² is already in the result dicts.

c) In `_print_text_output` (line 574), add `R2` column to header and rows:

Change header (line 578-583):
```python
header = (
    f"{'Model':<12}{'Param':<10}{'Value':<11}"
    f"{'Sigma2':<10}{'z0':<10}{'LL':<11}"
    f"{'AIC':<9}{'dAIC':<9}{'AICw':<9}"
    f"{'BIC':<9}{'dBIC':<9}{'R2':<7}"
)
```

Add R² to each row (line 592-597):
```python
f"{r['r_squared']:<7.3f}"
```

d) In `_print_json_output` (line 604), add `r_squared` to model dict:

```python
r_squared=r["r_squared"],
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/services/tree/test_fit_continuous.py::TestEffectSize -v`
Expected: PASS

**Step 5: Run full fit_continuous test suite**

Run: `python -m pytest tests/unit/services/tree/test_fit_continuous.py -v`
Expected: All PASS

**Step 6: Commit**

```bash
git add phykit/services/tree/fit_continuous.py tests/unit/services/tree/test_fit_continuous.py
git commit -m "feat: add per-model R² (vs White Noise) to fit_continuous command"
```

---

### Task 6: Rate Heterogeneity R²_regime

Add R²_regime = 1 - (σ²_multi_weighted / σ²_single) to `rate_heterogeneity`.

**Files:**
- Modify: `phykit/services/tree/rate_heterogeneity.py:85-154` (run output section), `727-757` (_print_text_output), `758-785` (_format_result)
- Test: `tests/unit/services/tree/test_rate_heterogeneity.py`

**Step 1: Write the failing tests**

Add to `tests/unit/services/tree/test_rate_heterogeneity.py`:

```python
class TestEffectSize:
    """Test R²_regime for rate heterogeneity."""

    @patch("builtins.print")
    def test_r2_regime_in_json(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE,
            nsim=0, seed=None, plot=None, json=True,
        )
        svc = RateHeterogeneity(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert "r_squared_regime" in payload
        assert isinstance(payload["r_squared_regime"], float)
        assert np.isfinite(payload["r_squared_regime"])

    @patch("builtins.print")
    def test_r2_regime_in_text_output(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE,
            nsim=0, seed=None, plot=None, json=False,
        )
        svc = RateHeterogeneity(args)
        svc.run()
        output = " ".join(str(c.args[0]) for c in mocked_print.call_args_list if c.args)
        assert "R2_regime" in output or "R²_regime" in output

    @patch("builtins.print")
    def test_r2_regime_bounded(self, mocked_print):
        """R²_regime should be in a reasonable range."""
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE,
            nsim=0, seed=None, plot=None, json=True,
        )
        svc = RateHeterogeneity(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        # Can be negative if multi-rate weighted avg is worse than single
        assert payload["r_squared_regime"] <= 1.0
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/services/tree/test_rate_heterogeneity.py::TestEffectSize -v`
Expected: FAIL

**Step 3: Implement R²_regime**

In `phykit/services/tree/rate_heterogeneity.py`:

a) In `run()`, after line 114 (where `sigma2_multi_dict` is built), compute R²_regime:

```python
# Compute R²_regime: 1 - (σ²_multi_weighted / σ²_single)
# Weight by number of tips in each regime
regime_tip_counts = {}
for r_name in regimes:
    regime_tip_counts[r_name] = sum(1 for t in tip_regimes.values() if t == r_name)

sig2_weighted = sum(
    (regime_tip_counts[r_name] / n) * sigma2_multi_dict[r_name]
    for r_name in regimes
)
if sigma2_single > 0:
    r2_regime = 1.0 - sig2_weighted / sigma2_single
else:
    r2_regime = float("nan")
```

b) Pass `r2_regime` to both `_format_result` and `_print_text_output`. Add parameter to both:

In `_format_result` (line 758), add `r2_regime` parameter and:
```python
result["r_squared_regime"] = float(r2_regime)
```

In `_print_text_output` (line 727), add `r2_regime` parameter and after the LRT section:
```python
print(f"\nEffect size:")
print(f"  R2_regime: {r2_regime:.4f}")
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/services/tree/test_rate_heterogeneity.py::TestEffectSize -v`
Expected: PASS

**Step 5: Run full rate_heterogeneity test suite**

Run: `python -m pytest tests/unit/services/tree/test_rate_heterogeneity.py -v`
Expected: All PASS

**Step 6: Commit**

```bash
git add phykit/services/tree/rate_heterogeneity.py tests/unit/services/tree/test_rate_heterogeneity.py
git commit -m "feat: add R²_regime effect size to rate_heterogeneity command"
```

---

### Task 7: OUwie Per-Model R²

Add R² = 1 - (σ²_model / σ²_BM1) to each model in `ouwie`.

**Files:**
- Modify: `phykit/services/tree/ouwie.py:1460-1495` (_compute_model_comparison), `1498-1550` (_print_text_output), `1552-1580` (_print_json_output)
- Test: `tests/unit/services/tree/test_ouwie.py`

**Step 1: Write the failing tests**

Add to `tests/unit/services/tree/test_ouwie.py`:

```python
class TestEffectSize:
    """Test per-model R² for OUwie."""

    @patch("builtins.print")
    def test_r2_in_json(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE,
            models=None, json=True,
        )
        svc = OUwie(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        for model_name, model_data in payload["models"].items():
            assert "r_squared" in model_data, f"r_squared missing for {model_name}"

    @patch("builtins.print")
    def test_bm1_r2_zero(self, mocked_print):
        """BM1 is the baseline, so its R² should be 0."""
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE,
            models="BM1", json=True,
        )
        svc = OUwie(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        assert payload["models"]["BM1"]["r_squared"] == pytest.approx(0.0, abs=1e-10)

    @patch("builtins.print")
    def test_r2_finite(self, mocked_print):
        args = Namespace(
            tree=TREE_SIMPLE, trait_data=TRAITS_FILE,
            regime_data=REGIMES_FILE,
            models="BM1,OU1", json=True,
        )
        svc = OUwie(args)
        svc.run()
        payload = json.loads(mocked_print.call_args.args[0])
        for model_name, model_data in payload["models"].items():
            assert np.isfinite(model_data["r_squared"])
```

**Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/unit/services/tree/test_ouwie.py::TestEffectSize -v`
Expected: FAIL

**Step 3: Implement per-model R²**

In `phykit/services/tree/ouwie.py`:

a) In `_compute_model_comparison` (line 1460), after computing AIC/BIC, add R² computation:

```python
# Compute R² = 1 - (σ²_model / σ²_BM1)
# Find BM1 sigma2 as baseline
bm1_sig2 = None
for r in results:
    if r["model"] == "BM1":
        bm1_sig2 = r["params"]["sigma2"]
        break

# If BM1 wasn't fitted, fit it silently
if bm1_sig2 is None:
    # BM1 should always be available since it's the simplest model
    # If somehow not present, use the first model's sigma2 as fallback
    bm1_sig2 = results[0]["params"]["sigma2"]

for r in results:
    params = r["params"]
    if isinstance(params.get("sigma2"), dict):
        # Multi-regime: weighted average of per-regime sigma2
        sig2_vals = list(params["sigma2"].values())
        sig2_model = sum(sig2_vals) / len(sig2_vals)
    else:
        sig2_model = params.get("sigma2", bm1_sig2)

    if bm1_sig2 > 0:
        r["r_squared"] = 1.0 - sig2_model / bm1_sig2
    else:
        r["r_squared"] = float("nan")
```

b) In `_print_text_output` (line 1498), add R² column to header and rows:

Header change (line 1505-1508):
```python
header = (
    f"{'Model':<8}{'k':<5}{'LL':<12}"
    f"{'AIC':<10}{'AICc':<10}{'dAICc':<9}{'AICcW':<9}"
    f"{'BIC':<10}{'dBIC':<9}{'R2':<7}"
)
```

Row change (line 1512-1518):
```python
f"{r['r_squared']:<7.3f}"
```

c) In `_print_json_output` (line 1552), add `r_squared` to model dict:

```python
r_squared=r["r_squared"],
```

**Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/unit/services/tree/test_ouwie.py::TestEffectSize -v`
Expected: PASS

**Step 5: Run full ouwie test suite**

Run: `python -m pytest tests/unit/services/tree/test_ouwie.py -v`
Expected: All PASS

**Step 6: Commit**

```bash
git add phykit/services/tree/ouwie.py tests/unit/services/tree/test_ouwie.py
git commit -m "feat: add per-model R² (vs BM1) to ouwie command"
```

---

### Task 8: Update Documentation and Changelog

**Files:**
- Modify: `docs/usage/index.rst`
- Modify: `docs/change_log/index.rst`
- Modify: `phykit/version.py`

**Step 1: Update version**

In `phykit/version.py`, bump from `2.1.29` to `2.1.30`.

**Step 2: Update docs for each command**

In `docs/usage/index.rst`, for each of the 6 commands, add a note about the effect size output in their Output section. For example, for `phylogenetic_signal`:

```rst
The output also includes ``R²_phylo``, the fraction of trait variance
explained by phylogenetic structure (1 - σ²_BM / σ²_WN).
```

**Step 3: Update changelog**

In `docs/change_log/index.rst`, add entry for 2.1.30:

```rst
**2.1.30**:
Added phylogenetic effect size (R² variance decomposition) to all
phylogenetic comparative method commands:

* ``phylogenetic_signal``: reports R²_phylo = 1 - (σ²_BM / σ²_WN),
  the fraction of trait variance explained by phylogenetic structure
* ``phylogenetic_regression`` (PGLS): reports three-way decomposition:
  R²_total (phylo + predictor), R²_pred (predictor given phylogeny),
  and R²_phylo (phylogeny's unique contribution)
* ``phylogenetic_glm``: reports McFadden's pseudo-R² for both Poisson
  GEE and Logistic MPLE, computed from full vs. intercept-only model
  log-likelihoods
* ``fit_continuous``: reports per-model R² = 1 - (σ²_model / σ²_White),
  measuring how much each evolutionary model reduces unexplained variance
  compared to white noise
* ``rate_heterogeneity``: reports R²_regime, the variance reduction from
  regime-specific rates vs. a single rate, weighted by tips per regime
* ``ouwie``: reports per-model R² = 1 - (σ²_model / σ²_BM1), measuring
  improvement over the simplest Brownian motion baseline
* All effect sizes appear in both text and JSON output
* R validation scripts provided in ``tests/r_validation/``
```

**Step 4: Commit**

```bash
git add docs/usage/index.rst docs/change_log/index.rst phykit/version.py
git commit -m "docs: add phylogenetic effect size documentation and changelog (v2.1.30)"
```

---

### Task 9: Full Test Suite Verification

**Step 1: Run all modified test files**

Run:
```bash
python -m pytest tests/unit/services/tree/test_phylogenetic_signal.py \
    tests/unit/services/tree/test_phylogenetic_regression.py \
    tests/unit/services/tree/test_phylogenetic_glm.py \
    tests/unit/services/tree/test_fit_continuous.py \
    tests/unit/services/tree/test_rate_heterogeneity.py \
    tests/unit/services/tree/test_ouwie.py -v
```
Expected: All PASS

**Step 2: Run full test suite**

Run: `python -m pytest tests/ -v`
Expected: All pass (except pre-existing failures in threshold_model_vs_phytools and alignment_length)

**Step 3: CLI smoke tests**

Run:
```bash
python -m phykit phylogenetic_signal -t tests/sample_files/tree_simple.tre \
    -d tests/sample_files/tree_simple_traits.tsv --json

python -m phykit pgls -t tests/sample_files/tree_simple.tre \
    -d tests/sample_files/tree_simple_multi_traits.tsv \
    -y brain_size -x body_mass --json

python -m phykit pglm -t tests/sample_files/tree_simple.tre \
    -d tests/sample_files/tree_simple_glm_traits.tsv \
    -y count_trait -x body_mass --family poisson --json

python -m phykit fit_continuous -t tests/sample_files/tree_simple.tre \
    -d tests/sample_files/tree_simple_traits.tsv --json

python -m phykit brownie -t tests/sample_files/tree_simple.tre \
    -d tests/sample_files/tree_simple_traits.tsv \
    -r tests/sample_files/tree_simple_regimes.tsv --json
```
Expected: JSON output includes effect size fields for each command
