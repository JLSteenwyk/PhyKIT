#!/usr/bin/env Rscript
# validate_phylo_anova.R
#
# Cross-validate PhyKIT phylo_anova against manual Cholesky-transformed
# ANOVA/MANOVA, following the RRPP approach of Adams & Collyer (2018).
#
# We validate deterministic components (SS, MS, F, Pillai's trace) exactly
# and verify that permutation-based values (Z, p) are reasonable.
#
# Requires: ape
#
# Usage:
#   Rscript tests/r_validation/validate_phylo_anova.R

suppressPackageStartupMessages(library(ape))

cat("=== Phylogenetic ANOVA / MANOVA (RRPP) Validation ===\n\n")

# ---------------------------------------------------------------
# Load data
# ---------------------------------------------------------------
tree <- read.tree("../sample_files/tree_simple.tre")
C <- vcv(tree)
ord <- sort(tree$tip.label)
C <- C[ord, ord]

cat("Sorted taxa:", ord, "\n\n")

# Cholesky transform
L <- t(chol(C))
L_inv <- solve(L)

# ---------------------------------------------------------------
# 1. Univariate Phylogenetic ANOVA
# ---------------------------------------------------------------
cat("--- Univariate Phylogenetic ANOVA ---\n\n")

d <- read.delim("../sample_files/tree_simple_anova_single.tsv",
  header = TRUE, stringsAsFactors = FALSE)
rownames(d) <- d$taxon
d <- d[ord, ]

Y <- as.matrix(d$body_mass)
Y_star <- L_inv %*% Y

group <- d$diet
unique_groups <- sort(unique(group))
n <- nrow(Y)
g <- length(unique_groups)

# Design matrix: intercept + treatment dummies
X <- matrix(0, nrow = n, ncol = g)
X[, 1] <- 1  # intercept
for (i in 1:n) {
  idx <- which(unique_groups == group[i])
  if (idx > 1) X[i, idx] <- 1
}
X_star <- L_inv %*% X

# Reduced model (intercept only)
X0 <- matrix(1, nrow = n, ncol = 1)
X0_star <- L_inv %*% X0

# Fit models
beta0 <- solve(t(X0_star) %*% X0_star) %*% t(X0_star) %*% Y_star
hat0 <- X0_star %*% beta0
resid0 <- Y_star - hat0

beta1 <- solve(t(X_star) %*% X_star) %*% t(X_star) %*% Y_star
hat1 <- X_star %*% beta1
resid1 <- Y_star - hat1

SS_total <- sum(resid0^2)
SS_resid <- sum(resid1^2)
SS_model <- SS_total - SS_resid

df_model <- g - 1
df_resid <- n - g
MS_model <- SS_model / df_model
MS_resid <- SS_resid / df_resid
F_obs <- MS_model / MS_resid

cat("Deterministic values (should match PhyKIT exactly):\n")
cat(sprintf("  SS_model: %.6f\n", SS_model))
cat(sprintf("  SS_resid: %.6f\n", SS_resid))
cat(sprintf("  SS_total: %.6f\n", SS_total))
cat(sprintf("  MS_model: %.6f\n", MS_model))
cat(sprintf("  MS_resid: %.6f\n", MS_resid))
cat(sprintf("  F_stat:   %.6f\n", F_obs))
cat(sprintf("  Df_model: %d\n", df_model))
cat(sprintf("  Df_resid: %d\n", df_resid))
cat("\n")

# RRPP permutations
set.seed(42)
n_perm <- 999
F_perms <- numeric(n_perm)
for (p in 1:n_perm) {
  perm_idx <- sample(n)
  resid_perm <- resid0[perm_idx]
  Y_perm <- hat0 + resid_perm

  beta1_p <- solve(t(X_star) %*% X_star) %*% t(X_star) %*% Y_perm
  hat1_p <- X_star %*% beta1_p
  resid1_p <- Y_perm - hat1_p

  SS_resid_p <- sum(resid1_p^2)
  SS_model_p <- SS_total - SS_resid_p
  MS_model_p <- SS_model_p / df_model
  MS_resid_p <- SS_resid_p / df_resid
  F_perms[p] <- MS_model_p / MS_resid_p
}

p_value <- mean(F_perms >= F_obs)
z_score <- (F_obs - mean(F_perms)) / sd(F_perms)

cat("Permutation-based values (seed=42, 999 permutations):\n")
cat(sprintf("  Z_score: %.4f\n", z_score))
cat(sprintf("  p_value: %.4f\n", p_value))
cat("  Note: permutation values may differ slightly from PhyKIT due\n")
cat("  to RNG implementation differences between R and Python.\n")
cat("  PhyKIT values should be in a similar range.\n\n")

# ---------------------------------------------------------------
# 2. Multivariate Phylogenetic MANOVA
# ---------------------------------------------------------------
cat("--- Multivariate Phylogenetic MANOVA ---\n\n")

d2 <- read.delim("../sample_files/tree_simple_anova_traits.tsv",
  header = TRUE, stringsAsFactors = FALSE)
rownames(d2) <- d2$taxon
d2 <- d2[ord, ]

Y_multi <- as.matrix(d2[, c("body_mass", "brain_size")])
Y_multi_star <- L_inv %*% Y_multi

group2 <- d2$diet

# Design matrix (same groups)
X2 <- matrix(0, nrow = n, ncol = g)
X2[, 1] <- 1
for (i in 1:n) {
  idx <- which(unique_groups == group2[i])
  if (idx > 1) X2[i, idx] <- 1
}
X2_star <- L_inv %*% X2
X02_star <- L_inv %*% X0

# Fit models
beta0_m <- solve(t(X02_star) %*% X02_star) %*% t(X02_star) %*% Y_multi_star
hat0_m <- X02_star %*% beta0_m
resid0_m <- Y_multi_star - hat0_m

beta1_m <- solve(t(X2_star) %*% X2_star) %*% t(X2_star) %*% Y_multi_star
hat1_m <- X2_star %*% beta1_m
resid1_m <- Y_multi_star - hat1_m

# SSCP matrices
SS_total_m <- t(resid0_m) %*% resid0_m
SS_resid_m <- t(resid1_m) %*% resid1_m
SS_model_m <- SS_total_m - SS_resid_m

# Scalar SS (trace)
ss_model_trace <- sum(diag(SS_model_m))
ss_resid_trace <- sum(diag(SS_resid_m))
ss_total_trace <- sum(diag(SS_total_m))

# Pillai's trace
H_E_inv <- SS_model_m %*% solve(SS_resid_m)
eig <- eigen(H_E_inv)$values
pillai <- sum(Re(eig) / (1 + Re(eig)))

cat("Deterministic values (should match PhyKIT exactly):\n")
cat(sprintf("  SS_model (trace): %.6f\n", ss_model_trace))
cat(sprintf("  SS_resid (trace): %.6f\n", ss_resid_trace))
cat(sprintf("  SS_total (trace): %.6f\n", ss_total_trace))
cat(sprintf("  MS_model: %.6f\n", ss_model_trace / df_model))
cat(sprintf("  MS_resid: %.6f\n", ss_resid_trace / df_resid))
cat(sprintf("  Pillai's trace: %.6f\n", pillai))
cat(sprintf("  Df_model: %d\n", df_model))
cat(sprintf("  Df_resid: %d\n", df_resid))
cat("\n")

# RRPP permutations for MANOVA
set.seed(42)
pillai_perms <- numeric(n_perm)
for (p in 1:n_perm) {
  perm_idx <- sample(n)
  resid_perm_m <- resid0_m[perm_idx, ]
  Y_perm_m <- hat0_m + resid_perm_m

  beta1_mp <- solve(t(X2_star) %*% X2_star) %*% t(X2_star) %*% Y_perm_m
  hat1_mp <- X2_star %*% beta1_mp
  resid1_mp <- Y_perm_m - hat1_mp

  SS_resid_mp <- t(resid1_mp) %*% resid1_mp
  SS_model_mp <- SS_total_m - SS_resid_mp

  tryCatch({
    H_E_inv_p <- SS_model_mp %*% solve(SS_resid_mp)
    eig_p <- Re(eigen(H_E_inv_p)$values)
    pillai_perms[p] <- sum(eig_p / (1 + eig_p))
  }, error = function(e) {
    pillai_perms[p] <<- 0
  })
}

p_value_m <- mean(pillai_perms >= pillai)
z_score_m <- (pillai - mean(pillai_perms)) / sd(pillai_perms)

cat("Permutation-based values (seed=42, 999 permutations):\n")
cat(sprintf("  Z_score: %.4f\n", z_score_m))
cat(sprintf("  p_value: %.4f\n", p_value_m))
cat("\n")

cat("=== Validation complete ===\n")
