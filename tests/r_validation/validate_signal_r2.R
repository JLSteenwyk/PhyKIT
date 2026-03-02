#!/usr/bin/env Rscript
# validate_signal_r2.R
#
# Compute R^2_phylo for phylogenetic signal via variance decomposition.
#
# Formula:
#   R^2_phylo = 1 - (sigma^2_BM / sigma^2_WN)
#
# where:
#   sigma^2_BM  = ML variance under Brownian Motion (phylogenetic covariance)
#   sigma^2_WN  = ML variance under White Noise (star tree / no phylogeny)
#
# The BM model assumes trait evolution follows a random walk along the
# phylogeny, so the covariance between species is proportional to shared
# evolutionary history. The WN model assumes all species are independent
# (equivalent to a star tree). When BM fits much better than WN,
# sigma^2_BM << sigma^2_WN and R^2_phylo approaches 1, indicating strong
# phylogenetic signal. When BM offers no improvement, R^2_phylo ~ 0.
#
# We also compute Pagel's lambda and Blomberg's K as reference values,
# and verify our manual GLS computation against geiger::fitContinuous().
#
# Requires: ape, geiger, phytools
#
# Usage:
#   Rscript tests/r_validation/validate_signal_r2.R

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
  library(phytools)
})

cat("=== Phylogenetic Signal R^2 Validation ===\n\n")

# ---------------------------------------------------------------
# Load data
# ---------------------------------------------------------------
tree <- read.tree("../sample_files/tree_simple.tre")

trait_data <- read.delim("../sample_files/tree_simple_traits.tsv",
  header = FALSE, comment.char = "#",
  col.names = c("taxon", "trait"), stringsAsFactors = FALSE)

# Named vector matching tree tip order
y <- setNames(trait_data$trait, trait_data$taxon)
y <- y[tree$tip.label]
n <- length(y)

cat("Tree tips:", paste(tree$tip.label, collapse = ", "), "\n")
cat("Trait values:", paste(sprintf("%.2f", y), collapse = ", "), "\n")
cat("n =", n, "\n\n")

# ---------------------------------------------------------------
# Method 1: Manual GLS computation
# ---------------------------------------------------------------
cat("--- Manual GLS computation ---\n\n")

# BM covariance matrix (proportional to shared branch lengths)
C_bm <- vcv(tree, corr = FALSE)
C_bm <- C_bm[names(y), names(y)]  # ensure matching order

# WN covariance = identity (all species independent)
C_wn <- diag(n)

# GLS ML estimate of sigma^2 under BM:
#   beta_hat = (1' C^{-1} 1)^{-1} (1' C^{-1} y)
#   sigma^2_ML = (1/n) (y - 1*beta_hat)' C^{-1} (y - 1*beta_hat)
ones <- rep(1, n)

C_bm_inv <- solve(C_bm)
beta_bm <- as.numeric(solve(t(ones) %*% C_bm_inv %*% ones) %*% (t(ones) %*% C_bm_inv %*% y))
resid_bm <- y - beta_bm
sigma2_bm <- as.numeric(t(resid_bm) %*% C_bm_inv %*% resid_bm) / n

C_wn_inv <- solve(C_wn)
beta_wn <- as.numeric(solve(t(ones) %*% C_wn_inv %*% ones) %*% (t(ones) %*% C_wn_inv %*% y))
resid_wn <- y - beta_wn
sigma2_wn <- as.numeric(t(resid_wn) %*% C_wn_inv %*% resid_wn) / n

# Compute log-likelihoods for reference
# logL = -n/2 * log(2*pi) - 1/2 * log(det(sigma^2 * C)) - n/2
logL_bm <- -n/2 * log(2 * pi) - 0.5 * determinant(sigma2_bm * C_bm, logarithm = TRUE)$modulus - n/2
logL_wn <- -n/2 * log(2 * pi) - 0.5 * determinant(sigma2_wn * C_wn, logarithm = TRUE)$modulus - n/2

r2_phylo <- 1 - (sigma2_bm / sigma2_wn)

cat(sprintf("BM: beta_hat = %.6f, sigma^2_ML = %.6f\n", beta_bm, sigma2_bm))
cat(sprintf("WN: beta_hat = %.6f, sigma^2_ML = %.6f\n", beta_wn, sigma2_wn))
cat(sprintf("BM logLik = %.6f\n", logL_bm))
cat(sprintf("WN logLik = %.6f\n", logL_wn))
cat(sprintf("\nR^2_phylo = 1 - (%.6f / %.6f) = %.6f\n\n", sigma2_bm, sigma2_wn, r2_phylo))

# ---------------------------------------------------------------
# Method 2: geiger::fitContinuous cross-check
# ---------------------------------------------------------------
cat("--- geiger::fitContinuous cross-check ---\n\n")

fit_bm <- fitContinuous(tree, y, model = "BM")
fit_wn <- fitContinuous(tree, y, model = "white")

cat(sprintf("geiger BM:  sigma^2 = %.6f, z0 = %.6f, logLik = %.6f\n",
            fit_bm$opt$sigsq, fit_bm$opt$z0, fit_bm$opt$lnL))
cat(sprintf("geiger WN:  sigma^2 = %.6f, z0 = %.6f, logLik = %.6f\n",
            fit_wn$opt$sigsq, fit_wn$opt$z0, fit_wn$opt$lnL))

r2_geiger <- 1 - (fit_bm$opt$sigsq / fit_wn$opt$sigsq)
cat(sprintf("\nR^2_phylo (geiger) = 1 - (%.6f / %.6f) = %.6f\n\n",
            fit_bm$opt$sigsq, fit_wn$opt$sigsq, r2_geiger))

# Note: geiger's sigma^2 is the rate parameter (sigsq), not the ML residual
# variance. The ratio may differ slightly from the manual GLS approach
# because geiger internally scales by the VCV. We need to be precise about
# what sigma^2 means in each context.
#
# For the R^2 we want:
#   sigma^2_BM = ML residual variance = (1/n)(y-mu)' C_BM^{-1} (y-mu)
#   sigma^2_WN = ML residual variance = (1/n)(y-mu)' I^{-1} (y-mu) = var_ML(y)
# The geiger sigsq is the BM rate parameter. The residual variance in the
# scaled-C formulation is: sigma^2_residual = sigsq (when C is not
# standardized). We use the manual computation as the canonical value.

# ---------------------------------------------------------------
# Reference: Pagel's lambda and Blomberg's K
# ---------------------------------------------------------------
cat("--- Reference signal measures ---\n\n")

lambda_fit <- fitContinuous(tree, y, model = "lambda")
cat(sprintf("Pagel's lambda = %.6f (logLik = %.6f)\n",
            lambda_fit$opt$lambda, lambda_fit$opt$lnL))

K_result <- phylosig(tree, y, method = "K", test = TRUE)
cat(sprintf("Blomberg's K = %.6f (p = %.4f)\n", K_result$K, K_result$P))

lambda_result <- phylosig(tree, y, method = "lambda", test = TRUE)
cat(sprintf("Pagel's lambda (phylosig) = %.6f (p = %.4f, logLik = %.6f)\n",
            lambda_result$lambda, lambda_result$P, lambda_result$logL))

# ---------------------------------------------------------------
# Summary of key values for Python unit tests
# ---------------------------------------------------------------
cat("\n=== KEY VALUES FOR PYTHON TESTS ===\n")
cat(sprintf("sigma2_bm_manual  = %.10f\n", sigma2_bm))
cat(sprintf("sigma2_wn_manual  = %.10f\n", sigma2_wn))
cat(sprintf("r2_phylo_manual   = %.10f\n", r2_phylo))
cat(sprintf("beta_bm           = %.10f\n", beta_bm))
cat(sprintf("beta_wn           = %.10f\n", beta_wn))
cat(sprintf("logL_bm_manual    = %.10f\n", logL_bm))
cat(sprintf("logL_wn_manual    = %.10f\n", logL_wn))
cat(sprintf("geiger_bm_sigsq   = %.10f\n", fit_bm$opt$sigsq))
cat(sprintf("geiger_wn_sigsq   = %.10f\n", fit_wn$opt$sigsq))
cat(sprintf("blomberg_K        = %.10f\n", K_result$K))
cat(sprintf("pagel_lambda      = %.10f\n", lambda_result$lambda))

cat("\nDone.\n")
