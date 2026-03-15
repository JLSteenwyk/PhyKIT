#!/usr/bin/env Rscript
# validate_fit_discrete.R
#
# Compare ER, SYM, and ARD models of discrete trait evolution using
# geiger::fitDiscrete() and print expected values for PhyKIT cross-validation.
#
# geiger::fitDiscrete() fits Mk models (Lewis 2001) to discrete character
# data on a phylogeny. ER has 1 rate parameter, SYM has k(k-1)/2,
# ARD has k(k-1) where k is the number of states.
#
# Note: geiger requires bifurcating trees (multi2di), while PhyKIT handles
# multifurcations natively. This means log-likelihoods may differ by ~0.5-1.0.
# The model ranking should agree.
#
# Requires: ape, geiger
#
# Usage:
#   cd tests/r_validation
#   Rscript validate_fit_discrete.R

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
})

tree <- read.tree("../sample_files/tree_simple.tre")
tree <- multi2di(tree)  # Resolve multifurcations for geiger
traits <- read.delim("../sample_files/tree_simple_discrete_traits.tsv")
trait_vec <- setNames(traits$diet, traits$taxon)
n <- length(trait_vec)

cat("=== KEY VALUES FOR PYTHON TESTS ===\n")
cat(sprintf("Number of taxa: %d\n", n))
cat(sprintf("Number of states: %d\n", length(unique(trait_vec))))
cat(sprintf("States: %s\n", paste(sort(unique(trait_vec)), collapse=", ")))
cat("\n")

for (model in c("ER", "SYM", "ARD")) {
  # Use equal root prior to match PhyKIT's fixed flat prior (pi = 1/k)
  fit <- fitDiscrete(tree, trait_vec, model = model, type = "equal")
  bic <- -2 * fit$opt$lnL + fit$opt$k * log(n)
  aic <- fit$opt$aic
  delta_aic <- aic  # Will be adjusted below
  cat(sprintf("Model: %s  lnL: %.4f  AIC: %.4f  BIC: %.4f  n_params: %d\n",
              model, fit$opt$lnL, aic, bic, fit$opt$k))
}

# Compute delta-AIC
er_fit <- fitDiscrete(tree, trait_vec, model = "ER", type = "equal")
sym_fit <- fitDiscrete(tree, trait_vec, model = "SYM", type = "equal")
ard_fit <- fitDiscrete(tree, trait_vec, model = "ARD", type = "equal")
aics <- c(er_fit$opt$aic, sym_fit$opt$aic, ard_fit$opt$aic)
cat(sprintf("\nBest model by AIC: %s\n", c("ER", "SYM", "ARD")[which.min(aics)]))
cat(sprintf("Delta-AIC (ER): %.4f\n", aics[1] - min(aics)))
cat(sprintf("Delta-AIC (SYM): %.4f\n", aics[2] - min(aics)))
cat(sprintf("Delta-AIC (ARD): %.4f\n", aics[3] - min(aics)))
