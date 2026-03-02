#!/usr/bin/env Rscript
# validate_brownie_r2.R
#
# Compute R^2_regime for rate heterogeneity across phylogenetic regimes
# using phytools::brownie.lite().
#
# Formula:
#   R^2_regime = 1 - (sigma^2_multi_weighted / sigma^2_single)
#
# where:
#   sigma^2_single = BM rate parameter from a single-rate model (BM1)
#   sigma^2_multi_weighted = weighted average of per-regime rates from
#                            multi-rate model (BMS), weighted by the
#                            proportion of total tree length in each regime
#
# Interpretation: What proportion of the total evolutionary rate is
# explained by regime-specific rates? When regimes have very different
# rates, the weighted multi-rate will be less than the single rate,
# and R^2_regime will be positive.
#
# phytools::brownie.lite() implements O'Meara et al. (2006) and
# Thomas et al. (2006) for testing rate heterogeneity.
#
# Requires: ape, phytools
#
# Usage:
#   Rscript tests/r_validation/validate_brownie_r2.R

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

cat("=== Brownie R^2_regime Validation ===\n\n")

# ---------------------------------------------------------------
# Load data
# ---------------------------------------------------------------
tree <- read.tree("../sample_files/tree_simple.tre")

trait_data <- read.delim("../sample_files/tree_simple_traits.tsv",
  header = FALSE, comment.char = "#",
  col.names = c("taxon", "trait"), stringsAsFactors = FALSE)

regime_data <- read.delim("../sample_files/tree_simple_regimes.tsv",
  header = FALSE, comment.char = "#",
  col.names = c("taxon", "regime"), stringsAsFactors = FALSE)

# Named trait vector
y <- setNames(trait_data$trait, trait_data$taxon)

# Named regime vector
regimes <- setNames(regime_data$regime, regime_data$taxon)

cat("Tree tips:", paste(tree$tip.label, collapse = ", "), "\n")
cat("Trait values:", paste(sprintf("%s=%.2f", names(y), y), collapse = ", "), "\n")
cat("Regimes:", paste(sprintf("%s=%s", names(regimes), regimes), collapse = ", "), "\n\n")

# ---------------------------------------------------------------
# Paint regimes onto tree using phytools
# ---------------------------------------------------------------
cat("--- Painting regimes onto tree ---\n\n")

# Use make.simmap or paintSubTree to assign regimes
# For brownie.lite, we need a simmap tree with mapped regimes
# The simplest approach: use make.simmap with the regime assignments

# Create regime factor for tips
tip_regimes <- regimes[tree$tip.label]

# Use make.simmap to reconstruct ancestral states and create mapped tree
# nsim=1 for a single stochastic map; set Q="empirical" for data-driven rates
set.seed(42)  # reproducibility
simmap_tree <- make.simmap(tree, tip_regimes, model = "ER", nsim = 1,
                           message = FALSE)

cat("Simmap tree mapped states:\n")
cat("  Total tree length:", sum(tree$edge.length), "\n")

# Summarize regime branch lengths
regime_lengths <- sapply(simmap_tree$maps, function(m) {
  tapply(m, names(m), sum)
})
# Aggregate across edges
all_regimes <- unique(unlist(lapply(simmap_tree$maps, names)))
total_per_regime <- sapply(all_regimes, function(r) {
  sum(sapply(simmap_tree$maps, function(m) {
    if (r %in% names(m)) m[r] else 0
  }))
})
total_tree_length <- sum(tree$edge.length)
cat("  Per-regime branch lengths:\n")
for (r in names(total_per_regime)) {
  cat(sprintf("    %s: %.4f (%.1f%%)\n", r,
              total_per_regime[r],
              100 * total_per_regime[r] / total_tree_length))
}
cat("\n")

# ---------------------------------------------------------------
# Run brownie.lite
# ---------------------------------------------------------------
cat("--- brownie.lite results ---\n\n")

brownie_result <- brownie.lite(simmap_tree, y)

cat("Single-rate model (BM1):\n")
cat(sprintf("  sigma^2_single = %.10f\n", brownie_result$sig2.single))
cat(sprintf("  a (ancestral state) = %.10f\n", brownie_result$a.single))
cat(sprintf("  logLik = %.10f\n", brownie_result$logL1))
cat(sprintf("  k = %d\n", brownie_result$k1))

cat("\nMulti-rate model (BMS):\n")
for (i in seq_along(brownie_result$sig2.multiple)) {
  regime_name <- names(brownie_result$sig2.multiple)[i]
  cat(sprintf("  sigma^2_%s = %.10f\n", regime_name,
              brownie_result$sig2.multiple[i]))
}
cat(sprintf("  a (ancestral state) = %.10f\n", brownie_result$a.multiple))
cat(sprintf("  logLik = %.10f\n", brownie_result$logL.multiple))
cat(sprintf("  k = %d\n", brownie_result$k2))

cat(sprintf("\nP-value (chi-sq test): %.6f\n", brownie_result$P.chisq))
cat(sprintf("Convergence: %s\n\n", brownie_result$convergence))

# ---------------------------------------------------------------
# Compute R^2_regime
# ---------------------------------------------------------------
cat("--- R^2_regime computation ---\n\n")

sigma2_single <- brownie_result$sig2.single
sigma2_multi <- brownie_result$sig2.multiple

# Compute weighted average of multi-regime rates
# Weights = proportion of total tree length in each regime
weights <- total_per_regime / total_tree_length
# Ensure ordering matches
weights_ordered <- weights[names(sigma2_multi)]

sigma2_multi_weighted <- sum(sigma2_multi * weights_ordered)

r2_regime <- 1 - (sigma2_multi_weighted / sigma2_single)

cat("Regime weights (proportion of tree length):\n")
for (r in names(weights_ordered)) {
  cat(sprintf("  %s: %.10f\n", r, weights_ordered[r]))
}
cat(sprintf("\nsigma^2_single         = %.10f\n", sigma2_single))
cat(sprintf("sigma^2_multi_weighted = %.10f\n", sigma2_multi_weighted))
cat(sprintf("  = sum(sigma^2_i * w_i) = ", ""))
terms <- sprintf("%.6f * %.6f", sigma2_multi, weights_ordered)
cat(paste(terms, collapse = " + "), "\n")
cat(sprintf("\nR^2_regime = 1 - (%.10f / %.10f) = %.10f\n\n",
            sigma2_multi_weighted, sigma2_single, r2_regime))

# ---------------------------------------------------------------
# Alternative: Likelihood-based comparison
# ---------------------------------------------------------------
cat("--- Likelihood comparison ---\n")
cat(sprintf("logLik_single    = %.10f\n", brownie_result$logL1))
cat(sprintf("logLik_multi     = %.10f\n", brownie_result$logL.multiple))
cat(sprintf("Delta_logLik     = %.10f\n",
            brownie_result$logL.multiple - brownie_result$logL1))
cat(sprintf("LRT chi-sq       = %.10f (df = %d)\n",
            2 * (brownie_result$logL.multiple - brownie_result$logL1),
            brownie_result$k2 - brownie_result$k1))

# ---------------------------------------------------------------
# Cross-check with OUwie BM1/BMS
# ---------------------------------------------------------------
cat("\n--- OUwie BM1/BMS cross-check ---\n\n")

tryCatch({
  suppressPackageStartupMessages(library(OUwie))

  # Prepare data for OUwie
  data_ouwie <- merge(regime_data, trait_data, by = "taxon")
  colnames(data_ouwie) <- c("species", "regime", "trait")

  # Root the tree and assign node labels
  tree_rooted <- root(tree, "dog", resolve.root = TRUE)
  tree_rooted$edge.length[tree_rooted$edge.length == 0] <- 1e-6

  # Assign regime labels to nodes via ML ACE
  tip_regime_factor <- as.factor(setNames(data_ouwie$regime, data_ouwie$species)[tree_rooted$tip.label])
  names(tip_regime_factor) <- tree_rooted$tip.label
  anc <- ace(tip_regime_factor, tree_rooted, type = "discrete", method = "ML")
  tree_rooted$node.label <- apply(anc$lik.anc, 1, function(x)
    colnames(anc$lik.anc)[which.max(x)])

  bm1 <- OUwie(tree_rooted, data_ouwie, model = "BM1", root.station = FALSE,
               algorithm = "invert", quiet = TRUE, check.identify = FALSE)
  bms <- OUwie(tree_rooted, data_ouwie, model = "BMS", root.station = FALSE,
               algorithm = "invert", quiet = TRUE, check.identify = FALSE)

  cat(sprintf("OUwie BM1: sigma^2 = %.10f, logLik = %.6f\n",
              bm1$solution["sigma.sq", 1], bm1$loglik))
  cat(sprintf("OUwie BMS:\n"))
  for (col in colnames(bms$solution)) {
    cat(sprintf("  sigma^2_%s = %.10f\n", col, bms$solution["sigma.sq", col]))
  }
  cat(sprintf("  logLik = %.6f\n", bms$loglik))
}, error = function(e) {
  cat(sprintf("OUwie cross-check skipped: %s\n", e$message))
})

# ---------------------------------------------------------------
# Summary
# ---------------------------------------------------------------
cat("\n=== KEY VALUES FOR PYTHON TESTS ===\n")
cat(sprintf("sigma2_single          = %.10f\n", sigma2_single))
for (r in names(sigma2_multi)) {
  cat(sprintf("sigma2_%s      = %.10f\n", r, sigma2_multi[r]))
}
cat(sprintf("sigma2_multi_weighted  = %.10f\n", sigma2_multi_weighted))
cat(sprintf("r2_regime              = %.10f\n", r2_regime))
cat(sprintf("logLik_single          = %.10f\n", brownie_result$logL1))
cat(sprintf("logLik_multi           = %.10f\n", brownie_result$logL.multiple))
cat(sprintf("ancestral_single       = %.10f\n", brownie_result$a.single))
cat(sprintf("ancestral_multi        = %.10f\n", brownie_result$a.multiple))
for (r in names(weights_ordered)) {
  cat(sprintf("weight_%s      = %.10f\n", r, weights_ordered[r]))
}

cat("\nDone.\n")
