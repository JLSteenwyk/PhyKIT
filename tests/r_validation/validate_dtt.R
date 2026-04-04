#!/usr/bin/env Rscript
# validate_dtt.R
#
# Cross-validate PhyKIT dtt against geiger::dtt().
#
# Validates deterministic values:
#   - DTT curve (relative disparity at each branching time)
#   - Branching times
#   - MDI statistic
#
# Requires: ape, geiger
#
# Usage:
#   Rscript tests/r_validation/validate_dtt.R

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
})

cat("=== Disparity Through Time (DTT) Validation ===\n\n")

# ---------------------------------------------------------------
# Load data (same as PhyKIT test data)
# ---------------------------------------------------------------
tree <- read.tree("../sample_files/ultrametric_tree.tre")
traits <- read.delim("../sample_files/tree_simple_multi_traits.tsv",
  header = TRUE, stringsAsFactors = FALSE)
rownames(traits) <- traits$taxon
traits$taxon <- NULL

cat("Tree tips:", tree$tip.label, "\n")
cat("Trait columns:", names(traits), "\n\n")

# ---------------------------------------------------------------
# Single trait DTT (body_mass)
# ---------------------------------------------------------------
cat("--- DTT: body_mass (single trait, avg.sq) ---\n\n")

dat <- traits[, "body_mass"]
names(dat) <- rownames(traits)

result <- dtt(tree, dat, index = "avg.sq", nsim = 0, plot = FALSE,
              calculateMDIp = FALSE)

cat("Branching times (relative):\n")
cat(sprintf("  %s\n", paste(round(result$times, 6), collapse = ", ")))
cat("\n")

cat("DTT values (relative disparity):\n")
cat(sprintf("  %s\n", paste(round(result$dtt, 6), collapse = ", ")))
cat("\n")

cat("MDI:", result$MDI, "\n\n")

# Print in a format easy to paste into Python tests
cat("--- Python-friendly format ---\n")
cat("times = [", paste(sprintf("%.6f", result$times), collapse=", "), "]\n")
cat("dtt_values = [", paste(sprintf("%.6f", result$dtt), collapse=", "), "]\n")
cat("mdi =", sprintf("%.6f", result$MDI), "\n\n")

# ---------------------------------------------------------------
# Multi-trait DTT (all traits)
# ---------------------------------------------------------------
cat("--- DTT: all traits (multivariate, avg.sq) ---\n\n")

dat_multi <- as.matrix(traits)

result_multi <- dtt(tree, dat_multi, index = "avg.sq", nsim = 0,
                    plot = FALSE, calculateMDIp = FALSE)

cat("Branching times (relative):\n")
cat(sprintf("  %s\n", paste(round(result_multi$times, 6), collapse = ", ")))
cat("\n")

cat("DTT values:\n")
cat(sprintf("  %s\n", paste(round(result_multi$dtt, 6), collapse = ", ")))
cat("\n")

cat("MDI:", result_multi$MDI, "\n\n")

cat("--- Python-friendly format (multi) ---\n")
cat("times_multi = [", paste(sprintf("%.6f", result_multi$times), collapse=", "), "]\n")
cat("dtt_multi = [", paste(sprintf("%.6f", result_multi$dtt), collapse=", "), "]\n")
cat("mdi_multi =", sprintf("%.6f", result_multi$MDI), "\n\n")

# ---------------------------------------------------------------
# Disparity function (subclade disparities)
# ---------------------------------------------------------------
cat("--- Subclade disparities (avg.sq) ---\n\n")

disp <- disparity(phy = tree, data = dat)
cat("Subclade disparities:\n")
print(round(disp, 6))
cat("\n")

# Total disparity
total_disp <- disparity(data = dat)
cat("Total disparity (no tree):", round(total_disp, 6), "\n\n")

# ---------------------------------------------------------------
# With simulations (for MDI p-value range check)
# ---------------------------------------------------------------
cat("--- DTT with simulations (nsim=100) ---\n\n")

set.seed(42)
result_sim <- dtt(tree, dat, index = "avg.sq", nsim = 100,
                  plot = FALSE, calculateMDIp = TRUE)

cat("MDI:", result_sim$MDI, "\n")
cat("MDI p-value:", result_sim$MDIp, "\n")
cat("Note: p-value is simulation-dependent; PhyKIT should be in similar range\n\n")

cat("=== Validation complete ===\n")
