#!/usr/bin/env Rscript
# validate_simmap_summary.R
#
# Cross-validate PhyKIT simmap_summary against phytools::make.simmap()
# and phytools::describe.simmap().
#
# Validates:
#   - Q matrix (fitted rates) against phytools::fitMk()
#   - Per-branch dwelling proportions against describe.simmap()$ace
#   - Expected transitions against describe.simmap()$count
#   - Node posteriors against describe.simmap()$ace
#
# Requires: ape, phytools
#
# Usage:
#   Rscript tests/r_validation/validate_simmap_summary.R

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

cat("=== SIMMAP Summary Validation ===\n\n")

# ---------------------------------------------------------------
# Load data
# ---------------------------------------------------------------
tree <- read.tree("../sample_files/tree_simple.tre")
traits_df <- read.delim("../sample_files/tree_simple_discrete_traits.tsv",
  header = TRUE, stringsAsFactors = FALSE)

# Build named vector
diet <- traits_df$diet
names(diet) <- traits_df$taxon

cat("Tree tips:", tree$tip.label, "\n")
cat("Trait states:", sort(unique(diet)), "\n\n")

# ---------------------------------------------------------------
# Fit Q matrix using fitMk (ER model)
# ---------------------------------------------------------------
cat("--- Q matrix (ER model, fitMk) ---\n\n")

fit <- fitMk(tree, diet, model = "ER")
cat("Fitted Q matrix:\n")
print(fit$Q)
cat("\nLog-likelihood:", logLik(fit), "\n\n")

# ---------------------------------------------------------------
# Run stochastic maps
# ---------------------------------------------------------------
cat("--- Running 100 stochastic maps ---\n\n")

set.seed(42)
simmaps <- make.simmap(tree, diet, model = "ER", nsim = 100,
                       message = FALSE)

# ---------------------------------------------------------------
# Describe simmap (per-branch summary)
# ---------------------------------------------------------------
cat("--- describe.simmap() output ---\n\n")

desc <- describe.simmap(simmaps, message = FALSE)

cat("Expected transitions (count):\n")
print(desc$count)
cat("\n")

cat("Node posterior probabilities (ace):\n")
print(desc$ace)
cat("\n")

cat("Per-branch dwelling time proportions (tips):\n")
print(desc$tips)
cat("\n")

# Total expected transitions
cat("Total expected transitions:", sum(desc$count[, "N"]), "\n\n")

# ---------------------------------------------------------------
# Summary of key values to compare with PhyKIT
# ---------------------------------------------------------------
cat("--- Key values for PhyKIT comparison ---\n\n")

# Q matrix rate (ER has one rate)
rate <- fit$Q[1, 2]
cat(sprintf("Q rate (off-diagonal): %.6f\n", rate))
cat(sprintf("Q diagonal: %.6f\n", fit$Q[1, 1]))
cat(sprintf("Log-likelihood: %.4f\n", logLik(fit)))
cat("\n")

# Total transitions by type
cat("Expected transition counts:\n")
for (col in colnames(desc$count)) {
  cat(sprintf("  %s: %.4f\n", col, mean(desc$count[, col])))
}
cat("\n")

# Node posteriors for root
cat("Root posterior probabilities:\n")
root_idx <- which(rownames(desc$ace) == as.character(Ntip(tree) + 1))
if (length(root_idx) > 0) {
  print(desc$ace[root_idx, ])
}
cat("\n")

cat("=== Validation complete ===\n")
cat("\nNote: PhyKIT and phytools use different RNG implementations,\n")
cat("so permutation-based values (posteriors, transitions) will differ\n")
cat("slightly. Q matrix and log-likelihood should match closely.\n")
