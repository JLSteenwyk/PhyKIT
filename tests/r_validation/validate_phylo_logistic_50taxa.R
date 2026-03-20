#!/usr/bin/env Rscript
# validate_phylo_logistic_50taxa.R
#
# Cross-validate phylo_logistic against phylolm::phyloglm
# using a 50-taxon dataset for more stable alpha estimation.
#
# Requires: ape, phylolm
#
# Usage:
#   cd tests/r_validation
#   Rscript validate_phylo_logistic_50taxa.R

suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
})

tree <- read.tree("../sample_files/tree_50taxa_logistic.tre")
data <- read.table("../sample_files/logistic_50taxa_traits.tsv",
                    header = TRUE, sep = "\t", row.names = 1)

cat("=== 50-taxon Phylogenetic Logistic Regression ===\n\n")

fit <- phyloglm(binary_y ~ continuous_x, data = data, phy = tree,
                method = "logistic_MPLE")
cat(capture.output(summary(fit)), sep = "\n")

cat("\n\n=== KEY VALUES FOR PYTHON COMPARISON ===\n")
cat(sprintf("Alpha: %.6f\n", fit$alpha))
cat(sprintf("Intercept estimate: %.6f\n", coef(fit)[1]))
cat(sprintf("Intercept SE: %.6f\n", fit$sd[1]))
cat(sprintf("continuous_x estimate: %.6f\n", coef(fit)[2]))
cat(sprintf("continuous_x SE: %.6f\n", fit$sd[2]))
cat(sprintf("Log-likelihood: %.6f\n", fit$logLik))
cat(sprintf("Penalized log-likelihood: %.6f\n", fit$penlogLik))
cat(sprintf("AIC: %.6f\n", fit$aic))
cat(sprintf("Convergence: %d\n", fit$convergence))
