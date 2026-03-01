#!/usr/bin/env Rscript
# Verification script for PhyKIT phylogenetic GLM
# Generates reference values from R's phylolm package
#
# Install: install.packages("phylolm")
# Run:     Rscript tests/r_verification/verify_phyloglm.R

library(phylolm)
library(ape)

cat("=== PhyKIT phylogenetic GLM R verification ===\n\n")

# Read tree and trait data
tree <- read.tree("tests/sample_files/tree_simple.tre")
traits <- read.delim("tests/sample_files/tree_simple_glm_traits.tsv", row.names = 1)

cat("Tree tips:", tree$tip.label, "\n")
cat("Trait data:\n")
print(traits)
cat("\n")

# --- Logistic MPLE (binomial) ---
cat("=== Logistic MPLE: binary_trait ~ body_mass ===\n")
fit_log <- phyloglm(binary_trait ~ body_mass, phy = tree, data = traits,
                    method = "logistic_MPLE", btol = 10, log.alpha.bound = 4)
cat("\nSummary:\n")
print(summary(fit_log))
cat("\nCoefficients:\n")
print(coef(fit_log))
cat("\nAlpha:", fit_log$alpha, "\n")
cat("Log-likelihood:", fit_log$logLik, "\n")
cat("AIC:", fit_log$aic, "\n")
cat("Std errors:", fit_log$sd, "\n")
cat("z-values:", coef(fit_log) / fit_log$sd, "\n")
cat("p-values:", 2 * pnorm(-abs(coef(fit_log) / fit_log$sd)), "\n")
cat("Penalized log-likelihood:", fit_log$penlogLik, "\n")

# --- Poisson GEE ---
cat("\n=== Poisson GEE: count_trait ~ body_mass ===\n")
fit_poi <- phyloglm(count_trait ~ body_mass, phy = tree, data = traits,
                    method = "poisson_GEE")
cat("\nSummary:\n")
print(summary(fit_poi))
cat("\nCoefficients:\n")
print(coef(fit_poi))
cat("\nLog-likelihood:", fit_poi$logLik, "\n")
cat("AIC:", fit_poi$aic, "\n")
cat("Std errors:", fit_poi$sd, "\n")
cat("z-values:", coef(fit_poi) / fit_poi$sd, "\n")
cat("p-values:", 2 * pnorm(-abs(coef(fit_poi) / fit_poi$sd)), "\n")
cat("Overdispersion (phi):", fit_poi$phi, "\n")
cat("n.iter:", fit_poi$n.iter, "\n")

cat("\n=== Done ===\n")
