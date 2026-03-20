#!/usr/bin/env Rscript
# Validation script for PhyKIT phylo_logistic command
# Compare output against R's phylolm package (Ives & Garland 2010)
#
# Usage: Rscript validate_phylo_logistic.R
#
# Requires: ape, phylolm

library(ape)
library(phylolm)

tree <- read.tree("../sample_files/tree_simple.tre")

# --- Test 1: Using tree_simple_glm_traits.tsv (binary_trait ~ body_mass) ---

data_glm <- read.table(
    "../sample_files/tree_simple_glm_traits.tsv",
    header = TRUE, sep = "\t", row.names = 1
)

fit_glm <- phyloglm(
    binary_trait ~ body_mass,
    data = data_glm,
    phy = tree,
    method = "logistic_MPLE"
)

cat("=== Test 1: binary_trait ~ body_mass (tree_simple_glm_traits.tsv) ===\n")
summary(fit_glm)

cat("\n=== KEY VALUES FOR PYTHON TESTS (Test 1) ===\n")
cat(sprintf("Alpha: %.6f\n", fit_glm$alpha))
cat(sprintf("Intercept: %.6f\n", coef(fit_glm)[1]))
cat(sprintf("body_mass coef: %.6f\n", coef(fit_glm)[2]))
cat(sprintf("Log-likelihood: %.6f\n", fit_glm$logLik))
cat(sprintf("AIC: %.6f\n", fit_glm$aic))
cat(sprintf("Penalized log-lik: %.6f\n", fit_glm$penlogLik))

# --- Test 2: Using tree_simple_binary_traits.tsv (has_wings ~ body_mass) ---

data_bin <- read.table(
    "../sample_files/tree_simple_binary_traits.tsv",
    header = TRUE, sep = "\t", row.names = 1
)

fit_bin <- phyloglm(
    has_wings ~ body_mass,
    data = data_bin,
    phy = tree,
    method = "logistic_MPLE"
)

cat("\n=== Test 2: has_wings ~ body_mass (tree_simple_binary_traits.tsv) ===\n")
summary(fit_bin)

cat("\n=== KEY VALUES FOR PYTHON TESTS (Test 2) ===\n")
cat(sprintf("Alpha: %.6f\n", fit_bin$alpha))
cat(sprintf("Intercept: %.6f\n", coef(fit_bin)[1]))
cat(sprintf("body_mass coef: %.6f\n", coef(fit_bin)[2]))
cat(sprintf("Log-likelihood: %.6f\n", fit_bin$logLik))
cat(sprintf("AIC: %.6f\n", fit_bin$aic))
cat(sprintf("Penalized log-lik: %.6f\n", fit_bin$penlogLik))
