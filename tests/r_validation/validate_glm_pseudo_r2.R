#!/usr/bin/env Rscript
# validate_glm_pseudo_r2.R
#
# Compute McFadden's pseudo-R^2 for phylogenetic GLMs.
#
# Formula:
#   McFadden's pseudo-R^2 = 1 - (LL_full / LL_null)
#
# where:
#   LL_full = log-likelihood of the full model (with predictors)
#   LL_null = log-likelihood of the null/intercept-only model
#
# McFadden's pseudo-R^2 is analogous to the coefficient of determination
# for linear models. Values between 0.2-0.4 are considered excellent fit
# in discrete choice models (McFadden, 1979).
#
# We compute this for two models:
#   1. Binomial (logistic): binary_trait ~ body_mass
#   2. Poisson: count_trait ~ body_mass
#
# For each, we fit both the full model and the intercept-only (null) model
# using phylolm::phyloglm().
#
# Requires: ape, phylolm
#
# Usage:
#   Rscript tests/r_validation/validate_glm_pseudo_r2.R

suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
})

cat("=== GLM McFadden's Pseudo-R^2 Validation ===\n\n")

# ---------------------------------------------------------------
# Load data
# ---------------------------------------------------------------
tree <- read.tree("../sample_files/tree_simple.tre")

traits <- read.delim("../sample_files/tree_simple_glm_traits.tsv",
  header = TRUE, row.names = 1, stringsAsFactors = FALSE)

cat("Tree tips:", paste(tree$tip.label, collapse = ", "), "\n")
cat("Trait data:\n")
print(traits)
cat("\n")

# Verify trait data alignment with tree
cat("Matching order:", all(rownames(traits) %in% tree$tip.label), "\n\n")

# ---------------------------------------------------------------
# 1. Binomial (logistic): binary_trait ~ body_mass
# ---------------------------------------------------------------
cat("=== Binomial (logistic): binary_trait ~ body_mass ===\n\n")

# Full model
fit_bin_full <- phyloglm(binary_trait ~ body_mass, data = traits, phy = tree,
  method = "logistic_MPLE", btol = 10, log.alpha.bound = 4)

cat("Full model:\n")
print(summary(fit_bin_full))
cat(sprintf("  logLik = %.10f\n", fit_bin_full$logLik))
cat(sprintf("  penalized logLik = %.10f\n", fit_bin_full$penlogLik))
cat(sprintf("  alpha = %.10f\n", fit_bin_full$alpha))
cat(sprintf("  coefficients: intercept = %.10f, body_mass = %.10f\n",
            coef(fit_bin_full)[1], coef(fit_bin_full)[2]))

# Null model (intercept only)
fit_bin_null <- phyloglm(binary_trait ~ 1, data = traits, phy = tree,
  method = "logistic_MPLE", btol = 10, log.alpha.bound = 4)

cat("\nNull model:\n")
print(summary(fit_bin_null))
cat(sprintf("  logLik = %.10f\n", fit_bin_null$logLik))
cat(sprintf("  penalized logLik = %.10f\n", fit_bin_null$penlogLik))
cat(sprintf("  alpha = %.10f\n", fit_bin_null$alpha))

# McFadden's pseudo-R^2 using standard log-likelihood
mcfadden_bin <- 1 - (fit_bin_full$logLik / fit_bin_null$logLik)

# McFadden's pseudo-R^2 using penalized log-likelihood
mcfadden_bin_pen <- 1 - (fit_bin_full$penlogLik / fit_bin_null$penlogLik)

cat(sprintf("\nMcFadden's pseudo-R^2 (logLik):    %.10f\n", mcfadden_bin))
cat(sprintf("McFadden's pseudo-R^2 (penlogLik): %.10f\n", mcfadden_bin_pen))
cat(sprintf("  = 1 - (%.6f / %.6f)\n\n", fit_bin_full$logLik, fit_bin_null$logLik))

# ---------------------------------------------------------------
# 2. Poisson: count_trait ~ body_mass
# ---------------------------------------------------------------
cat("=== Poisson: count_trait ~ body_mass ===\n\n")

# Full model
fit_poi_full <- phyloglm(count_trait ~ body_mass, data = traits, phy = tree,
  method = "poisson_GEE")

cat("Full model:\n")
print(summary(fit_poi_full))
cat(sprintf("  logLik = %.10f\n", fit_poi_full$logLik))
cat(sprintf("  coefficients: intercept = %.10f, body_mass = %.10f\n",
            coef(fit_poi_full)[1], coef(fit_poi_full)[2]))
cat(sprintf("  overdispersion (phi) = %.10f\n", fit_poi_full$phi))

# Null model (intercept only)
fit_poi_null <- phyloglm(count_trait ~ 1, data = traits, phy = tree,
  method = "poisson_GEE")

cat("\nNull model:\n")
print(summary(fit_poi_null))
cat(sprintf("  logLik = %.10f\n", fit_poi_null$logLik))
cat(sprintf("  overdispersion (phi) = %.10f\n", fit_poi_null$phi))

# McFadden's pseudo-R^2
mcfadden_poi <- 1 - (fit_poi_full$logLik / fit_poi_null$logLik)

cat(sprintf("\nMcFadden's pseudo-R^2: %.10f\n", mcfadden_poi))
cat(sprintf("  = 1 - (%.6f / %.6f)\n\n", fit_poi_full$logLik, fit_poi_null$logLik))

# ---------------------------------------------------------------
# Non-phylogenetic GLM cross-check
# ---------------------------------------------------------------
cat("=== Non-phylogenetic GLM cross-check ===\n\n")

# Standard GLM for reference
glm_bin_full <- glm(binary_trait ~ body_mass, data = traits, family = binomial)
glm_bin_null <- glm(binary_trait ~ 1, data = traits, family = binomial)
mcfadden_bin_std <- 1 - (logLik(glm_bin_full) / logLik(glm_bin_null))
cat(sprintf("Standard binomial GLM McFadden R^2: %.10f\n", mcfadden_bin_std))
cat(sprintf("  Full logLik = %.10f, Null logLik = %.10f\n",
            logLik(glm_bin_full), logLik(glm_bin_null)))

glm_poi_full <- glm(count_trait ~ body_mass, data = traits, family = poisson)
glm_poi_null <- glm(count_trait ~ 1, data = traits, family = poisson)
mcfadden_poi_std <- 1 - (logLik(glm_poi_full) / logLik(glm_poi_null))
cat(sprintf("Standard Poisson GLM McFadden R^2:  %.10f\n", mcfadden_poi_std))
cat(sprintf("  Full logLik = %.10f, Null logLik = %.10f\n\n",
            logLik(glm_poi_full), logLik(glm_poi_null)))

# ---------------------------------------------------------------
# Summary
# ---------------------------------------------------------------
cat("=== KEY VALUES FOR PYTHON TESTS ===\n")
cat("--- Binomial (logistic_MPLE) ---\n")
cat(sprintf("bin_full_logLik     = %.10f\n", fit_bin_full$logLik))
cat(sprintf("bin_null_logLik     = %.10f\n", fit_bin_null$logLik))
cat(sprintf("bin_full_penlogLik  = %.10f\n", fit_bin_full$penlogLik))
cat(sprintf("bin_null_penlogLik  = %.10f\n", fit_bin_null$penlogLik))
cat(sprintf("bin_mcfadden_r2     = %.10f\n", mcfadden_bin))
cat(sprintf("bin_mcfadden_r2_pen = %.10f\n", mcfadden_bin_pen))
cat(sprintf("bin_alpha           = %.10f\n", fit_bin_full$alpha))
cat(sprintf("bin_intercept       = %.10f\n", coef(fit_bin_full)[1]))
cat(sprintf("bin_slope           = %.10f\n", coef(fit_bin_full)[2]))
cat("--- Poisson (poisson_GEE) ---\n")
cat(sprintf("poi_full_logLik     = %.10f\n", fit_poi_full$logLik))
cat(sprintf("poi_null_logLik     = %.10f\n", fit_poi_null$logLik))
cat(sprintf("poi_mcfadden_r2     = %.10f\n", mcfadden_poi))
cat(sprintf("poi_intercept       = %.10f\n", coef(fit_poi_full)[1]))
cat(sprintf("poi_slope           = %.10f\n", coef(fit_poi_full)[2]))
cat(sprintf("poi_phi             = %.10f\n", fit_poi_full$phi))

cat("\nDone.\n")
