#!/usr/bin/env Rscript
# validate_pgls_r2.R
#
# Compute three-way R^2 variance decomposition for phylogenetic GLS (PGLS).
#
# Decomposition:
#   sigma^2_OLS = sigma^2_pred + sigma^2_phylo + sigma^2_resid
#
# Components:
#   R^2_pred  = 1 - (sigma^2_GLS / sigma^2_OLS)
#             = proportion of total variance explained by the predictor
#   R^2_phylo = (sigma^2_GLS - sigma^2_resid) / sigma^2_OLS
#             = proportion explained by phylogenetic covariance alone
#   R^2_resid = sigma^2_resid / sigma^2_OLS
#             = unexplained (residual) proportion
#
# where:
#   sigma^2_OLS   = var(y) with ML denominator = sum((y-mean(y))^2)/n
#   sigma^2_GLS   = ML residual variance from intercept-only PGLS
#   sigma^2_resid = ML residual variance from full PGLS model (y ~ x)
#
# These three components sum to 1:
#   R^2_pred + R^2_phylo + R^2_resid = 1
#
# We use caper::pgls() with lambda="ML" and cross-check with manual GLS.
#
# Requires: ape, caper, nlme
#
# Usage:
#   Rscript tests/r_validation/validate_pgls_r2.R

suppressPackageStartupMessages({
  library(ape)
  library(caper)
  library(nlme)
})

cat("=== PGLS Three-Way R^2 Decomposition Validation ===\n\n")

# ---------------------------------------------------------------
# Load data
# ---------------------------------------------------------------
tree <- read.tree("../sample_files/tree_simple.tre")

traits <- read.delim("../sample_files/tree_simple_multi_traits.tsv",
  header = TRUE, stringsAsFactors = FALSE)

cat("Trait data:\n")
print(traits)
cat("\n")

# ---------------------------------------------------------------
# caper PGLS: brain_size ~ body_mass
# ---------------------------------------------------------------
cat("--- caper PGLS: brain_size ~ body_mass ---\n\n")

comp_data <- comparative.data(tree, traits, names.col = "taxon",
  vcv = TRUE, vcv.dim = 3)

# Full model
pgls_full <- pgls(brain_size ~ body_mass, data = comp_data, lambda = "ML")
cat("Full model summary:\n")
print(summary(pgls_full))

# Intercept-only model (for sigma^2_GLS baseline)
pgls_null <- pgls(brain_size ~ 1, data = comp_data, lambda = "ML")
cat("\nNull (intercept-only) model summary:\n")
print(summary(pgls_null))

# Extract ML sigma^2 values
# caper pgls stores the ML estimated sigma^2 of residuals
sigma2_resid_full <- summary(pgls_full)$sigma2  # from full model
sigma2_gls_null <- summary(pgls_null)$sigma2     # from intercept-only model
lambda_full <- pgls_full$param["lambda"]
lambda_null <- pgls_null$param["lambda"]

cat(sprintf("\nFull model:  sigma^2 = %.10f, lambda = %.6f\n",
            sigma2_resid_full, lambda_full))
cat(sprintf("Null model:  sigma^2 = %.10f, lambda = %.6f\n",
            sigma2_gls_null, lambda_null))

# OLS ML variance (no phylogenetic correction)
y <- traits$brain_size[match(tree$tip.label, traits$taxon)]
n <- length(y)
sigma2_ols <- sum((y - mean(y))^2) / n

cat(sprintf("OLS var(y):  sigma^2 = %.10f\n\n", sigma2_ols))

# Three-way decomposition
r2_pred  <- 1 - (sigma2_gls_null / sigma2_ols)
# Note: For the phylo component, we need to be careful about what
# "sigma^2_GLS" and "sigma^2_resid" mean. The full model residual
# variance may use a different lambda than the null.
#
# Strict decomposition: using null model sigma^2 as baseline
r2_phylo_component <- (sigma2_gls_null - sigma2_resid_full) / sigma2_ols
r2_resid <- sigma2_resid_full / sigma2_ols

cat("--- Three-Way Decomposition: brain_size ~ body_mass ---\n")
cat(sprintf("R^2_pred   = 1 - (%.6f / %.6f) = %.10f\n",
            sigma2_gls_null, sigma2_ols, r2_pred))
cat(sprintf("R^2_phylo  = (%.6f - %.6f) / %.6f = %.10f\n",
            sigma2_gls_null, sigma2_resid_full, sigma2_ols, r2_phylo_component))
cat(sprintf("R^2_resid  = %.6f / %.6f = %.10f\n",
            sigma2_resid_full, sigma2_ols, r2_resid))
cat(sprintf("Sum check  = %.10f (should be ~1.0)\n\n",
            r2_pred + r2_phylo_component + r2_resid))

# ---------------------------------------------------------------
# Manual GLS cross-check
# ---------------------------------------------------------------
cat("--- Manual GLS cross-check ---\n\n")

C <- vcv(tree, corr = FALSE)
# Reorder to match tree tip labels
taxa_order <- tree$tip.label
y_vec <- setNames(traits$brain_size, traits$taxon)[taxa_order]
x_vec <- setNames(traits$body_mass, traits$taxon)[taxa_order]
C_ordered <- C[taxa_order, taxa_order]

ones <- rep(1, n)
C_inv <- solve(C_ordered)

# Intercept-only GLS: y ~ 1
beta0_gls <- as.numeric(solve(t(ones) %*% C_inv %*% ones) %*% (t(ones) %*% C_inv %*% y_vec))
resid0 <- y_vec - beta0_gls
sigma2_gls_manual <- as.numeric(t(resid0) %*% C_inv %*% resid0) / n

# Full GLS: y ~ 1 + x
X <- cbind(ones, x_vec)
beta_gls <- as.numeric(solve(t(X) %*% C_inv %*% X) %*% (t(X) %*% C_inv %*% y_vec))
resid_full <- y_vec - X %*% beta_gls
sigma2_full_manual <- as.numeric(t(resid_full) %*% C_inv %*% resid_full) / n

# OLS variance
sigma2_ols_manual <- sum((y_vec - mean(y_vec))^2) / n

cat(sprintf("Manual GLS null sigma^2 = %.10f\n", sigma2_gls_manual))
cat(sprintf("Manual GLS full sigma^2 = %.10f\n", sigma2_full_manual))
cat(sprintf("Manual OLS sigma^2      = %.10f\n", sigma2_ols_manual))

r2_pred_manual  <- 1 - (sigma2_gls_manual / sigma2_ols_manual)
r2_phylo_manual <- (sigma2_gls_manual - sigma2_full_manual) / sigma2_ols_manual
r2_resid_manual <- sigma2_full_manual / sigma2_ols_manual

cat(sprintf("\nManual R^2_pred  = %.10f\n", r2_pred_manual))
cat(sprintf("Manual R^2_phylo = %.10f\n", r2_phylo_manual))
cat(sprintf("Manual R^2_resid = %.10f\n", r2_resid_manual))
cat(sprintf("Manual sum check = %.10f\n\n", r2_pred_manual + r2_phylo_manual + r2_resid_manual))

# ---------------------------------------------------------------
# Second regression: longevity ~ body_mass
# ---------------------------------------------------------------
cat("--- caper PGLS: longevity ~ body_mass ---\n\n")

pgls_full2 <- pgls(longevity ~ body_mass, data = comp_data, lambda = "ML")
pgls_null2 <- pgls(longevity ~ 1, data = comp_data, lambda = "ML")

sigma2_resid_full2 <- summary(pgls_full2)$sigma2
sigma2_gls_null2 <- summary(pgls_null2)$sigma2

y2 <- traits$longevity[match(tree$tip.label, traits$taxon)]
sigma2_ols2 <- sum((y2 - mean(y2))^2) / n

r2_pred2  <- 1 - (sigma2_gls_null2 / sigma2_ols2)
r2_phylo2 <- (sigma2_gls_null2 - sigma2_resid_full2) / sigma2_ols2
r2_resid2 <- sigma2_resid_full2 / sigma2_ols2

cat(sprintf("Full sigma^2  = %.10f, lambda = %.6f\n", sigma2_resid_full2, pgls_full2$param["lambda"]))
cat(sprintf("Null sigma^2  = %.10f, lambda = %.6f\n", sigma2_gls_null2, pgls_null2$param["lambda"]))
cat(sprintf("OLS sigma^2   = %.10f\n", sigma2_ols2))
cat(sprintf("R^2_pred      = %.10f\n", r2_pred2))
cat(sprintf("R^2_phylo     = %.10f\n", r2_phylo2))
cat(sprintf("R^2_resid     = %.10f\n", r2_resid2))
cat(sprintf("Sum check     = %.10f\n\n", r2_pred2 + r2_phylo2 + r2_resid2))

# ---------------------------------------------------------------
# nlme::gls cross-check (with corBrownian)
# ---------------------------------------------------------------
cat("--- nlme::gls cross-check ---\n\n")

df <- data.frame(
  taxon = taxa_order,
  brain_size = y_vec,
  body_mass = x_vec,
  row.names = taxa_order
)

gls_full_nlme <- gls(brain_size ~ body_mass, data = df,
  correlation = corBrownian(1, tree, form = ~taxon),
  method = "ML")
gls_null_nlme <- gls(brain_size ~ 1, data = df,
  correlation = corBrownian(1, tree, form = ~taxon),
  method = "ML")

cat(sprintf("nlme GLS full: sigma = %.10f, logLik = %.6f\n",
            gls_full_nlme$sigma, logLik(gls_full_nlme)))
cat(sprintf("nlme GLS null: sigma = %.10f, logLik = %.6f\n",
            gls_null_nlme$sigma, logLik(gls_null_nlme)))
# Note: nlme sigma is the standard deviation, not variance
cat(sprintf("nlme GLS full sigma^2 = %.10f\n", gls_full_nlme$sigma^2))
cat(sprintf("nlme GLS null sigma^2 = %.10f\n\n", gls_null_nlme$sigma^2))

# ---------------------------------------------------------------
# Summary
# ---------------------------------------------------------------
cat("=== KEY VALUES FOR PYTHON TESTS ===\n")
cat("--- brain_size ~ body_mass (manual GLS, lambda=1/BM) ---\n")
cat(sprintf("sigma2_ols          = %.10f\n", sigma2_ols_manual))
cat(sprintf("sigma2_gls_null     = %.10f\n", sigma2_gls_manual))
cat(sprintf("sigma2_gls_full     = %.10f\n", sigma2_full_manual))
cat(sprintf("r2_pred             = %.10f\n", r2_pred_manual))
cat(sprintf("r2_phylo            = %.10f\n", r2_phylo_manual))
cat(sprintf("r2_resid            = %.10f\n", r2_resid_manual))
cat(sprintf("beta_intercept      = %.10f\n", beta_gls[1]))
cat(sprintf("beta_slope          = %.10f\n", beta_gls[2]))
cat("--- brain_size ~ body_mass (caper, lambda=ML) ---\n")
cat(sprintf("caper_sigma2_full   = %.10f\n", sigma2_resid_full))
cat(sprintf("caper_sigma2_null   = %.10f\n", sigma2_gls_null))
cat(sprintf("caper_lambda_full   = %.10f\n", lambda_full))
cat(sprintf("caper_lambda_null   = %.10f\n", lambda_null))
cat(sprintf("caper_r2_pred       = %.10f\n", r2_pred))
cat(sprintf("caper_r2_phylo      = %.10f\n", r2_phylo_component))
cat(sprintf("caper_r2_resid      = %.10f\n", r2_resid))

cat("\nDone.\n")
