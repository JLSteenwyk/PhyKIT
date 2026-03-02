#!/usr/bin/env Rscript
# validate_fit_continuous_r2.R
#
# Compute R^2_model for continuous trait evolution models via
# geiger::fitContinuous().
#
# Formula:
#   R^2_model = 1 - (sigma^2_model / sigma^2_White)
#
# where:
#   sigma^2_model = ML rate parameter (sigsq) from the fitted model
#   sigma^2_White = ML rate parameter from white noise (no phylogeny)
#
# Interpretation: How much better does this evolutionary model explain
# the trait variance compared to no evolutionary model at all (white noise)?
# For BM, this is equivalent to phylogenetic signal R^2. For OU, EB, etc.,
# it captures the additional structure imposed by the model.
#
# Models fitted:
#   - BM (Brownian motion)
#   - OU (Ornstein-Uhlenbeck)
#   - EB (Early Burst / ACDC)
#   - White (White Noise / star tree)
#   - Lambda (Pagel's lambda transformation)
#   - Kappa (Pagel's kappa transformation)
#   - Delta (Pagel's delta transformation)
#
# Requires: ape, geiger
#
# Usage:
#   Rscript tests/r_validation/validate_fit_continuous_r2.R

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
})

cat("=== fit_continuous R^2 Validation ===\n\n")

# ---------------------------------------------------------------
# Load data
# ---------------------------------------------------------------
tree <- read.tree("../sample_files/tree_simple.tre")

trait_data <- read.delim("../sample_files/tree_simple_traits.tsv",
  header = FALSE, comment.char = "#",
  col.names = c("taxon", "trait"), stringsAsFactors = FALSE)

y <- setNames(trait_data$trait, trait_data$taxon)
n <- length(y)

cat("Tree tips:", paste(tree$tip.label, collapse = ", "), "\n")
cat("Trait values:", paste(sprintf("%.2f", y), collapse = ", "), "\n")
cat("n =", n, "\n\n")

# ---------------------------------------------------------------
# Fit all models
# ---------------------------------------------------------------
cat("--- Fitting all models via geiger::fitContinuous ---\n\n")

models <- c("BM", "OU", "EB", "white", "lambda", "kappa", "delta")
results <- list()

for (model_name in models) {
  tryCatch({
    fit <- fitContinuous(tree, y, model = model_name)
    results[[model_name]] <- fit
    cat(sprintf("%-8s: sigsq = %.10f, z0 = %.10f, logLik = %.6f, AIC = %.4f, AICc = %.4f, k = %d\n",
                model_name,
                fit$opt$sigsq,
                fit$opt$z0,
                fit$opt$lnL,
                fit$opt$aic,
                fit$opt$aicc,
                fit$opt$k))

    # Print model-specific parameters
    if (model_name == "OU") {
      cat(sprintf("          alpha = %.10f\n", fit$opt$alpha))
    } else if (model_name == "EB") {
      cat(sprintf("          a (rate change) = %.10f\n", fit$opt$a))
    } else if (model_name == "lambda") {
      cat(sprintf("          lambda = %.10f\n", fit$opt$lambda))
    } else if (model_name == "kappa") {
      cat(sprintf("          kappa = %.10f\n", fit$opt$kappa))
    } else if (model_name == "delta") {
      cat(sprintf("          delta = %.10f\n", fit$opt$delta))
    }
  }, error = function(e) {
    cat(sprintf("%-8s: ERROR - %s\n", model_name, e$message))
  })
}

# ---------------------------------------------------------------
# Compute R^2 for each model relative to white noise
# ---------------------------------------------------------------
cat("\n--- R^2_model = 1 - (sigma^2_model / sigma^2_White) ---\n\n")

if (!is.null(results[["white"]])) {
  sigma2_white <- results[["white"]]$opt$sigsq
  cat(sprintf("sigma^2_White = %.10f\n\n", sigma2_white))

  cat(sprintf("%-8s  %15s  %15s  %10s\n", "Model", "sigma^2", "R^2_model", "Delta_logLik"))
  cat(paste(rep("-", 58), collapse = ""), "\n")

  for (model_name in models) {
    if (!is.null(results[[model_name]])) {
      sigma2_model <- results[[model_name]]$opt$sigsq
      r2 <- 1 - (sigma2_model / sigma2_white)
      delta_ll <- results[[model_name]]$opt$lnL - results[["white"]]$opt$lnL
      cat(sprintf("%-8s  %15.10f  %15.10f  %10.4f\n",
                  model_name, sigma2_model, r2, delta_ll))
    }
  }
} else {
  cat("ERROR: White noise model failed to fit.\n")
}

# ---------------------------------------------------------------
# Alternative R^2 using log-likelihoods (likelihood-based R^2)
# ---------------------------------------------------------------
cat("\n--- Alternative: Likelihood-ratio based pseudo-R^2 ---\n")
cat("--- R^2_LR = 1 - exp(-(2/n)(logL_model - logL_white)) ---\n\n")

if (!is.null(results[["white"]])) {
  logL_white <- results[["white"]]$opt$lnL

  for (model_name in models) {
    if (!is.null(results[[model_name]])) {
      logL_model <- results[[model_name]]$opt$lnL
      r2_lr <- 1 - exp(-(2/n) * (logL_model - logL_white))
      cat(sprintf("%-8s  R^2_LR = %.10f  (logL = %.6f)\n",
                  model_name, r2_lr, logL_model))
    }
  }
}

# ---------------------------------------------------------------
# Multi-trait: brain_size and longevity
# ---------------------------------------------------------------
cat("\n--- Multi-trait validation ---\n\n")

multi_traits <- read.delim("../sample_files/tree_simple_multi_traits.tsv",
  header = TRUE, stringsAsFactors = FALSE)

for (trait_col in c("body_mass", "brain_size", "longevity")) {
  trait_vec <- setNames(multi_traits[[trait_col]], multi_traits$taxon)

  fit_bm <- fitContinuous(tree, trait_vec, model = "BM")
  fit_wn <- fitContinuous(tree, trait_vec, model = "white")

  r2 <- 1 - (fit_bm$opt$sigsq / fit_wn$opt$sigsq)
  cat(sprintf("%-12s: BM sigsq = %.10f, WN sigsq = %.10f, R^2 = %.10f\n",
              trait_col, fit_bm$opt$sigsq, fit_wn$opt$sigsq, r2))
}

# ---------------------------------------------------------------
# Summary
# ---------------------------------------------------------------
cat("\n=== KEY VALUES FOR PYTHON TESTS ===\n")
cat("--- body_mass trait ---\n")
for (model_name in models) {
  if (!is.null(results[[model_name]])) {
    r2 <- 1 - (results[[model_name]]$opt$sigsq / results[["white"]]$opt$sigsq)
    cat(sprintf("%s_sigsq   = %.10f\n", model_name, results[[model_name]]$opt$sigsq))
    cat(sprintf("%s_z0      = %.10f\n", model_name, results[[model_name]]$opt$z0))
    cat(sprintf("%s_logLik  = %.10f\n", model_name, results[[model_name]]$opt$lnL))
    cat(sprintf("%s_r2      = %.10f\n", model_name, r2))

    # Model-specific parameters
    if (model_name == "OU") {
      cat(sprintf("%s_alpha   = %.10f\n", model_name, results[[model_name]]$opt$alpha))
    } else if (model_name == "EB") {
      cat(sprintf("%s_a       = %.10f\n", model_name, results[[model_name]]$opt$a))
    } else if (model_name == "lambda") {
      cat(sprintf("%s_lambda  = %.10f\n", model_name, results[[model_name]]$opt$lambda))
    } else if (model_name == "kappa") {
      cat(sprintf("%s_kappa   = %.10f\n", model_name, results[[model_name]]$opt$kappa))
    } else if (model_name == "delta") {
      cat(sprintf("%s_delta   = %.10f\n", model_name, results[[model_name]]$opt$delta))
    }
  }
}

cat("\nDone.\n")
