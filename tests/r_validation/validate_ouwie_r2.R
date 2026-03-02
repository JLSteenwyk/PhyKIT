#!/usr/bin/env Rscript
# validate_ouwie_r2.R
#
# Compute R^2_OU for each OUwie model relative to a BM1 baseline.
#
# Formula:
#   R^2_OU = 1 - (sigma^2_model / sigma^2_BM1)
#
# where:
#   sigma^2_model = estimated rate parameter from the fitted OU/BM model
#   sigma^2_BM1   = rate parameter from single-rate Brownian motion
#
# For multi-regime models (BMS, OUM, OUMV, OUMA, OUMVA), sigma^2_model
# is the weighted average of per-regime sigma^2 values, weighted by the
# proportion of total tree length assigned to each regime.
#
# Interpretation: How much of the BM-baseline evolutionary rate is
# "absorbed" by the more complex model structure (OU attraction,
# regime-specific parameters)?
#
# Models:
#   BM1   - single-rate Brownian motion (baseline; R^2 = 0 by definition)
#   BMS   - multi-rate Brownian motion (different sigma^2 per regime)
#   OU1   - single-optimum Ornstein-Uhlenbeck
#   OUM   - multi-optimum OU (different theta per regime)
#   OUMV  - multi-optimum OU with variable sigma^2
#   OUMA  - multi-optimum OU with variable alpha
#   OUMVA - multi-optimum OU with variable alpha and sigma^2
#
# Requires: ape, OUwie, phytools
#
# Usage:
#   Rscript tests/r_validation/validate_ouwie_r2.R

suppressPackageStartupMessages({
  library(ape)
  library(OUwie)
  library(phytools)
})

cat("=== OUwie R^2 Validation ===\n\n")

# ---------------------------------------------------------------
# Load and prepare data
# ---------------------------------------------------------------
tree <- read.tree("../sample_files/tree_simple.tre")

trait_data <- read.delim("../sample_files/tree_simple_traits.tsv",
  header = FALSE, comment.char = "#",
  col.names = c("taxon", "trait"), stringsAsFactors = FALSE)

regime_data <- read.delim("../sample_files/tree_simple_regimes.tsv",
  header = FALSE, comment.char = "#",
  col.names = c("taxon", "regime"), stringsAsFactors = FALSE)

# Merge into OUwie format: species, regime, trait
data <- merge(regime_data, trait_data, by = "taxon")
colnames(data) <- c("species", "regime", "trait")

cat("Data:\n")
print(data)
cat("\n")

# ---------------------------------------------------------------
# Root tree and assign regime labels to internal nodes
# ---------------------------------------------------------------

# Root at dog (to resolve the trifurcation, matching PhyKIT)
tree <- root(tree, "dog", resolve.root = TRUE)
tree$edge.length[tree$edge.length == 0] <- 1e-6

cat("Rooted tree tip labels:", paste(tree$tip.label, collapse = ", "), "\n")
cat("Number of internal nodes:", tree$Nnode, "\n")

# Assign regime labels to internal nodes via ML ACE
tip_regimes <- setNames(data$regime, data$species)
regime_vec <- tip_regimes[tree$tip.label]
regime_factor <- as.factor(regime_vec)
names(regime_factor) <- tree$tip.label

anc <- ace(regime_factor, tree, type = "discrete", method = "ML")
node_regimes <- apply(anc$lik.anc, 1, function(x)
  colnames(anc$lik.anc)[which.max(x)])
tree$node.label <- node_regimes

cat("Node labels:", tree$node.label, "\n\n")

# ---------------------------------------------------------------
# Compute regime-specific tree lengths (for weighted averages)
# ---------------------------------------------------------------

# Use simmap to get precise per-regime branch lengths
set.seed(42)
simmap_tree <- make.simmap(tree, regime_factor, model = "ER", nsim = 1,
                           message = FALSE)

all_regimes <- unique(c(as.character(regime_factor), node_regimes))
total_per_regime <- sapply(all_regimes, function(r) {
  sum(sapply(simmap_tree$maps, function(m) {
    if (r %in% names(m)) m[r] else 0
  }))
})
total_tree_length <- sum(tree$edge.length)
regime_weights <- total_per_regime / total_tree_length

cat("Regime-specific tree lengths:\n")
for (r in names(total_per_regime)) {
  cat(sprintf("  %s: %.4f (weight = %.6f)\n", r,
              total_per_regime[r], regime_weights[r]))
}
cat(sprintf("Total tree length: %.4f\n\n", total_tree_length))

# ---------------------------------------------------------------
# Fit all OUwie models
# ---------------------------------------------------------------
cat("--- Fitting OUwie models (root.station=FALSE) ---\n\n")

model_names <- c("BM1", "BMS", "OU1", "OUM", "OUMV", "OUMA", "OUMVA")
results <- list()

for (model_name in model_names) {
  cat(sprintf("Fitting %s...\n", model_name))
  tryCatch({
    fit <- OUwie(tree, data, model = model_name, root.station = FALSE,
                 algorithm = "invert", quiet = TRUE,
                 check.identify = FALSE)
    results[[model_name]] <- fit

    cat(sprintf("  logLik = %.10f\n", fit$loglik))
    cat(sprintf("  AIC    = %.4f\n", fit$AIC))
    cat(sprintf("  AICc   = %.4f\n", fit$AICc))
    cat(sprintf("  k      = %d\n", fit$param.count))

    cat("  Solution matrix:\n")
    print(fit$solution)

    # Extract sigma^2 values
    sigma2_vals <- fit$solution["sigma.sq", ]
    cat("  sigma.sq:", paste(sprintf("%.10f", sigma2_vals), collapse = ", "), "\n")

    # Extract alpha values if OU model
    if (model_name %in% c("OU1", "OUM", "OUMV", "OUMA", "OUMVA")) {
      alpha_vals <- fit$solution["alpha", ]
      cat("  alpha:", paste(sprintf("%.10f", alpha_vals), collapse = ", "), "\n")
    }

    # Extract theta (optima) if available
    if (!is.null(fit$theta)) {
      cat("  theta:", paste(sprintf("%.10f", fit$theta[, 1]), collapse = ", "), "\n")
    }

    cat("\n")
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n\n", e$message))
  })
}

# ---------------------------------------------------------------
# Compute R^2_OU for each model relative to BM1
# ---------------------------------------------------------------
cat("--- R^2_OU = 1 - (sigma^2_model / sigma^2_BM1) ---\n\n")

if (!is.null(results[["BM1"]])) {
  sigma2_bm1 <- results[["BM1"]]$solution["sigma.sq", 1]
  cat(sprintf("BM1 sigma^2 (baseline) = %.10f\n\n", sigma2_bm1))

  cat(sprintf("%-8s  %15s  %15s  %10s  %10s\n",
              "Model", "sigma^2_eff", "R^2_OU", "logLik", "AICc"))
  cat(paste(rep("-", 68), collapse = ""), "\n")

  for (model_name in model_names) {
    if (!is.null(results[[model_name]])) {
      fit <- results[[model_name]]
      sigma2_vals <- fit$solution["sigma.sq", ]

      # For single-regime models, use the single sigma^2
      # For multi-regime models, compute weighted average
      if (length(sigma2_vals) == 1 || all(sigma2_vals == sigma2_vals[1])) {
        sigma2_eff <- sigma2_vals[1]
      } else {
        # Match regime names to weights
        regime_names <- colnames(fit$solution)
        w <- regime_weights[regime_names]
        # Normalize weights in case they don't sum to 1
        w <- w / sum(w)
        sigma2_eff <- sum(sigma2_vals * w)
      }

      r2_ou <- 1 - (sigma2_eff / sigma2_bm1)

      cat(sprintf("%-8s  %15.10f  %15.10f  %10.4f  %10.4f\n",
                  model_name, sigma2_eff, r2_ou,
                  fit$loglik, fit$AICc))
    }
  }
} else {
  cat("ERROR: BM1 model failed; cannot compute R^2.\n")
}

# ---------------------------------------------------------------
# Alternative R^2 using only fitted sigma^2 (no weighting)
# ---------------------------------------------------------------
cat("\n--- Alternative: Mean sigma^2 (unweighted) ---\n\n")

if (!is.null(results[["BM1"]])) {
  for (model_name in model_names) {
    if (!is.null(results[[model_name]])) {
      sigma2_vals <- results[[model_name]]$solution["sigma.sq", ]
      sigma2_mean <- mean(sigma2_vals)
      r2_mean <- 1 - (sigma2_mean / sigma2_bm1)
      cat(sprintf("%-8s  mean_sigma^2 = %.10f  R^2 = %.10f\n",
                  model_name, sigma2_mean, r2_mean))
    }
  }
}

# ---------------------------------------------------------------
# Detailed solution matrices
# ---------------------------------------------------------------
cat("\n--- Full solution matrices ---\n\n")

for (model_name in model_names) {
  if (!is.null(results[[model_name]])) {
    cat(sprintf("=== %s ===\n", model_name))
    cat("Solution:\n")
    print(results[[model_name]]$solution)
    if (!is.null(results[[model_name]]$theta)) {
      cat("Theta:\n")
      print(results[[model_name]]$theta)
    }
    cat("\n")
  }
}

# ---------------------------------------------------------------
# Summary
# ---------------------------------------------------------------
cat("=== KEY VALUES FOR PYTHON TESTS ===\n")

if (!is.null(results[["BM1"]])) {
  cat(sprintf("bm1_sigma2       = %.10f\n", sigma2_bm1))
  cat(sprintf("bm1_logLik       = %.10f\n", results[["BM1"]]$loglik))
  cat(sprintf("bm1_AICc         = %.10f\n", results[["BM1"]]$AICc))

  for (model_name in model_names) {
    if (!is.null(results[[model_name]])) {
      fit <- results[[model_name]]
      sigma2_vals <- fit$solution["sigma.sq", ]

      cat(sprintf("\n--- %s ---\n", model_name))
      cat(sprintf("%s_logLik   = %.10f\n", model_name, fit$loglik))
      cat(sprintf("%s_AICc     = %.10f\n", model_name, fit$AICc))
      cat(sprintf("%s_k        = %d\n", model_name, fit$param.count))

      for (i in seq_along(sigma2_vals)) {
        col_name <- colnames(fit$solution)[i]
        cat(sprintf("%s_sigma2_%s = %.10f\n", model_name, col_name, sigma2_vals[i]))
      }

      if (model_name %in% c("OU1", "OUM", "OUMV", "OUMA", "OUMVA")) {
        alpha_vals <- fit$solution["alpha", ]
        for (i in seq_along(alpha_vals)) {
          col_name <- colnames(fit$solution)[i]
          cat(sprintf("%s_alpha_%s  = %.10f\n", model_name, col_name, alpha_vals[i]))
        }
      }

      if (!is.null(fit$theta)) {
        theta_vals <- fit$theta[, 1]
        for (i in seq_along(theta_vals)) {
          cat(sprintf("%s_theta_%d  = %.10f\n", model_name, i, theta_vals[i]))
        }
      }

      # Effective sigma^2 and R^2
      if (length(sigma2_vals) == 1 || all(sigma2_vals == sigma2_vals[1])) {
        sigma2_eff <- sigma2_vals[1]
      } else {
        regime_names <- colnames(fit$solution)
        w <- regime_weights[regime_names]
        w <- w / sum(w)
        sigma2_eff <- sum(sigma2_vals * w)
      }
      r2_ou <- 1 - (sigma2_eff / sigma2_bm1)
      cat(sprintf("%s_sigma2_eff = %.10f\n", model_name, sigma2_eff))
      cat(sprintf("%s_r2_ou     = %.10f\n", model_name, r2_ou))
    }
  }
}

cat("\nDone.\n")
