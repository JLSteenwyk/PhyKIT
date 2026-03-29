#!/usr/bin/env Rscript
# validate_phylo_path.R
#
# Cross-validate PhyKIT phylo_path against R phylopath package.
#
# Validates:
#   - Fisher's C statistic per model
#   - CICc, delta CICc, and model weights
#   - Path coefficients for best model
#
# Requires: ape, phylopath
#
# Usage:
#   Rscript tests/r_validation/validate_phylo_path.R

suppressPackageStartupMessages({
  library(ape)
  library(phylopath)
})

cat("=== Phylogenetic Path Analysis Validation ===\n\n")

# ---------------------------------------------------------------
# Load data
# ---------------------------------------------------------------
tree <- read.tree("../sample_files/tree_simple.tre")
traits <- read.delim("../sample_files/tree_simple_multi_traits.tsv",
  header = TRUE, stringsAsFactors = FALSE)
rownames(traits) <- traits$taxon
traits$taxon <- NULL

cat("Tree tips:", tree$tip.label, "\n")
cat("Trait columns:", names(traits), "\n\n")

# ---------------------------------------------------------------
# Define candidate DAGs (matching PhyKIT models file)
# ---------------------------------------------------------------
direct <- define_model_set(
  direct   = c(body_mass %->% brain_size, body_mass %->% longevity),
  mediated = c(body_mass %->% brain_size, brain_size %->% longevity),
  full     = c(body_mass %->% brain_size, body_mass %->% longevity, brain_size %->% longevity)
)

cat("--- Running phylo_path ---\n\n")
result <- phylo_path(direct, data = traits, tree = tree, model = "lambda")

cat("Model comparison:\n")
s <- summary(result)
print(s)
cat("\n")

cat("Detailed model results:\n")
for (nm in names(result$d_sep)) {
  cat(sprintf("  %s: C = %.4f, p = %.4f, CICc = %.4f\n",
    nm, result$d_sep[[nm]]$C, result$d_sep[[nm]]$p,
    result$d_sep[[nm]]$CICc))
}
cat("\n")

# ---------------------------------------------------------------
# Best model coefficients
# ---------------------------------------------------------------
cat("--- Best model ---\n\n")
best <- best(result)
cat("Best model coefficients:\n")
print(coef_plot(best, error_bar = "se"))
cat("\n")

# ---------------------------------------------------------------
# Model-averaged coefficients
# ---------------------------------------------------------------
cat("--- Model-averaged coefficients ---\n\n")
avg <- average(result, avg_method = "conditional")
cat("Averaged coefficients:\n")
print(coef_plot(avg, error_bar = "se"))
cat("\n")

cat("=== Validation complete ===\n")
