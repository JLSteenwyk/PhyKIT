#!/usr/bin/env Rscript
# validate_phylo_impute.R
#
# Cross-validate phylo_impute against Rphylopars::phylopars()
#
# Requires: ape, Rphylopars
#
# Usage:
#   cd tests/r_validation
#   Rscript validate_phylo_impute.R

suppressPackageStartupMessages({
  library(ape)
  library(Rphylopars)
})

tree <- read.tree("../sample_files/tree_simple.tre")

# Read trait data with missing values
data <- read.table("../sample_files/tree_simple_traits_missing.tsv",
                    header = TRUE, sep = "\t")

# Rphylopars requires a 'species' column
colnames(data)[1] <- "species"

# Replace ? with NA
for (col in 2:ncol(data)) {
  data[[col]] <- as.numeric(as.character(data[[col]]))
}

cat("=== Input data ===\n")
print(data)

# Run phylopars imputation (BM model)
result <- phylopars(trait_data = data, tree = tree)

cat("\n=== Imputed values (anc_recon for tips) ===\n")
print(result$anc_recon[1:length(tree$tip.label), ])

cat("\n=== KEY VALUES FOR PYTHON COMPARISON ===\n")
tip_names <- tree$tip.label
recon <- result$anc_recon[1:length(tip_names), ]
rownames(recon) <- tip_names

# Print imputed values for missing entries
missing_entries <- list(
  c("bear", "brain_size"),
  c("sea_lion", "longevity"),
  c("monkey", "longevity")
)

for (entry in missing_entries) {
  taxon <- entry[1]
  trait <- entry[2]
  val <- recon[taxon, trait]
  cat(sprintf("%s %s: %.6f\n", taxon, trait, val))
}
