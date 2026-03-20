#!/usr/bin/env Rscript
# validate_kmult.R
# Cross-validate K_mult against geomorph::physignal()
# Requires: ape, geomorph

suppressPackageStartupMessages({
  library(ape)
  library(geomorph)
})

tree <- read.tree("../sample_files/tree_simple.tre")
data <- read.table("../sample_files/tree_simple_multi_traits.tsv",
                    header=TRUE, sep="\t", row.names=1)

# Subset to shared taxa
shared <- intersect(tree$tip.label, rownames(data))
tree <- keep.tip(tree, shared)
data <- data[shared, , drop=FALSE]

# Convert to matrix
Y <- as.matrix(data)

# Run physignal (K_mult)
set.seed(42)
result <- physignal(Y, tree, iter=999)

cat("=== K_mult Validation ===\n")
cat(sprintf("K_mult: %.6f\n", result$phy.signal))
cat(sprintf("p-value: %.6f\n", result$pvalue))
cat(sprintf("Effect size (Z): %.6f\n", result$Z))
