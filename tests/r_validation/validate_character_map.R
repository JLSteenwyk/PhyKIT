#!/usr/bin/env Rscript
# validate_character_map.R
#
# Cross-validate character_map CI/RI and parsimony score against phangorn.
#
# Requires: ape, phangorn
#
# Usage:
#   cd tests/r_validation
#   Rscript validate_character_map.R

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})

tree <- read.tree("../sample_files/tree_character_map.tre")
tree <- multi2di(tree)

# Read character matrix as phyDat
data_raw <- read.table(
  "../sample_files/character_matrix_simple.tsv",
  header = TRUE, sep = "\t", row.names = 1,
  colClasses = "character"
)

# Convert to phyDat with USER type
data_matrix <- as.matrix(data_raw)
levels <- sort(unique(as.vector(data_matrix)))
levels <- levels[levels != "?" & levels != "-"]
phydat <- phyDat(data_matrix, type = "USER", levels = levels)

# Compute parsimony score (tree length)
fitch_score <- parsimony(tree, phydat, method = "fitch")

# Compute CI and RI
ci <- CI(tree, phydat)
ri <- RI(tree, phydat)

cat("=== KEY VALUES FOR PYTHON TESTS ===\n")
cat(sprintf("Number of taxa: %d\n", length(tree$tip.label)))
cat(sprintf("Number of characters: %d\n", ncol(data_matrix)))
cat(sprintf("Fitch parsimony score (tree length): %d\n", fitch_score))
cat(sprintf("Consistency index (CI): %.6f\n", ci))
cat(sprintf("Retention index (RI): %.6f\n", ri))

# Per-character parsimony scores
cat("\n=== PER-CHARACTER PARSIMONY SCORES ===\n")
for (i in seq_len(ncol(data_matrix))) {
  single_data <- phyDat(data_matrix[, i, drop = FALSE], type = "USER", levels = levels)
  score_i <- parsimony(tree, single_data, method = "fitch")
  cat(sprintf("char%d: %d steps\n", i - 1, score_i))
}

# ACCTRAN ancestral reconstruction
cat("\n=== ANCESTRAL STATES (ACCTRAN) ===\n")
anc <- ancestral.pars(tree, phydat, type = "ACCTRAN")
n_tips <- length(tree$tip.label)
n_internal <- tree$Nnode
for (i in seq_len(n_internal)) {
  node_idx <- n_tips + i
  cat(sprintf("Node %d:", node_idx))
  node_data <- anc[[node_idx]]
  if (is.matrix(node_data)) {
    for (ch in seq_len(nrow(node_data))) {
      state <- levels[which.max(node_data[ch, ])]
      cat(sprintf(" %s", state))
    }
  }
  cat("\n")
}
