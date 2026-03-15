#!/usr/bin/env Rscript
# validate_pic.R
#
# Compute phylogenetically independent contrasts using ape::pic()
# and compare against PhyKIT's independent_contrasts command.
#
# PIC (Felsenstein 1985) computes n-1 standardized contrasts for
# n tips by postorder traversal:
#   contrast_i = (x_left - x_right) / sqrt(v_left + v_right)
#
# Note: ape::pic() requires a fully bifurcating tree. We resolve
# multifurcations with multi2di(). Different resolution orders
# may produce different individual contrast values at polytomy
# nodes, but the sum of squared contrasts is invariant.
#
# Requires: ape
#
# Usage:
#   cd tests/r_validation
#   Rscript validate_pic.R

suppressPackageStartupMessages({
  library(ape)
})

tree <- read.tree("../sample_files/tree_simple.tre")
tree <- multi2di(tree)
traits <- read.delim("../sample_files/tree_simple_traits.tsv",
                     comment.char = "#", header = FALSE)
trait_vec <- setNames(traits$V2, traits$V1)

contrasts <- pic(trait_vec, tree)

cat("=== KEY VALUES FOR PYTHON TESTS ===\n")
cat(sprintf("Number of taxa: %d\n", length(trait_vec)))
cat(sprintf("Number of contrasts: %d\n", length(contrasts)))
cat(sprintf("Sum of squared contrasts: %.6f\n", sum(contrasts^2)))
cat(sprintf("Mean absolute contrast: %.6f\n", mean(abs(contrasts))))
cat(sprintf("Variance of contrasts: %.6f\n", var(contrasts)))
cat("\nIndividual contrasts (sorted by node):\n")
for (i in seq_along(contrasts)) {
    cat(sprintf("  Node %d: %.6f\n", as.integer(names(contrasts)[i]),
                contrasts[i]))
}

# Contrasts that are invariant regardless of polytomy resolution:
# (from fully resolved parts of the tree)
cat("\n=== INVARIANT CONTRASTS (fully resolved clades) ===\n")
cat(sprintf("bear-raccoon contrast: %.6f\n",
            contrasts[names(contrasts) == "9"]))
cat(sprintf("sea_lion-seal contrast: %.6f\n",
            contrasts[names(contrasts) == "10" | names(contrasts) == "11"][1]))
