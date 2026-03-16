#!/usr/bin/env Rscript
# validate_parsimony.R
#
# Compute Fitch parsimony score using phangorn::parsimony()
# and compare against PhyKIT's parsimony_score command.
#
# Fitch parsimony (Fitch 1971) counts the minimum number of
# character state changes on a tree. The downpass algorithm
# assigns state sets at internal nodes: intersection of children
# if shared, union otherwise (adding 1 to score).
#
# Requires: ape, phangorn
#
# Usage:
#   cd tests/r_validation
#   Rscript validate_parsimony.R

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})

tree <- read.tree("../sample_files/tree_simple.tre")
tree <- multi2di(tree)
aln <- read.phyDat("../sample_files/tree_simple_alignment.fa", format = "fasta")

fitch_score <- parsimony(tree, aln, method = "fitch")

cat("=== KEY VALUES FOR PYTHON TESTS ===\n")
cat(sprintf("Number of taxa: %d\n", length(tree$tip.label)))
cat(sprintf("Alignment length: %d\n", ncol(as.character(aln))))
cat(sprintf("Fitch parsimony score: %d\n", fitch_score))
