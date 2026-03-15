#!/usr/bin/env Rscript
# validate_kf_distance.R
#
# Compute KF (branch score) distance using phangorn::KF.dist()
# and compare against PhyKIT's kf_distance command.
#
# The KF distance (Kuhner & Felsenstein 1994) is:
#   KF = sqrt( sum_over_all_splits( (b1_i - b2_i)^2 ) )
# where b1_i and b2_i are the branch lengths for split i in
# each tree, and 0 is used for splits absent from a tree.
# Both internal and terminal (external) branches are included.
#
# Requires: ape, phangorn
#
# Usage:
#   cd tests/r_validation
#   Rscript validate_kf_distance.R

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})

# --- Test case 1: trees with same taxa, different topology ---
t0 <- read.tree("../sample_files/tree_simple.tre")
t1 <- read.tree("../sample_files/tree_simple_other_topology.tre")

shared <- intersect(t0$tip.label, t1$tip.label)
t0_pruned <- keep.tip(t0, shared)
t1_pruned <- keep.tip(t1, shared)

kf <- KF.dist(t0_pruned, t1_pruned)
total_bl <- sum(t0_pruned$edge.length) + sum(t1_pruned$edge.length)
normalized_kf <- as.numeric(kf) / total_bl

cat("=== TEST CASE 1: Same taxa, different topology ===\n")
cat(sprintf("Number of shared taxa: %d\n", length(shared)))
cat(sprintf("Total branch length tree 0: %.6f\n", sum(t0_pruned$edge.length)))
cat(sprintf("Total branch length tree 1: %.6f\n", sum(t1_pruned$edge.length)))
cat(sprintf("KF distance (plain): %.4f\n", kf))
cat(sprintf("KF distance (normalized): %.4f\n", normalized_kf))

# --- Test case 2: identical trees ---
kf_identical <- KF.dist(t0_pruned, t0_pruned)
cat("\n=== TEST CASE 2: Identical trees ===\n")
cat(sprintf("KF distance (plain): %.4f\n", kf_identical))

# --- Test case 3: trees with different taxa (pruning) ---
t2 <- read.tree("../sample_files/tree_simple_other_topology_incomplete_taxon_representation.tre")
shared2 <- intersect(t0$tip.label, t2$tip.label)
t0_p2 <- keep.tip(t0, shared2)
t2_p2 <- keep.tip(t2, shared2)
kf2 <- KF.dist(t0_p2, t2_p2)
cat("\n=== TEST CASE 3: Different taxa (pruned to shared) ===\n")
cat(sprintf("Number of shared taxa: %d\n", length(shared2)))
cat(sprintf("KF distance (plain): %.4f\n", kf2))

cat("\n=== KEY VALUES FOR PYTHON TESTS ===\n")
cat(sprintf("Test 1 plain KF: %.4f\n", kf))
cat(sprintf("Test 2 plain KF (identical): %.4f\n", kf_identical))
cat(sprintf("Test 3 plain KF (different taxa): %.4f\n", kf2))
